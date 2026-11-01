"""Calvo (1983) across the CPI basket: do goods that reprice often behave as the model says?

Calvo's slope of the Phillips curve depends on one number, the quarterly probability theta
that a price is left alone: kappa = (1 - theta)(1 - beta*theta) / theta. Bryan and Meyer
(2010, Cleveland Fed Economic Commentary 2010-2, Table 1) list the 45 published CPI components
with the average months between price changes that Bils and Klenow (2004) measured in the
BLS microdata. Under Calvo with a constant hazard, theta_i = exp(-3 / D_i) for a category
whose prices last D_i months on average.

Two tests, both on the BLS CPI component indexes as republished on FRED:
  1. Slope test. For each category, estimate the Phillips curve
         pi_it - beta*pi_i,t+1 = a_i + kappa_i * gap_t + e_it
     by two-stage least squares (instruments: lagged gap and lagged own inflation) and ask
     whether the estimated kappa_i line up with the Calvo kappa_i implied by D_i.
  2. Timing test. In the 2021-23 inflation, Calvo implies that categories with long price
     spells should have peaked later. Compare each category's peak quarter with D_i.

Usage (from the repository root):
    python src/calvo_categories.py            # uses committed data/raw files, downloads any missing
    python src/calvo_categories.py --download # refresh every series first
"""

import argparse
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
FRED_DIR = RAW / "fred"
OUT = ROOT / "results"
FRED_URL = "https://fred.stlouisfed.org/graph/fredgraph.csv?id="
CPI_START_YEAR = 1987  # first year kept from each CPI component series

BETA = 0.99                       # as in the aggregate script and Eichenbaum-Fisher (2003)
SLOPE_SAMPLE = ("1990Q1", "2019Q4")  # after the 1980s disinflation, before the pandemic; fixed in advance
MIN_QUARTERS = 60                 # a category needs 15 years inside the window to enter the slope test
HAC_LAGS = 4                      # Newey-West lags for the category regressions
PEAK_WINDOW = ("2021Q1", "2024Q4")  # search window for each category's peak four-quarter inflation
NAMED = ["Motor fuel", "Gas (piped) and electricity", "Used cars and trucks", "New vehicles",
         "Meats, poultry, fish, and eggs", "Food away from home", "Rent of primary residence",
         "Motor vehicle insurance", "Medical care services", "Personal care services"]


def fetch_cpi(series_ids, force=False):
    """Download the BLS CPI component indexes as republished on FRED (same BLS identifiers).

    Result: data/raw/cpi_components_fred.csv with columns series_id, year, period, value."""
    path = RAW / "cpi_components_fred.csv"
    if path.exists() and not force:
        return path
    rows = []
    for sid in series_ids:
        r = requests.get(FRED_URL + sid, timeout=60)
        r.raise_for_status()
        for line in r.text.splitlines()[1:]:
            date, value = line.split(",")
            if value != "." and int(date[:4]) >= CPI_START_YEAR:
                rows.append((sid, int(date[:4]), "M" + date[5:7], value))
        time.sleep(0.5)
    d = pd.DataFrame(rows, columns=["series_id", "year", "period", "value"])
    d.to_csv(path, index=False)
    print(f"  downloaded {d.series_id.nunique()} CPI series from FRED, {len(d)} monthly observations")
    return path


def bls_quarterly(bls, sid):
    """Quarterly mean of a monthly BLS index; quarters with fewer than three months are dropped."""
    d = bls[bls.series_id == sid].copy()
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d = d.dropna(subset=["value"])
    d["date"] = pd.to_datetime(d["year"].astype(str) + "-" + d["period"].str[1:] + "-01")
    s = d.set_index("date")["value"].sort_index()
    q = s.index.to_period("Q")
    counts = s.groupby(q).count()
    return s.groupby(q).mean()[counts == 3]


def fetch_plain(sid, force=False):
    path = FRED_DIR / f"{sid}.csv"
    if not path.exists() or force:
        r = requests.get(FRED_URL + sid, timeout=60)
        r.raise_for_status()
        path.write_text(r.text)
        print(f"  downloaded {sid}")
    return path


def quarterly(sid):
    d = pd.read_csv(FRED_DIR / f"{sid}.csv", na_values=".")
    d.columns = ["date", "value"]
    d["date"] = pd.to_datetime(d["date"])
    s = d.dropna().set_index("date")["value"].astype(float)
    counts = s.groupby(s.index.to_period("Q")).count()
    qs = s.groupby(s.index.to_period("Q")).mean()
    return qs[counts == counts.max()] if counts.max() == 3 else qs  # drop partial quarters


def calvo_theta(months):
    """Quarterly probability of no price change for an average spell of `months` months."""
    return np.exp(-3.0 / months)


def calvo_kappa(theta, beta=BETA):
    return (1 - theta) * (1 - beta * theta) / theta


def slope_2sls(pi, gap):
    """2SLS of pi_t - beta*pi_{t+1} on gap_t; instruments gap_{t-1}, gap_{t-2}, pi_{t-1}, pi_{t-2}."""
    d = pd.DataFrame({"pi": pi, "gap": gap})
    d["y"] = d["pi"] - BETA * d["pi"].shift(-1)
    for j in (1, 2):
        d[f"gap_l{j}"] = d["gap"].shift(j)
        d[f"pi_l{j}"] = d["pi"].shift(j)
    d = d.dropna()
    Z = sm.add_constant(d[["gap_l1", "gap_l2", "pi_l1", "pi_l2"]])
    first = sm.OLS(d["gap"], Z).fit()
    X_hat = sm.add_constant(first.fittedvalues.rename("gap"))
    X = sm.add_constant(d["gap"])
    b = np.linalg.lstsq(X_hat.values, d["y"].values, rcond=None)[0]
    u = d["y"].values - X.values @ b                          # structural residuals
    T = len(d)
    A = np.linalg.inv(X_hat.values.T @ X_hat.values)
    m = X_hat.values * u[:, None]
    S = m.T @ m
    for lag in range(1, HAC_LAGS + 1):
        w = 1 - lag / (HAC_LAGS + 1)
        G = m[lag:].T @ m[:-lag]
        S += w * (G + G.T)
    V = A @ S @ A
    first_F = first.f_test("gap_l1 = 0, gap_l2 = 0, pi_l1 = 0, pi_l2 = 0").fvalue
    return b[1], np.sqrt(V[1, 1]), T, float(np.squeeze(first_F))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true", help="re-download every series")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    FRED_DIR.mkdir(parents=True, exist_ok=True)

    bm = pd.read_csv(RAW / "bryan_meyer_2010_table1.csv")
    print(f"Loaded bryan_meyer_2010_table1.csv: {len(bm)} categories "
          f"({(bm.group == 'flexible').sum()} flexible, {(bm.group == 'sticky').sum()} sticky)")

    # Seven categories have no CPI series on FRED; they are dropped and listed here.
    unmatched = bm[bm.fred_series_id.isna()]
    print(f"Dropped {len(unmatched)} categories with no FRED series "
          f"({unmatched.relative_importance.sum():.1f}% of CPI weight): {', '.join(unmatched.category)}")
    bm = bm.dropna(subset=["fred_series_id"]).reset_index(drop=True)
    bls = pd.read_csv(fetch_cpi(list(bm.fred_series_id), force=args.download), dtype={"series_id": str})
    print(f"Loaded cpi_components_fred.csv: {bls.series_id.nunique()} series, "
          f"{bls.year.min()}-{bls.year.max()}, {len(bls)} monthly observations")

    for sid in ("GDPC1", "GDPPOT"):
        fetch_plain(sid, force=args.download)
    gap = 100 * (quarterly("GDPC1") / quarterly("GDPPOT") - 1)   # CBO output gap, percent
    print(f"Loaded output gap (GDPC1 vs. CBO GDPPOT): {gap.index.min()}-{gap.index.max()}")

    infl = {}
    for r in bm.itertuples():
        lvl = bls_quarterly(bls, r.fred_series_id)
        pi = 100 * np.log(lvl).diff()                      # percent per quarter
        if r.fred_series_id.startswith("CUUR"):           # not seasonally adjusted by BLS
            # remove quarter-of-year means, estimated over the whole history of the series
            pi = pi - pi.groupby(pi.index.quarter).transform("mean") + pi.mean()
        infl[r.category] = pi
    infl = pd.DataFrame(infl)
    infl.to_csv(OUT / "category_quarterly_inflation.csv", index_label="quarter")

    # ---- Test 1: slopes ----
    lo, hi = (pd.Period(x, "Q") for x in SLOPE_SAMPLE)
    out = []
    for r in bm.itertuples():
        pi = infl[r.category]
        win = pi[(pi.index >= lo - 2) & (pi.index <= hi + 1)]
        g = gap.reindex(win.index)
        n_avail = int(pi[(pi.index >= lo) & (pi.index <= hi)].notna().sum())
        theta = calvo_theta(r.months_between_changes)
        row = dict(category=r.category, group=r.group, months_between_changes=r.months_between_changes,
                   relative_importance=r.relative_importance, calvo_theta=theta,
                   calvo_reset_prob_quarter=1 - theta, calvo_kappa=calvo_kappa(theta),
                   quarters_in_window=n_avail)
        if n_avail >= MIN_QUARTERS:
            k, se, T, F = slope_2sls(win, g)
            row.update(kappa_hat=k, kappa_se=se, n_obs=T, first_stage_F=F)
        out.append(row)
    slopes = pd.DataFrame(out)
    kept = slopes.dropna(subset=["kappa_hat"])
    dropped = slopes[slopes.kappa_hat.isna()] if "kappa_hat" in slopes else slopes
    print(f"\nSlope test, {SLOPE_SAMPLE[0]}-{SLOPE_SAMPLE[1]}: {len(kept)} categories kept, "
          f"{len(dropped)} dropped for fewer than {MIN_QUARTERS} quarters: {', '.join(dropped.category)}")
    slopes.to_csv(OUT / "category_slopes.csv", index=False)

    # Cross-category comparison: does kappa_hat rise with Calvo's kappa?
    rho, rho_p = stats.spearmanr(kept["calvo_kappa"], kept["kappa_hat"])
    rho_f, rho_fp = stats.spearmanr(kept["months_between_changes"], kept["kappa_hat"])
    X = sm.add_constant(np.log(kept["calvo_kappa"]))
    wls = sm.WLS(kept["kappa_hat"], X, weights=kept["relative_importance"]).fit(cov_type="HC1")
    # Group indexes: CPI-weighted averages of category inflation, and their own slopes
    grp_rows = []
    for gname in ("flexible", "sticky"):
        sub = kept[kept.group == gname]
        w = sub.set_index("category")["relative_importance"]
        pi_g = (infl[w.index] * w).sum(axis=1, min_count=len(w)) / w.sum()
        win = pi_g[(pi_g.index >= lo - 2) & (pi_g.index <= hi + 1)]
        k, se, T, F = slope_2sls(win, gap.reindex(win.index))
        grp_rows.append(dict(group=gname, categories=len(sub), weight=w.sum(),
                             mean_months=float(np.average(sub.months_between_changes, weights=w)),
                             calvo_kappa_weighted=float(np.average(sub.calvo_kappa, weights=w)),
                             kappa_hat=k, kappa_se=se, n_obs=T, first_stage_F=F))
    grp = pd.DataFrame(grp_rows)
    ratio_pred = grp.loc[0, "calvo_kappa_weighted"] / grp.loc[1, "calvo_kappa_weighted"]
    ratio_est = grp.loc[0, "kappa_hat"] / grp.loc[1, "kappa_hat"]
    diff = grp.loc[0, "kappa_hat"] - grp.loc[1, "kappa_hat"]
    diff_se = np.hypot(grp.loc[0, "kappa_se"], grp.loc[1, "kappa_se"])   # ignores covariance
    grp.to_csv(OUT / "group_slopes.csv", index=False)

    summary = pd.DataFrame([
        dict(statistic="categories in slope test", value=len(kept), p_value=np.nan),
        dict(statistic="Spearman rho, kappa_hat vs Calvo kappa", value=rho, p_value=rho_p),
        dict(statistic="Spearman rho, kappa_hat vs months between changes", value=rho_f, p_value=rho_fp),
        dict(statistic="WLS slope of kappa_hat on log Calvo kappa", value=wls.params.iloc[1],
             p_value=wls.pvalues.iloc[1]),
        dict(statistic="WLS slope s.e. (HC1)", value=wls.bse.iloc[1], p_value=np.nan),
        dict(statistic="categories with kappa_hat > 0", value=int((kept.kappa_hat > 0).sum()), p_value=np.nan),
        dict(statistic="flexible/sticky ratio of Calvo kappa (CPI-weighted)", value=ratio_pred, p_value=np.nan),
        dict(statistic="flexible/sticky ratio of estimated kappa", value=ratio_est, p_value=np.nan),
        dict(statistic="flexible minus sticky estimated kappa", value=diff,
             p_value=2 * (1 - stats.norm.cdf(abs(diff / diff_se)))),
        dict(statistic="flexible minus sticky s.e. (independence)", value=diff_se, p_value=np.nan),
    ])
    summary.to_csv(OUT / "slope_test_summary.csv", index=False)
    print(grp.to_string(index=False))
    print(summary.to_string(index=False))

    # ---- Test 2: timing of the 2021-23 peak ----
    yoy = 100 * (np.exp(infl.rolling(4).sum() / 100) - 1)    # four-quarter inflation, percent
    plo, phi = (pd.Period(x, "Q") for x in PEAK_WINDOW)
    pk = []
    for r in bm.itertuples():
        s = yoy.loc[(yoy.index >= plo) & (yoy.index <= phi), r.category].dropna()
        pre = yoy.loc[(yoy.index >= pd.Period("2017Q1", "Q")) & (yoy.index <= pd.Period("2019Q4", "Q")),
                      r.category].mean()
        if len(s) == 0:
            continue
        q = s.idxmax()
        pk.append(dict(category=r.category, group=r.group, months_between_changes=r.months_between_changes,
                       relative_importance=r.relative_importance, peak_quarter=str(q),
                       quarters_after_2021Q1=(q - plo).n, peak_yoy_pct=s.max(),
                       avg_yoy_2017_2019_pct=pre))
    peaks = pd.DataFrame(pk)
    peaks.to_csv(OUT / "peak_timing.csv", index=False)
    rt, rt_p = stats.spearmanr(peaks.months_between_changes, peaks.quarters_after_2021Q1)
    Xp = sm.add_constant(np.log(peaks.months_between_changes))
    tw = sm.WLS(peaks.quarters_after_2021Q1, Xp, weights=peaks.relative_importance).fit(cov_type="HC1")
    by_group = peaks.groupby("group").apply(
        lambda d: pd.Series({"categories": len(d),
                             "weighted_mean_quarters_after_2021Q1": np.average(d.quarters_after_2021Q1,
                                                                               weights=d.relative_importance),
                             "weighted_mean_peak_yoy_pct": np.average(d.peak_yoy_pct,
                                                                      weights=d.relative_importance)}),
        include_groups=False).reset_index()
    timing = pd.DataFrame([
        dict(statistic="categories in timing test", value=len(peaks), p_value=np.nan),
        dict(statistic="Spearman rho, peak quarter vs months between changes", value=rt, p_value=rt_p),
        dict(statistic="WLS quarters later per log-month of stickiness", value=tw.params.iloc[1],
             p_value=tw.pvalues.iloc[1]),
        dict(statistic="WLS s.e. (HC1)", value=tw.bse.iloc[1], p_value=np.nan),
        dict(statistic="implied delay, 2-month vs 11-month category (quarters)",
             value=tw.params.iloc[1] * np.log(11 / 2), p_value=np.nan),
    ])
    timing.to_csv(OUT / "peak_timing_summary.csv", index=False)
    by_group.to_csv(OUT / "peak_timing_by_group.csv", index=False)
    print(f"\nTiming test: {len(peaks)} categories")
    print(timing.to_string(index=False))
    print(by_group.to_string(index=False))

    # ---- Named table ----
    named = slopes.merge(peaks[["category", "peak_quarter", "peak_yoy_pct", "avg_yoy_2017_2019_pct"]],
                         on="category", how="left")
    named = named[named.category.isin(NAMED)].sort_values("months_between_changes")
    named.to_csv(OUT / "named_categories.csv", index=False)
    print("\n" + named[["category", "months_between_changes", "calvo_reset_prob_quarter", "calvo_kappa",
                        "kappa_hat", "kappa_se", "peak_quarter", "peak_yoy_pct"]].to_string(index=False))

    # ---- Figures ----
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for gname, col in (("flexible", "#B24A3B"), ("sticky", "#2F6690")):
        d = kept[kept.group == gname]
        ax.errorbar(d.calvo_kappa, d.kappa_hat, yerr=1.96 * d.kappa_se, fmt="o", color=col, ms=5,
                    elinewidth=0.8, alpha=0.8, label=f"{gname.capitalize()}-price categories")
    slope_labels = {"Motor fuel": (-30, 8), "Fuel oil and other fuels": (6, 0),
                    "Public transportation": (6, 4), "Medical care services": (-40, -14)}
    for r in kept[kept.category.isin(slope_labels)].itertuples():
        ax.annotate(r.category, (r.calvo_kappa, r.kappa_hat), fontsize=7.5,
                    xytext=slope_labels[r.category], textcoords="offset points")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Phillips-curve slope Calvo predicts from price-change frequency (log scale)")
    ax.set_ylabel("Estimated slope on the output gap")
    ax.set_title(f"Category slopes, {SLOPE_SAMPLE[0]}-{SLOPE_SAMPLE[1]} (bars: 95% interval)")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "category_slopes.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for gname, col in (("flexible", "#B24A3B"), ("sticky", "#2F6690")):
        d = peaks[peaks.group == gname]
        ax.scatter(d.months_between_changes, d.quarters_after_2021Q1, s=20 + 12 * d.relative_importance,
                   color=col, alpha=0.75, label=f"{gname.capitalize()}-price categories")
    peak_labels = {"Motor fuel": (5, -2), "Used cars and trucks": (5, -10), "New vehicles": (-20, 8),
                   "Gas (piped) and electricity": (-60, 8), "Food away from home": (-95, -4),
                   "Rent of primary residence": (2, -16), "Medical care services": (5, -3),
                   "Personal care services": (-30, 7), "Education": (6, -3)}
    for r in peaks[peaks.category.isin(peak_labels)].itertuples():
        ax.annotate(r.category, (r.months_between_changes, r.quarters_after_2021Q1), fontsize=7.5,
                    xytext=peak_labels[r.category], textcoords="offset points")
    ax.set_xscale("log")
    ticks = [plo + i for i in range(0, (phi - plo).n + 1, 2)]
    ax.set_yticks([(t - plo).n for t in ticks])
    ax.set_yticklabels([str(t) for t in ticks])
    ax.set_xlabel("Average months between price changes (Bils-Klenow; log scale)")
    ax.set_ylabel("Quarter of peak four-quarter inflation")
    ax.set_title("When each CPI category peaked, 2021-2024")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "peak_timing.png", dpi=150)
    plt.close(fig)
    print("\nWrote results/category_slopes.png and results/peak_timing.png")


if __name__ == "__main__":
    main()
