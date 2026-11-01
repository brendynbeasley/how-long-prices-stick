"""Calvo (1983) on aggregate US data: how long do prices stay fixed, judged from inflation?

Calvo's firms get a chance to reset their price with a fixed probability 1 - theta each
quarter. Log-linearized, that gives the New Keynesian Phillips curve, which Eichenbaum and
Fisher (2003, Chicago Fed Economic Perspectives, Table 1, panel A) estimate by GMM from the
orthogonality condition

    E[ (theta*pi_t - beta*theta*pi_{t+1} - (1 - beta*theta)(1 - theta)*s_t) * X_t ] = 0

with beta = 0.99, s_t = log labor share minus its mean, and X_t = {1, Z_t, ..., Z_{t-3}}.
This script rebuilds their table on today's FRED vintage (1959Q1-2001Q4, four inflation
measures, L = 0 and L = 12 lags in the weighting matrix), then reruns it on later samples,
and compares the implied price duration with the micro evidence that Bryan and Meyer (2010)
compiled from Bils and Klenow (2004).

Usage (from the repository root):
    python src/calvo_aggregate.py            # uses data/raw/fred/*.csv, downloads any missing
    python src/calvo_aggregate.py --download # refresh every FRED series first
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
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
FRED_DIR = RAW / "fred"
OUT = ROOT / "results"
FRED_URL = "https://fred.stlouisfed.org/graph/fredgraph.csv?id="

BETA = 0.99            # discount factor fixed by Eichenbaum and Fisher (2003)
N_LAGS_Z = 4           # X_t = {1, Z_t, Z_{t-1}, Z_{t-2}, Z_{t-3}}, their instrument set A
HAC_LAGS = (0, 12)     # L = 0 imposes the model's serial-correlation structure; L = 12 is Gali-Gertler's

# Eichenbaum and Fisher (2003), Table 1, panel A: theta (s.e.), J statistic [p-value], 1959Q1-2001Q4
EF_TABLE1A = {
    ("GDP deflator", 0): (0.89, 0.03, 49.4, 0.001), ("GDP deflator", 12): (0.91, 0.02, 13.2, 0.95),
    ("NFB deflator", 0): (0.86, 0.03, 41.1, 0.01), ("NFB deflator", 12): (0.86, 0.02, 12.8, 0.96),
    ("CPI", 0): (0.88, 0.05, 38.9, 0.02), ("CPI", 12): (0.86, 0.02, 12.6, 0.96),
    ("PCE deflator", 0): (0.87, 0.02, 44.8, 0.004), ("PCE deflator", 12): (0.88, 0.02, 12.8, 0.96),
}

SERIES = {
    "GDPDEF": "GDP implicit price deflator (quarterly)",
    "IPDNBS": "Nonfarm business sector implicit price deflator (quarterly)",
    "CPIAUCSL": "CPI-U all items, seasonally adjusted (monthly)",
    "PCECTPI": "PCE chain-type price index (quarterly)",
    "PRS85006173": "Nonfarm business sector labor share (quarterly index)",
    "COMPNFB": "Nonfarm business sector hourly compensation (quarterly index)",
    "GDPC1": "Real GDP (quarterly)",
    "GS10": "10-year Treasury constant maturity yield (monthly)",
    "TB3MS": "3-month Treasury bill secondary market rate (monthly)",
    "PPIACO": "PPI all commodities (monthly), the commodity-price instrument",
}
INFLATION = {"GDP deflator": "GDPDEF", "NFB deflator": "IPDNBS", "CPI": "CPIAUCSL",
             "PCE deflator": "PCECTPI"}

SAMPLES = {
    "1959Q1-2001Q4 (Eichenbaum-Fisher)": ("1959Q1", "2001Q4"),
    "1959Q1-2019Q4": ("1959Q1", "2019Q4"),
    "1985Q1-2019Q4": ("1985Q1", "2019Q4"),
    "1959Q1-latest": ("1959Q1", None),
}


def fetch(sid, force=False):
    """Download one FRED series to data/raw/fred/<sid>.csv unless it is already there."""
    path = FRED_DIR / f"{sid}.csv"
    if path.exists() and not force:
        return path
    FRED_DIR.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        r = requests.get(FRED_URL + sid, timeout=60)
        if r.status_code == 200 and r.text.startswith("observation_date"):
            path.write_text(r.text)
            print(f"  downloaded {sid}: {len(r.text.splitlines()) - 1} observations")
            time.sleep(0.5)
            return path
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"FRED download failed for {sid} (status {r.status_code})")


def load_quarterly(sid):
    """Read a FRED csv and return a quarterly series (quarterly mean of monthly data)."""
    d = pd.read_csv(FRED_DIR / f"{sid}.csv", na_values=".")
    d.columns = ["date", "value"]
    d["date"] = pd.to_datetime(d["date"])
    s = d.dropna().set_index("date")["value"].astype(float)
    return s.groupby(s.index.to_period("Q")).mean()


def build_panel():
    """Quarterly panel of inflation measures, real marginal cost and the six instruments."""
    q = {sid: load_quarterly(sid) for sid in SERIES}
    for sid, s in q.items():
        print(f"Loaded {sid}: {s.index.min()}-{s.index.max()}, {len(s)} quarters ({SERIES[sid]})")
    p = pd.DataFrame(index=pd.period_range("1947Q1", max(s.index.max() for s in q.values()), freq="Q"))
    for name, sid in INFLATION.items():
        p[name] = np.log(q[sid]).diff()            # quarterly inflation, log points
    p["log_labor_share"] = np.log(q["PRS85006173"])
    p["wage_inflation"] = np.log(q["COMPNFB"]).diff()
    p["commodity_inflation"] = np.log(q["PPIACO"]).diff()
    p["spread"] = (q["GS10"] - q["TB3MS"]) / 100
    p["log_gdp"] = np.log(q["GDPC1"])
    return p


def sample_data(p, infl, start, end):
    """Assemble pi_t, pi_{t+1}, s_t and the instrument matrix for one sample and inflation measure."""
    d = p.copy()
    end = end or str(d.dropna(subset=[infl, "log_labor_share"]).index[-2])  # need pi_{t+1}
    win = (d.index >= pd.Period(start, "Q")) & (d.index <= pd.Period(end, "Q"))
    # Demean inflation and labor share, and quadratically detrend real GDP, within the sample
    # (Eichenbaum and Fisher: pi-hat is inflation less its sample average; s-hat is log share less its mean).
    mean_pi = d.loc[win, infl].mean()
    d["pi"] = d[infl] - mean_pi
    d["s"] = d["log_labor_share"] - d.loc[win, "log_labor_share"].mean()
    tt = np.arange(len(d), dtype=float)
    ok = win & d["log_gdp"].notna()
    coef = np.polyfit(tt[ok], d.loc[ok, "log_gdp"], 2)
    d["gdp_gap"] = d["log_gdp"] - np.polyval(coef, tt)
    d["pi_lead"] = d["pi"].shift(-1)
    z_cols = ["s", "gdp_gap", "pi", "commodity_inflation", "spread", "wage_inflation"]
    inst = [pd.Series(1.0, index=d.index, name="const")]
    for j in range(N_LAGS_Z):
        for c in z_cols:
            inst.append(d[c].shift(j).rename(f"{c}_l{j}"))
    X = pd.concat(inst, axis=1)
    keep = win & X.notna().all(axis=1) & d[["pi", "pi_lead", "s"]].notna().all(axis=1)
    return d.loc[keep, "pi"].values, d.loc[keep, "pi_lead"].values, d.loc[keep, "s"].values, \
        X.loc[keep].values, int(win.sum()), int(keep.sum()), mean_pi


def psi(theta, pi, pi_lead, s, beta=BETA):
    return theta * pi - beta * theta * pi_lead - (1 - beta * theta) * (1 - theta) * s


def long_run_cov(u, L):
    """Newey-West (Bartlett) long-run covariance of the moment contributions u (T x k).

    Uncentered, as in Hansen (1982): under the null the moments have mean zero. This is the
    choice that reproduces the size of Eichenbaum and Fisher's J statistics."""
    T = u.shape[0]
    S = u.T @ u / T
    for lag in range(1, L + 1):
        w = 1 - lag / (L + 1)
        G = u[lag:].T @ u[:-lag] / T
        S += w * (G + G.T)
    return S


def gmm_theta(pi, pi_lead, s, X, L):
    """Two-step (iterated to convergence) GMM estimate of theta with beta fixed."""
    T, k = X.shape

    def gbar(th):
        return (X * psi(th, pi, pi_lead, s)[:, None]).mean(axis=0)

    W = np.linalg.inv(X.T @ X / T)
    theta = 0.8
    for _ in range(50):
        res = optimize.minimize_scalar(lambda th: gbar(th) @ W @ gbar(th), bounds=(0.01, 0.999),
                                       method="bounded", options={"xatol": 1e-8})
        new = res.x
        S = long_run_cov(X * psi(new, pi, pi_lead, s)[:, None], L)
        W = np.linalg.pinv(S)
        if abs(new - theta) < 1e-7:
            break
        theta = new
    res = optimize.minimize_scalar(lambda th: gbar(th) @ W @ gbar(th), bounds=(0.01, 0.999),
                                   method="bounded", options={"xatol": 1e-9})
    theta = res.x
    g = gbar(theta)
    J = T * g @ W @ g
    h = 1e-6
    G = (gbar(theta + h) - gbar(theta - h)) / (2 * h)
    se = np.sqrt(1 / (T * G @ W @ G))
    df = k - 1
    return theta, se, J, df, 1 - stats.chi2.cdf(J, df)


def calvo_slope(theta, beta=BETA):
    return (1 - theta) * (1 - beta * theta) / theta


def survival_figure(theta_macro, bm):
    """Share of prices still unchanged k months after being set: macro theta vs. the CPI basket."""
    k = np.arange(0, 37)
    w = bm["relative_importance"] / bm["relative_importance"].sum()
    micro = np.array([(w * np.exp(-m / bm["months_between_changes"])).sum() for m in k])
    macro = theta_macro ** (k / 3)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(k, macro * 100, color="#B24A3B", lw=2.2,
            label=f"Calvo, aggregate estimate (\u03b8 = {theta_macro:.2f} a quarter)")
    ax.plot(k, micro * 100, color="#2F6690", lw=2.2,
            label="CPI basket, Bils-Klenow frequencies (Bryan-Meyer weights)")
    ax.set_xlabel("Months since the price was set")
    ax.set_ylabel("Prices not yet changed (%)")
    ax.set_xlim(0, 36)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("How long prices last: what inflation implies vs. what stores do")
    fig.tight_layout()
    fig.savefig(OUT / "price_survival.png", dpi=150)
    plt.close(fig)
    return pd.DataFrame({"months": k, "unchanged_macro_pct": macro * 100, "unchanged_micro_pct": micro * 100})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true", help="re-download every FRED series")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    for sid in SERIES:
        fetch(sid, force=args.download)
    p = build_panel()

    rows = []
    for sample, (start, end) in SAMPLES.items():
        for infl in INFLATION:
            pi, pi_lead, s, X, n_win, n_keep, mean_pi = sample_data(p, infl, start, end)
            for L in HAC_LAGS:
                th, se, J, df, pj = gmm_theta(pi, pi_lead, s, X, L)
                dur_q = 1 / (1 - th)
                dur_se = se / (1 - th) ** 2                    # delta method
                rows.append(dict(sample=sample, inflation=infl, L=L, n_quarters=n_keep,
                                 dropped_quarters=n_win - n_keep, theta=th, theta_se=se,
                                 J=J, J_df=df, J_p=pj, duration_quarters=dur_q,
                                 duration_quarters_se=dur_se, duration_months=3 * dur_q,
                                 slope_kappa=calvo_slope(th), mean_inflation_annual_pct=400 * mean_pi))
    est = pd.DataFrame(rows)
    est.to_csv(OUT / "nkpc_gmm_estimates.csv", index=False)
    print(f"\nGMM estimates: {len(est)} rows written to results/nkpc_gmm_estimates.csv")
    print(f"  quarters used per sample (GDP deflator, L=0): "
          + ", ".join(f"{r.sample}: {r.n_quarters} kept, {r.dropped_quarters} dropped"
                      for r in est[(est.inflation == "GDP deflator") & (est.L == 0)].itertuples()))

    # Side-by-side with Eichenbaum and Fisher's published Table 1, panel A
    base = est[est["sample"].str.startswith("1959Q1-2001Q4")]
    comp = []
    for (infl, L), (th, se, J, pj) in EF_TABLE1A.items():
        r = base[(base.inflation == infl) & (base.L == L)].iloc[0]
        comp.append(dict(inflation=infl, L=L, ef2003_theta=th, ef2003_se=se, ef2003_J=J, ef2003_J_p=pj,
                         rerun_theta=round(r.theta, 3), rerun_se=round(r.theta_se, 3),
                         rerun_J=round(r.J, 1), rerun_J_p=round(r.J_p, 4), rerun_J_df=r.J_df,
                         rerun_duration_quarters=round(r.duration_quarters, 1)))
    comp = pd.DataFrame(comp)
    comp.to_csv(OUT / "ef2003_table1_comparison.csv", index=False)
    print("\nEichenbaum-Fisher (2003) Table 1A vs. rerun on today's data:")
    print(comp.to_string(index=False))

    # Micro benchmark: the CPI basket's frequencies (Bils-Klenow 2004 via Bryan-Meyer 2010, Table 1)
    bm = pd.read_csv(RAW / "bryan_meyer_2010_table1.csv")
    print(f"\nLoaded bryan_meyer_2010_table1.csv: {len(bm)} CPI categories, "
          f"weights sum to {bm.relative_importance.sum():.1f}")
    w = bm["relative_importance"] / bm["relative_importance"].sum()
    srt = bm.sort_values("months_between_changes")
    cum = (srt["relative_importance"] / srt["relative_importance"].sum()).cumsum()
    median_months = srt.loc[cum >= 0.5, "months_between_changes"].iloc[0]
    # Quarterly probability that a price survives, averaged over the basket with CPI weights
    theta_micro = float((w * np.exp(-3 / bm["months_between_changes"])).sum())
    gdp = est[(est.inflation == "GDP deflator") & (est.L == 0)].set_index("sample")
    micro = pd.DataFrame([
        dict(measure="Weighted mean months between price changes, CPI basket",
             value=float((w * bm["months_between_changes"]).sum())),
        dict(measure="Weighted median months between price changes, CPI basket", value=float(median_months)),
        dict(measure="Basket-average quarterly probability a price is unchanged (micro theta)",
             value=theta_micro),
        dict(measure="Months implied by micro theta, 3/(1 - theta)", value=3 / (1 - theta_micro)),
        dict(measure="Calvo slope implied by micro theta", value=calvo_slope(theta_micro)),
        dict(measure="Share of CPI weight in categories changing price at least every 4.3 months (%)",
             value=float(100 * w[bm.months_between_changes <= 4.3].sum())),
        dict(measure="Aggregate theta, GDP deflator, 1959Q1-2019Q4, L=0",
             value=float(gdp.loc["1959Q1-2019Q4", "theta"])),
        dict(measure="Months implied by aggregate theta, 1959Q1-2019Q4",
             value=float(gdp.loc["1959Q1-2019Q4", "duration_months"])),
        dict(measure="Calvo slope implied by aggregate theta, 1959Q1-2019Q4",
             value=float(gdp.loc["1959Q1-2019Q4", "slope_kappa"])),
        dict(measure="Ratio of micro to aggregate Calvo slope",
             value=calvo_slope(theta_micro) / float(gdp.loc["1959Q1-2019Q4", "slope_kappa"])),
    ])
    micro.to_csv(OUT / "micro_vs_macro.csv", index=False)
    print(micro.to_string(index=False))

    # Worked examples of the two Calvo formulas at durations quoted in the write-up
    ex = pd.DataFrame({"months_between_changes": [0.7, 4.3, 10.7, 23.7]})
    ex["theta_quarterly"] = np.exp(-3 / ex["months_between_changes"])
    ex["expected_duration_quarters"] = 1 / (1 - ex["theta_quarterly"])
    ex["kappa"] = calvo_slope(ex["theta_quarterly"])
    ex.to_csv(OUT / "calvo_formula_examples.csv", index=False)
    print(ex.to_string(index=False))

    surv = survival_figure(float(gdp.loc["1959Q1-2019Q4", "theta"]), bm)
    surv.to_csv(OUT / "price_survival.csv", index=False)
    print("\nWrote results/price_survival.png and price_survival.csv")


if __name__ == "__main__":
    main()
