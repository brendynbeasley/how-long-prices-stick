# Technical notes

Companion to the [README](README.md). Every table the scripts produce is in
[results/FINDINGS.md](results/FINDINGS.md); a review of the original paper's methods is in
[METHODS_REVIEW.md](METHODS_REVIEW.md).

## The original study

Guillermo A. Calvo, "Staggered Prices in a Utility-Maximizing Framework," *Journal of Monetary Economics* 12(3),
1983, pp. 383–398. A continuous-time theory paper with no data. Firms may change price only when a random signal
arrives, at a constant hazard independent of the price's age, so the expected life of a price is the inverse of
the hazard. In the discrete-time form used since, a price survives a quarter with probability θ and the New
Keynesian Phillips curve slope is κ = (1 − θ)(1 − βθ)/θ.

The paper has no headline number. The benchmark empirical test is Eichenbaum and Fisher (2003, Table 1, panel A),
who estimate θ by GMM on 1959Q1–2001Q4 data with β = 0.99 and find θ between 0.84 and 0.91, "between six and 11
quarters" before reoptimizing. With the GDP deflator they report θ = 0.89 (0.03), J = 49.4 [0.001] with L = 0 and
θ = 0.91 (0.02), J = 13.2 [0.95] with L = 12.

## Headline results

| Test | Estimate | Uncertainty | N |
|---|---|---|---|
| θ, GDP deflator, 1959Q1–2001Q4, L = 0 (Eichenbaum-Fisher rebuilt) | 0.925 | s.e. 0.039; J = 48.3, p = 0.002 | 172 |
| θ, GDP deflator, 1959Q1–2019Q4, L = 0 | 0.952 | s.e. 0.027; J = 52.9, p = 0.001 | 244 |
| Implied price duration from that θ | 62 months | | |
| CPI-weighted mean duration, Bils-Klenow via Bryan-Meyer | 7.9 months | | 45 categories |
| Calvo κ, micro θ (0.59) vs. aggregate θ | 0.29 vs. 0.0030 | ratio 99 | |
| Category slope test, WLS of κ̂ on log Calvo κ | 0.039 | s.e. 0.022, p = 0.07 | 38 |
| Spearman ρ, κ̂ vs. Calvo κ | 0.26 | p = 0.12 | 38 |
| Timing test, quarters later per log-month of duration | 1.90 | s.e. 0.35, p < 0.001 | 38 |
| Spearman ρ, peak quarter vs. duration | 0.63 | p < 0.001 | 38 |

## Data

All files were downloaded on November 1, 2026 and are committed under `data/raw/`.

- `data/raw/fred/*.csv`, from FRED (`fredgraph.csv?id=`). GDPDEF (GDP deflator), IPDNBS (nonfarm business
  deflator), CPIAUCSL (CPI-U), PCECTPI (PCE price index), PRS85006173 (nonfarm business labor share), COMPNFB
  (nonfarm business hourly compensation), GDPC1 (real GDP), GDPPOT (CBO potential GDP), GS10 and TB3MS (Treasury
  yields) and PPIACO (PPI, all commodities).
- `data/raw/cpi_components_fred.csv`, 38 BLS CPI component indexes as republished on FRED (seasonally adjusted
  `CUSR...` series where FRED carries them, otherwise not adjusted `CUUR...`), monthly from 1987.
- `data/raw/bryan_meyer_2010_table1.csv`, transcribed from Bryan and Meyer (2010), Table 1: 45 CPI categories, the
  average months between price changes (Bils and Klenow 2004) and CPI relative importance, plus the FRED series
  used for each.

## Method

Rules fixed before seeing results:

- **Aggregate GMM.** Moment condition E[(θπ̂ₜ − βθπ̂ₜ₊₁ − (1 − θ)(1 − βθ)ŝₜ)Xₜ] = 0 with β = 0.99; π̂ is
  quarterly log inflation minus its sample mean; ŝ is the log labor share minus its sample mean. Instruments are
  Xₜ = {1, Zₜ, ..., Zₜ₋₃} with Z = (ŝ, quadratically detrended log real GDP, π, commodity-price inflation,
  10-year minus 3-month spread, hourly-compensation inflation), 25 in all. Monthly series enter as quarterly
  means. The quadratic trend and means are fitted within each sample.
- **Substitutes for their series.** Eichenbaum and Fisher used Haver data. FRED lacks several of their series,
  so commodity prices are the PPI for all commodities (PPIACO) in place of the Commodity Research Bureau index;
  the long rate is the 10-year yield (GS10) in place of a composite of bonds of ten years and longer; wage
  growth uses hourly compensation (COMPNFB) in place of total nonfarm business compensation; marginal cost is
  BLS's nonfarm business labor share (PRS85006173) in place of unit labor cost over the deflator; and the
  quadratic trend is fitted to log rather than level real GDP.
- **Weighting.** Two-step GMM iterated to convergence, with the uncentered Newey-West (Bartlett) long-run
  covariance at L = 0 and L = 12 lags. Standard errors from the same matrix; Hansen's J with 24 degrees of freedom.
  A centered version gave J statistics far larger than Eichenbaum and Fisher's at L = 12. I switched to the
  uncentered form (Hansen 1982) after seeing that comparison, and use it for every sample.
- **Samples.** 1959Q1–2001Q4 (their sample), 1959Q1–2019Q4, 1985Q1–2019Q4, and 1959Q1 to the latest quarter.
- **Micro θ.** θᵢ = exp(−3/Dᵢ) for a category whose prices last Dᵢ months, which assumes a constant hazard within
  the category. Basket figures use Bryan and Meyer's CPI weights.
- **Category slopes.** Quarterly inflation is 100 × Δlog of the quarterly-average index; non-adjusted series have
  quarter-of-year means removed. For each category, 2SLS of πᵢₜ − 0.99πᵢ,ₜ₊₁ on a constant and the CBO output gap
  xₜ, instrumented by xₜ₋₁, xₜ₋₂, πᵢ,ₜ₋₁, πᵢ,ₜ₋₂, over 1990Q1–2019Q4 (a category needs 60 quarters), Newey-West
  standard errors with 4 lags. Cross-category WLS of κ̂ᵢ on ln κᵢ^Calvo with CPI weights and HC1 errors. Group
  indexes are CPI-weighted averages of the categories in each group, defined only in quarters where every member
  has data; four series start in 1998Q2, so the two group regressions run over 1998Q4–2019Q4 (85 quarters).
  Eight categories with later-starting series enter the category regressions with 85 to 117 quarters.
- **Timing.** Four-quarter inflation from the quarterly series; the peak is the maximum over 2021Q1–2024Q4.
  Education peaks in 2024Q4, the last quarter of the window, so its peak is censored.
  WLS of the peak quarter (counted from 2021Q1) on ln Dᵢ with CPI weights and HC1 errors.
- **Dropped categories.** Car and truck rental, leased cars and trucks, processed fruits and vegetables, jewelry
  and watches, motor vehicle insurance, miscellaneous personal goods and motor vehicle fees have no index on FRED
  (4.1% of the basket). Fresh fruits and vegetables uses the all fruits and vegetables index (CUSR0000SAF113);
  recreation uses FRED's CPIRECSL; personal care products uses the non-adjusted index.

## Full test statistics

Eichenbaum and Fisher's table rebuilt, 1959Q1–2001Q4 (θ, s.e., J, p):

| Inflation | L | Theirs | Mine |
|---|---|---|---|
| GDP deflator | 0 | 0.89 (0.03), 49.4 [0.001] | 0.925 (0.039), 48.3 [0.002] |
| GDP deflator | 12 | 0.91 (0.02), 13.2 [0.95] | 0.977 (0.049), 11.0 [0.989] |
| NFB deflator | 0 | 0.86 (0.03), 41.1 [0.01] | 0.850 (0.030), 45.4 [0.005] |
| NFB deflator | 12 | 0.86 (0.02), 12.8 [0.96] | 0.695 (0.018), 12.8 [0.969] |
| CPI | 0 | 0.88 (0.05), 38.9 [0.02] | 0.873 (0.042), 36.2 [0.052] |
| CPI | 12 | 0.86 (0.02), 12.6 [0.96] | 0.921 (0.024), 12.4 [0.975] |
| PCE deflator | 0 | 0.87 (0.02), 44.8 [0.004] | 0.860 (0.025), 39.6 [0.024] |
| PCE deflator | 12 | 0.88 (0.02), 12.8 [0.96] | 0.918 (0.019), 12.2 [0.978] |

The L = 0 estimates are within about one of their standard errors of the published ones. The L = 12 estimates
are less stable, and the nonfarm business deflator at L = 12 lands at 0.695, well below theirs. All 32
estimates, including the later samples, are in `results/nkpc_gmm_estimates.csv`. With L = 0, θ ranges from 0.929
to 0.952 over 1959–2019 and from 0.931 to 0.968 over 1985–2019.

Group slopes, 1990Q1–2019Q4: flexible κ̂ = 0.093 (0.086), sticky κ̂ = −0.005 (0.006); CPI-weighted Calvo κ 10.18
vs. 0.112 (ratio 91); difference in κ̂ 0.098, p = 0.26 treating the groups as independent. First-stage F
statistics exceed 200 in every category.

The 2021–23 peaks: flexible group mean 5.0 quarters after 2021Q1, sticky group 8.2; peak four-quarter inflation
20.9% vs. 7.6% (CPI-weighted). The slope of 1.90 implies a category repricing every 11 months peaks 3.2 quarters
after one repricing every 2 months.

The estimating equations in display form:

```math
\begin{aligned}
&\mathbb{E}\Big[\big(\theta\hat{\pi}_t - \beta\theta\hat{\pi}_{t+1}
 - (1-\theta)(1-\beta\theta)\hat{s}_t\big)X_t\Big] = 0 \\
&\pi_{it} - \beta\pi_{i,t+1} = a_i + \kappa_i x_t + e_{it} \\
&\text{peak}_i = d_0 + d_1 \ln D_i + v_i
\end{aligned}
```

## Reproducing the results

```bash
pip install -r requirements.txt
python src/calvo_aggregate.py
python src/calvo_categories.py
```

Both scripts read the committed files in `data/raw/`; `--download` refreshes them from FRED. Outputs go to
`results/`.

## References

- Bils, M. & Klenow, P. J. (2004). Some evidence on the importance of sticky prices. *Journal of Political
  Economy* 112(5).
- Bryan, M. F. & Meyer, B. (2010). Are some prices in the CPI more forward looking than others? We think so.
  Federal Reserve Bank of Cleveland, *Economic Commentary* 2010-2.
- Calvo, G. A. (1983). Staggered prices in a utility-maximizing framework. *Journal of Monetary Economics* 12(3).
- Eichenbaum, M. & Fisher, J. D. M. (2003). Testing the Calvo model of sticky prices. *Economic Perspectives*,
  Federal Reserve Bank of Chicago, 2Q/2003.
- Galí, J. & Gertler, M. (1999). Inflation dynamics: a structural econometric analysis. *Journal of Monetary
  Economics*.
- Hansen, L. P. (1982). Large sample properties of generalized method of moments estimators. *Econometrica*.
- Nakamura, E. & Steinsson, J. (2008). Five facts about prices: a reevaluation of menu cost models. *Quarterly
  Journal of Economics* 123(4).
