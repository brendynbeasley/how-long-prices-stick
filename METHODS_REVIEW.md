# Methods review: Calvo (1983)

**Paper.** Guillermo A. Calvo, "Staggered Prices in a Utility-Maximizing Framework," *Journal of Monetary
Economics* 12(3), September 1983, pp. 383–398. doi:10.1016/0304-3932(83)90060-0.

**Statistical focus.** The paper's one load-bearing assumption and the formal result built on it: that each
price setter receives a chance to change its price at random moments, with a constant probability per unit of
time that does not depend on how long the current price has stood, so that the age of prices follows an
exponential (in discrete time, geometric) distribution with mean 1/δ; and the implication, developed in the
later literature into the New Keynesian Phillips curve, that inflation depends on expected future inflation and
current real marginal cost with a slope fixed by that probability.

**Review standard.** Internal consistency of the formal results, explicitness of assumptions, and the
testability of the predictions, as is appropriate for a theory paper. Where the paper has later been taken to
data, the review notes what those tests can and cannot establish about the paper itself.

**What was read.** The full text of the article sits behind the publisher's paywall and was not available to
me. I examined the Crossref bibliographic record (citation, issue date, pages and the paper's 17-item reference
list); the OpenAlex record (author affiliation, Columbia University); the encyclopedia entries on Calvo contracts
and on Guillermo Calvo; and two secondary sources that restate the model in detail, Eichenbaum and Fisher,
"Testing the Calvo model of sticky prices," *Economic Perspectives* (Federal Reserve Bank of Chicago), 2Q/2003,
pp. 40–53, and Bryan and Meyer, "Are Some Prices in the CPI More Forward Looking than Others? We Think So,"
Federal Reserve Bank of Cleveland *Economic Commentary* 2010-2. Everything below about the paper's content is
drawn from those accounts and is flagged where it depends on them.

---

## Methods summary

The research design is a formal macroeconomic model, not an empirical study. According to the encyclopedia
entry and the later literature, the original model is written in continuous time. Price setters are firms with
some power over their own price (the reference list includes Barro's 1972 "A theory of monopolistic price
adjustment" and Sheshinski and Weiss's 1977 "Inflation and costs of price adjustment," the two menu-cost papers
of the time). Each firm may revise its price only when a random signal arrives; signals arrive with a constant
hazard δ, independently across firms and over time. By the law of large numbers a fraction δ of all prices is
revised per unit of time, and the expected life of a price is 1/δ. A firm that gets the signal sets its price
with an eye to the whole expected period the price will stand, so the reset price is a weighted average of
expected future conditions with exponentially declining weights. The paper places this pricing rule inside a
utility-maximizing (optimizing household) monetary economy, following Phelps (1978) and Taylor (1979, 1980),
both of which appear in its reference list, and studies the dynamics of the price level and output under
alternative monetary policies. The reference list also includes Coddington's text on differential equations and
Gantmacher's on matrices, consistent with an analysis of a linear differential-equation system and its
stability.

Most of the standard empirical headings do not apply, and I say so item by item.

Sampling strategy, sample size, inclusion and exclusion criteria: not applicable. There is no sample.

Independent and dependent variables: in the model, the money supply rule (and the interest rate) is the
exogenous policy variable; the price level, inflation and output are the endogenous outcomes. The key structural
parameter is the hazard δ.

Statistical methods: not applicable. The "statistics" of a theory paper are its formal results: the derivation
of the aggregate price level as an exponentially weighted average of past reset prices, the reset-price
condition, and the resulting dynamic system.

Reported effect sizes, confidence intervals, p-values: none. The paper reports no estimates. Whether it contains
numerical examples or a calibration could not be checked without the full text.

Missing data treatment: not applicable.

Stated limitations: could not be extracted, because I did not read the full text. The literature's own reading
is that the constant hazard is an analytical convenience rather than a description of firm behavior; Eichenbaum
and Fisher describe the Calvo mechanism as capturing "firms' response to various costs of changing prices" in a
reduced-form way.

## Statistics extraction table

| Statistic | Value and unit | Source location | What it measures | What it does not establish |
|---|---|---|---|---|
| Hazard of a price change | δ per unit of time, constant, independent of a price's age | Model setup, as restated in the encyclopedia entry and in Eichenbaum and Fisher (2003), p. 42 | The probability that a given firm may change its price in a short interval | That real firms behave this way; the hazard is assumed, not derived from costs or estimated |
| Expected life of a price | 1/δ units of time | Implied by the exponential distribution; Eichenbaum and Fisher (2003), p. 42 give the discrete form (1 − θ)⁻¹ | The mean waiting time to the next price change | The value of δ; the paper does not measure it |
| Aggregate price level | Exponentially weighted average of past reset prices, weights δe^{−δ(t−s)} | Continuous-time setup, as restated in secondary accounts | How staggering makes the price level sticky even though each reset price is forward looking | That the price level actually adjusts at that speed |
| Discrete-time inflation equation (later form) | θπ̂_t = βθE π̂_{t+1} + (1 − βθ)(1 − θ)ŝ_t | Eichenbaum and Fisher (2003), equation 5, p. 43; derived by later authors, not by Calvo in this form | Inflation as forward-looking, driven by real marginal cost, with slope set by θ | That this equation is in Calvo (1983); the discrete-time, marginal-cost form came later |
| Example duration in the later literature | θ = 0.75 a quarter implies reoptimizing once a year | Eichenbaum and Fisher (2003), p. 42 | An illustrative calibration | Anything about the data |
| Estimated θ, GDP deflator, 1959Q1–2001Q4 | 0.89 (s.e. 0.03) with L = 0; 0.91 (0.02) with L = 12 | Eichenbaum and Fisher (2003), Table 1, panel A | Average reoptimization every 9 to 11 quarters, if the Calvo restrictions hold | A test of the paper's internal logic; the J test rejects the model with L = 0 (J = 49.4, p = 0.001) |

## Interpretation notes

The paper is famous for a modeling device rather than a finding. Its contribution is a price-setting rule that
is easier to work with than the fixed-length contracts of Taylor and Phelps. The tractability comes from memorylessness. Because the hazard does
not depend on a price's age, the whole cross-section of past prices collapses into one state variable, the
aggregate price level, and its law of motion is a single first-order equation. That is why the device spread.

The convenience is also the testable part. A constant hazard predicts that the probability a price changes this
month is unrelated to how long it has been in place, and that the distribution of price spells is geometric. It
also predicts that the slope of the inflation equation depends on the hazard alone (given the discount factor
and, in later versions, the degree of strategic complementarity). Both are falsifiable with micro price data and
with aggregate inflation data, and both have been tested heavily since.

Neither the paper's logic nor its derivations are in dispute in anything I read. The open questions concern
whether the assumption is a good approximation, which is an empirical matter outside the paper.

## Unanswered questions

Whether the full text contains a numerical illustration, a calibration of δ, or a stated view on its empirical
magnitude could not be determined from the sources I examined.

The paper does not, on the secondary accounts, give the firm a reason to wait. Why the hazard should be constant,
rather than rising with the gap between a firm's price and its desired price (as in the menu-cost models it
cites), is left open, and the later state-dependent pricing literature (for example Dotsey, King and Wolman 1999,
cited by Eichenbaum and Fisher) exists to answer it.

How large δ is, and whether it is the same across goods and over time, is not something the model says. The
aggregate estimates in Eichenbaum and Fisher (θ between 0.84 and 0.91 a quarter, so prices fixed for six to 11
quarters) and the micro evidence summarized by Bryan and Meyer (half of 350 CPI categories changing price at
least every 4.3 months, from Bils and Klenow 2004) point to very different values, which the model alone cannot
reconcile.

## Items requiring expert verification

The statements above about the paper's continuous-time structure, its demand specification and its policy
experiments rest on secondary sources; a reader with access to the full text should confirm them, in particular
whether Calvo's firms face constant-elasticity demand of the kind later formalized with Dixit and Stiglitz
preferences, or a simpler linear demand.

Whether the discrete-time New Keynesian Phillips curve can fairly be attributed to Calvo (1983), as it often is,
or should be credited to later authors (the encyclopedia entry credits Roberts 1995 for the curve) is a question
of intellectual history that a specialist should settle.

The Eichenbaum and Fisher estimates reported here come from their 2003 summary article, which in turn draws on
their longer 2003 working paper; the instrument set and the treatment of the weighting matrix differ in the
details across those versions.

This repository's own rerun of the Eichenbaum and Fisher table on today's data uses revised national accounts,
a different labor-share vintage and my own GMM code, so its numbers are a fresh estimate on the same design and
should not be read as a reproduction of the 2003 results or as evidence about Calvo's derivations.
