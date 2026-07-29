# nmSurvMed

`nmSurvMed` implements causal mediation analysis for non-Markov
illness-death and multi-state survival models.

## Reference

Li-Sheng Zhuang, Jih-Chang Yu, and Yen-Tsung Huang (2026).
“Non-Parametric Mediation Analysis of Non-Markov Illness-Death Model.”
Accepted for publication in *Statistics in Medicine*.
[doi:10.1002/sim.70685](https://doi.org/10.1002/sim.70685).

## Installation

Install the package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "Clisten0225/nmSurvMed",
  subdir = "nmSurvMed"
)
```

Then load the package:

```r
library(nmSurvMed)
```

The package dependencies are installed automatically. The
`install.packages("remotes")` command only needs to be run if `remotes` is not
already installed.

## Platform support

Parallel execution is cross-platform. Linux and macOS use forked processes;
Windows uses a PSOCK cluster. Choose a value that does not exceed the detected
number of logical cores:

```r
parallel = TRUE
cores = 10
```

## Main functions

- `nmSurvMed()` runs bandwidth selection, effect estimation, and bootstrap
  inference in a single call.
- `estimate_cf_hazard()` estimates the counterfactual hazard-increment curve
  for a joint intervention `(z_a, z_b)`.
- `estimate_direct_effect()` estimates a time-varying natural direct effect.
- `estimate_indirect_effect()` estimates a time-varying natural indirect
  effect.
- `select_bandwidth()` selects scenario-specific kernel bandwidths by
  stratified K-fold cross-validation.
- `bootstrap_nm()` performs an event-stratified bootstrap and returns
  pointwise percentile confidence intervals.
- `bandwidth_sensitivity()` compares every combination of the minimum,
  optimal, and maximum bandwidths across the three counterfactual scenarios
  using fixed-bandwidth bootstrap intervals.
- `plot_bandwidth_sensitivity()` draws the four sensitivity effect panels.
- `bandwidth_decision()` classifies each time interval as `Positive`,
  `Negative`, `Not significant`, or `Mixed` across all bandwidth settings.

## Minimal workflow

```r
library(nmSurvMed)

data <- read.csv("../examples/nonmarkov_simulated_data.csv")
time <- seq(0, 2.5, by = 0.01)

fit <- nmSurvMed(
  data,
  T1 = "T1_hat",
  D1 = "D1",
  T2 = "T2_hat",
  D2 = "D2",
  Z = "Z",
  time = time,
  bandwidth = "paper_recommand",
  folds = 5,
  bootstrap = 1000L,
  parallel = TRUE,
  cores = 10L,
  seed = 123
)

fit$effects
fit$effects$survival
fit$effects$cumulative_hazard
fit$bandwidth
fit$bandwidth_values
fit$bandwidth_scores
fit$diagnostics

plot_nmSurvMed(
  fit,
  file = "nmSurvMed_effects.png"
)
```

For a fixed-bandwidth bootstrap sensitivity analysis:

```r
sensitivity_fit <- bandwidth_sensitivity(
  data = data,
  T1 = "T1_hat", D1 = "D1",
  T2 = "T2_hat", D2 = "D2", Z = "Z",
  time = time,
  bandwidth = "paper_recommand",
  folds = 5,
  bootstrap = 1000L,
  parallel = TRUE,
  cores = 10L,
  seed = 123
)

plot_bandwidth_sensitivity(sensitivity_fit)
decision <- bandwidth_decision(sensitivity_fit)
decision$summary
decision$intervals
```

For `bandwidth_decision()`:

- `Positive`: all bandwidth-specific confidence intervals are above zero.
- `Negative`: all bandwidth-specific confidence intervals are below zero.
- `Not significant`: all bandwidth-specific confidence intervals contain
  zero.
- `Mixed`: conclusions differ across bandwidth settings.
- `summary`: reports only consistently significant `Positive` and `Negative`
  intervals.
- `intervals` and `pointwise`: retain the complete `Not significant` and
  `Mixed` results.

Complete teaching examples and their data are kept outside the package in the
repository's top-level `examples/` directory.

Low-level functions use `T1_hat`, `T2_hat`, `D1`, `D2`, and `Z` by default.
Alternative names can be supplied through their `columns` argument.

With `bandwidth = "paper_recommand"`, the candidate set is

```math
\mathcal{H}
=
\left\{
C m^{-1/K}
:
C \in \{0.2, 0.5, 0.7\},
\quad
K \in \{2.5, 3, 3.5, 4, 4.5\}
\right\}.
```

where $m$ is the sample size. Selected bandwidths are held fixed during
bootstrap resampling.
