# nmSurvMed

`nmSurvMed` provides non-parametric causal mediation analysis for
illness-death data without requiring the Markov assumption. It is intended
for researchers studying how a binary exposure affects a terminal event
directly and indirectly through an intermediate event in multi-state
survival data.

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
  pointwise percentile confidence intervals. It is normally called
  automatically by `nmSurvMed()`.
- `bandwidth_sensitivity()` compares every combination of the minimum,
  optimal, and maximum bandwidths across the three counterfactual scenarios
  using fixed-bandwidth bootstrap intervals.
- `plot_bandwidth_sensitivity()` draws the four sensitivity effect panels.
- `bandwidth_decision()` classifies each time interval as `Positive`,
  `Negative`, `Not significant`, or `Mixed` across all bandwidth settings.

## Data format

The input must be a data frame with one row per participant. Column names are
specified in the function call and do not need to match the example names.

| Argument | Example column | Meaning |
|---|---|---|
| `T1` | `T1_hat` | Observed time of the intermediate event or its censoring time |
| `D1` | `D1` | Intermediate-event indicator: 1 = observed, 0 = censored |
| `T2` | `T2_hat` | Observed terminal-event or follow-up time |
| `D2` | `D2` | Terminal-event indicator: 1 = observed, 0 = censored |
| `Z` | `Z` | Binary exposure or treatment indicator |

The event indicators must contain only 0 and 1, times must be finite and
non-negative, and `T1` should not exceed `T2`.

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
  bandwidth = "paper_recommend",
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

The main outputs are:

| Output | Contents |
|---|---|
| `fit$effects$survival` | Direct and indirect survival effects, standard errors, and pointwise confidence intervals |
| `fit$effects$cumulative_hazard` | Direct and indirect cumulative-hazard effects with the same inference summaries |
| `fit$bandwidth` | Selected scenario-specific bandwidths |
| `fit$bandwidth_values` | The same selected bandwidths as a full-precision named numeric vector for programmatic use |
| `fit$bandwidth_scores` | Cross-validation scores for the candidate bandwidths |
| `fit$data_roles` | Mapping between Exposure, Mediator, Outcome, and the input columns |
| `fit$diagnostics` | Sample counts, event counts, risk-set checks, and analysis warnings |
| `fit$runtime` | Elapsed time for bandwidth selection, estimation, bootstrap, and the full analysis |

## Bandwidth sensitivity analysis

> **Computational warning:** This analysis can be extremely time-consuming.
> It runs a separate bootstrap analysis for every unique combination of the
> minimum, optimal, and maximum bandwidths across three counterfactual
> scenarios—up to 27 settings. Consider testing the workflow with a small
> `bootstrap` value before running the final analysis with 1,000 replicates.

For a fixed-bandwidth bootstrap sensitivity analysis:

```r
sensitivity_fit <- bandwidth_sensitivity(
  data = data,
  T1 = "T1_hat", D1 = "D1",
  T2 = "T2_hat", D2 = "D2", Z = "Z",
  time = time,
  bandwidth = "paper_recommend",
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

With `bandwidth = "paper_recommend"`, the candidate set is

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
