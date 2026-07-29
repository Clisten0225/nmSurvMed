library(nmSurvMed)

# Parallel computation is optional.
# Use FALSE for sequential computation.
use_parallel <- TRUE

# Used only for parallel computation and must not exceed detectCores().
n_cores <- 10L

demo_data <- read.csv("nonmarkov_simulated_data.csv")
time <- seq(0, 2.5, by = 0.01)

# Select bandwidths, estimate direct and indirect effects, and calculate
# bootstrap confidence intervals.
fit <- nmSurvMed(
  data = demo_data,
  T1 = "T1_hat",
  D1 = "D1",
  T2 = "T2_hat",
  D2 = "D2",
  Z = "Z",
  time = time,
  bandwidth = "paper_recommend",
  scale = "both",
  folds = 5,
  bootstrap = 1000L,
  parallel = use_parallel,
  cores = n_cores,
  seed = 123,
  keep_bootstrap = FALSE,
  verbose = TRUE
)

print(fit)
plot_nmSurvMed(
  fit,
  file = "nonmarkov_simulated_effects.png"
)

# Compare all minimum, optimal, and maximum bandwidth combinations.
sensitivity_fit <- bandwidth_sensitivity(
  data = demo_data,
  T1 = "T1_hat",
  D1 = "D1",
  T2 = "T2_hat",
  D2 = "D2",
  Z = "Z",
  time = time,
  bandwidth = "paper_recommend",
  scale = "both",
  folds = 5,
  bootstrap = 1000L,
  parallel = use_parallel,
  cores = n_cores,
  seed = 123,
  keep_bootstrap = FALSE,
  verbose = TRUE
)

print(sensitivity_fit$bandwidth_settings)
plot_bandwidth_sensitivity(
  sensitivity_fit,
  file = "nonmarkov_simulated_bandwidth_sensitivity.png"
)

# Report consistently significant positive and negative intervals.
decision <- bandwidth_decision(sensitivity_fit)
print(decision$summary)
