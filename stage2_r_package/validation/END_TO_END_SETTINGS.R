# Reproduce the selected formal Module B KDR-CC G1 cohort.
reproduction_settings <- list(
  replicates = 1:100,
  workers = 32L,
  resume = TRUE,
  retry_failed = TRUE,
  prepare_results_after_fit = TRUE,
  require_complete_cohort = TRUE,
  progress_seconds = 20,
  output_directory = NULL,
  report = list(font = "Times New Roman", dpi = 300L, base_size = 14),
  agreement = list(log_or = 2e-4, standard_error = 2e-4,
    coefficient_relative = 1e-3, covariance_relative = 1e-3,
    penalty_relative = 1e-3, contribution_percentage_points = 0.02,
    p_value = 5e-4)
)
# N=20,000; L=20,000; M=144; kappa=0.0001; lambda=1.
# These statistical settings, all random streams and evaluation coordinates
# are read from reference/stage_one_G1 and remain fixed across reproduction.
