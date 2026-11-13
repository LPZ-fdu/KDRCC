# Immutable study defaults. Use 05_submit.R for run-specific settings.
# Counts N and L refer to complete matched sets, not individual rows.
simulation_config <- function(root) {
  list(
    schema = "kdrcc-formal-v1", seed = 27182818L,
    modules = "A",                       # After configuration freezing: c("B", "C")
    baseline = simulation_mechanism(), profile = "eight",
    noise = list(group_fraction = 0.20),
    configuration = list(
      fixed_file = file.path("config", "fixed_parameters.R"),
      # Named entries override the fixed configuration. Leave list() to inherit.
      # Example: B = list(N = 20000L, M = 96L, lambda = 10)
      overrides = list(B = list(), C = list(), D = list(), E = list())
    ),
    A = list(
      repetitions = 100L,
      reference = list(N = 20000L, L = 20000L, M = 144L, kappa = 1e-4, lambda = 10),
      grid = list(
        N = c(1000L,2000L,3000L,5000L,7500L,10000L,15000L,20000L,30000L,40000L,50000L),
        L = c(2000L,5000L,7500L,10000L,15000L,20000L,30000L,40000L,50000L),
        M = c(24L,48L,64L,80L,96L,112L,128L,144L,160L,192L),
        kappa = c(1e-6,1e-5,3e-5,1e-4,3e-4,1e-3),
        lambda = c(0,1,3,5,7,10,15,20),
        crossed_kappa = c(1e-5,3e-5,1e-4), crossed_lambda = c(3,5,7,10),
        crossed_L = c(10000L,20000L,30000L), crossed_M = c(96L,144L,192L)
      ),
      views = c("N", "L", "M", "kappa", "lambda", "kappa_lambda", "L_M"),
      # Optional data.frame with columns N,L,M,kappa,lambda replaces the grid.
      manual_design = NULL,
      recommendation = list(min_success_fraction = .95, min_repetitions = 80L,
        weights = c(joint_nrmse=.20,individual_nrmse=.20,
                    joint_nrms_bias=.10,individual_nrms_bias=.10,
                    joint_or_error=.05,individual_or_error=.05,
                    contribution_share_error=.10,interaction_miss=.10,
                    direction_miss=.05,direction_wrong=.05))
    ),
    B = list(repetitions = 100L),
    E = sim_E_defaults(),
    D=list(repetitions=100L, baseline_B_run=NULL,
      correlation_levels=data.frame(correlation_id=c("DM","DH"),
        correlation_label=c("Moderate","High"),target_mean=c(.60,.85),
        min_abs=c(.45,.75),mean_tolerance=c(.03,.03),stringsAsFactors=FALSE),
      calibration=list(reference_groups=10000L,calibration_groups=10000L,
        validation_groups=10000L,maximum_mix=.9999,root_tolerance=1e-5)),
    C = list(repetitions = 100L, null_counts = c(8L,16L,40L,80L),
             reuse_B_baseline = FALSE),
    kernel = list(coverage_fraction=.25, kmeans_starts=2L, kmeans_maxit=25L,
      initial_scale=.0625, start_factors=1, rho_bounds=c(.01,100),
      maxit=80L, continuations=1L, projected_gradient_tol=1e-5, penalty_rows=256L),
    nuisance = list(folds=3L, extra_degree=2L),
    comparator = list(knots=c(.25,.5,.75), boundary=c(.005,.995)),
    solver = list(maxit=70L, tolerance=1e-7, step_tolerance=1e-6, max_backtrack=20L),
    evaluation = list(groups=5000L, function_points=1000L, coverage_points=1000L,
      quantiles=sort(unique(c(seq(.1,.9,length.out=21),.25,.5,.75))),
      confidence=.95, zero_tolerance=1e-12),
    interpretation = list(contribution_points=1000L, chunk_points=128L,
      single_points=41L, interaction_points=9L,
      alpha=.05, tail_accuracy=1e-8, tail_max_terms=100000L,
      tail_retry_terms=1000000L, joint_OR_quantiles=c(.25,.75,.90)),
    selection = list(alpha=.05, adjust="holm", evaluation_points=1000L,
      accuracy=1e-8, max_terms=100000L, retry_terms=1000000L,
      # Fallback is used only if numerical integration fails or is ambiguous.
      fallback=TRUE, fallback_min_draws=199999L, expected_tail_exceedances=200L,
      fallback_batch=10000L),
    execution = list(workers=32L, resume=TRUE, retry_failed=FALSE,
      raw_memory_MB=256, feature_memory_MB=256,
      block_groups=2000L, poll_seconds=1, progress_seconds=30,
      worker_status_seconds=5, keep_work_cache=FALSE,
      # NULL uses an owned directory inside this project's results directory.
      cache_parent=NULL, clear_stale_lock=FALSE),
    resources = list(enabled=TRUE, sample_seconds=1, disk_seconds=5),
    report = list(font="Times New Roman", dpi=300L, width_inches=12,
      base_cex=1.12,
      save_pointwise_tables=TRUE, scatter_points=300L,
      C_null_counts=NULL, C_panels_per_page=6L)
  )
}
