list(action = "report", modules = c("C", "E"), seed = 27182818L, 
    workers = 32L, repetitions = list(A = 100L, B = 100L, C = 100L, 
        D = 100L, E = 100L), null_counts = c(8L, 16L, 40L, 80L
    ), resume = TRUE, retry_failed = FALSE, clear_stale_lock = FALSE, 
    common_fixed = list(), B_overrides = list(), C_overrides = list(), 
    D_overrides = list(), E_overrides = list(), E_settings = list(
        null_counts = c(8L, 16L, 40L, 80L), dependence = c("Original", 
        "DM", "DH"), include_baseline = TRUE, tuning_sets = 30000L, 
        validation_folds = 3L, kernel_candidates = structure(list(
            M = c(96L, 144L, 144L, 144L, 192L, 288L), kappa = c(1e-04, 
            3e-05, 1e-04, 3e-04, 1e-04, 1e-04)), class = "data.frame", row.names = c(NA, 
        -6L)), lambda = c(0.3, 1, 3), manual_candidates = NULL), 
    D_baseline_run = NULL, D_correlation_levels = structure(list(
        correlation_id = c("DM", "DH"), correlation_label = c("Moderate", 
        "High"), target_mean = c(0.6, 0.85), min_abs = c(0.45, 
        0.75), mean_tolerance = c(0.03, 0.03)), class = "data.frame", row.names = c(NA, 
    -2L)), A_overrides = list(), advanced_overrides = list(), 
    freeze = list(source_run = NULL, config_id = NULL, manual_overrides = list()), 
    report_run_directories = NULL, report_overrides = list(font = "Times New Roman", 
        dpi = 300L, width_inches = 12, base_cex = 1.12, scatter_points = 300L, 
        save_pointwise_tables = TRUE, C_null_counts = NULL, C_panels_per_page = 6L))
