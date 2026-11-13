# Installation settings. Run with source() or Rscript --vanilla.
settings <- list(
  library = NULL,                 # NULL: use R's user library
  archive = "auto",              # auto, source, or binary
  install_dependencies = TRUE,
  repository = "https://cloud.r-project.org",
  output = "installation_test"
)

install_and_test_kdrcc <- function(settings) {
  script <- NULL
  for (fr in sys.frames()) if (!is.null(fr$ofile)) script <- fr$ofile
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (is.null(script) && length(args)) script <- sub("^--file=", "", args[1])
  root <- if (is.null(script)) getwd() else dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
  if (getRversion() < "3.6.3") stop("Install R 3.6.3 or a newer release.")
  if (is.null(settings$library)) {
    lib <- Sys.getenv("R_LIBS_USER")
    if (!nzchar(lib)) lib <- file.path(root, "R_library")
    lib <- path.expand(strsplit(lib, .Platform$path.sep, fixed = TRUE)[[1]][1])
  } else lib <- path.expand(settings$library)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  out <- file.path(root, settings$output)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
  dep_ok <- function() requireNamespace("CompQuadForm", quietly = TRUE) && packageVersion("CompQuadForm") >= "1.4.4"
  win44 <- .Platform$OS.type == "windows" && getRversion() >= "4.4" && getRversion() < "4.5"
  archives <- file.path(root, "installers")
  if (!dep_ok()) {
    if (!settings$install_dependencies) stop("Install CompQuadForm >= 1.4.4 before continuing.")
    bin <- file.path(archives, "CompQuadForm_1.4.4.zip")
    if (win44 && file.exists(bin)) {
      msg("Installing bundled CompQuadForm binary")
      try(install.packages(bin, repos = NULL, type = "win.binary", lib = lib), silent = TRUE)
    }
    if (!dep_ok()) {
      msg("Installing CompQuadForm from", settings$repository)
      try(install.packages("CompQuadForm", repos = settings$repository, lib = lib), silent = TRUE)
    }
    if (!dep_ok()) {
      msg("Installing bundled CompQuadForm source; platform compilation tools are used")
      install.packages(file.path(archives, "CompQuadForm_1.4.4.tar.gz"), repos = NULL, type = "source", lib = lib)
    }
    if (!dep_ok()) stop("CompQuadForm installation requires the platform's R development tools or a compatible repository binary.")
  }
  if (!settings$archive %in% c("auto", "source", "binary")) stop("archive must be auto, source or binary.")
  binary <- settings$archive == "binary" || (settings$archive == "auto" && win44)
  if (binary && !win44) stop("The bundled binary targets Windows R 4.4; select archive='source' for this installation.")
  archive <- file.path(archives, if (binary) "kdrcc_0.1.0.zip" else "kdrcc_0.1.0.tar.gz")
  if (!file.exists(archive)) stop("Copy the complete release, including installers, to this computer.")
  if ("kdrcc" %in% loadedNamespaces()) stop("Restart R before installing kdrcc in this session.")
  msg("Installing", basename(archive), "to", lib)
  install.packages(archive, repos = NULL, type = if (binary) "win.binary" else "source", lib = lib)
  library("kdrcc", character.only = TRUE, lib.loc = lib)
  checks <- list()
  check <- function(name, pass) {
    checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(pass), stringsAsFactors = FALSE)
    write.csv(do.call(rbind, checks), file.path(out, "installation_checks.csv"), row.names = FALSE)
    if (!isTRUE(pass)) stop("Installation test failed: ", name)
    msg("PASS", name)
  }
  check("Installed package version", as.character(packageVersion("kdrcc")) == "0.1.0")
  d <- kdrcc::kdrcc_example(800L, seed = 702L)
  parts <- kdrcc::kdrcc_split(d, "set_id", learning = 250L, seed = 703L)
  control <- kdrcc::kdrcc_control(centers = 10L, kappa = 1e-3, lambda = 1, kernel_tolerance = 1e-4)
  fit <- kdrcc::kdrcc(parts$estimation, exposures = c("X1", "X2", "X3", "X4"),
    adjust = ~ W + I(W^2), learning_data = parts$learning,
    nuisance = kdrcc::kdrcc_nuisance(g_df = 2L, kernel_rank = 4L, f_df = 2L), control = control)
  check("Kernel learning and final estimation", fit$diagnostics$converged)
  check("Full coefficient covariance", fit$covariance$status == "available" && all(is.finite(vcov(fit))))
  p <- predict(fit, parts$estimation[1:12, ], type = "or")
  check("Effect prediction and intervals", all(is.finite(as.matrix(p))))
  vars <- kdrcc::kdrcc_variables(fit)
  check("Variable participation tests", nrow(vars$data) == 4L && all(is.finite(vars$data$p_adjusted)))
  interaction <- kdrcc::kdrcc_interaction(fit, matrix(c("X1", "X3"), 1, 2), points = 5L)
  check("Four-point interaction tests", all(is.finite(interaction$overall$p_adjusted)))
  contribution <- kdrcc::kdrcc_contribution(fit, variables = c("X1", "X2", "X3"))
  check("Exact contribution allocation", contribution$identity_error < 1e-8)
  saveRDS(fit, file.path(out, "example_model.rds"), version = 2)
  check("Saved model prediction", identical(p, predict(readRDS(file.path(out, "example_model.rds")), parts$estimation[1:12, ], type = "or")))
  write.csv(p, file.path(out, "predictions.csv"), row.names = FALSE)
  write.csv(vars$data, file.path(out, "variable_tests.csv"), row.names = FALSE)
  write.csv(contribution$data, file.path(out, "contribution_shares.csv"), row.names = FALSE)
  write.csv(interaction$overall, file.path(out, "interaction_tests.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()), capture.output(print(fit))), file.path(out, "session_info.txt"))
  png(file.path(out, "joint_effect.png"), width = 1800, height = 1200, res = 180)
  plot(kdrcc::kdrcc_joint(fit));dev.off()
  png(file.path(out, "single_effects.png"), width = 2400, height = 1600, res = 180)
  plot(kdrcc::kdrcc_single(fit));dev.off()
  msg("Installation and workflow tests completed. Results:", out)
  invisible(list(library = lib, results = out, fit = fit))
}

installation_result <- install_and_test_kdrcc(settings)
