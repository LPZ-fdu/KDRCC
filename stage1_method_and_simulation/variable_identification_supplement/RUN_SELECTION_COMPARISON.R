selection_arguments <- function(args) {
  ans <- list(workers = 32L, report_only = FALSE)
  while (length(args)) {
    flag <- args[1L]
    args <- args[-1L]
    if (flag == "--report-only") {
      ans$report_only <- TRUE
      next
    }
    fields <- c("--original-work" = "original", "--output" = "output",
                "--workers" = "workers", "--python" = "python")
    if (!flag %in% names(fields) || !length(args) || startsWith(args[1L], "--")) {
      stop("Invalid submission option: ", flag, call. = FALSE)
    }
    ans[[fields[[flag]]]] <- args[1L]
    args <- args[-1L]
  }
  if (is.null(ans$original) || is.null(ans$output)) {
    stop("Required arguments: --original-work and --output", call. = FALSE)
  }
  workers <- suppressWarnings(as.numeric(ans$workers))
  if (!is.finite(workers) || workers < 1 || workers != floor(workers) ||
      workers > .Machine$integer.max) stop("Invalid worker count", call. = FALSE)
  ans$workers <- as.integer(workers)
  ans
}

selection_main <- function(root, args) {
  original <- normalizePath(args$original, winslash = "/", mustWork = TRUE)
  dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
  work <- normalizePath(args$output, winslash = "/", mustWork = TRUE)
  overlap <- function(a, b) {
    if (.Platform$OS.type == "windows") {
      a <- tolower(a)
      b <- tolower(b)
    }
    a == b || startsWith(a, paste0(b, "/")) || startsWith(b, paste0(a, "/"))
  }
  if (overlap(work, original) || overlap(work, root)) {
    stop("The input, supplement source, and working directories must be separate", call. = FALSE)
  }
  available <- if (is.null(args$python)) unname(Sys.which(c("python3", "python"))) else args$python
  suitable <- function(executable) {
    status <- suppressWarnings(tryCatch(
      system2(executable, c("-c", shQuote("import sys; sys.exit(sys.version_info < (3, 6))")),
              stdout = TRUE, stderr = TRUE), error = function(e) NULL))
    !is.null(status) && is.null(attr(status, "status"))
  }
  python <- Filter(suitable, unique(available[nzchar(available)]))[1L]
  if (is.na(python)) stop("Python 3.6 or later is required for the delivery archive")
  if (!requireNamespace("CompQuadForm", quietly = TRUE)) stop("CompQuadForm is required")
  source(file.path(root, "R", "selection.R"), encoding = "UTF-8")
  source(file.path(root, "R", "report.R"), encoding = "UTF-8")
  source_root <- file.path(original, "results", selection_run_ids[["C"]], "source_snapshot")
  source(file.path(source_root, "load.R"), encoding = "UTF-8")
  load_kdrcc_simulation(source_root)
  runtime <- file.path(work, "runtime")
  dir.create(runtime, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(runtime, "run.log")
  note <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""), "\n")
    cat(line)
    cat(line, file = logfile, append = TRUE)
    flush.console()
  }
  context <- selection_context(original, runtime)
  context$source_root <- source_root
  context$supplement_root <- root
  runs <- unlist(context$runs, use.names = FALSE)
  plans <- file.path(runs, "plan.rds")
  configs <- file.path(runs, "config.rds")
  eval_files <- unlist(Map(function(run, entries) {
    file.path(run, "evaluation", paste0(names(entries), ".rds"))
  }, runs, context$evaluation), use.names = FALSE)
  engine_files <- sort(list.files(file.path(source_root, "R"), pattern = "[.]R$", full.names = TRUE))
  hashes <- c(tools::md5sum(c(plans, configs, eval_files)), tools::md5sum(engine_files),
              tools::md5sum(file.path(root, "R", "selection.R")))
  identity <- list(plan = context$plan, hashes = unname(hashes),
                   R = as.character(getRversion()), platform = R.version$platform,
                   CompQuadForm = as.character(utils::packageVersion("CompQuadForm")))
  identity_file <- file.path(runtime, "checkpoint_identity.rds")
  previous <- sim_read(identity_file)
  if (!is.null(previous) && !identical(previous, identity)) {
    stop("The working directory belongs to another source, input, or runtime configuration")
  }
  if (is.null(previous)) kdr_save(identity, identity_file)
  kdr_csv(context$plan, file.path(runtime, "task_plan.csv"))
  kdr_csv(data.frame(module = names(selection_run_ids), input_run = unname(selection_run_ids),
    plan_md5 = unname(tools::md5sum(plans)), configuration_md5 = unname(tools::md5sum(configs))),
    file.path(runtime, "original_run_references.csv"))
  kdr_csv(data.frame(field = c("R", "Platform", "CompQuadForm", "PSOCK workers", "Selection alpha", "Correction"),
    value = c(identity$R, identity$platform, identity$CompQuadForm, args$workers, "0.05", "Holm")),
    file.path(runtime, "environment.csv"))
  note("Variable identification comparison | modules B, C and D")
  if (!args$report_only) selection_execute(context, args$workers, note)
  paths <- file.path(runtime, "checkpoints", paste0(context$plan$id, ".rds"))
  if (!all(file.exists(paths))) stop("The comparison is incomplete; run the estimation stage to continue")
  note("Preparing publication figures, tables, and captions.")
  selection_report(context, file.path(work, "analysis"))
  note("Analysis complete. Packaging the publication report and detailed records.")
  code <- system2(python, c(shQuote(file.path(root, "tools", "package_results.py")), shQuote(work)))
  if (!identical(as.integer(code), 0L)) stop("Delivery packaging failed")
  cat("Delivery ready.\n")
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Submit this file with Rscript")
  root <- dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
  selection_main(root, selection_arguments(commandArgs(TRUE)))
}
