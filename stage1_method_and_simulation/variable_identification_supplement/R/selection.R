selection_run_ids <- c(B = "module_B_6fab5de7950a", C = "module_C_09bb6e6143c2",
                       D = "module_D_4f38f6e38982")
selection_methods <- c("KDR", "SplineCLR", "LinearCLR")

selection_plan <- function(runs) {
  out <- list()
  for (module in names(runs)) {
    tasks <- readRDS(file.path(runs[[module]], "plan.rds"))$tasks
    fields <- c("key", "kernel_key", "config_id", "method", "scenario_gf", "null_count",
                "p", "replicate", "N", "L", "M", "kappa", "lambda", "batch")
    d <- tasks[, fields, drop = FALSE]
    d$module <- module
    d$correlation <- if (module == "D") tasks$correlation_id else "D0"
    d$original_key <- d$key
    d$origin <- ifelse(d$method == "KDR" & d$scenario_gf == "G1",
                       "Archived test", "Archived fit")
    records <- file.path(runs[[module]], "fits", paste0(d$key, ".rds"))
    if (!all(file.exists(records))) stop("Incomplete input fit archive for module ", module)
    if (module == "C") {
      if (!identical(sort(unique(as.integer(d$null_count))), c(8L, 16L, 40L, 80L)) ||
          nrow(d) != 400L || any(d$method != "KDR" | d$scenario_gf != "G1")) {
        stop("Module C tasks do not match the frozen design")
      }
      add <- lapply(c("LinearCLR", "SplineCLR"), function(method) {
        extra <- d
        extra$key <- sub("KDR_G1$", paste0(method, "_correct"), extra$key)
        extra$method <- method
        extra$scenario_gf <- "correct"
        extra$origin <- "New fit"
        extra
      })
      d <- do.call(rbind, c(list(d), add))
    }
    d$id <- paste(module, d$key, sep = "_")
    out[[module]] <- d
  }
  plan <- do.call(rbind, out)
  rownames(plan) <- NULL
  stopifnot(!anyDuplicated(plan$id), nrow(plan) == 3600L, sum(plan$origin == "New fit") == 800L)
  plan
}

selection_context <- function(original_work, runtime) {
  runs <- as.list(setNames(file.path(original_work, "results", selection_run_ids), names(selection_run_ids)))
  for (module in names(runs)) {
    required <- c("config.rds", "plan.rds", "manifest.rds")
    if (!all(file.exists(file.path(runs[[module]], required)))) {
      stop("Incomplete input metadata for ", selection_run_ids[[module]])
    }
  }
  configs <- lapply(runs, function(run) readRDS(file.path(run, "config.rds")))
  if (!all(vapply(configs, function(cfg) identical(cfg$selection$adjust, "holm") &&
                 isTRUE(all.equal(cfg$selection$alpha, .05)), logical(1L)))) {
    stop("The frozen Holm testing configuration could not be verified")
  }
  evaluation <- list(B = list(), C = list(), D = list())
  evaluation$B[["p8"]] <- readRDS(file.path(runs$B, "evaluation", "p8.rds"))
  for (p in c(16L, 24L, 48L, 88L)) {
    name <- paste0("p", p)
    evaluation$C[[name]] <- readRDS(file.path(runs$C, "evaluation", paste0(name, ".rds")))
  }
  for (correlation in c("DM", "DH")) {
    evaluation$D[[paste0(correlation, "_p8")]] <- readRDS(
      file.path(runs$D, "evaluation", paste0(correlation, "_p8.rds")))
  }
  list(runs = runs, configs = configs, evaluation = evaluation, runtime = runtime,
       plan = selection_plan(runs))
}

selection_wald <- function(fit, specification, p, cfg) {
  if (is.null(fit$V)) stop("Comparator coefficient covariance is unavailable")
  # In the original additive basis theta = B beta. Blocks in theta correspond to exposures.
  B <- specification$B
  coefficients <- drop(B %*% fit$beta)
  covariance <- kdr_sym(B %*% fit$V %*% t(B))
  width <- length(coefficients) %/% p
  if (length(coefficients) != width * p || !width ||
      (specification$kind == "linear" && width != 1L)) stop("Invalid comparator basis blocks")
  scores <- vapply(seq_len(p), function(j) {
    idx <- (j - 1L) * width + seq_len(width)
    block <- covariance[idx, idx, drop = FALSE]
    if (any(!is.finite(block)) || any(diag(block) <= 0)) return(NA_real_)
    sd <- sqrt(diag(block))
    correlation <- block / outer(sd, sd)
    U <- tryCatch(chol(correlation), error = function(e) NULL)
    if (is.null(U)) return(NA_real_)
    sum(forwardsolve(t(U), coefficients[idx] / sd)^2)
  }, numeric(1L))
  p_value <- stats::pchisq(scores, df = width, lower.tail = FALSE)
  adjusted <- stats::p.adjust(p_value, method = cfg$selection$adjust, n = p)
  table <- data.frame(variable = kdr_names(p), true_active = seq_len(p) <= 8L,
    score = scores, p_value = p_value, algorithm = if (width == 1L) "Wald" else "Joint Wald",
    absolute_accuracy = NA_real_, ifault = NA_integer_, draws = NA_integer_,
    p_adjusted = adjusted, selected = adjusted < cfg$selection$alpha, df = width,
    stringsAsFactors = FALSE)
  valid <- all(is.finite(p_value))
  list(table = table, valid = valid, error = if (valid) "" else "Unavailable Wald covariance block")
}

selection_c_blocks <- function(run, cfg, job, role) {
  total <- if (role == "learning") job$L else job$N
  block_size <- cfg$execution$block_groups
  source <- file.path(run, "cache")
  lapply(seq_len(ceiling(total / block_size)), function(b) {
    file <- sim_pool_file(source, role, job$replicate, b)
    seed <- kdr_seed(cfg, "data_pool", role, job$replicate, b)
    archived <- file.exists(file)
    block <- if (archived) readRDS(file) else kdr_generate(block_size, 8L, cfg, seed)
    if (block$n != block_size || ncol(block$A) != 8L ||
        any(drop(group_sum(block$event, block$group)) != 1L)) stop("Invalid input data block")
    block <- sim_add_null(block, job$p, cfg, seed)
    block <- sim_subset(block, min(block_size, total - (b - 1L) * block_size))
    list(A = block$A, W = block$W, group = block$group, event = block$event,
         n = block$n, source = if (archived) "Archived block" else "Regenerated block")
  })
}

selection_c_fit <- function(job, run, cfg, input = NULL) {
  if (is.null(input)) {
    learning <- selection_c_blocks(run, cfg, job, "learning")
    estimation <- selection_c_blocks(run, cfg, job, "estimation")
  } else {
    learning <- input$learning
    estimation <- input$estimation
  }
  offsets <- c(0, head(cumsum(vapply(learning, function(x) x$n, integer(1L))), -1L))
  training <- list(A = do.call(rbind, lapply(learning, "[[", "A")),
    group = unlist(Map(function(x, offset) x$group + offset, learning, offsets), use.names = FALSE))
  kind <- if (job$method == "LinearCLR") "linear" else "spline"
  specification <- kdr_comparator_setup(training, cfg, kind)
  features <- lapply(estimation, function(x) {
    list(X = kdr_basis(x$A, specification), W = x$W, group = x$group, event = x$event)
  })
  each <- function(FUN) invisible(lapply(features, FUN))
  fit <- kdr_cl_model(each, specification, cfg, wrong = FALSE)
  list(fit = fit, specification = specification,
       block_sources = table(vapply(c(learning, estimation), "[[", character(1L), "source")))
}

selection_sanitize <- function(value, context) {
  for (path in c(context$runs, dirname(dirname(context$runs[[1L]])), context$runtime,
                 dirname(context$runtime), context$source_root, context$supplement_root)) {
    if (!is.null(path) && nzchar(path)) {
      value <- gsub(path, "<workspace>", value, fixed = TRUE)
      value <- gsub(chartr("/", "\\", path), "<workspace>", value, fixed = TRUE)
    }
  }
  # File diagnostics may include system-generated temporary names.
  gsub("([A-Za-z]:[/\\\\]|/)([^[:space:]'\"]+[/\\\\])+[^[:space:]'\"]*", "<path>", value)
}

selection_record <- function(job, context, input = NULL) {
  path <- file.path(context$runtime, "checkpoints", paste0(job$id, ".rds"))
  existing <- sim_read(path)
  if (!is.null(existing)) return(invisible(existing$valid))
  start <- proc.time()[3L]
  run <- context$runs[[job$module]]
  cfg <- sim_job_config(context$configs[[job$module]], job)
  warnings <- character()
  fit_ok <- FALSE
  new_fit <- NULL
  attempt <- withCallingHandlers(tryCatch({
    if (job$origin == "New fit") {
      new_fit <- selection_c_fit(job, run, cfg, input)
      fit_ok <- TRUE
      selection_wald(new_fit$fit, new_fit$specification, job$p, cfg)
    } else {
      original <- readRDS(file.path(run, "fits", paste0(job$original_key, ".rds")))
      fit_ok <- isTRUE(original$ok)
      if (!fit_ok) stop(original$error)
      if (job$origin == "Archived test") {
        if (is.null(original$selection)) stop("Archived variable test unavailable")
        original$selection
      } else if (job$method == "KDR") {
        evaluation_key <- if (job$module == "D") paste0(job$correlation, "_p", job$p) else paste0("p", job$p)
        ev <- context$evaluation[[job$module]][[evaluation_key]]
        spec <- readRDS(file.path(run, "kernels", paste0(job$kernel_key, ".rds")))
        kdr_selection(original$fit, spec, ev, cfg, original$inference$V, job)
      } else {
        spec <- readRDS(file.path(run, "comparators", paste0(job$kernel_key, "_", job$method, ".rds")))
        selection_wald(original$fit, spec, job$p, cfg)
      }
    }
  }, error = function(e) list(table = data.frame(), valid = FALSE, error = conditionMessage(e))),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (nrow(attempt$table) && !"df" %in% names(attempt$table)) attempt$table$df <- NA_integer_
  valid <- isTRUE(attempt$valid) && identical(as.character(attempt$table$variable), kdr_names(job$p)) &&
    all(!is.na(attempt$table$selected))
  record <- list(job = job, fit_ok = fit_ok, valid = valid, selection = attempt$table,
    error = selection_sanitize(attempt$error, context),
    warnings = unique(selection_sanitize(warnings, context)), new_fit = new_fit,
    elapsed_seconds = unname(proc.time()[3L] - start))
  kdr_save(record, path)
  invisible(valid)
}

selection_task <- function(group, context) {
  rows <- context$plan[group, , drop = FALSE]
  input <- NULL
  if (all(rows$origin == "New fit") && nrow(rows) == 2L) {
    missing <- !file.exists(file.path(context$runtime, "checkpoints", paste0(rows$id, ".rds")))
    if (any(missing)) {
      job <- as.list(rows[1L, , drop = FALSE])
      cfg <- sim_job_config(context$configs$C, job)
      input <- tryCatch(list(
        learning = selection_c_blocks(context$runs$C, cfg, job, "learning"),
        estimation = selection_c_blocks(context$runs$C, cfg, job, "estimation")),
        error = function(e) NULL)
    }
  }
  for (i in seq_len(nrow(rows))) selection_record(as.list(rows[i, , drop = FALSE]), context, input)
  invisible(TRUE)
}

selection_execute <- function(context, workers, note) {
  plan <- context$plan
  pending <- !file.exists(file.path(context$runtime, "checkpoints", paste0(plan$id, ".rds")))
  archived <- as.list(which(pending & plan$origin != "New fit"))
  fresh <- which(pending & plan$origin == "New fit")
  c_groups <- split(fresh, paste(plan$replicate[fresh], plan$p[fresh], sep = ":"))
  phases <- list("Archived records" = archived, "New comparator fits" = unname(c_groups))
  n_groups <- sum(lengths(phases))
  completed <- sum(!pending)
  note("Checkpoints available: ", completed, "/", nrow(plan), ". Workers: ", workers, ".")
  if (!n_groups) return(invisible(NULL))
  workers <- min(workers, n_groups)
  cl <- NULL
  if (workers > 1L) {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
               VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(source_root, supplement_root) {
      source(file.path(source_root, "load.R"), encoding = "UTF-8")
      load_kdrcc_simulation(source_root)
      source(file.path(supplement_root, "R", "selection.R"), encoding = "UTF-8")
      NULL
    }, context$source_root, context$supplement_root)
  }
  for (phase_name in names(phases)) {
    phase <- phases[[phase_name]]
    if (!length(phase)) next
    phase_total <- sum(lengths(phase))
    start <- proc.time()[3L]
    initial <- completed
    note(phase_name, " | ", phase_total, " pending records")
    width <- max(1L, workers * 2L)
    for (offset in seq.int(1L, length(phase), by = width)) {
      batch <- phase[seq.int(offset, min(length(phase), offset + width - 1L))]
      if (is.null(cl)) lapply(batch, selection_task, context = context)
      else parallel::parLapplyLB(cl, batch, selection_task, context = context)
      completed <- completed + sum(lengths(batch))
      elapsed <- proc.time()[3L] - start
      remaining <- elapsed * (phase_total - completed + initial) / (completed - initial)
      note("Completed ", completed, "/", nrow(plan), " | elapsed ", sim_seconds(elapsed),
           " | phase time remaining ", sim_seconds(remaining))
    }
  }
  invisible(NULL)
}
