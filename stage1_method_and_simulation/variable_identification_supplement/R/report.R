selection_collect <- function(context) {
  plan <- context$plan
  records <- lapply(plan$id, function(id) {
    record <- readRDS(file.path(context$runtime, "checkpoints", paste0(id, ".rds")))
    record$block_sources <- record$new_fit$block_sources
    record$new_fit <- NULL
    record
  })
  status <- plan
  status$fit_ok <- vapply(records, function(record) record$fit_ok, logical(1L))
  status$valid <- vapply(records, function(record) record$valid, logical(1L))
  status$elapsed_seconds <- vapply(records, function(record) record$elapsed_seconds, numeric(1L))
  status$warnings <- vapply(records, function(record) paste(record$warnings, collapse = " | "), character(1L))
  status$error <- vapply(records, function(record) record$error, character(1L))
  status$data_blocks <- vapply(records, function(record) {
    sources <- record$block_sources
    paste(paste(names(sources), as.integer(sources), sep = "="), collapse = "; ")
  }, character(1L))
  tests <- kdr_bind(lapply(seq_along(records), function(i) {
    data <- records[[i]]$selection
    if (!nrow(data)) return(NULL)
    cbind(plan[rep(i, nrow(data)), c("id", "method", "scenario_gf", "module", "correlation",
      "null_count", "p", "replicate", "origin")], data, stringsAsFactors = FALSE)
  }))
  per_run <- plan[, c("id", "module", "correlation", "null_count", "p", "replicate", "method", "scenario_gf")]
  per_run$valid <- status$valid
  for (exposure in exposure_names) per_run[[exposure]] <- NA_real_
  per_run$TPR <- per_run$FPR <- NA_real_
  for (i in which(status$valid)) {
    data <- records[[i]]$selection
    active <- data$true_active
    per_run[i, exposure_names] <- as.list(as.numeric(data$selected[match(exposure_names, data$variable)]))
    per_run$TPR[i] <- mean(data$selected[active])
    if (any(!active)) per_run$FPR[i] <- mean(data$selected[!active])
  }
  list(status = status, tests = tests, per_run = per_run)
}

selection_scenarios <- function(d) {
  parts <- list()
  for (module in c("B", "C", "D")) {
    for (scenario in if (module == "C") "G1" else paste0("G", 1:4)) {
      adjustment <- if (scenario %in% c("G1", "G2")) "correct" else "wrong"
      target <- ifelse(d$method == "KDR", scenario, adjustment)
      z <- d[d$module == module & d$scenario_gf == target, , drop = FALSE]
      z$scenario <- rep(scenario, nrow(z))
      parts[[length(parts) + 1L]] <- z
    }
  }
  out <- kdr_bind(parts)
  # The low-correlation D reference uses the same B repetitions and results.
  baseline <- out[out$module == "B", , drop = FALSE]
  baseline$module <- rep("D", nrow(baseline))
  out <- rbind(out, baseline)
  out$method <- factor(out$method, levels = selection_methods)
  out <- out[order(out$module, out$scenario, match(out$correlation, c("D0", "DM", "DH")),
                   out$null_count, out$method, out$replicate), , drop = FALSE]
  out$method <- as.character(out$method)
  out
}

selection_interval <- function(values, lower_bound = 0, upper_bound = 1, wilson = FALSE) {
  values <- values[is.finite(values)]
  n <- length(values)
  if (!n) return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, mcse = NA_real_))
  rate <- mean(values)
  mcse <- if (n > 1L) stats::sd(values) / sqrt(n) else NA_real_
  limits <- rate + c(-1, 1) * if (n > 1L) stats::qt(.975, n - 1L) * mcse else NA_real_
  if (wilson) {
    z <- stats::qnorm(.975)
    denom <- 1 + z^2 / n
    center <- (rate + z^2 / (2 * n)) / denom
    margin <- z * sqrt(rate * (1 - rate) / n + z^2 / (4 * n^2)) / denom
    limits <- center + c(-margin, margin)
  }
  100 * c(mean = rate, lower = max(lower_bound, limits[1L]),
          upper = min(upper_bound, limits[2L]), mcse = mcse)
}

selection_summary <- function(d) {
  fields <- c("module", "scenario", "correlation", "null_count", "p")
  overall <- specific <- paired <- list()
  keys <- do.call(paste, c(d[fields], sep = "\r"))
  for (i in split(seq_len(nrow(d)), factor(keys, levels = unique(keys)))) {
    z <- d[i, , drop = FALSE]
    setting <- z[1L, fields, drop = FALSE]
    for (method in selection_methods) {
      m <- z[z$method == method, , drop = FALSE]
      if (!nrow(m)) next
      tpr <- selection_interval(m$TPR)
      fpr <- selection_interval(m$FPR)
      overall[[length(overall) + 1L]] <- cbind(setting, data.frame(method = method,
        planned_repetitions = nrow(m), usable_repetitions = sum(m$valid),
        TPR_percent = unname(tpr["mean"]), TPR_lower = unname(tpr["lower"]),
        TPR_upper = unname(tpr["upper"]), TPR_MCSE = unname(tpr["mcse"]),
        FPR_percent = unname(fpr["mean"]), FPR_lower = unname(fpr["lower"]),
        FPR_upper = unname(fpr["upper"]), FPR_MCSE = unname(fpr["mcse"])))
      for (v in exposure_names) {
        a <- selection_interval(m[[v]], wilson = TRUE)
        specific[[length(specific) + 1L]] <- cbind(setting, data.frame(method = method,
          exposure = sim_label(v), exposure_id = v, usable_repetitions = sum(is.finite(m[[v]])),
          detected_repetitions = sum(m[[v]] == 1, na.rm = TRUE), identification_percent = unname(a["mean"]),
          lower = unname(a["lower"]), upper = unname(a["upper"])))
      }
    }
    metrics <- c("TPR", if (setting$null_count > 0L) "FPR", exposure_names)
    for (method in setdiff(selection_methods, "KDR")) {
      a <- z[z$method == "KDR", c("replicate", metrics), drop = FALSE]
      b <- z[z$method == method, c("replicate", metrics), drop = FALSE]
      m <- merge(a, b, by = "replicate", suffixes = c(".KDR", ".CLR"))
      for (metric in metrics) {
        diff <- m[[paste0(metric, ".KDR")]] - m[[paste0(metric, ".CLR")]]
        est <- selection_interval(diff, -1, 1)
        paired[[length(paired) + 1L]] <- cbind(setting, data.frame(
          comparator = sim_label(method), target = sim_label(metric), paired_repetitions = sum(is.finite(diff)),
          difference_percentage_points = unname(est["mean"]), lower = unname(est["lower"]),
          upper = unname(est["upper"]), MCSE = unname(est["mcse"])))
      }
    }
  }
  list(overall = kdr_bind(overall), specific = kdr_bind(specific), paired = kdr_bind(paired))
}

selection_plot <- function(path, width, height, draw) {
  font <- sim_report_font(list(font = "Times New Roman"))$font
  if (.Platform$OS.type == "windows" && !capabilities("cairo")) {
    grDevices::windowsFonts(SelectionTimes = grDevices::windowsFont("TT Times New Roman"))
    grDevices::png(path, width = width, height = height, units = "in", res = 300,
                   type = "windows", family = "SelectionTimes")
    font <- "SelectionTimes"
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = 300, type = "cairo")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(family = font)
  draw()
  invisible(font)
}

selection_figure_c <- function(d, path, by_exposure = FALSE) {
  levels <- c(8, 16, 40, 80)
  colors <- setNames(c("#0072B2", "#D55E00", "#009E73"), selection_methods)
  symbols <- setNames(c(16, 17, 15), selection_methods)
  panel <- function(data, rate, low, high, label, ylim, tag, ylab = label) {
    graphics::plot(c(.65, length(levels) + .35), ylim, type = "n", xaxt = "n",
                   xlab = "", ylab = ylab, cex.lab = 1.06)
    graphics::axis(1, at = seq_along(levels), labels = levels)
    graphics::abline(h = pretty(ylim), col = "#E4E7EA", lwd = .8)
    for (i in seq_along(selection_methods)) {
      method <- selection_methods[i]
      z <- data[data$method == method, , drop = FALSE]
      z <- z[order(z$null_count), , drop = FALSE]
      x <- match(z$null_count, levels) + (i - 2L) * .08
      good <- is.finite(z[[low]]) & is.finite(z[[high]])
      graphics::segments(x[good], z[[low]][good], x[good], z[[high]][good], col = colors[method], lwd = 1.1)
      graphics::lines(x, z[[rate]], col = colors[method], lwd = 1.8)
      graphics::points(x, z[[rate]], col = colors[method], pch = symbols[method], cex = .85)
    }
    graphics::mtext(tag, side = 3, adj = 0, line = .25, font = 2)
    if (by_exposure) graphics::mtext(label, side = 3, line = .25, font = 2)
  }
  selection_plot(path, if (by_exposure) 13.6 else 10.5, if (by_exposure) 7.7 else 5.0, function() {
    graphics::par(mfrow = if (by_exposure) c(2, 4) else c(1, 2), mar = c(2.9, 4.7, 1.65, 1.05),
                  oma = c(2.4, .35, 2.65, .35), mgp = c(2.8, .7, 0), tcl = -.23, bty = "l", las = 1)
    graphics::par(cex = if (by_exposure) .91 else .98)
    if (by_exposure) {
      for (i in seq_along(exposure_names)) {
        exposure <- exposure_names[i]
        panel(d[d$exposure_id == exposure, , drop = FALSE], "identification_percent", "lower", "upper",
              sim_label(exposure), c(0, 100), paste0("(", LETTERS[i], ")"),
              ylab = if (i %% 4L == 1L) "Identification rate (%)" else "")
      }
    } else {
      panel(d, "TPR_percent", "TPR_lower", "TPR_upper", "Active exposures identified (%)", c(0, 100), "(A)")
      max_fpr <- max(5, d$FPR_upper, d$FPR_percent, na.rm = TRUE)
      panel(d, "FPR_percent", "FPR_lower", "FPR_upper", "Null exposures selected (%)",
            c(0, min(100, 1.15 * max_fpr)), "(B)")
    }
    graphics::mtext("Number of null exposures", side = 1, outer = TRUE, line = .95, cex = 1.06)
    graphics::par(xpd = NA)
    m <- graphics::grconvertX(.5, from = "ndc", to = "user")
    y <- graphics::grconvertY(.982, from = "ndc", to = "user")
    graphics::legend(m, y, legend = sim_label(selection_methods), xjust = .5, horiz = TRUE,
      bty = "n", col = colors, pch = symbols, lty = 1, lwd = 1.8, cex = 1.0)
  })
}

selection_report <- function(context, directory) {
  data <- selection_collect(context)
  kdr_csv(data$status, file.path(context$runtime, "fit_and_test_status.csv"))
  kdr_csv(data$tests, file.path(context$runtime, "variable_tests_by_repetition.csv"))
  rates <- selection_scenarios(data$per_run)
  kdr_csv(rates, file.path(context$runtime, "selection_rates_by_repetition.csv"))
  report <- selection_summary(rates)
  tables <- file.path(directory, "tables")
  figures <- file.path(directory, "figures")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  registry <- list()
  table <- function(data, id, title, note) {
    kdr_csv(data, file.path(tables, paste0(id, ".csv")))
    sim_html_table(data, title, note, file.path(tables, paste0(id, ".html")))
    registry[[id]] <<- data.frame(file = paste0("tables/", id, ".html"), title = title, note = note)
  }
  rate_note <- paste("Holm correction at 0.05 covers all candidate exposures within each fit.",
    "Rates average complete repetition-level decisions. Intervals are 95% Monte Carlo intervals.",
    "The Holm family includes every null exposure in Module C.")
  scenario_note <- paste("The correctly adjusted CLR fit supplies G1/G2 and the misspecified CLR fit supplies G3/G4.",
    "All eight candidate exposures are active in these settings.")
  table(subset(report$overall, module == "B"), "Table_SV1_B_identification",
    "Variable identification under nuisance-model specification", paste(rate_note, scenario_note))
  table(subset(report$overall, module == "C"), "Table_SV2_C_identification",
    "Variable identification as the number of null exposures increases",
    paste(rate_note, "The settings use KDR-CC G1 and correctly adjusted CLR fits."))
  table(subset(report$overall, module == "D"), "Table_SV3_D_identification",
    "Variable identification under increased exposure correlation",
    paste(rate_note, scenario_note, "D0 reuses the original Module B results; DM and DH are the moderate and high settings."))
  table(report$specific, "Table_SV4_exposure_specific", "Identification rates for individual active exposures",
    paste("Intervals are 95% Wilson Monte Carlo intervals.",
      "Relative humidity is active through interaction with temperature in the generating response."))
  table(report$paired, "Table_SV5_paired_differences", "Paired differences in variable identification",
    paste("Differences are KDR-CC minus the comparator in percentage points, on repetitions with both tests available.",
      "Intervals are 95% Monte Carlo intervals for the mean paired differences."))

  figure <- function(d, id, title, note, by_exposure = FALSE) {
    selection_figure_c(d, file.path(figures, paste0(id, ".png")), by_exposure)
    registry[[id]] <<- data.frame(file = paste0("figures/", id, ".png"), title = title, note = note)
  }
  figure(subset(report$overall, module == "C"), "Figure_SV1_C_selection_rates",
    "Active- and null-exposure identification across dimensions",
    "Points are mean rates and vertical lines are 95% Monte Carlo intervals.")
  figure(subset(report$specific, module == "C"), "Figure_SV2_C_active_exposures",
    "Active-exposure identification across dimensions", "Points are identification rates and lines are 95% Wilson Monte Carlo intervals.", TRUE)
  captions <- kdr_bind(registry)
  kdr_csv(captions, file.path(directory, "Figure_and_table_captions.csv"))
  entries <- vapply(seq_len(nrow(captions)), function(i) {
    paste0("<li><a href=\"", captions$file[i], "\">", sim_escape(captions$title[i]), "</a><p>",
           sim_escape(captions$note[i]), "</p></li>")
  }, character(1L))
  writeLines(c("<!doctype html><html lang='en'><meta charset='utf-8'>",
    "<title>Variable identification in matched case-crossover simulations</title>",
    "<style>body{font-family:'Times New Roman',serif;max-width:920px;margin:36px auto;line-height:1.45}p{margin-top:4px}</style>",
    "<h1>Variable identification in matched case-crossover simulations</h1><ol>", entries, "</ol></html>"),
    file.path(directory, "index.html"), useBytes = TRUE)
  invisible(report)
}
