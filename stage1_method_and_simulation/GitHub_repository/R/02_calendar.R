# KDR-CC simulation implementation.
calendar_windows <- function(year) {
    dates <- seq(as.Date(sprintf("%d-01-01", year)), as.Date(sprintf("%d-12-31", year)), by = "day")
    split(dates, paste(format(dates, "%m"), as.POSIXlt(dates)$wday, sep = "_"))
}

generate_candidates <- function(n_groups, cfg, seed) with_seed(seed, {
    windows <- calendar_windows(cfg$year)
    ans <- vector("list", n_groups)
    buffer <- if (is.null(cfg$lag_windows)) 
        max(cfg$pollutant_lag, cfg$weather_lag)
    else max(unlist(cfg$lag_windows))
    for (s in seq_len(n_groups)) {
        candidate <- windows[[sample.int(length(windows), 1L)]]
        month_start <- as.Date(format(candidate[1], "%Y-%m-01"))
        dates <- seq(month_start - buffer, max(candidate), by = "day")
        nt <- length(dates)
        W <- integer(nt)
        W[1] <- sample(c(-1L, 0L, 1L), 1L)
        for (t in 2:nt) {
            probs <- rep((1 - cfg$p_stay)/2, 3)
            probs[W[t - 1] + 2L] <- cfg$p_stay
            W[t] <- sample(c(-1L, 0L, 1L), 1L, prob = probs)
        }
        factors <- matrix(rnorm(nt * 9L), nt, 9L)
        for (t in 2:nt) factors[t, ] <- cfg$phi * factors[t - 1, ] + sqrt(1 - cfg$phi^2) * factors[t, 
            ]
        independent_multiplier <- if (is.null(cfg$independent_multiplier)) 
            rep(1, 8)
        else cfg$independent_multiplier
        latent <- sqrt(1 - cfg$common_fraction) * sweep(factors[, -1, drop = FALSE], 2, independent_multiplier, 
            "*") + sqrt(cfg$common_fraction) * outer(factors[, 1], cfg$direction) + cfg$confounding * 
            W
        if (!is.null(cfg$confounding_quadratic)) 
            latent <- latent + cfg$confounding_quadratic * (W^2 - 2/3)
        compressed <- 2.5 * tanh(latent/2.5)
        phase <- 2 * base::pi * as.integer(format(dates, "%j"))/365
        daily <- matrix(0, nt, 8L)
        daily[, 1:6] <- sweep(exp(0.35 * compressed[, 1:6] + 0.25 * outer(cos(phase), cfg$direction[1:6])), 
            2, cfg$base_concentration, "*")
        daily[, 7] <- 18 + 12 * sin(phase) + 5 * compressed[, 7]
        daily[, 8] <- 100 * plogis(qlogis(0.6) + 0.45 * compressed[, 8] + 0.25 * cos(phase))
        positions <- match(candidate, dates)
        A <- matrix(0, length(candidate), 8L)
        for (q in 1:8) {
            lags <- if (is.null(cfg$lag_windows)) 
                0:(if (q <= 6) 
                  cfg$pollutant_lag
                else cfg$weather_lag)
            else cfg$lag_windows[[q]]
            A[, q] <- vapply(positions, function(i) mean(daily[i - lags, q]), numeric(1))
        }
        colnames(A) <- exposure_names
        ans[[s]] <- data.frame(group = s, individual_id = s, date = candidate, W = W[positions], A, check.names = FALSE)
    }
    do.call(rbind, ans)
})

outcome_confounding <- function(W, cfg) {
    linear <- if (is.null(cfg$outcome_linear)) 
        0.35
    else cfg$outcome_linear
    quadratic <- if (is.null(cfg$outcome_quadratic)) 
        0.2
    else cfg$outcome_quadratic
    linear * W + quadratic * (W^2 - 2/3)
}

