# KDR-CC simulation implementation.
kdr_nuisance_raw <- function(A, W, cfg) {
    if (!is.null(cfg$nuisance$raw_spec))
        return(sim_E_nuisance_raw(A, W, cfg$nuisance$raw_spec))
    k <- if (identical(cfg$profile, "pollutants")) 
        6L
    else 8L
    z <- sweep(sweep(A[, seq_len(k), drop = FALSE], 2, cfg$baseline$reference$center[seq_len(k)]), 2, 
        cfg$baseline$reference$scale[seq_len(k)], "/")
    main <- sapply(1:6, function(j) kdr_anchor(A[, j], cfg$baseline$truth$anchors[[j]]))
    weather <- NULL
    if (k == 8L) {
        t <- (A[, 7] - cfg$baseline$truth$temperature$optimum)/10
        h <- (A[, 8] - 60)/20
        weather <- cbind(pmin(t, 0)^2, pmax(t, 0)^2, h^2)
    }
    pp <- combn(seq_len(k), 2)
    ix <- sapply(seq_len(ncol(pp)), function(i) z[, pp[1, i]] * z[, pp[2, i]])
    T <- cbind(W = W, W2 = W^2, main, weather, z, ix)
    if (ncol(A) > k) {
        E <- A[, -seq_len(k), drop = FALSE]
        T <- cbind(T, E)
        if (cfg$nuisance$extra_degree >= 2) 
            T <- cbind(T, E^2)
    }
    T
}

kdr_pairs <- function(d) {
    cases <- which(d$event == 1L)
    controls <- which(d$event == 0L)
    gc <- d$group[controls]
    ev <- cases[gc]
    cnt <- tabulate(d$group) - 1L
    plus <- c(ev, controls)
    minus <- c(controls, ev)
    g <- c(gc, gc)
    list(plus = plus, minus = minus, group = g, D = c(rep(1, length(ev)), rep(0, length(ev))), w = 1/(2 * 
        cnt[g]), cell = (d$W[plus] + 1L) * 3L + d$W[minus] + 2L, DW = cbind(d$W[plus] - d$W[minus], d$W[plus]^2 - 
        d$W[minus]^2))
}

