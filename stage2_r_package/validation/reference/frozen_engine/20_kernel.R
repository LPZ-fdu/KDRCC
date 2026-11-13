# KDR-CC simulation implementation.
kdr_kernel_raw <- function(A, s) {
    Z <- sweep(sweep(as.matrix(A), 2, s$center), 2, s$scale, "/")
    K <- exp(-kdr_distance(sweep(Z, 2, sqrt(s$rho), "*"), sweep(s$centers, 2, sqrt(s$rho), "*")))
    sweep(K, 2, exp(-drop(s$centers^2 %*% s$rho)))
}

kdr_basis <- function(A, s) {
    if (s$kind == "kernel") 
        return(kdr_kernel_raw(A, s) %*% s$B)
    if (s$kind == "linear") 
        raw <- sweep(sweep(A, 2, s$center), 2, s$scale, "/")
    else raw <- do.call(cbind, lapply(seq_len(ncol(A)), function(j) splines::ns(A[, j], knots = s$spec[[j]]$knots, 
        Boundary.knots = s$spec[[j]]$boundary, intercept = FALSE)))
    raw %*% s$B
}

kdr_derivative <- function(A, s, j) {
    Z <- sweep(sweep(A, 2, s$center), 2, s$scale, "/")
    K <- exp(-kdr_distance(sweep(Z, 2, sqrt(s$rho), "*"), sweep(s$centers, 2, sqrt(s$rho), "*")))
    (-2 * s$rho[j] * outer(Z[, j], s$centers[, j], "-") * K) %*% s$B
}

kdr_weighted_centers <- function(Z, w, M, cfg, seed) with_seed(seed, {
    nc <- max(2L, M - round(M * cfg$kernel$coverage_fraction))
    best <- NULL
    for (st in seq_len(cfg$kernel$kmeans_starts)) {
        C <- matrix(0, nc, ncol(Z))
        C[1, ] <- Z[sample.int(nrow(Z), 1, prob = w), ]
        dist <- drop(kdr_distance(Z, C[1, , drop = FALSE]))
        for (j in 2:nc) {
            C[j, ] <- Z[sample.int(nrow(Z), 1, prob = w * pmax(dist, 1e-30)), ]
            dist <- pmin(dist, drop(kdr_distance(Z, C[j, , drop = FALSE])))
        }
        for (it in seq_len(cfg$kernel$kmeans_maxit)) {
            D <- kdr_distance(Z, C)
            cl <- max.col(-D, ties.method = "first")
            old <- C
            for (j in seq_len(nc)) {
                ix <- which(cl == j)
                C[j, ] <- if (length(ix)) 
                  colSums(Z[ix, , drop = FALSE] * w[ix])/sum(w[ix])
                else Z[which.max(w * apply(D, 1, min)), ]
            }
            if (max(abs(C - old)) < 1e-05) 
                break
        }
        loss <- sum(w * apply(kdr_distance(Z, C), 1, min))
        if (is.null(best) || loss < best$loss) 
            best <- list(C = C, loss = loss)
    }
    C <- best$C
    lo <- apply(Z, 2, wquantile, w = w, p = 0.005)
    hi <- apply(Z, 2, wquantile, w = w, p = 0.995)
    eligible <- which(rowSums(sweep(Z, 2, lo, ">=") & sweep(Z, 2, hi, "<=")) == ncol(Z))
    if (length(eligible) < M - nc) 
        eligible <- seq_len(nrow(Z))
    for (j in seq_len(M - nc)) {
        dist <- apply(kdr_distance(Z[eligible, , drop = FALSE], C), 1, min)
        ii <- which.max(dist)
        C <- rbind(C, Z[eligible[ii], ])
        eligible <- eligible[-ii]
    }
    C
})

kdr_geometry <- function(d, cfg, seed) {
    A <- d$A
    w <- 1/(d$n * tabulate(d$group)[d$group])
    center <- apply(A, 2, wquantile, w = w, p = 0.5)
    mu <- colSums(A * w)
    scale <- sqrt(colSums(sweep(A, 2, mu)^2 * w))
    if (any(scale <= 0)) 
        stop("An exposure has zero scale")
    Z <- sweep(sweep(A, 2, center), 2, scale, "/")
    C <- kdr_weighted_centers(Z, w, cfg$M, cfg, seed)
    D <- as.vector(kdr_distance(Z, C))
    base <- 1/median(D[D > 0])
    list(kind = "kernel", cols = seq_len(ncol(A)), center = center, scale = scale, centers = C, base_rho = rep(base, 
        ncol(A)), rho = rep(base * cfg$kernel$initial_scale, ncol(A)), M = cfg$M)
}

kdr_rho_gradient <- function(Z, centers, K, rho, kr, residual, coefficients, N) {
    WK <- sweep(K * residual, 2, coefficients, "*")
    first <- colSums(Z^2 * rowSums(WK))
    second <- colSums(centers^2 * colSums(WK))
    cross <- colSums(centers * crossprod(WK, Z))
    reference <- sum(residual) * colSums(centers^2 * (coefficients * kr))
    rho * (-first - second + 2 * cross + reference)/N
}

kdr_learning_fit <- function(d, s, cfg) {
    ctrl <- cfg$kernel
    p <- ncol(d$A)
    M <- s$M
    kap <- ctrl$kappa
    Z <- sweep(sweep(d$A, 2, s$center), 2, s$scale, "/")
    starts <- list()
    best <- NULL
    lower <- log(s$base_rho * ctrl$rho_bounds[1])
    upper <- log(s$base_rho * ctrl$rho_bounds[2])
    for (st in seq_along(ctrl$start_factors)) {
        kdr_bar("Kernel initializations", st - 1, length(ctrl$start_factors), "Profile conditional likelihood with learned distance parameters")
        env <- new.env(parent = emptyenv())
        env$last <- NULL
        evals <- 0L
        evaluate <- function(logrho) {
            if (!is.null(env$last) && identical(logrho, env$last)) 
                return(env$value)
            ss <- s
            ss$rho <- exp(logrho)
            Phi <- kdr_kernel_raw(d$A, ss)
            X <- cbind(Phi, d$W, d$W^2)
            b <- if (is.null(env$beta)) 
                numeric(ncol(X))
            else env$beta
            converged <- FALSE
            for (it in seq_len(cfg$solver$maxit)) {
                a <- kdr_cl_block(b, X, d$event, d$group)
                grad <- -a$U/d$n + kap * b
                H <- a$H/d$n + diag(kap, length(b))
                step <- drop(solve(H, -grad))
                obj <- a$loss/d$n + kap * sum(b^2)/2
                if (max(abs(grad)) < min(cfg$solver$tolerance, 1e-09) && max(abs(step))/(1 + max(abs(b))) < 
                  cfg$solver$step_tolerance) {
                  converged <- TRUE
                  break
                }
                ok <- FALSE
                rate <- 1
                for (bt in 0:cfg$solver$max_backtrack) {
                  bb <- b + rate * step
                  val <- kdr_cl_block(bb, X, d$event, d$group, FALSE)$loss/d$n + kap * sum(bb^2)/2
                  if (is.finite(val) && val <= obj + 1e-04 * rate * sum(grad * step) + 1e-15) {
                    ok <- TRUE
                    break
                  }
                  rate <- rate/2
                }
                if (!ok) 
                  stop("Kernel inner Newton line search failed")
                b <- bb
            }
            if (!converged) 
                stop("Kernel inner coefficients did not converge")
            rho <- ss$rho
            K <- exp(-kdr_distance(sweep(Z, 2, sqrt(rho), "*"), sweep(s$centers, 2, sqrt(rho), "*")))
            kr <- exp(-drop(s$centers^2 %*% rho))
            res <- a$probability - d$event
            gr <- kdr_rho_gradient(Z, s$centers, K, rho, kr, res, b[seq_len(M)], d$n)
            evals <<- evals + 1L
            kdr_heartbeat(paste("Kernel objective evaluations", evals, ""))
            env$last <- logrho
            env$beta <- b
            env$value <- list(value = obj, gradient = gr, beta = b)
            env$value
        }
        initial <- pmax(lower, pmin(upper, log(s$rho * ctrl$start_factors[st])))
        for (attempt in seq_len(1L + ctrl$continuations)) {
            opt <- optim(initial, function(x) evaluate(x)$value, function(x) evaluate(x)$gradient, method = "L-BFGS-B", 
                lower = lower, upper = upper, control = list(maxit = ctrl$maxit, pgtol = ctrl$projected_gradient_tol, 
                  factr = 1000))
            fin <- evaluate(opt$par)
            pg <- max(abs(opt$par - pmax(lower, pmin(upper, opt$par - fin$gradient))))
            if (pg <= ctrl$projected_gradient_tol || !opt$convergence %in% c(0L, 1L)) 
                break
            initial <- opt$par
            if (attempt <= ctrl$continuations) 
                kdr_log("Continuing kernel optimization from the current solution")
        }
        good <- is.finite(pg) && pg <= ctrl$projected_gradient_tol
        starts[[st]] <- data.frame(start = st, objective = fin$value, projected_gradient = pg, converged = good, 
            optimizer_code = opt$convergence, optimization_blocks = attempt, evaluations = evals, boundaries = sum(abs(opt$par - 
                lower) < 1e-06 | abs(opt$par - upper) < 1e-06))
        if (good && (is.null(best) || fin$value < best$objective)) 
            best <- list(rho = exp(opt$par), beta = fin$beta, objective = fin$value)
        kdr_bar("Kernel initializations", st, length(ctrl$start_factors), paste(if (good) 
            "Converged"
        else "Projected-gradient tolerance not met", " | projected gradient", format(pg, digits = 3), " | optimizer status", 
            opt$convergence))
    }
    if (is.null(best)) 
        stop("No kernel initialization met the convergence criterion; ", paste(vapply(starts, function(z) paste0("initialization", 
            z$start, ": gradient=", signif(z$projected_gradient, 4), ", status=", z$optimizer_code), character(1)), 
            collapse = "; "))
    best$diagnostics <- kdr_bind(starts)
    best
}

kdr_nuisance_spec <- function(des, cfg) {
    T <- kdr_nuisance_raw(des$A, des$W, cfg)
    out <- lapply(c("correct", "wrong"), function(mode) {
        keep <- if (mode == "correct") 
            seq_len(ncol(T))
        else setdiff(seq_len(ncol(T)), 2L)
        B <- kdr_whiten(T[, keep, drop = FALSE], des$group)
        C <- matrix(0, 2, ncol(B))
        C[1, ] <- B[match(1L, keep), ]
        if (2L %in% keep) 
            C[2, ] <- B[match(2L, keep), ]
        list(keep = keep, B = B, C = C)
    })
    names(out) <- c("correct", "wrong")
    out
}

kdr_kernel_penalty <- function(des, s, cfg, seed) {
    ids <- with_seed(seed, sample.int(nrow(des$A), min(cfg$kernel$penalty_rows, nrow(des$A))))
    Z <- sweep(sweep(des$A[ids, , drop = FALSE], 2, s$center), 2, s$scale, "/")
    P <- matrix(0, s$M, s$M)
    for (i in seq_len(nrow(Z))) {
        D <- sweep(s$centers, 2, Z[i, ], "-")
        K <- exp(-drop(D^2 %*% s$rho))
        E <- sweep(D, 2, s$rho, "*")
        Q <- (tcrossprod(E)^2 - tcrossprod(E^2))/2
        P <- P + 16 * outer(K, K) * Q/nrow(Z)
    }
    kdr_psd(crossprod(s$B, P %*% s$B))
}

kdr_comparator_setup <- function(des, cfg, kind) {
    A <- des$A
    s <- list(kind = kind, center = apply(A, 2, median), scale = apply(A, 2, sd), cols = seq_len(ncol(A)))
    if (kind == "spline") 
        s$spec <- lapply(seq_len(ncol(A)), function(j) list(knots = quantile(A[, j], cfg$comparator$knots), 
            boundary = quantile(A[, j], cfg$comparator$boundary)))
    s$B <- diag(ncol(A) * if (kind == "spline") 
        length(cfg$comparator$knots) + 1L
    else 1L)
    X <- kdr_basis(A, s)
    s$B <- kdr_whiten(X, des$group)
    s$dimension <- ncol(s$B)
    s
}

kdr_feature_iterator <- function(raw, s, folder, cfg) {
    memo <- new.env(parent=emptyenv()); used <- 0
    limit <- cfg$execution$feature_memory_MB * 1024^2
    if (!length(limit)) limit <- 0
    function(FUN, progress = FALSE) {
        b <- 0L
        raw(function(z) {
            b <<- b + 1L
            path <- file.path(folder, paste0("X_", b, "_", z$n, ".rds"))
            key <- paste0(b,"_",z$n)
            if (exists(key,memo,inherits=FALSE)) X <- memo[[key]] else if (file.exists(path))
                X <- readRDS(path)
            else {
                X <- kdr_basis(z$A, s)
                kdr_save(X, path)
            }
            if (!exists(key,memo,inherits=FALSE)) {
                bytes <- as.numeric(object.size(X))
                if (used+bytes<=limit) {memo[[key]]<-X;used<<-used+bytes}
            }
            z$X <- X
            z$pairs$Z <- X[z$pairs$plus, , drop = FALSE] - X[z$pairs$minus, , drop = FALSE]
            FUN(z)
        }, progress)
    }
}

