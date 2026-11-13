# KDR-CC simulation implementation.
kdr_objective <- function(beta, each, aux, s, lambda, N, derivatives = TRUE) {
    q <- kdr_dr_state(beta, each, aux, rep(1, 9), N, derivatives)
    W <- crossprod(s$T)
    u <- drop(s$T %*% q$U)
    P <- lambda * s$penalty
    list(value = sum(u^2)/2 + drop(crossprod(beta, P %*% beta))/2, state = q, gradient = if (derivatives) drop(crossprod(q$J, 
        W %*% q$U) + P %*% beta) else NULL, gauss_newton = if (derivatives) crossprod(s$T %*% q$J) + 
        P else NULL)
}

kdr_fit <- function(each, aux, s, cfg, N) {
    lambda <- cfg$penalty$a * cfg$penalty$reference_N/N
    b <- numeric(s$dimension)
    converged <- FALSE
    history <- list()
    for (it in seq_len(cfg$solver$maxit)) {
        q <- kdr_objective(b, each, aux, s, lambda, N)
        ok <- FALSE
        for (damp in c(0, 1e-10, 1e-08, 1e-06, 1e-04, 0.01)) {
            H <- q$gauss_newton + diag(damp * max(1, mean(diag(q$gauss_newton))), length(b))
            step <- tryCatch(drop(solve(H, -q$gradient)), error = function(e) NULL)
            if (is.null(step)) 
                next
            if (max(abs(q$gradient)) < cfg$solver$tolerance && max(abs(step))/(1 + max(abs(b))) < cfg$solver$step_tolerance) {
                converged <- TRUE
                break
            }
            slope <- sum(q$gradient * step)
            if (!is.finite(slope) || slope >= 0) 
                next
            rate <- 1
            for (bt in 0:cfg$solver$max_backtrack) {
                a <- kdr_attempt(kdr_objective(b + rate * step, each, aux, s, lambda, N, FALSE))
                if (a$ok && is.finite(a$value$value) && a$value$value <= q$value + 1e-04 * rate * slope + 
                  1e-15) {
                  ok <- TRUE
                  break
                }
                rate <- rate/2
            }
            if (ok) 
                break
        }
        history[[it]] <- data.frame(iteration=it,objective=q$value,
            gradient=max(abs(q$gradient)),residual=max(abs(s$T %*% q$state$U)),
            step_size=if(converged)0 else if(ok)rate else NA_real_,damping=damp)
        if (converged)
            break
        if (!ok) 
            stop("Regularized estimating-equation line search failed")
        b <- b + rate * step
        kdr_heartbeat(paste("Regularized optimization iteration", it, ""))
    }
    if (!converged) 
        stop("Regularized estimation did not converge within the iteration limit")
    list(history=kdr_bind(history),converged=converged,beta = b, lambda = lambda, objective = q$value, gradient = max(abs(q$gradient)), iterations = it, 
        residual = max(abs(s$T %*% q$state$U)), jacobian_condition = kappa(q$state$J, exact = TRUE), 
        gamma = Reduce(`+`, lapply(aux$g, `[[`, "gamma"))/length(aux$g))
}

kdr_curvature <- function(fit, each, aux, s, N) {
    b <- fit$beta
    p <- length(b)
    K <- length(aux$g)
    C <- ncol(aux$f$mass)
    st <- kdr_dr_state(b, each, aux, rep(1, 9), N)
    W <- crossprod(s$T)
    r <- drop(W %*% st$U)
    JW <- t(st$J) %*% W
    E <- matrix(0, p, p)
    Ag <- lapply(seq_len(K), function(k) JW %*% st$G[[k]])
    Af <- lapply(seq_len(K), function(k) lapply(seq_len(C), function(c) st$F[k, c] * JW))
    each(function(z) {
        pa <- z$pairs
        ff <- z$fold[pa$group]
        cc <- if (aux$f$wrong) 
            rep(1L, length(ff))
        else pa$cell
        for (k in seq_len(K)) {
            i <- which(ff == k)
            if (!length(i)) 
                next
            Z <- pa$Z[i, , drop = FALSE]
            R <- Z - aux$f$means[[k]][cc[i], , drop = FALSE]
            DW <- pa$DW[i, , drop = FALSE]
            eta <- drop(Z %*% b + DW %*% aux$g[[k]]$gamma)
            ee <- numeric(length(i))
            yes <- pa$D[i] == 1
            if (any(-eta[yes] > 650)) 
                stop("Exponential overflow in exact curvature")
            ee[yes] <- exp(-eta[yes])
            a <- pa$w[i] * ee/N
            rr <- drop(R %*% r)
            E <<- E + crossprod(Z, Z * (a * rr))
            Ag[[k]] <<- Ag[[k]] + crossprod(Z, DW * (a * rr))
            for (c in seq_len(C)) {
                ix <- which(cc[i] == c)
                if (length(ix)) 
                  Af[[k]][[c]] <<- Af[[k]][[c]] + outer(colSums(Z[ix, , drop = FALSE] * a[ix]), r)
            }
        }
    })
    H <- kdr_sym(crossprod(s$T %*% st$J) + E + fit$lambda * s$penalty)
    list(H = H, Ag = Ag, Af = Af, state = st, r = r, JW = JW)
}

kdr_dr_covariance <- function(fit, each, aux, s, cfg, N) {
    a <- kdr_curvature(fit, each, aux, s, N)
    H <- a$H
    if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0) 
        stop("The exact Hessian is not positive definite; covariance is unavailable")
    Hi <- kdr_inverse(H)
    p <- length(fit$beta)
    K <- length(aux$g)
    C <- ncol(aux$f$mass)
    Q <- matrix(0, p, p)
    sumS <- numeric(p)
    seen <- 0L
    each(function(z) {
        pa <- z$pairs
        ff <- z$fold[pa$group]
        cc <- if (aux$f$wrong) 
            rep(1L, length(ff))
        else pa$cell
        rows <- matrix(0, length(ff), p)
        for (k in seq_len(K)) {
            i <- which(ff == k)
            if (!length(i)) 
                next
            Z <- pa$Z[i, , drop = FALSE]
            R <- Z - aux$f$means[[k]][cc[i], , drop = FALSE]
            eta <- drop(Z %*% fit$beta + pa$DW[i, , drop = FALSE] %*% aux$g[[k]]$gamma)
            ee <- numeric(length(i))
            yes <- pa$D[i] == 1
            if (any(-eta[yes] > 650)) 
                stop("Exponential overflow in covariance scores")
            ee[yes] <- exp(-eta[yes])
            res <- ee - (1 - pa$D[i])
            w <- pa$w[i]
            rows[i, ] <- (R * (w * res)) %*% t(a$JW) - Z * (w * ee * drop(R %*% a$r))
        }
        S <- group_sum(rows, pa$group)
        Traw <- kdr_nuisance_raw(z$A, z$W, cfg)
        for (k in seq_len(K)) {
            gm <- aux$g[[k]]
            X <- Traw[, gm$keep, drop = FALSE]
            eta <- drop(X %*% gm$raw_beta)
            mx <- as.numeric(tapply(eta, z$group, max))
            ee <- exp(eta - mx[z$group])
            pr <- ee/drop(group_sum(ee, z$group))[z$group]
            IG <- group_sum((X %*% gm$influence_map) * (z$event - pr), z$group)
            IG[z$fold == k, ] <- 0
            S <- S + (N/gm$ntrain) * IG %*% t(a$Ag[[k]])
            for (c in seq_len(C)) {
                ix <- which(pa$D == 0 & ff != k & cc == c)
                if (!length(ix)) 
                  next
                R <- sweep(pa$Z[ix, , drop = FALSE], 2, aux$f$means[[k]][c, ])
                IF <- group_sum(R * (pa$w[ix]/aux$f$mass[k, c]), pa$group[ix])
                ids <- sort(unique(pa$group[ix]))
                S[ids, ] <- S[ids, , drop = FALSE] + (N/aux$f$ntrain[k]) * IF %*% t(a$Af[[k]][[c]])
            }
        }
        sumS <<- sumS + colSums(S)
        Q <<- Q + crossprod(S)
        seen <<- seen + z$n
        kdr_heartbeat("Matched-set covariance: direct and nuisance-estimation contributions")
    })
    if (seen != N) 
        stop("Matched-set count mismatch during covariance estimation")
    meat <- (Q - outer(sumS, sumS)/N)/(N * (N - 1))
    V <- kdr_sym(Hi %*% meat %*% t(Hi))
    kdr_vroot(V)
    list(V = V, hessian_condition = kappa(H, exact = TRUE), minimum_curvature = min(eigen(H, symmetric = TRUE, 
        only.values = TRUE)$values), scope = "Conditional on the independently learned kernel; matched-set linearization of the exact regularized objective, including nuisance-estimation contributions; empirical bias and coverage are evaluated against the generating response")
}

