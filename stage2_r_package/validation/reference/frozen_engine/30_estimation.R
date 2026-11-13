# KDR-CC simulation implementation.
kdr_cl_block <- function(beta, V, y, g, hessian = TRUE) {
    eta <- drop(V %*% beta)
    mx <- as.numeric(tapply(eta, g, max))
    a <- exp(eta - mx[g])
    den <- drop(group_sum(a, g))
    pr <- a/den[g]
    score <- group_sum(V * (y - pr), g)
    H <- NULL
    if (hessian) {
        mu <- group_sum(V * pr, g)
        H <- crossprod(V, V * pr) - crossprod(mu)
    }
    list(loss = sum(mx + log(den)) - sum(eta * y), U = colSums(score), H = H, scores = score, probability = pr, 
        group_loss = mx + log(den) - drop(group_sum(eta * y, g)))
}

kdr_cl_state <- function(beta, each, hessian = TRUE, moments = FALSE) {
    p <- length(beta)
    U <- numeric(p)
    H <- matrix(0, p, p)
    Q <- matrix(0, p, p)
    loss <- 0
    n <- 0L
    each(function(z) {
        a <- kdr_cl_block(beta, z$V, z$event, z$group, hessian)
        U <<- U + a$U
        if (hessian) 
            H <<- H + a$H
        if (moments) 
            Q <<- Q + crossprod(a$scores)
        loss <<- loss + a$loss
        n <<- n + max(z$group)
    })
    list(loss = loss/n, U = U/n, H = H/n, Q = Q, n = n)
}

kdr_cl_fit <- function(each, p, cfg, start = NULL) {
    b <- if (is.null(start)) 
        numeric(p)
    else start
    ctrl <- cfg$solver
    converged <- FALSE
    history <- list()
    for (it in seq_len(ctrl$maxit)) {
        s <- kdr_cl_state(b, each)
        delta <- drop(kdr_inverse(s$H) %*% s$U)
        history[[it]] <- data.frame(iteration=it,objective=s$loss,
          gradient=max(abs(s$U)),residual=max(abs(s$U)),step_size=NA_real_,damping=0)
        if (max(abs(s$U)) < ctrl$tolerance && max(abs(delta))/(1 + max(abs(b))) < ctrl$step_tolerance) {
            converged <- TRUE
            break
        }
        step <- 1
        ok <- FALSE
        for (bt in 0:ctrl$max_backtrack) {
            q <- kdr_cl_state(b + step * delta, each, FALSE)
            if (is.finite(q$loss) && q$loss <= s$loss - 1e-04 * step * sum(s$U * delta) + 1e-14) {
                ok <- TRUE
                break
            }
            step <- step/2
        }
        if (!ok) 
            stop("Conditional logistic line search failed")
        history[[it]]$step_size <- step
        b <- b + step * delta
    }
    if (!converged) 
        stop("Conditional logistic estimation did not converge")
    list(history=kdr_bind(history),converged=converged,beta = b, H = s$H, n = s$n, iterations = it, loss = s$loss)
}

kdr_g_each <- function(raw, nu, cfg, k = 0L, cache = NULL) {
    function(FUN) {
        ii <- 0L
        raw(function(z) {
            ii <<- ii + 1L
            path <- if (!is.null(cache)) 
                file.path(cache, paste0("nu_", ii, "_", z$n, ".rds"))
            else NULL
            if (!is.null(path) && file.exists(path)) 
                V <- readRDS(path)
            else {
                T <- kdr_nuisance_raw(z$A, z$W, cfg)
                V <- T[, nu$keep, drop = FALSE] %*% nu$B
                if (!is.null(path)) 
                  kdr_save(V, path)
            }
            take <- if (k == 0L) 
                rep(TRUE, nrow(V))
            else z$fold[z$group] != k
            if (any(take)) {
                g <- match(z$group[take], unique(z$group[take]))
                FUN(list(V = V[take, , drop = FALSE], event = z$event[take], group = g))
            }
        })
    }
}

kdr_g_fit <- function(raw, s, cfg, mode, cache) {
    path <- file.path(cache, paste0("g_", mode, ".rds"))
    if (file.exists(path)) 
        return(readRDS(path))
    nu <- s$nuisance[[mode]]
    K <- cfg$nuisance$folds
    ans <- vector("list", K)
    for (k in seq_len(K)) {
        kdr_bar("Nuisance g folds", k - 1L, K, if (mode == "correct") 
            "Correct working family"
        else "Quadratic confounder term omitted")
        f <- kdr_cl_fit(kdr_g_each(raw, nu, cfg, k, file.path(cache, paste0("g_design_", mode))), ncol(nu$B), 
            cfg)
        ans[[k]] <- list(gamma = drop(nu$C %*% f$beta), raw_beta = drop(nu$B %*% f$beta), keep = nu$keep, 
            influence_map = nu$B %*% kdr_inverse(f$H) %*% t(nu$C), ntrain = f$n, iterations = f$iterations)
        kdr_bar("Nuisance g folds", k, K, "Complete")
    }
    kdr_save(ans, path)
    ans
}

kdr_f_fit <- function(each, cfg, N, wrong = FALSE) {
    K <- cfg$nuisance$folds
    mass <- NULL
    num <- NULL
    ng <- integer(K)
    each(function(z) {
        pa <- z$pairs
        p <- ncol(pa$Z)
        C <- if (wrong) 
            1L
        else 9L
        if (is.null(mass)) {
            mass <<- matrix(0, K, C)
            num <<- array(0, c(K, C, p))
        }
        cell <- if (wrong) 
            rep(1L, length(pa$D))
        else pa$cell
        for (k in seq_len(K)) {
            ng[k] <<- ng[k] + sum(z$fold != k)
            for (c in seq_len(C)) {
                i <- which(pa$D == 0 & z$fold[pa$group] != k & cell == c)
                mass[k, c] <<- mass[k, c] + sum(pa$w[i])
                if (length(i)) 
                  num[k, c, ] <<- num[k, c, ] + colSums(pa$Z[i, , drop = FALSE] * pa$w[i])
            }
        }
    })
    if (any(mass <= 0)) 
        stop("An auxiliary-distribution cell is empty")
    list(means = lapply(seq_len(K), function(k) matrix(num[k, , ], nrow = ncol(mass))/mass[k, ]), mass = mass/ng, 
        ntrain = ng, wrong = wrong, N = N)
}

kdr_dr_state <- function(beta, each, aux, weights, N, derivatives = TRUE) {
    p <- length(beta)
    K <- length(aux$g)
    C <- ncol(aux$f$mass)
    U <- numeric(p)
    J <- matrix(0, p, p)
    G <- lapply(seq_len(K), function(k) matrix(0, p, 2))
    F <- matrix(0, K, C)
    each(function(z) {
        pa <- z$pairs
        ff <- z$fold[pa$group]
        cell <- if (aux$f$wrong) 
            rep(1L, length(ff))
        else pa$cell
        for (k in seq_len(K)) {
            i <- which(ff == k)
            if (!length(i)) 
                next
            Z <- pa$Z[i, , drop = FALSE]
            R <- Z - aux$f$means[[k]][cell[i], , drop = FALSE]
            eta <- drop(Z %*% beta + pa$DW[i, , drop = FALSE] %*% aux$g[[k]]$gamma)
            ee <- numeric(length(i))
            yes <- pa$D[i] == 1
            if (any(-eta[yes] > 650)) 
                stop("Exponential overflow in the estimating equation")
            ee[yes] <- exp(-eta[yes])
            res <- ee - (1 - pa$D[i])
            a <- pa$w[i] * weights[pa$cell[i]]
            U <<- U + colSums(R * (a * res))
            if (derivatives) {
                J <<- J - crossprod(R, Z * (a * ee))
                G[[k]] <<- G[[k]] - crossprod(R, pa$DW[i, , drop = FALSE] * (a * ee))
                F[k, ] <<- F[k, ] - vapply(seq_len(C), function(c) sum(a[cell[i] == c] * res[cell[i] == 
                  c]), numeric(1))
            }
        }
    })
    list(U = U/N, J = J/N, G = lapply(G, function(x) x/N), F = F/N)
}

kdr_cl_model <- function(each, s, cfg, wrong = FALSE) {
    qdim <- if (wrong) 
        1L
    else 2L
    ce <- function(FUN) each(function(z) FUN(list(V = cbind(z$X, if (wrong) matrix(z$W, ncol = 1) else cbind(z$W, 
        z$W^2)), event = z$event, group = z$group)))
    f <- kdr_cl_fit(ce, s$dimension + qdim, cfg)
    st <- kdr_cl_state(f$beta, ce, TRUE, TRUE)
    Hi <- kdr_inverse(st$H)
    meat <- (st$Q - outer(st$U, st$U) * st$n)/(st$n * (st$n - 1))
    V <- Hi %*% meat %*% t(Hi)
    gam <- c(tail(f$beta, qdim), if (wrong) 0)
    list(history=f$history,converged=f$converged,objective=f$loss,gradient=max(abs(st$U)),beta = head(f$beta, s$dimension), inference_beta = head(f$beta, s$dimension), V = V[seq_len(s$dimension), 
        seq_len(s$dimension), drop = FALSE], gamma = gam, iterations = f$iterations, a = NA_real_, lambda = 0, 
        residual = max(abs(st$U)), jacobian_condition = kappa(st$H, exact = TRUE))
}

