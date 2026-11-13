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
