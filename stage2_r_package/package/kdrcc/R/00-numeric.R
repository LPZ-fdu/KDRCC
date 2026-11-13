with_seed <- function(seed, expr) {
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had) 
        old <- get(".Random.seed", envir = .GlobalEnv)
    kind <- RNGkind()
    on.exit({
        do.call(RNGkind, as.list(kind))
        if (had) assign(".Random.seed", old, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, 
            inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    })
    RNGversion("3.6.0")
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(seed)
    force(expr)
}

wquantile <- function(x, w, p) {
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    vapply(p, function(t) x[which(cumsum(w) >= t)[1L]], numeric(1))
}

group_sum <- function(x, group) {
    if (is.null(dim(x))) 
        x <- matrix(x, ncol = 1L)
    unname(rowsum(x, group, reorder = TRUE))
}

kdr_bind <- function(x) {
    x <- Filter(function(z) !is.null(z) && NROW(z) > 0, x)
    if (!length(x)) 
        return(data.frame())
    do.call(rbind, x)
}

kdr_center <- function(X, g) {
    mu <- group_sum(X, g)/tabulate(g)
    X - mu[g, , drop = FALSE]
}

kdr_whiten <- function(X, g, tol = 1e-12, reduce = FALSE) {
    Xc <- kdr_center(X, g)
    sc <- sqrt(colMeans(Xc^2))
    if (any(!is.finite(sc)))
        stop("Nonfinite within-set variation in the design")
    active <- which(sc >= 1e-14)
    if (!reduce && length(active) < ncol(X))
        stop("A design column has no within-set variation")
    if (!length(active))
        stop("No design column has usable within-set variation")
    if (length(active) < ncol(X))
        Xc <- Xc[, active, drop = FALSE]
    q <- qr(sweep(Xc, 2, sc[active], "/")/sqrt(max(g)), tol = tol)
    if (!reduce && q$rank < ncol(X))
        stop("The within-set design is rank deficient")
    if (!q$rank) stop("The within-set design has rank zero")
    r <- seq_len(q$rank)
    keep <- active[q$pivot[r]]
    P <- diag(ncol(X))[, keep, drop = FALSE]
    B <- P %*% diag(1/sc[keep], q$rank) %*% backsolve(qr.R(q)[r, r, drop = FALSE], diag(q$rank))
    if (reduce) {
        status <- rep("Linearly dependent", ncol(X))
        status[sc < 1e-14] <- "Below variation threshold"
        status[keep] <- "Retained"
        attr(B, "kernel_basis") <- data.frame(center = seq_len(ncol(X)),
            within_set_rms = sc, status = status, stringsAsFactors = FALSE)
    }
    B
}

kdr_inverse <- function(A, tol = 1e-13) {
    if (any(!is.finite(A)) || rcond(A) < tol) 
        stop("The information matrix is nearly singular")
    solve(A)
}

kdr_psd <- function(P) {
    P <- (P + t(P))/2
    e <- eigen(P, symmetric = TRUE)
    if (min(e$values) < -1e-07 * max(1, max(e$values))) 
        stop("The penalty matrix is not positive semidefinite")
    P <- tcrossprod(sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*"))
    P/max(mean(diag(P)), 1e-30)
}

kdr_distance <- function(A, B) pmax(outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B), 0)

kdr_sym <- function(A) (A + t(A))/2

kdr_vroot <- function(V) {
    e <- eigen(kdr_sym(V), symmetric = TRUE)
    if (any(!is.finite(e$values)) || min(e$values) < -1e-08 * max(1, max(e$values))) 
        stop("The covariance matrix is not positive semidefinite")
    sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*")
}

kdr_log <- function(...) {
  if (isTRUE(getOption("kdrcc.trace", FALSE))) {
    cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), paste(..., collapse=" "), "\n", sep="")
    flush.console()
  }
}
kdr_bar <- function(stage, i, n, detail="") kdr_log(stage, paste0(i,"/",n), detail)
kdr_heartbeat <- function(...) invisible(NULL)
.kdr_scalar <- function(x, name, lower=0, integer=FALSE, strict=FALSE) {
  if (length(x)!=1L || !is.numeric(x) || !is.finite(x) ||
      (if(strict)x<=lower else x<lower) || (integer && x!=floor(x)))
    stop(name, " has an invalid value.", call.=FALSE)
  x
}
.kdr_columns <- function(data, cols) {
  if (!is.data.frame(data)) stop("data must be a data.frame.",call.=FALSE)
  if (!all(cols %in% names(data))) stop("Missing columns: ",paste(setdiff(cols,names(data)),collapse=", "),call.=FALSE)
}
