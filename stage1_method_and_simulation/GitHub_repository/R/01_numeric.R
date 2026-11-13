# KDR-CC simulation implementation.
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

canonical_bytes <- function(x) {
    if (is.environment(x) || is.function(x) || is.language(x)) 
        stop("Unsupported object in stable identifier")
    if (is.null(x)) 
        return(charToRaw("NULL;"))
    attributes_x <- attributes(x)
    if (is.list(x) && !is.null(names(x))) {
        if (anyNA(names(x)) || any(names(x) == "") || anyDuplicated(names(x))) 
            stop("Named identifier lists require unique nonempty keys")
        o <- order(enc2utf8(names(x)), method = "radix")
        x <- x[o]
        attributes_x <- attributes(x)
    }
    names_x <- attributes_x$names
    attributes_x$names <- NULL
    if (!length(attributes_x)) 
        attributes_x <- NULL
    attributes(x) <- NULL
    tag <- if (is.numeric(x)) 
        "number"
    else typeof(x)
    header <- charToRaw(paste0(tag, ":", length(x), ":"))
    payload <- switch(tag, number = {
        v <- as.double(x)
        v[!is.na(v) & v == 0] <- 0
        missing <- ifelse(is.nan(v), 2L, ifelse(is.na(v), 1L, 0L))
        v[missing > 0] <- 0
        c(as.raw(missing), writeBin(v, base::raw(), size = 8, endian = "big"))
    }, logical = as.raw(ifelse(is.na(x), 2L, as.integer(x))), character = do.call(base::c, lapply(x, 
        function(s) {
            if (is.na(s)) return(charToRaw("NA;"))
            bytes <- charToRaw(enc2utf8(s))
            c(charToRaw(paste0(length(bytes), ":")), bytes)
        })), raw = x, list = do.call(base::c, lapply(x, canonical_bytes)), stop("Unsupported stable identifier type: ", 
        tag))
    c(header, payload, charToRaw("|names|"), canonical_bytes(names_x), charToRaw("|attributes|"), canonical_bytes(attributes_x))
}

stage_seed <- function(master, ...) {
    bytes <- as.integer(canonical_bytes(list(policy = "kdrcc-rng-v2", keys = list(...))))
    h <- as.double(master)
    for (b in bytes) h <- (h * 131 + b)%%2147483646
    as.integer(h + 1)
}

object_id <- function(x) {
    f <- tempfile()
    on.exit(unlink(f))
    writeBin(canonical_bytes(x), f)
    unname(tools::md5sum(f))
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

trap_weights <- function(x) {
    if (length(x) < 2L || any(diff(x) <= 0)) 
        stop("Integration grid must increase")
    dx <- diff(x)
    w <- c(dx[1], head(dx, -1) + tail(dx, -1), tail(dx, 1))/2
    w/sum(w)
}

kdr_save <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile("pending_", dirname(path))
    saveRDS(x, tmp, version = 3)
    if (file.exists(path)) {
        old <- paste0(path, ".previous")
        if (file.exists(old)) 
            unlink(old)
        if (!file.rename(path, old)) 
            stop("Cannot back up checkpoint: ", path)
    }
    else old <- NULL
    if (!file.rename(tmp, path)) {
        if (!is.null(old)) 
            file.rename(old, path)
        stop("Cannot commit checkpoint")
    }
    if (!is.null(old) && file.exists(old)) 
        unlink(old)
    invisible(path)
}

kdr_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    con <- file(path, "wb")
    writeBin(as.raw(c(239, 187, 191)), con)
    text <- capture.output(write.table(x, row.names = FALSE, sep = ",", na = "", qmethod = "double"))
    writeBin(charToRaw(enc2utf8(paste0(paste(text, collapse = "\r\n"), "\r\n"))), con)
    close(con)
}

kdr_bind <- function(x) {
    x <- Filter(function(z) !is.null(z) && NROW(z) > 0, x)
    if (!length(x)) 
        return(data.frame())
    do.call(rbind, x)
}

kdr_attempt <- function(expr) {
    t <- proc.time()[3]
    tryCatch(list(ok = TRUE, value = force(expr), error = "", seconds = unname(proc.time()[3] - t)), 
        error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e), seconds = unname(proc.time()[3] - 
            t)))
}

kdr_center <- function(X, g) {
    mu <- group_sum(X, g)/tabulate(g)
    X - mu[g, , drop = FALSE]
}

kdr_whiten <- function(X, g, tol = 1e-12) {
    Xc <- kdr_center(X, g)
    sc <- sqrt(colMeans(Xc^2))
    if (any(sc < 1e-14)) 
        stop("A design column has no within-set variation")
    q <- qr(sweep(Xc, 2, sc, "/")/sqrt(max(g)), tol = tol)
    if (q$rank < ncol(X)) 
        stop("The within-set design is rank deficient")
    P <- diag(ncol(X))[, q$pivot, drop = FALSE]
    P %*% diag(1/sc[q$pivot], ncol(X)) %*% backsolve(qr.R(q), diag(ncol(X)))
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

