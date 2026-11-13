library(kdrcc)

set.seed(30)
group <- rep(seq_len(5L), each = 3L)
X <- matrix(rnorm(45), 15L)
X_bad <- cbind(X, group, X[, 1L], 1e-17 * X[, 2L])
B <- kdrcc:::kdr_kernel_whiten(X_bad, group)
status <- attr(B, "kernel_basis")$status
stopifnot(ncol(B) == 3L, sum(status == "Below variation threshold") == 2L,
          sum(status == "Linearly dependent") == 1L,
          max(abs(crossprod(kdrcc:::kdr_center(X_bad, group) %*% B) / 5 - diag(3))) < 1e-12)
stopifnot(max(abs(kdrcc:::kdr_kernel_whiten(X, group) - kdrcc:::kdr_whiten(X, group))) < 1e-12)

data <- kdrcc_example(800L, discrete = TRUE)
parts <- kdrcc_split(data, "set_id", learning = .25, tuning = .20)
control <- kdrcc_control(centers = 12L, penalty_rows = 24L, evaluation_points = 50L)
grid <- data.frame(M = 12L, kappa = 1e-4, lambda = c(.3, 1))
ns <- asNamespace("kdrcc")
learner <- get(".kdr_learn", ns)
geometry <- get("kdr_geometry", ns)
with_degenerate_centers <- function(...) {
  s <- geometry(...)
  s$centers[11L, ] <- 1e8
  s$centers[12L, ] <- s$centers[1L, ]
  s
}
environment(learner) <- list2env(list(kdr_geometry = with_degenerate_centers), parent = ns)
task <- get(".kdr_cv_task", ns)
covariance <- get(".kdr_covariance", ns)
check_covariance <- function(beta, z, aux, s, lambda) {
  if (lambda == .3) stop("Test unavailable covariance")
  covariance(beta, z, aux, s, lambda)
}
environment(task) <- list2env(list(.kdr_learn = learner, .kdr_covariance = check_covariance), parent = ns)
tune <- kdrcc_tune
environment(tune) <- list2env(list(.kdr_cv_task = task), parent = ns)

selection <- tune(parts$tuning, parts$learning, exposures = paste0("X", 1:4),
                  adjust = ~ W, grid = grid, control = control, fallback = FALSE)
stopifnot(all(selection$grid$score_complete),
          identical(selection$grid$eligible, c(FALSE, TRUE)),
          identical(selection$grid$covariance_folds, c(0L, 3L)),
          all(selection$grid$kernel_dimension == 10L),
          all(selection$scores$covariance_error[selection$scores$candidate == 1L] ==
              "Test unavailable covariance"))
fit <- kdrcc(parts$estimation, exposures = paste0("X", 1:4), adjust = ~ W,
             kernel = selection$kernel, control = selection$control)
diagnostic <- kdrcc_diagnostics(fit)
stopifnot(nrow(vcov(fit)) == 10L, diagnostic$data$candidate_centers == 12L,
          diagnostic$data$kernel_dimension == 10L,
          diagnostic$covariance_status == "available",
          all(is.finite(predict(fit, head(data))$or_lower)))
