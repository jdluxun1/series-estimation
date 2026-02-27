## ============================================================
## empirics_cay.R
## Time-varying coefficient regression on the `cay` sheet:
##   c_t = a + b_w(t) w_t + b_y(t) y_t + u_t
## Intercept is time-invariant; slope functions are estimated by Chebyshev series.
## ============================================================

source("aux_functions.R")

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package `readxl` is required. Install it with install.packages('readxl').")
}

raw <- readxl::read_excel("data.xlsx", sheet = "cay")
needed <- c("date", "c", "w", "y")

missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in `cay` sheet: %s",
               paste(missing_cols, collapse = ", ")))
}

dat <- raw[, needed]
dat <- dat[stats::complete.cases(dat), ]

n <- nrow(dat)
tau <- (1:n) / n

y_dep <- as.numeric(dat$c)
X <- as.matrix(dat[, c("w", "y")])
colnames(X) <- c("w", "y")

k_grid <- 1:12
k_hat <- gcv_select_k_multivar(
  y = y_dep,
  X = X,
  tau = tau,
  k_grid = k_grid,
  include_intercept = TRUE,
  intercept_type = "constant"
)

fit <- fit_tv_multivar_given_k(
  y = y_dep,
  X = X,
  tau = tau,
  k = k_hat,
  include_intercept = TRUE,
  intercept_type = "constant"
)

coef_paths <- coef_paths_from_fit(
  coef_vec = fit$coef,
  tau = tau,
  k = k_hat,
  block_names = fit$block_names,
  block_sizes = fit$block_sizes
)

coef_out <- cbind(
  data.frame(date = dat$date, tau = tau),
  coef_paths
)

fit_out <- data.frame(
  date = dat$date,
  c = y_dep,
  fitted = fit$fitted,
  residual = fit$resid
)

write.csv(coef_out, "cay_time_varying_coefficients.csv", row.names = FALSE)
write.csv(fit_out, "cay_fitted_values.csv", row.names = FALSE)

cat("Selected series order k:", k_hat, "\n")
cat("Estimated constant intercept:", fit$coef[1], "\n")
cat("Average time-varying coefficient on w:", mean(coef_paths$w), "\n")
cat("Average time-varying coefficient on y:", mean(coef_paths$y), "\n")
