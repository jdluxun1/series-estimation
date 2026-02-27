## ============================================================
## aux_functions.R
## Helpers for Section 4 Monte Carlo:
##  - DGP: y_t = A_t X_t + u0_t,  X_t = X_{t-1} + ux_t
##  - (u0_t, ux_t)' is bivariate MA(1): u_t = eps_t + Theta eps_{t-1}
##  - Series estimator with Chebyshev time polynomials
##  - Optional augmentation with leads/lags of ΔX_t
##  - GCV selection of k
## ============================================================

cheb_basis <- function(tau, k) {
  # tau in [0,1], length n
  # psi0 = 1, psij = sqrt(2) cos(j*pi*tau), j>=1
  n <- length(tau)
  Psi <- matrix(1, n, k)
  if (k >= 2) {
    j <- 1:(k - 1)
    Psi[, 2:k] <- sqrt(2) * cos(outer(tau, j, function(a, b) b * pi * a))
  }
  Psi
}

A_true <- function(tau, model = c("M1", "M2", "M3")) {
  model <- match.arg(model)
  if (model == "M1") return(1 + tau)
  if (model == "M2") return(cos(2 * pi * tau))
  if (model == "M3") return(exp(-tau))
  stop("Unknown model.")
}

make_ma1_bivar <- function(n, Theta, Sigma_e, burn = 200) {
  # Simulate u_t = eps_t + Theta eps_{t-1}
  # eps_t ~ iid N(0, Sigma_e)
  # Return u0, ux of length n
  stopifnot(is.matrix(Theta), all(dim(Theta) == c(2, 2)))
  stopifnot(is.matrix(Sigma_e), all(dim(Sigma_e) == c(2, 2)))
  
  n2 <- n + burn
  L <- chol(Sigma_e)
  
  eps <- matrix(rnorm(2 * n2), n2, 2) %*% L
  u <- matrix(0, n2, 2)
  eps_lag <- rbind(rep(0, 2), eps[-n2, , drop = FALSE])
  u <- eps + eps_lag %*% t(Theta)
  
  u <- u[(burn + 1):n2, , drop = FALSE]
  list(u0 = u[, 1], ux = u[, 2])
}

build_series_regressors <- function(X, tau, k) {
  # For scalar X: X^(k)_t = (X_t, psi1(tau)*X_t, ..., psi_{k-1}(tau)*X_t)
  # Design matrix: n x k
  Psi <- cheb_basis(tau, k)
  Xk <- Psi * as.numeric(X)
  Xk
}

build_aug_matrix <- function(dx, Ka = 0L, Ke = 0L) {
  # Build Z_t = (ΔX_{t-Ka},...,ΔX_{t-1}, ΔX_t, ΔX_{t+1},...,ΔX_{t+Ke})
  # Returned as a matrix aligned on "central" t, with rows t = Ka+1,...,n-Ke
  dx <- as.numeric(dx)
  n <- length(dx)
  if (Ka < 0 || Ke < 0) stop("Ka, Ke must be nonnegative integers.")
  if (Ka + Ke >= n - 1) stop("Ka+Ke too large for sample size.")
  
  idx_center <- (Ka + 1):(n - Ke)
  cols <- list()
  
  # lags: ΔX_{t-1},...,ΔX_{t-Ka}
  if (Ka > 0) {
    for (j in 1:Ka) cols[[length(cols) + 1]] <- dx[idx_center - j]
  }
  # contemporaneous ΔX_t
  cols[[length(cols) + 1]] <- dx[idx_center]
  # leads: ΔX_{t+1},...,ΔX_{t+Ke}
  if (Ke > 0) {
    for (j in 1:Ke) cols[[length(cols) + 1]] <- dx[idx_center + j]
  }
  
  Z <- do.call(cbind, cols)
  colnames(Z) <- c(
    if (Ka > 0) paste0("dX_lag", 1:Ka) else NULL,
    "dX_0",
    if (Ke > 0) paste0("dX_lead", 1:Ke) else NULL
  )
  list(Z = Z, idx = idx_center)
}

ols_fit_fast <- function(y, X) {
  # fast OLS with intercept excluded (we do NOT include intercept here
  # since psi0=1 already handles level in A(t), and y = A(t)X + u0)
  fit <- lm.fit(x = X, y = y)
  list(coef = fit$coefficients, resid = fit$residuals, fitted = fit$fitted.values)
}

fit_series_given_k <- function(y, X, tau, k) {
  Xk <- build_series_regressors(X, tau, k)
  fit <- ols_fit_fast(y, Xk)
  list(k = k, coef = fit$coef, resid = fit$resid)
}

fit_series_aug_given_k <- function(y, X, dx, tau, k, Ka = 0L, Ke = 0L) {
  # y_t ~ X^(k)_t + Z_t  where Z_t are leads/lags of ΔX_t
  aug <- build_aug_matrix(dx, Ka, Ke)
  idx <- aug$idx
  
  y2   <- y[idx]
  X2   <- X[idx]
  tau2 <- tau[idx]
  
  Xk2 <- build_series_regressors(X2, tau2, k)
  Xreg <- cbind(Xk2, aug$Z)
  
  fit <- ols_fit_fast(y2, Xreg)
  list(k = k, coef = fit$coef, resid = fit$resid, idx = idx,
       p_series = ncol(Xk2), p_aug = ncol(aug$Z))
}

gcv_select_k <- function(y, X, tau, k_grid) {
  # GCV(k) = (1 - df(k)/n)^(-2) * (1/n) sum uhat^2
  n <- length(y)
  out <- data.frame(k = k_grid, gcv = NA_real_)
  
  for (ii in seq_along(k_grid)) {
    k <- k_grid[ii]
    fit <- fit_series_given_k(y, X, tau, k)
    rss <- mean(fit$resid^2)
    df  <- k  # number of parameters (series only)
    out$gcv[ii] <- (1 - df / n)^(-2) * rss
  }
  out[which.min(out$gcv), "k"]
}

gcv_select_k_aug <- function(y, X, dx, tau, k_grid, Ka = 0L, Ke = 0L) {
  aug <- build_aug_matrix(dx, Ka, Ke)
  idx <- aug$idx
  y2 <- y[idx]
  X2 <- X[idx]
  tau2 <- tau[idx]
  n2 <- length(y2)
  
  out <- data.frame(k = k_grid, gcv = NA_real_)
  for (ii in seq_along(k_grid)) {
    k <- k_grid[ii]
    fit <- fit_series_aug_given_k(y, X, dx, tau, k, Ka, Ke)
    rss <- mean(fit$resid^2)
    df  <- k + ncol(aug$Z)  # series + augmented ΔX terms
    out$gcv[ii] <- (1 - df / n2)^(-2) * rss
  }
  out[which.min(out$gcv), "k"]
}

Ahat_from_coef <- function(beta_series, tau, k) {
  # beta_series: length k, corresponding to (B0,...,B_{k-1}) in scalar case
  Psi <- cheb_basis(tau, k)  # n x k
  as.numeric(Psi %*% beta_series)
}

one_rep <- function(n, model, theta21, sigma21,
                    Theta_fixed = matrix(c(0.3, -0.4, NA, 0.6), 2, 2, byrow = TRUE),
                    k_grid = NULL,
                    Ka = 0L, Ke = 0L,
                    use_aug = FALSE) {
  
  if (is.null(k_grid)) {
    # modest grid; you can widen if you want
    k_grid <- unique(pmax(3, round(seq(5, min(80, floor(n/3)), length.out = 20))))
  }
  
  Theta <- Theta_fixed
  Theta[2,1] <- theta21
  
  Sigma_e <- matrix(c(1, sigma21, sigma21, 1), 2, 2)
  
  # DGP
  tau <- (1:n) / n
  A   <- A_true(tau, model)
  
  u <- make_ma1_bivar(n, Theta, Sigma_e)
  u0 <- u$u0
  ux <- u$ux
  
  X <- cumsum(ux)  # X0=0; if you want X0=Op(1), add rnorm(1) constant shift
  dx <- ux         # ΔX_t = ux_t
  
  y <- A * X + u0
  
  # Estimation
  if (!use_aug) {
    k_hat <- gcv_select_k(y, X, tau, k_grid)
    fit <- fit_series_given_k(y, X, tau, k_hat)
    Ahat <- Ahat_from_coef(fit$coef, tau, k_hat)
    mse <- mean((Ahat - A)^2)
    return(list(mse = mse, k = k_hat))
  } else {
    k_hat <- gcv_select_k_aug(y, X, dx, tau, k_grid, Ka, Ke)
    fit <- fit_series_aug_given_k(y, X, dx, tau, k_hat, Ka, Ke)
    
    # series coefficients are the first k_hat entries
    beta_series <- fit$coef[1:k_hat]
    # Ahat only defined on idx (trimmed due to leads/lags)
    tau2 <- tau[fit$idx]
    A2   <- A[fit$idx]
    Ahat2 <- Ahat_from_coef(beta_series, tau2, k_hat)
    mse <- mean((Ahat2 - A2)^2)
    
    return(list(mse = mse, k = k_hat))
  }
}
