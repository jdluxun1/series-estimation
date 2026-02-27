## ============================================================
## aux_functions.R
## Helpers for empirical multivariate time-varying coefficient
## regressions using Chebyshev series.
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

build_tv_design <- function(X, tau, k, include_intercept = TRUE,
                            intercept_type = c("constant", "time_varying")) {
  # Build a series design matrix for:
  # y_t = a + b1(t) X1_t + ... + bp(t) Xp_t + u_t  (constant intercept)
  # or y_t = a(t) + ... when intercept_type = "time_varying".
  X <- as.matrix(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }

  n <- nrow(X)
  if (length(tau) != n) stop("tau and X must have the same number of rows.")

  intercept_type <- match.arg(intercept_type)
  Psi <- cheb_basis(tau, k)
  blocks <- list()
  block_names <- character(0)
  block_sizes <- integer(0)

  if (include_intercept) {
    if (intercept_type == "constant") {
      blocks[[length(blocks) + 1L]] <- matrix(1, nrow = n, ncol = 1)
      block_names <- c(block_names, "intercept")
      block_sizes <- c(block_sizes, 1L)
    } else {
      blocks[[length(blocks) + 1L]] <- Psi
      block_names <- c(block_names, "intercept")
      block_sizes <- c(block_sizes, k)
    }
  }

  for (j in seq_len(ncol(X))) {
    blocks[[length(blocks) + 1L]] <- Psi * as.numeric(X[, j])
    block_names <- c(block_names, colnames(X)[j])
    block_sizes <- c(block_sizes, k)
  }

  Xreg <- do.call(cbind, blocks)
  x_names <- unlist(Map(function(vn, bs) {
    if (bs == 1L) {
      vn
    } else {
      paste0(vn, "_psi", seq_len(bs) - 1L)
    }
  }, block_names, block_sizes))
  colnames(Xreg) <- x_names

  list(Xreg = Xreg, block_names = block_names, block_sizes = block_sizes, k = k,
       include_intercept = include_intercept, intercept_type = intercept_type)
}

fit_tv_multivar_given_k <- function(y, X, tau, k, include_intercept = TRUE,
                                    intercept_type = c("constant", "time_varying")) {
  design <- build_tv_design(X, tau, k, include_intercept = include_intercept,
                            intercept_type = intercept_type)
  fit <- lm.fit(x = design$Xreg, y = y)

  list(
    k = k,
    coef = fit$coefficients,
    resid = fit$residuals,
    fitted = fit$fitted.values,
    block_names = design$block_names,
    block_sizes = design$block_sizes,
    include_intercept = design$include_intercept,
    intercept_type = design$intercept_type
  )
}

gcv_select_k_multivar <- function(y, X, tau, k_grid, include_intercept = TRUE,
                                  intercept_type = c("constant", "time_varying")) {
  n <- length(y)
  p <- ncol(as.matrix(X))
  intercept_type <- match.arg(intercept_type)

  df_intercept <- if (!isTRUE(include_intercept)) 0L else if (intercept_type == "constant") 1L else k_grid

  out <- data.frame(k = k_grid, gcv = NA_real_)
  for (ii in seq_along(k_grid)) {
    k <- k_grid[ii]
    fit <- fit_tv_multivar_given_k(y, X, tau, k,
                                   include_intercept = include_intercept,
                                   intercept_type = intercept_type)
    rss <- mean(fit$resid^2)
    df <- p * k + if (length(df_intercept) == 1L) df_intercept else df_intercept[ii]
    out$gcv[ii] <- (1 - df / n)^(-2) * rss
  }

  out[which.min(out$gcv), "k"]
}

coef_paths_from_fit <- function(coef_vec, tau, k, block_names, block_sizes) {
  if (length(block_names) != length(block_sizes)) {
    stop("block_names and block_sizes must have the same length.")
  }

  expected <- sum(block_sizes)
  if (length(coef_vec) != expected) {
    stop("Coefficient length does not match sum(block_sizes).")
  }

  Psi <- cheb_basis(tau, k)
  out <- list()
  pos <- 1L
  for (jj in seq_along(block_names)) {
    bs <- block_sizes[jj]
    bj <- coef_vec[pos:(pos + bs - 1L)]
    if (bs == 1L) {
      out[[block_names[jj]]] <- rep(as.numeric(bj), length(tau))
    } else {
      out[[block_names[jj]]] <- as.numeric(Psi %*% bj)
    }
    pos <- pos + bs
  }

  as.data.frame(out)
}
