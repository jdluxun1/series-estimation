## ============================================================
## main_simulation.R
## Runs Section 4 Monte Carlo:
##  - Models: M1, M2, M3 for A_t
##  - theta21 in {0.8,0.4,0,-0.8}, sigma21 in {-0.85,-0.5,0.5}
##  - Compare baseline (no augmentation) vs augmentation with (Ka,Ke)
##    following Section 4 combinations:
##      (0,0), (2,0), (4,0), (2,1), (4,2)
## ============================================================



set.seed(123)

# Settings
n_vec   <- c(200, 500,1000)        # adjust as you like
R       <- 2000               # reps per design point (increase for final)
models  <- c("M1", "M2", "M3")

theta21_grid <- c(0.8, 0.4, 0.0, -0.8)
sigma21_grid <- c(-0.85, -0.5, 0.5)

# (Ka, Ke) combos from Section 4
aug_cases <- data.frame(
  Ka = c(0, 2, 4, 2, 4,1,4,8),
  Ke = c(0, 0, 0, 1, 2,1,4,8)
)

# k grid (you can make n-dependent inside loop if preferred)
make_k_grid <- function(n) {
  unique(pmax(3, round(seq(5, min(100, floor(n/2)), length.out = 25))))
}

# Storage
results <- list()

counter <- 0L
for (n in n_vec) {
  k_grid <- 1:12
  
  for (model in models) {
    for (theta21 in theta21_grid) {
      for (sigma21 in sigma21_grid) {
        
        # Baseline (no augmentation)
        mse0 <- numeric(R)
        k0   <- numeric(R)
        
        for (r in 1:R) {
          out <- one_rep(
            n = n, model = model,
            theta21 = theta21, sigma21 = sigma21,
            k_grid = k_grid,
            use_aug = FALSE
          )
          mse0[r] <- out$mse
          k0[r]   <- out$k
        }
        
        # Augmented cases
        mse_aug <- matrix(NA_real_, nrow = R, ncol = nrow(aug_cases))
        k_aug   <- matrix(NA_real_, nrow = R, ncol = nrow(aug_cases))
        
        for (cc in 1:nrow(aug_cases)) {
          Ka <- aug_cases$Ka[cc]
          Ke <- aug_cases$Ke[cc]
          
          for (r in 1:R) {
            out <- one_rep(
              n = n, model = model,
              theta21 = theta21, sigma21 = sigma21,
              k_grid = k_grid,
              Ka = Ka, Ke = Ke,
              use_aug = TRUE
            )
            mse_aug[r, cc] <- out$mse
            k_aug[r, cc]   <- out$k
          }
        }
        
        counter <- counter + 1L
        results[[counter]] <- list(
          n = n, model = model, theta21 = theta21, sigma21 = sigma21,
          baseline = list(mse_mean = mean(mse0), mse_sd = sd(mse0), k_mean = mean(k0)),
          aug = data.frame(
            Ka = aug_cases$Ka,
            Ke = aug_cases$Ke,
            mse_mean = colMeans(mse_aug),
            mse_sd   = apply(mse_aug, 2, sd),
            k_mean   = colMeans(k_aug)
          )
        )
        
        cat(sprintf("Done: n=%d, %s, theta21=%.2f, sigma21=%.2f\n",
                    n, model, theta21, sigma21))
      }
    }
  }
}

# Flatten to a single data.frame summary
summary_df <- do.call(rbind, lapply(results, function(x) {
  base <- x$baseline
  out_base <- data.frame(
    n = x$n, model = x$model,
    theta21 = x$theta21, sigma21 = x$sigma21,
    Ka = NA_integer_, Ke = NA_integer_,
    variant = "baseline",
    mse_mean = base$mse_mean, mse_sd = base$mse_sd, k_mean = base$k_mean
  )
  out_aug <- transform(x$aug,
                       n = x$n, model = x$model, theta21 = x$theta21, sigma21 = x$sigma21,
                       variant = "augmented"
  )
  out_aug <- out_aug[, c("n","model","theta21","sigma21","Ka","Ke","variant","mse_mean","mse_sd","k_mean")]
  rbind(out_base, out_aug)
}))

print(head(summary_df, 20))

# Example: for each design point, find the best (smallest) MSE among augmented cases
best_aug <- do.call(rbind, lapply(results, function(x) {
  idx <- which.min(x$aug$mse_mean)
  data.frame(
    n = x$n, model = x$model, theta21 = x$theta21, sigma21 = x$sigma21,
    best_Ka = x$aug$Ka[idx], best_Ke = x$aug$Ke[idx],
    best_mse = x$aug$mse_mean[idx],
    baseline_mse = x$baseline$mse_mean
  )
}))
print(head(best_aug, 20))

# Save outputs
write.csv(summary_df, "mc_summary_all.csv", row.names = FALSE)
write.csv(best_aug,   "mc_best_aug.csv", row.names = FALSE)
