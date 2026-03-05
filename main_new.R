
df <- read_excel("data.xlsx", sheet = "cay")

Dates <- df$Dates
X <- cbind(w = df$w, y = df$y)

n   <- nrow(df)
tau <- seq(0, 1, length.out = n)

K      <- 1
k_grid <- 1:12

idx <- build_dols_diff_matrix(X, K)$idx

## consumption measure 1: c_pce
cayf_c_pce   <- get_ecm_dols(df$c_pce, X, K)
caytvp_c_pce <- get_ecm_tvp_gcv(df$c_pce, X, tau, k_grid, K)

## consumption measure 2: c_2
cayf_c_2   <- get_ecm_dols(df$c_2, X, K)
caytvp_c_2 <- get_ecm_tvp_gcv(df$c_2, X, tau, k_grid, K)

out <- data.frame(
  Dates        = Dates[idx],
  cayf_c_pce   = cayf_c_pce,
  caytvp_c_pce = caytvp_c_pce,
  cayf_c_2     = cayf_c_2,
  caytvp_c_2   = caytvp_c_2
)

df_eq <- read_excel("data.xlsx", sheet = "equity")

rx_out <- merge(
  out["Dates"],
  data.frame(Dates = df_eq$Dates, rx = df_eq$r),
  by = "Dates",
  all.x = FALSE,
  all.y = FALSE
)

dat <- merge(out, rx_out, by = "Dates")
dat <- dat[order(dat$Dates), ]
dat <- dat[!is.na(dat$rx), ]

## ---------- Panel A: mean, sd, ac(1) ----------
ac1 <- function(x) as.numeric(stats::acf(x, plot = FALSE, lag.max = 1)$acf[2])

panelA <- data.frame(
  series = c("rx", names(out)[names(out) != "Dates"]),
  mean   = c(mean(dat$rx, na.rm = TRUE), sapply(dat[names(out)[names(out) != "Dates"]], mean, na.rm = TRUE)),
  sd     = c(sd(dat$rx, na.rm = TRUE),   sapply(dat[names(out)[names(out) != "Dates"]], sd,   na.rm = TRUE)),
  ac1    = c(ac1(dat$rx),                sapply(dat[names(out)[names(out) != "Dates"]], ac1))
)

## ---------- Panel B: correlation matrix ----------
vars <- c("rx", names(out)[names(out) != "Dates"])
panelB <- cor(dat[, vars], use = "pairwise.complete.obs")

## ---------- Some plot ----------

## ---- read recession indicator from xlsx ----
rec_df <- read_excel("data.xlsx", sheet = "recession")  # <- edit sheet name if needed
rec_df <- data.frame(Dates = rec_df$Dates, rec = rec_df$USRECQ)  # <- edit column names if needed
rec_df <- rec_df[order(rec_df$Dates), ]

## ---- convert 0/1 recession series into shading intervals ----
rec_on <- rec_df$rec == 1
sw <- diff(c(FALSE, rec_on, FALSE))
starts <- rec_df$Dates[which(sw ==  1)]
ends   <- rec_df$Dates[which(sw == -1) - 1]
recessions <- data.frame(xmin = starts, xmax = ends)

## ---- scaling: "units of standardized deviations" (divide by SD, no de-mean) ----
std_units <- function(x) x / sd(x, na.rm = TRUE)

## merge everything (handles time differences automatically)
dat <- merge(out, rx_out, by = "Dates", all = FALSE)
dat <- merge(dat, rec_df, by = "Dates", all = FALSE)
dat <- dat[order(dat$Dates), ]

## ---------- helper for plotting ----------
make_plot <- function(dat, var1, var2, title_txt) {
  
  df <- data.frame(
    Dates = dat$Dates,
    rx  = std_units(dat$rx),
    f   = std_units(dat[[var1]]),
    tvp = std_units(dat[[var2]])
  )
  
  ggplot(df, aes(x = Dates)) +
    geom_rect(
      data = recessions,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey70", alpha = 0.35
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    
    ## excess returns
    geom_line(aes(y = rx, color = "rx"), linewidth = 1.5) +
    
    ## cay series
    geom_line(aes(y = f, color = "f"), linewidth = 2) +
    geom_line(aes(y = tvp, color = "tvp"), linewidth = 2) +
    
    scale_color_manual(
      name = NULL,
      values = c(
        rx  = "grey60",
        f   = "black",
        tvp = "blue"
      ),
      breaks = c("rx", "f", "tvp"),
      labels = c(
        rx  = expression(r[t] - r[f*","*t]),
        f   = expression(widehat(cay)[f*","*t]),
        tvp = expression(widehat(cay)[tvp*","*t])
      )
    ) +
    
    labs(
      x = NULL,
      y = "Standardized Units"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.05, 0.08),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

## ---------- Figure for PCE ----------
p_pce <- make_plot(dat, "cayf_c_pce", "caytvp_c_pce",
                   "cay (PCE) vs Excess Returns")

## ---------- Figure for other consumption ----------
p_c2 <- make_plot(dat, "cayf_c_2", "caytvp_c_2",
                  "cay (Alt Consumption) vs Excess Returns")

print(p_pce)
print(p_c2)

## In-sample predictive regression--------------------------

## horizons exactly as in Table 3
H_set <- c(1, 2, 4, 8, 12, 16, 20)

## merge cay + returns (keeps common dates only; handles differing sample periods)
dat <- merge(out, rx_out, by = "Dates", all = FALSE)
dat <- dat[order(dat$Dates), ]

## forward H-quarter cumulative excess return: sum_{j=1}^H rx_{t+j}
fwd_sum <- function(rx, H) {
  n <- length(rx)
  y <- rep(NA_real_, n)
  for (t in 1:(n - H)) y[t] <- sum(rx[(t + 1):(t + H)], na.rm = FALSE)
  y
}

## run predictive regression with NW (default lag selection)
one_reg <- function(x, rx, H) {
  y <- fwd_sum(rx, H)
  ok <- is.finite(x) & is.finite(y)
  fit <- lm(y[ok] ~ x[ok])
  
  Vnw <- NeweyWest(fit, lag = H-1, prewhite = FALSE, adjust = TRUE)
  ct  <- coeftest(fit, vcov. = Vnw)
  
  beta  <- unname(coef(fit)[2])
  tstat <- unname(ct[2, "t value"])
  adjR2 <- summary(fit)$r.squared
  
  c(beta = beta, t = tstat, adjR2 = adjR2)
}

## predictors in your out
preds <- c("cayf_c_pce", "caytvp_c_pce", "cayf_c_2", "caytvp_c_2")

## compute results
res_list <- lapply(preds, function(p) {
  mat <- sapply(H_set, function(H) one_reg(dat[[p]], dat$rx, H))
  t(mat)  ## rows = H, cols = beta/t/adjR2
})
names(res_list) <- preds

## format like Table 3: beta \n (t) \n [adjR2]
fmt_cell <- function(beta, t, adjR2, d_beta = 3, d_t = 3, d_r2 = 3) {
  paste0(
    formatC(beta,  format = "f", digits = d_beta), "\n",
    "(", formatC(t, format = "f", digits = d_t), ")\n",
    "[", formatC(adjR2, format = "f", digits = d_r2), "]"
  )
}

table3_like <- data.frame(Predictor = preds, stringsAsFactors = FALSE)
for (j in seq_along(H_set)) {
  H <- H_set[j]
  colname <- paste0("H=", H)
  table3_like[[colname]] <- sapply(preds, function(p) {
    beta  <- res_list[[p]][j, "beta"]
    tstat <- res_list[[p]][j, "t"]
    adjR2 <- res_list[[p]][j, "adjR2"]
    fmt_cell(beta, tstat, adjR2)
  })
}

## write to Excel (single sheet)
write_xlsx(list(Table3 = table3_like), "Table3_my_cay.xlsx")



