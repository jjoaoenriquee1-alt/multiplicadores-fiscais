# ============================================================
# Fiscal multipliers in Brazil (1997-2019) - paper figures
# 1) Monthly tvVAR  2) Quarterly tvVAR  3) Quarterly TVAR (regimes)
# 4) Correlation of the social multiplier with social variables
#
# Run with the working directory set to the repository folder
# (the folder that contains this script and the data files).
# ============================================================

library(tvReg)
library(vars)
library(tsDyn)      # from https://matthieustigler.r-universe.dev
library(readxl)
library(seasonal)
library(showtext)
library(tseries)
library(writexl)

# In RStudio, move to the script's folder automatically
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  p_script <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  if (nzchar(p_script)) setwd(dirname(p_script))
}
stopifnot("Data files not found: set the working directory to the repository folder" =
            file.exists("dadosBRFDesagreg.xlsx"))

# Font: Mena Grotesk if its .otf files are in ./fonts, otherwise the default sans
font_dir   <- "fonts"
font_files <- file.path(font_dir, c("Mena Grotesk Book.otf", "Mena Grotesk Bold.otf",
                                    "Mena Grotesk Book Italic.otf", "Mena Grotesk Bold Italic.otf"))
if (all(file.exists(font_files))) {
  FONT <- "MenaGrotesk"
  font_add(FONT, regular = font_files[1], bold = font_files[2],
           italic = font_files[3], bolditalic = font_files[4])
} else {
  FONT <- "sans"
  message("Mena Grotesk not found in ./fonts - using the default sans font.")
}
showtext_auto()

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
sa_all  <- function(x) sapply(colnames(x), function(v) final(seas(x[, v])), simplify = FALSE)
dlog    <- function(g, rec, y) diff(log(cbind(dg = g, drec = rec, dpib = y)))
aic_lag <- function(d, lmax) VARselect(d, lag.max = lmax, type = "const")$selection["AIC(n)"]
irf_h   <- function(irf, h, resp = "dpib") irf$irf[["dg"]][, resp, h + 1]
dates_m <- function(n) tail(seq(as.Date("1997-02-01"), by = "month",   length.out = 275), n)
dates_q <- function(n) tail(seq(as.Date("1997-04-01"), by = "quarter", length.out = 91),  n)

# Trailing moving average (expanding mean for the first k-1 obs); k = 1 -> no smoothing
ma <- function(x, k) {
  if (k == 1) return(as.numeric(x))
  s <-as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))
  s[1:(k - 1)] <- cumsum(x[1:(k - 1)]) / seq_len(k - 1)
  s
}

# Multiplier = (cum. GDP response / cum. spending response) / smoothed G/Y
mult_ma <- function(irf, gy, h, k) {
  elast <- irf_h(irf, h) / irf_h(irf, h, "dg")
  as.numeric(elast / ma(tail(as.numeric(gy), length(elast)), k))
}

# ============================================================
# 1) MONTHLY tvVAR (Feb/1997 - Dec/2019)
# ============================================================
fis  <- read_excel("dadosBRFDesagreg.xlsx")
defl <- fis$deflator
pib  <- read.csv2("pib_mensal_4380.csv.csv", fileEncoding = "latin1",
                  stringsAsFactors = FALSE)[[3]][1:276]
pib  <- as.numeric(gsub(",", ".", gsub("\\.", "", pib))) * tail(defl, 1) / defl[1:276]

m <- sa_all(ts(cbind(g = fis$gastoReal[1:276], rec = fis$recreal[1:276],
                     gs = fis$gsocialreal[1:276], y = pib),
               start = c(1997, 1), frequency = 12))

ym  <- list(primary = dlog(m$g, m$rec, m$y), social = dlog(m$gs, m$rec, m$y))
gym <- list(primary = m$g / m$y, social = m$gs / m$y)
pm  <- lapply(ym, aic_lag, lmax = 12)

tvm   <- Map(function(d, p) tvVAR(d, p = p, type = "const", est = "lc"), ym, pm)
irm   <- lapply(tvm, tvIRF, cumulative = TRUE, n.ahead = 24)
irm36 <- lapply(tvm, tvIRF, cumulative = TRUE, n.ahead = 36)
mm    <- Map(mult_ma, irm, gym, h = 12, k = 12)

# Bootstrap CIs (68% and 90%), cached
f_ic <- "ic_bootstrap_68_90.rds"
if (file.exists(f_ic)) {
  ic <- readRDS(f_ic)
} else {
  ci <- function(x, lv) { set.seed(20261008); confint(x, runs = 1000, level = lv) }
  ic <- list(g68 = ci(irm$primary, 0.6827), g90 = ci(irm$primary, 0.90),
             s68 = ci(irm$social,  0.6827), s90 = ci(irm$social,  0.90))
  saveRDS(ic, f_ic)
}

# Robustness: exogenous controls (real exchange rate, 12m inflation)
ex   <- ts(read_excel("Exogenas.xlsx"), start = c(1997, 1), end = c(2019, 12), frequency = 12)
exog <- data.frame(dpi = diff(log(ex[, "icpca12"])), dca = diff(log(ex[, "caRe"])))
tvr  <- Map(function(d, p) tvVAR(d, p = p, type = "const", est = "lc", exogen = exog), ym, pm)
mr   <- Map(mult_ma, lapply(tvr, tvIRF, cumulative = TRUE, n.ahead = 24), gym, h = 12, k = 12)

# Robustness: bandwidth
bws    <- list(NULL, rep(1.96, 3), rep(0.5, 3), rep(0.3, 3))
cols_b <- c("#45ff66", "#EB52ff", "#3366ff", "#feff41")
bw_mult <- function(d, p, gy) lapply(bws, function(b) {
  mod <- if (is.null(b)) tvVAR(d, p = p, type = "const") else tvVAR(d, p = p, type = "const", bw = b)
  list(y   = mult_ma(tvIRF(mod, cumulative = TRUE, n.ahead = 24), gy, 12, 12),
       lab = paste0("bw = [", paste(format(round(mod$bw, 3), trim = TRUE), collapse = ", "), "]",
                    if (is.null(b)) " (cross-val)"))
})
mb <- Map(bw_mult, ym, pm, gym)

# ============================================================
# 2) QUARTERLY tvVAR (1997Q2 - 2019Q4)
#    Seasonal adjustment window 1997Q1-2019Q4 (pandemic left out).
#    The TVAR in section 3 reuses these same SA series.
# ============================================================
qraw <- read_excel("dadosBR_trimestral.xlsx")
q <- sa_all(ts(cbind(g  = qraw[["Despesa primária"]],  rec = qraw[["Receita primária"]],
                     gs = qraw[["Benefícios sociais"]], y  = qraw[["PIB"]]),
               start = c(1997, 1), end = c(2019, 4), frequency = 4))

yq  <- list(primary = dlog(q$g, q$rec, q$y), social = dlog(q$gs, q$rec, q$y))
gyq <- list(primary = q$g / q$y, social = q$gs / q$y)
pq  <- lapply(yq, aic_lag, lmax = 8)
tvq <- Map(function(d, p) tvVAR(d, p = p, type = "const", est = "lc"), yq, pq)
mq  <- Map(mult_ma, lapply(tvq, tvIRF, cumulative = TRUE, n.ahead = 8), gyq, h = 4, k = 4)

# ============================================================
# RESIDUAL DIAGNOSTICS - auxiliary constant-coefficient VAR
# (same data and AIC lags used for lag selection)
# ============================================================
test_row <- function(model, test, h)
  data.frame(model = model, test = test, stat = round(unname(h$statistic), 3),
             df = unname(h$parameter[1]), p_value = round(h$p.value, 4))

var_tests <- function(d, p, lags_pt, lags_arch, model) {
  v  <- vars::VAR(d, p = as.integer(p), type = "const")
  nt <- vars::normality.test(v, multivariate.only = FALSE)
  rbind(
    test_row(model, "Portmanteau (adjusted)",       # H0: no serial correlation
             vars::serial.test(v, lags.pt = lags_pt, type = "PT.adjusted")$serial),
    test_row(model, "ARCH-LM (multivariate)",       # H0: homoskedasticity
             vars::arch.test(v, lags.multi = lags_arch, multivariate.only = TRUE)$arch.mul),
    test_row(model, "Jarque-Bera (multivariate)", nt$jb.mul$JB),   # H0: normality
    do.call(rbind, Map(function(h, nm) test_row(model, paste("Jarque-Bera", nm), h),
                       nt$jb.uni, colnames(d))),
    data.frame(model = model, test = "Max root modulus (stable if < 1)",
               stat = round(max(vars::roots(v)), 4), df = NA, p_value = NA))
}

var_diag <- rbind(
  var_tests(ym$primary, pm$primary, 24, 12, "Monthly - primary"),
  var_tests(ym$social,  pm$social,  24, 12, "Monthly - social"),
  var_tests(yq$primary, pq$primary, 12, 4,  "Quarterly - primary"),
  var_tests(yq$social,  pq$social,  12, 4,  "Quarterly - social"))
rownames(var_diag) <- NULL
print(var_diag)
write.csv(var_diag, "var_residual_tests.csv", row.names = FALSE)

# ============================================================
# 3) QUARTERLY TVAR - social spending by regime (1998Q1 - 2019Q4)
#    VAR (dg, drec, dpib), dg = social spending
#    External transition variables (thVar): Gini, poverty, extreme poverty
#    Lag order: HQ (fewer parameters per regime than AIC)
#    Multiplier (Pires, 2014; Hall, 2009): dY/dX = elasticity / (X/Y)
# ============================================================

# Same SA series as the quarterly tvVAR (section 2)
dfgs <- window(dlog(q$gs, q$rec, q$y), start = c(1998, 1), end = c(2019, 4))
XY   <- window(q$gs / q$y,              start = c(1998, 1), end = c(2019, 4))

trans_ts <- function(tab, col) {
  yr <- as.integer(substr(tab$Periodo, 1, 4))
  qt <- as.integer(substr(tab$Periodo, 6, 6))
  window(ts(tab[[col]] / 100, start = c(yr[1], qt[1]), frequency = 4),
         start = c(1998, 1), end = c(2019, 4))
}
gini  <- read_excel("GiniBR_trimestral.xlsx")
pob   <- read_excel("Pobreza_trimestral.xlsx")
trans <- list(gini    = trans_ts(gini, "Gini_trim"),
              poverty = trans_ts(pob,  "Pobreza"),
              extreme = trans_ts(pob,  "Extrema_Pobreza"))
labs  <- list(gini    = c("Gini",            "Low Inequality",      "High Inequality"),
              poverty = c("Poverty",         "Low Poverty",         "High Poverty"),
              extreme = c("Extreme Poverty", "Low Extreme Poverty", "High Extreme Poverty"))

# thVar must match the VAR sample exactly
for (z in trans) stopifnot(length(z) == nrow(dfgs), all(time(z) == time(dfgs)))
stopifnot(all(time(XY) == time(dfgs)))
cat("TVAR sample:", nrow(dfgs), "quarters\n")

# Lag selection (HQ); report all criteria in the paper
sel_tv <- VARselect(dfgs, lag.max = 8, type = "const")
print(sel_tv$selection)
p_tv <- as.integer(sel_tv$selection["HQ(n)"])
cat("TVAR lag (HQ):", p_tv, "\n")

TH_DELAY <- 1       # regime given by transition variable at t-1
TRIM     <- 0.15
H_TV     <- 0:9
NBOOT    <- 500
CI_PLOT  <- 0.68    # band level in the TVAR figures

tvar <- lapply(trans, function(z) TVAR(dfgs, lag = p_tv, include = "const", nthresh = 1,
                                       thVar = z, thDelay = TH_DELAY, trim = TRIM,
                                       plot = FALSE, trace = FALSE))

for (v in names(tvar)) {
  cat("\n", labs[[v]][1], "- threshold:", round(getTh(tvar[[v]]), 4), "\n")
  print(table(regime(tvar[[v]]), useNA = "ifany"))
}

# Regime of each quarter (1 = low, 2 = high), aligned to dfgs
reg_al <- function(mod, n = nrow(dfgs)) {
  r <- as.numeric(regime(mod))
  if (length(r) < n) r <- c(rep(NA, n - length(r)), r)
  r
}
regs <- sapply(tvar, reg_al)

# Do the three transition variables classify quarters the same way?
print(table(gini = regs[, "gini"],       poverty = regs[, "poverty"]))
print(table(gini = regs[, "gini"],       extreme = regs[, "extreme"]))
print(table(poverty = regs[, "poverty"], extreme = regs[, "extreme"]))

# Cumulative IRF within regime r (Cholesky: dg, drec, dpib; dg +1% on impact).
# Linear IRF per regime: economy assumed to stay in regime r over the horizon.
fri_regime <- function(mod, vars, r, response, n.ahead = max(H_TV), shock = 0.01) {
  k <- length(vars)
  B <- coef(mod)[[r]]
  A <- B[, !grepl("Intercept|Trend", colnames(B)), drop = FALSE]
  p <- ncol(A) / k
  stopifnot(p == round(p))
  Alag <- lapply(seq_len(p), function(i) A[, ((i - 1) * k + 1):(i * k), drop = FALSE])
  Phi <- vector("list", n.ahead + 1)
  Phi[[1]] <- diag(k)
  for (h in seq_len(n.ahead)) {
    Phi[[h + 1]] <- matrix(0, k, k)
    for (j in seq_len(min(h, p))) Phi[[h + 1]] <- Phi[[h + 1]] + Alag[[j]] %*% Phi[[h + 1 - j]]
  }
  res <- residuals(mod)
  reg <- regime(mod)
  reg <- tail(reg[!is.na(reg)], nrow(res))
  Sig <- crossprod(res[reg == r, , drop = FALSE]) / sum(reg == r)
  P   <- t(chol(Sig))
  i   <- match("dg", vars)
  e   <- P[, i] * shock / P[i, i]
  cumsum(sapply(Phi, function(Ph) (Ph %*% e)[match(response, vars)]))
}

# Residual bootstrap with fixed regimes (transition variable is external).
# reest_th = TRUE re-estimates the threshold in each draw (robustness:
# adds threshold uncertainty, so bands get wider).
boot_tvar <- function(mod, dados, z, p, nboot = NBOOT, h = H_TV, reest_th = FALSE,
                      seed = 123) {
  set.seed(seed)
  vars <- colnames(dados); k <- length(vars)
  Y    <- as.matrix(dados); n <- nrow(Y)
  B    <- lapply(coef(mod), function(b) b[, !grepl("Trend", colnames(b)), drop = FALSE])
  res  <- residuals(mod)
  reg  <- regime(mod)
  reg  <- tail(as.numeric(reg[!is.na(reg)]), nrow(res))   # regime for t = p+1, ..., n
  th   <- getTh(mod)
  for (r in 1:2) {                                         # center within regime
    idx <- which(reg == r)
    res[idx, ] <- sweep(res[idx, , drop = FALSE], 2, colMeans(res[idx, , drop = FALSE]))
  }
  one <- function() {
    u <- matrix(0, nrow(res), k)
    for (r in 1:2) {
      idx <- which(reg == r)
      u[idx, ] <- res[idx[sample.int(length(idx), length(idx), replace = TRUE)], ]
    }
    ys <- Y
    for (t in (p + 1):n) {
      x <- c(1, as.vector(t(ys[(t - 1):(t - p), , drop = FALSE])))
      ys[t, ] <- B[[reg[t - p]]] %*% x + u[t - p, ]
    }
    colnames(ys) <- vars
    mb <- tryCatch(
      if (reest_th)
        TVAR(ys, lag = p, include = "const", nthresh = 1, thVar = z, thDelay = TH_DELAY,
             trim = TRIM, plot = FALSE, trace = FALSE)
      else
        TVAR(ys, lag = p, include = "const", nthresh = 1, thVar = z, thDelay = TH_DELAY,
             gamma = th, plot = FALSE, trace = FALSE),
      error = function(e) NULL)
    if (is.null(mb)) return(NULL)
    out <- tryCatch(lapply(1:2, function(r) list(y = fri_regime(mb, vars, r, "dpib", max(h)),
                                                 g = fri_regime(mb, vars, r, "dg",   max(h)))),
                    error = function(e) NULL)
    if (is.null(out) || any(!is.finite(unlist(out)))) NULL else out
  }
  draws <- Filter(Negate(is.null), replicate(nboot, one(), simplify = FALSE))
  cat("Valid bootstrap draws:", length(draws), "of", nboot, "\n")
  pick  <- function(r, s) do.call(rbind, lapply(draws, function(d) d[[r]][[s]]))
  list(y = list(pick(1, "y"), pick(2, "y")), g = list(pick(1, "g"), pick(2, "g")))
}

bands <- function(d, lv = 0.90) apply(d, 2, quantile, probs = c((1 - lv) / 2, (1 + lv) / 2), na.rm = TRUE)

fri  <- lapply(tvar, function(mod) lapply(1:2, function(r)
  list(y = fri_regime(mod, colnames(dfgs), r, "dpib"), g = fri_regime(mod, colnames(dfgs), r, "dg"))))
boot <- Map(boot_tvar, tvar, list(dfgs), trans, p = p_tv)

# TVAR multipliers at h = 4 and h = 8
#   mult_avgXY: elasticity / sample-mean X/Y
#   mult_regXY: elasticity / mean X/Y of the regime's quarters
#   68% and 90% percentile bands for both
XY_avg <- mean(XY)
k_h    <- match(c(4, 8), H_TV)

ci_cols <- function(d, pref) do.call(cbind, lapply(c(68, 90), function(lv) {
  b <- bands(d, lv / 100)[, k_h, drop = FALSE]
  setNames(data.frame(b[1, ], b[2, ]), paste0(pref, c("_lo", "_hi"), lv))
}))

tvar_mult <- do.call(rbind, lapply(names(tvar), function(v) {
  do.call(rbind, lapply(1:2, function(r) {
    xy_r  <- mean(XY[which(regs[, v] == r)])
    eps   <- fri[[v]][[r]]$y / fri[[v]][[r]]$g
    eps_b <- boot[[v]]$y[[r]] / boot[[v]]$g[[r]]
    cbind(data.frame(transition = labs[[v]][1], regime = labs[[v]][r + 1], h = c(4, 8),
                     n_quarters = sum(regs[, v] == r, na.rm = TRUE),
                     XY_regime  = xy_r,
                     elasticity = eps[k_h],
                     mult_avgXY = eps[k_h] / XY_avg),
          ci_cols(eps_b / XY_avg, "avg"),
          data.frame(mult_regXY = eps[k_h] / xy_r),
          ci_cols(eps_b / xy_r, "reg"))
  }))
}))
rownames(tvar_mult) <- NULL
print(tvar_mult, digits = 3)
write.csv(tvar_mult, "tvar_multipliers.csv", row.names = FALSE)

# ------------------------------------------------------------
# TVAR residual diagnostics (per equation)
#   Ljung-Box (H0: no autocorrelation), Jarque-Bera (H0: normality)
#   Ljung-Box within each regime as an extra check
# ------------------------------------------------------------
tvar_tests <- function(mod, model, p, lags = c(8, 12)) {
  res <- as.matrix(residuals(mod)); colnames(res) <- colnames(dfgs)
  reg <- regime(mod)
  reg <- tail(as.numeric(reg[!is.na(reg)]), nrow(res))
  do.call(rbind, lapply(colnames(res), function(nm) rbind(
    do.call(rbind, lapply(lags, function(l)
      test_row(model, paste0("Ljung-Box ", nm, " (lag ", l, ")"),
               Box.test(res[, nm], lag = l, type = "Ljung-Box", fitdf = p)))),
    do.call(rbind, lapply(1:2, function(r)
      test_row(model, paste0("Ljung-Box ", nm, " (lag ", lags[1], ", regime ", r, ")"),
               Box.test(res[reg == r, nm], lag = lags[1], type = "Ljung-Box", fitdf = p)))),
    test_row(model, paste("Jarque-Bera", nm), jarque.bera.test(res[, nm])))))
}

tvar_diag <- do.call(rbind, lapply(names(tvar), function(v)
  tvar_tests(tvar[[v]], paste("TVAR -", labs[[v]][1]), p_tv)))
rownames(tvar_diag) <- NULL
print(tvar_diag)
write.csv(tvar_diag, "tvar_residual_tests.csv", row.names = FALSE)

# ============================================================
# PLOT FUNCTIONS
# ============================================================
std_par <- function() par(mfrow = c(1, 1), mar = c(4, 4, 3, 1), las = 0, family = FONT)
qlim    <- function(...) range(sapply(list(...), quantile, c(0.01, 0.99), na.rm = TRUE), 0, 1)
srs     <- function(x, y, lab, col, lty = 1, lwd = 2, type = "l")
  list(x = x, y = y, lab = lab, col = col, lty = lty, lwd = lwd, type = type)

plot_ci <- function(c68, c90, main, zero = FALSE) {
  pick <- function(o, s) o[[s]][["dg"]][, "dpib", 13]
  ef <- pick(c68, "irf"); x <- dates_m(length(ef))
  lo68 <- pick(c68, "Lower"); hi68 <- pick(c68, "Upper")
  lo90 <- pick(c90, "Lower"); hi90 <- pick(c90, "Upper")
  std_par()
  plot(x, ef, type = "n", main = main, ylab = "Accumulated GDP change", xlab = "",
       ylim = range(lo90, hi90, if (zero) 0, na.rm = TRUE))
  polygon(c(x, rev(x)), c(hi90, rev(lo90)), col = adjustcolor("grey80", 0.3), border = NA)
  polygon(c(x, rev(x)), c(hi68, rev(lo68)), col = adjustcolor("grey50", 0.5), border = NA)
  abline(h = 0, lty = 2, col = "red", lwd = 1.5)
  lines(x, ef, lwd = 2)
}

plot_lines <- function(s, ylim = NULL, main = "", ylab = "Multiplier", refs = TRUE,
                       pos = "topright", ncol = 1, cex = 0.85) {
  if (is.null(ylim)) ylim <- range(unlist(lapply(s, `[[`, "y")), na.rm = TRUE)
  std_par()
  plot(s[[1]]$x, s[[1]]$y, type = "n", main = main, ylab = ylab, xlab = "", ylim = ylim)
  if (refs) {
    abline(h = 0, lty = 2, col = "red",  lwd = 1.2)
    abline(h = 1, lty = 3, col = "blue", lwd = 1.2)
  }
  for (i in s) lines(i$x, i$y, col = i$col, lty = i$lty, lwd = i$lwd, type = i$type, pch = 16, cex = 0.7)
  pick <- function(f) sapply(s, `[[`, f)
  legend(pos, legend = pick("lab"), col = pick("col"), lty = pick("lty"), lwd = pick("lwd"),
         pch = ifelse(pick("type") == "o", 16, NA), bty = "n", cex = cex, ncol = ncol)
}

plot_fri_ic <- function(low, high, h, ic_low, ic_high, title = "",
                        labels = c("Low", "High"), alpha = 0.18) {
  par(mfrow = c(1, 1), family = FONT, mar = c(7, 5, 3, 2), las = 1)
  yl <- range(c(0, low, high, ic_low, ic_high), na.rm = TRUE)
  plot(h, low, type = "n", ylim = yl, xaxt = "n", bty = "n", xlab = "", ylab = "", main = title)
  axis(1, at = h, tick = FALSE)
  title(xlab = "Quarters", line = 2.8)
  abline(h = pretty(yl), col = "grey85")
  abline(h = 0, col = "grey60")
  cols <- c("black", "grey55")
  polygon(c(h, rev(h)), c(ic_low[1, ],  rev(ic_low[2, ])),  col = adjustcolor(cols[1], alpha), border = NA)
  polygon(c(h, rev(h)), c(ic_high[1, ], rev(ic_high[2, ])), col = adjustcolor(cols[2], alpha), border = NA)
  for (r in 1:2) {
    lines(h, ic_low[r, ],  lty = 2, col = cols[1])
    lines(h, ic_high[r, ], lty = 2, col = cols[2])
  }
  lines(h, low,  lwd = 3, col = cols[1])
  lines(h, high, lwd = 3, col = cols[2])
  y_leg <- grconvertY(grconvertY(par("usr")[3], "user", "inches") - 5 * par("csi"), "inches", "user")
  legend(x = mean(range(h)), y = y_leg, xjust = 0.5, yjust = 0.5, xpd = NA, horiz = TRUE,
         bty = "n", lwd = 3, legend = labels, col = cols)
}

plot_threshold <- function(z, mod, lab) {
  std_par()
  plot(z, lwd = 2, ylab = lab, xlab = "",
       main = paste(lab, "and estimated threshold (TVAR - social spending)"))
  abline(h = getTh(mod), lty = 2, col = "red")
  legend("topright", bty = "n", lty = c(1, 2), col = c("black", "red"), lwd = c(2, 1),
         legend = c(lab, paste0("Threshold = ", round(getTh(mod), 4))))
}

# ============================================================
# FIGURES
# ============================================================

# Accumulated effect (h = 12) with 68% / 90% CI
plot_ci(ic$g68, ic$g90, "Accumulated Effect (h=12) - Primary Spending")
plot_ci(ic$s68, ic$s90, "Accumulated Effect (h=12) - Social Spending", zero = TRUE)

# Fig. 2 - Dynamic fiscal multiplier (G/Y MA12)
n2 <- min(lengths(mm))
plot_lines(list(srs(dates_m(n2), tail(mm$primary, n2), "Primary Spending", "#45ff66", lty = 5),
                srs(dates_m(n2), tail(mm$social,  n2), "Social Spending",  "#EB52ff", lty = 5)),
           ylim = qlim(tail(mm$primary, n2), tail(mm$social, n2)),
           main = "Dynamic Fiscal Multiplier", ylab = "R$ GDP per R$ 1 of spending")

# TVAR - transition variables and estimated thresholds
for (v in names(tvar)) plot_threshold(trans[[v]], tvar[[v]], labs[[v]][1])

# TVAR - cumulative GDP response by regime (CI_PLOT bands)
for (v in names(tvar))
  plot_fri_ic(fri[[v]][[1]]$y, fri[[v]][[2]]$y, H_TV,
              bands(boot[[v]]$y[[1]], CI_PLOT), bands(boot[[v]]$y[[2]], CI_PLOT),
              paste0("Cumulative response of GDP to a social spending shock (", labs[[v]][1], ")"),
              labels = labs[[v]][2:3])

# Figs. 7-8 - Multipliers with and without controls
ctrl <- function(v, col, ...) {
  n <- min(length(mm[[v]]), length(mr[[v]])); x <- dates_m(n)
  plot_lines(list(srs(x, tail(mm[[v]], n), "Without controls", "black"),
                  srs(x, tail(mr[[v]], n), "With controls (real exchange rate, inflation)", col, lty = 2)),
             ...)
}
ctrl("primary", "#45ff66", ylim = qlim(mm$primary, mr$primary), refs = FALSE, pos = "bottomright")
ctrl("social",  "#EB52ff")

# Figs. 9-10 - Monthly vs quarterly multipliers (smoothed G/Y)
mvq <- function(v)
  plot_lines(list(srs(dates_m(length(mm[[v]])), mm[[v]], "Monthly",   "black",   lwd = 1.3),
                  srs(dates_q(length(mq[[v]])), mq[[v]], "Quarterly", "#3366ff", lwd = 2.2, type = "o")),
             ylab = "R$ GDP per R$ 1 of spending")
mvq("primary")
mvq("social")

# Figs. 11-12 - Multipliers at different horizons
hz     <- c(6, 12, 18, 24, 36)
cols_h <- c("#45ff66", "#EB52ff", "#3366ff", "#feff41", "#FFAB4D")
by_h <- function(v, pos) {
  x <- dates_m(dim(irm36[[v]]$irf[["dg"]])[1])
  plot_lines(Map(function(h, col) srs(x, mult_ma(irm36[[v]], gym[[v]], h, 12), paste0("h=", h, "m"), col),
                 hz, cols_h),
             pos = pos, ncol = 2, cex = 0.75)
}
by_h("primary", "topleft")
by_h("social",  "topright")

# Figs. 13-14 - Multipliers with different bandwidths
by_bw <- function(v, pos) {
  x <- dates_m(length(mb[[v]][[1]]$y))
  plot_lines(Map(function(b, col) srs(x, b$y, b$lab, col, lwd = 2.5), mb[[v]], cols_b),
             pos = pos, cex = 1)
}
by_bw("primary", "bottomleft")
by_bw("social",  "bottomright")




# ============================================================
# 4) CORRELATION: SOCIAL MULTIPLIER vs SOCIAL VARIABLES (1997-2019)
#    Annual social-spending multipliers from the monthly tvVAR (h = 12)
#    Denominator: monthly G/Y (k = 1) or smoothed G/Y MA12 (k = 12)
#    Annual aggregation: mean of the year or December value
#    Social variables (annual) come from base_padronizada.xlsx
# ============================================================
MULT_FROM_MODEL <- TRUE   # FALSE -> use the mult_social_* columns stored in the base

mult_vars <- c("mult_social_mensal_media", "mult_social_mensal_dez",
               "mult_social_ma12_media",   "mult_social_ma12_dez")
soc_vars  <- c("top_1_percent", "top_10_percent", "pobreza", "extrema_pobreza", "gini")
soc_labs  <- c(top_1_percent = "Top 1%", top_10_percent = "Top 10%", pobreza = "Poverty",
               extrema_pobreza = "Extreme Poverty", gini = "Gini")

# Monthly social multipliers (same tvVAR / IRF as Fig. 2)
ms  <- list(mensal = mult_ma(irm$social, gym$social, h = 12, k = 1),
            ma12   = mm$social)
dts <- dates_m(length(ms$mensal))
ano <- as.integer(format(dts, "%Y"))
dez <- format(dts, "%m") == "12"

mult_anual <- data.frame(ano = sort(unique(ano)))
for (d in names(ms)) {
  mult_anual[[paste0("mult_social_", d, "_media")]] <- as.numeric(tapply(ms[[d]], ano, mean))
  mult_anual[[paste0("mult_social_", d, "_dez")]]   <- ms[[d]][dez][match(mult_anual$ano, ano[dez])]
}

base_soc <- as.data.frame(read_excel("base_padronizada.xlsx"))
if (!MULT_FROM_MODEL) mult_anual <- base_soc[, c("ano", mult_vars)]
dados_cor <- merge(mult_anual[, c("ano", mult_vars)], base_soc[, c("ano", soc_vars)], by = "ano")
dados_cor <- na.omit(dados_cor[dados_cor$ano >= 1997 & dados_cor$ano <= 2019, ])
cat("Correlation sample:", min(dados_cor$ano), "-", max(dados_cor$ano),
    "(", nrow(dados_cor), "years )\n")

# Correlation and p-value matrices (Pearson)
ct       <- function(m, s) cor.test(dados_cor[[m]], dados_cor[[s]], method = "pearson")
cor_sub  <- cor(dados_cor[, mult_vars], dados_cor[, soc_vars])
p_matrix <- sapply(soc_vars, function(s) sapply(mult_vars, function(m) ct(m, s)$p.value))

# Long table with significance stars
cor_res <- expand.grid(Multiplicador = mult_vars, Variavel_Social = soc_vars,
                       stringsAsFactors = FALSE)
idx <- cbind(cor_res$Multiplicador, cor_res$Variavel_Social)
cor_res$Correlacao    <- cor_sub[idx]
cor_res$P_valor       <- p_matrix[idx]
cor_res$N             <- nrow(dados_cor)
cor_res$Significancia <- as.character(cut(cor_res$P_valor, c(-Inf, 0.01, 0.05, 0.10, Inf),
                                          labels = c("***", "**", "*", "ns"), right = FALSE))
cor_res$Correlacao_Formatada <- paste0(round(cor_res$Correlacao, 3),
                                       ifelse(cor_res$Significancia == "ns", "", cor_res$Significancia))
cor_res <- cor_res[order(cor_res$Multiplicador, -abs(cor_res$Correlacao)), ]
rownames(cor_res) <- NULL
print(cor_res, digits = 3)

as_tab <- function(m) data.frame(Multiplicador = rownames(m), m, check.names = FALSE)
write_xlsx(list(Resultados_Completos = cor_res,
                Matriz_Correlacao    = as_tab(cor_sub),
                Matriz_Pvalor        = as_tab(p_matrix),
                Base_Anual           = dados_cor),
           "correlacoes_denominador_variando.xlsx")

# ------------------------------------------------------------
# Figure - social variables vs social multiplier (G/Y MA12, annual mean)
# Panels ordered by |r|; multiplier rescaled onto the left axis,
# right axis shows its original values
# ------------------------------------------------------------
mult_foco <- "mult_social_ma12_media"
ord <- names(sort(abs(cor_sub[mult_foco, ]), decreasing = TRUE))

par(family = FONT, mfrow = c(length(ord), 1), mar = c(3, 4.5, 2.5, 4),
    oma = rep(0.5, 4), las = 1)
for (s in ord) {
  d  <- dados_cor[, c("ano", mult_foco, s)]
  rs <- range(d[[s]]); rmu <- range(d[[mult_foco]])
  fator <- diff(rs) / diff(rmu)
  const <- rs[1] - rmu[1] * fator

  plot(d$ano, d[[s]], type = "l", col = "#FFAB4D", lwd = 1.5, ylab = soc_labs[s],
       xlab = "", ylim = rs, axes = FALSE, cex.lab = 0.9)
  axis(1, at = seq(1997, 2019, by = 2), cex.axis = 0.8)
  axis(2, cex.axis = 0.8)
  abline(h = 0, lty = 2, col = "red", lwd = 1.5)
  lines(d$ano, d[[mult_foco]] * fator + const, col = "#EB52ff", lwd = 1.5, lty = 2)
  tk <- pretty(rmu); tk <- tk[tk >= rmu[1] & tk <= rmu[2]]
  axis(4, at = tk * fator + const, labels = tk, cex.axis = 0.8)
  box()

  r_i <- cor_res[cor_res$Multiplicador == mult_foco & cor_res$Variavel_Social == s, ]
  mtext(paste0("r = ", r_i$Correlacao_Formatada), side = 3, adj = 0, line = 0.3, cex = 0.65)
  legend("topright", legend = c(soc_labs[s], "Social Multiplier MA12"),
         col = c("#FFAB4D", "#EB52ff"), lwd = 1.5, lty = c(1, 2), bty = "n",
         cex = 0.7, inset = c(0.01, -0.2), y.intersp = 0.5, xpd = NA)
}
par(mfrow = c(1, 1))
