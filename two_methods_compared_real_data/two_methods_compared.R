script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  script_path <- normalizePath("two_methods_compared.R", mustWork = FALSE)
}

setwd("/Users/nominora/Desktop/MRCE")
pkgload::load_all(".")

## MRCE vs PN comparison with per-mouse 10/5 splits, 5-fold CV (8/2 within the 10), 100 reps
suppressPackageStartupMessages({
  library(readxl)
  library(MRCE)
  library(ggplot2)
})


t_start <- Sys.time()

message("Reading data...")
data_path <- file.path(dirname(script_path), "Data_Cortex_Nuclear.xls")
df <- read_excel(data_path)
df <- as.data.frame(df)

# Factorize
for (nm in c("Genotype", "Treatment", "Behavior")) df[[nm]] <- factor(df[[nm]])

# Protein cols 2:78
stopifnot(ncol(df) >= 78)
y_cols <- 2:78

# Median impute Y
for (j in y_cols) {
  v <- as.numeric(df[[j]])
  v[!is.finite(v)] <- NA
  med <- median(v, na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  v[is.na(v)] <- med
  df[[j]] <- v
}

Y_full <- as.matrix(df[, y_cols, drop = FALSE])
X_full <- model.matrix(~ Genotype + Treatment + Behavior - 1, data = df)

stopifnot(all(is.finite(X_full)), all(is.finite(Y_full)))

# Mouse grouping: strip suffix after underscore
df$MouseID_base <- sub("_.*$", "", as.character(df$MouseID))

# Train/test split helper: for every mouse with >=15 rows, sample 10 for train (for CV) and 5 for test
draw_train_test <- function() {
  mouse_tab <- table(df$MouseID_base)
  mouse_ids <- names(mouse_tab)[mouse_tab >= 15]
  if (!length(mouse_ids)) stop("No MouseID_base with >=15 measurements; check data.")
  train_idx <- integer(0)
  test_idx <- integer(0)
  train_mouse <- character(0)
  for (id in mouse_ids) {
    rows <- which(df$MouseID_base == id)
    rows <- rows[1:15]
    tr <- sample(rows, 10, replace = FALSE) # 10 used for CV training
    te <- setdiff(rows, tr) # remaining 5 for held-out test
    train_idx <- c(train_idx, tr)
    test_idx <- c(test_idx, te)
    train_mouse <- c(train_mouse, rep(id, length(tr)))
  }
  list(train_idx = train_idx, test_idx = test_idx, train_mouse = train_mouse, mouse_ids = mouse_ids)
}

# Helper: scale using training stats
scale_with <- function(X, Y, rows) {
  centerY <- colMeans(Y[rows, , drop = FALSE])
  scaleY <- apply(Y[rows, , drop = FALSE], 2, sd)
  centerX <- colMeans(X[rows, , drop = FALSE])
  scaleX <- apply(X[rows, , drop = FALSE], 2, sd)
  scaleY[scaleY == 0] <- 1
  scaleX[scaleX == 0] <- 1
  list(
    X_tr = scale(X[rows, , drop = FALSE], center = centerX, scale = scaleX),
    Y_tr = scale(Y[rows, , drop = FALSE], center = centerY, scale = scaleY),
    centerX = centerX, scaleX = scaleX,
    centerY = centerY, scaleY = scaleY
  )
}

# Fold assignment on training mice: each mouse’s 10 rows -> 5 folds (2 rows per fold)
make_folds <- function(mouse_ids, k = 5) {
  folds <- vector("list", k)
  uniq_ids <- unique(mouse_ids)
  for (id in uniq_ids) {
    rows <- which(mouse_ids == id)
    if (length(rows) %% k != 0) stop("Cannot split mouse ", id, " rows evenly into ", k, " folds.")
    per_fold <- length(rows) / k # here 10 / 5 = 2
    assign_seq <- sample(rep(seq_len(k), each = per_fold))
    for (f in seq_len(k)) folds[[f]] <- c(folds[[f]], rows[assign_seq == f])
  }
  folds
}

# Lambda grids
# MRCE: lam1.vec, lam2.vec (decreasing, independent)
lam_mrce_vec <- rev(exp(seq(log(1e-2), log(1), length.out = 8)))  # decreasing as required by mrce(cv)

# PN/PG centered 8x8 grid factors (max at index 4)
center_factors <- c(1e-3, 1e-2, 1e-1, 1, 10^-0.5, 10^-1.5, 10^-2.5, 10^-3.5)
stopifnot(length(center_factors) == 8L, which.max(center_factors) == 4L)

# Load PN/PG functions only
load_pn_funs <- function(path) {
  if (!file.exists(path)) stop("missing ", path)
  env <- new.env(parent = baseenv())
  exprs <- parse(path)
  for (e in exprs) {
    if (is.call(e) && length(e) >= 3 && as.character(e[[1]]) %in% c("<-", "=")) {
      rhs <- e[[3]]
      if (is.call(rhs) && identical(rhs[[1]], as.symbol("function"))) {
        eval(e, env)
      }
    }
  }
  env
}
pn_path <- file.path(dirname(script_path), "11.7_need_to_change.R")
pn_env <- load_pn_funs(pn_path)
assign("pn_prox_ctrl", list(max_admm = 1000, tol = 1e-5, eta = 1e-8), envir = pn_env)
assign("prox_ctrl", list(max_admm = 1000, tol = 1e-5, eta = 1e-8), envir = pn_env)
assign("armijo_c", 0.1, envir = pn_env)
pn_sparse_mvreg <- get("pn_sparse_mvreg", envir = pn_env)
pg_sparse_mvreg <- get("pg_sparse_mvreg", envir = pn_env)

save_cv_heatmap <- function(cv_mse, x_vals, y_vals, best_idx, out_path, title, x_lab, y_lab) {
  heat_df <- expand.grid(x_idx = seq_along(x_vals), y_idx = seq_along(y_vals))
  heat_df$mse <- as.vector(cv_mse)
  heat_df$x <- log10(x_vals[heat_df$x_idx])
  heat_df$y <- log10(y_vals[heat_df$y_idx])
  min_pt <- data.frame(
    x = log10(x_vals[best_idx[1]]),
    y = log10(y_vals[best_idx[2]]),
    mse = cv_mse[best_idx[1], best_idx[2]]
  )
  heat_plot <- ggplot(heat_df, aes(x = x, y = y, fill = mse)) +
    geom_tile() +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey90") +
    geom_point(data = min_pt, aes(x = x, y = y), color = "#cc0000", size = 2) +
    geom_text(data = min_pt, aes(x = x, y = y, label = "min"), vjust = -0.8, color = "#cc0000", size = 3) +
    labs(title = title, x = x_lab, y = y_lab, fill = "CV MSE") +
    theme_minimal()
  ggsave(out_path, heat_plot, width = 6, height = 5, dpi = 150)
}

# Configuration
target_reps <- 1L
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1) {
  nproc_out <- tryCatch(system("getconf _NPROCESSORS_ONLN 2>/dev/null", intern = TRUE), error = function(e) "")
  if (!length(nproc_out) || !nzchar(nproc_out[1])) {
    nproc_out <- tryCatch(system("sysctl -n hw.logicalcpu 2>/dev/null", intern = TRUE), error = function(e) "")
  }
  detected_cores <- suppressWarnings(as.integer(nproc_out[1]))
}
if (!is.finite(detected_cores) || detected_cores < 1) detected_cores <- 1L
cv_cores <- as.integer(detected_cores)
cv_fit_max_iter <- 10L
cv_fit_tol <- 1e-3
refit_max_iter <- 30L
refit_tol <- 1e-4
results_df <- data.frame(
  rep = integer(target_reps),
  mse_mrce = numeric(target_reps),
  mse_pn_mse = numeric(target_reps),
  mse_pg_mse = numeric(target_reps),
  time_mrce = numeric(target_reps),
  time_pn_mse = numeric(target_reps),
  time_pg_mse = numeric(target_reps)
)

message(sprintf("Starting %d repetitions...", target_reps))
message(sprintf("PN/PG CV cores: %d", cv_cores))

for (rep_i in 1:target_reps) {
  set.seed(12345 + rep_i)
  
  split <- draw_train_test()
  train_idx <- split$train_idx
  test_idx <- split$test_idx
  train_mouse <- split$train_mouse
  
  X_tr_raw <- X_full[train_idx, , drop = FALSE]
  Y_tr_raw <- Y_full[train_idx, , drop = FALSE]
  X_te_raw <- X_full[test_idx, , drop = FALSE]
  Y_te_raw <- Y_full[test_idx, , drop = FALSE]
  mouse_tr <- train_mouse
  
  sdY <- apply(Y_tr_raw, 2, sd)
  keep_y <- which(is.finite(sdY) & sdY > 0)
  Y_tr_raw <- Y_tr_raw[, keep_y, drop = FALSE]
  Y_te_raw <- Y_te_raw[, keep_y, drop = FALSE]
  
  sdX <- apply(X_tr_raw, 2, sd)
  keep_x <- which(is.finite(sdX) & sdX > 0)
  X_tr_raw <- X_tr_raw[, keep_x, drop = FALSE]
  X_te_raw <- X_te_raw[, keep_x, drop = FALSE]
  
  folds <- make_folds(mouse_tr, k = 5)
  
  # MRCE CV
  mrce_cv_fit <- mrce(
    X = X_tr_raw, Y = Y_tr_raw,
    lam1.vec = lam_mrce_vec, lam2.vec = lam_mrce_vec,
    method = "cv", kfold = 5,folds = folds,
    standardize = FALSE, silent = TRUE
  )
  mrce_cv_err <- mrce_cv_fit$cv.err
  best_idx_mrce <- which(mrce_cv_err == min(mrce_cv_err), arr.ind = TRUE)[1, ]
  lam_best_mrce1 <- lam_mrce_vec[best_idx_mrce[1]]
  lam_best_mrce2 <- lam_mrce_vec[best_idx_mrce[2]]

  mrce_heatmap_path <- file.path(dirname(script_path), sprintf("mrce_cv_heatmap_rep%03d.png", rep_i))
  save_cv_heatmap(
    cv_mse = mrce_cv_err,
    x_vals = lam_mrce_vec,
    y_vals = lam_mrce_vec,
    best_idx = best_idx_mrce,
    out_path = mrce_heatmap_path,
    title = sprintf("MRCE CV Heatmap (Rep %d)", rep_i),
    x_lab = "log10(lambda_1)",
    y_lab = "log10(lambda_2)"
  )
  
  # PN/PG centered 8x8 lambda grids derived from training stats
  n_tr <- nrow(X_tr_raw)
  X_tr_pn_center <- scale(X_tr_raw, center = TRUE, scale = FALSE)
  X_tr_pn_full <- cbind(Intercept = 1, X_tr_pn_center)
  Sxy_tr <- crossprod(X_tr_pn_full, Y_tr_raw) / n_tr
  Syy_tr <- crossprod(Y_tr_raw) / n_tr
  lambda_gamma_0 <- max(abs(Sxy_tr[-1, , drop = FALSE]))
  tmp_off <- Syy_tr - diag(diag(Syy_tr))
  lambda_Omega_0 <- max(abs(tmp_off))
  if (!is.finite(lambda_gamma_0) || lambda_gamma_0 <= 0) lambda_gamma_0 <- 1
  if (!is.finite(lambda_Omega_0) || lambda_Omega_0 <= 0) lambda_Omega_0 <- 1
  lam_gamma_vec <- lambda_gamma_0 * center_factors
  lam_omega_vec <- lambda_Omega_0 * center_factors
  stopifnot(
    length(lam_gamma_vec) == 8L,
    length(lam_omega_vec) == 8L,
    which.max(lam_gamma_vec) == 4L,
    which.max(lam_omega_vec) == 4L
  )

  fold_data <- lapply(folds, function(val_rows) {
    tr_rows <- setdiff(seq_len(nrow(X_tr_raw)), val_rows)
    sc <- scale_with(X_tr_raw, Y_tr_raw, tr_rows)
    X_tr_f <- sc$X_tr
    Y_tr_f <- sc$Y_tr
    X_val_f <- scale(X_tr_raw[val_rows, , drop = FALSE], center = sc$centerX, scale = sc$scaleX)
    Y_val_f <- scale(Y_tr_raw[val_rows, , drop = FALSE], center = sc$centerY, scale = sc$scaleY)
    list(
      X_tr_aug = cbind(Intercept = rep(1, nrow(X_tr_f)), X_tr_f),
      Y_tr_f = Y_tr_f,
      X_val_aug = cbind(Intercept = rep(1, nrow(X_val_f)), X_val_f),
      Y_val_f = Y_val_f
    )
  })

  run_sparse_cv <- function(solver_fun) {
    cv_mse <- matrix(NA_real_, nrow = length(lam_gamma_vec), ncol = length(lam_omega_vec))
    grid_pairs <- expand.grid(gi = seq_along(lam_gamma_vec), oi = seq_along(lam_omega_vec))

    eval_one_pair <- function(pair_i) {
      gi <- grid_pairs$gi[pair_i]
      oi <- grid_pairs$oi[pair_i]
      lam_g <- lam_gamma_vec[gi]
      lam_o <- lam_omega_vec[oi]
      fold_mse <- numeric(length(fold_data))
      for (i in seq_along(fold_data)) {
        fd <- fold_data[[i]]
        fit <- solver_fun(
          Xc = fd$X_tr_aug,
          Y = fd$Y_tr_f,
          lambda_gamma = lam_g,
          lambda_Omega = lam_o,
          gamma_init = NULL,
          Omega_init = diag(ncol(fd$Y_tr_f)),
          max_iter = cv_fit_max_iter,
          tol = cv_fit_tol
        )
        beta <- fit$gamma %*% solve(fit$Omega)
        mu <- as.numeric(colMeans(fd$Y_tr_f) - t(beta) %*% colMeans(fd$X_tr_aug))
        Yhat <- tcrossprod(matrix(1, nrow = nrow(fd$X_val_aug), ncol = 1), mu) + fd$X_val_aug %*% beta
        fold_mse[i] <- mean((fd$Y_val_f - Yhat)^2)
      }
      list(gi = gi, oi = oi, mse = mean(fold_mse))
    }

    if (.Platform$OS.type != "windows" && cv_cores > 1L) {
      cv_res <- parallel::mclapply(
        seq_len(nrow(grid_pairs)),
        eval_one_pair,
        mc.cores = cv_cores,
        mc.preschedule = FALSE
      )
    } else {
      cv_res <- lapply(seq_len(nrow(grid_pairs)), eval_one_pair)
    }

    for (res in cv_res) {
      cv_mse[res$gi, res$oi] <- res$mse
    }
    best_idx <- which(cv_mse == min(cv_mse, na.rm = TRUE), arr.ind = TRUE)[1, ]
    list(
      cv_mse = cv_mse,
      best_idx = best_idx,
      lambda_gamma = lam_gamma_vec[best_idx[1]],
      lambda_omega = lam_omega_vec[best_idx[2]]
    )
  }

  pn_cv <- run_sparse_cv(pn_sparse_mvreg)
  pg_cv <- run_sparse_cv(pg_sparse_mvreg)

  pn_cv_mse <- pn_cv$cv_mse
  pg_cv_mse <- pg_cv$cv_mse
  best_idx_pn <- pn_cv$best_idx
  best_idx_pg <- pg_cv$best_idx
  lam_best_pn_mse <- pn_cv$lambda_gamma
  lam_best_pn_mse_Omega <- pn_cv$lambda_omega
  lam_best_pg_mse <- pg_cv$lambda_gamma
  lam_best_pg_mse_Omega <- pg_cv$lambda_omega

  stopifnot(
    all(dim(pn_cv_mse) == c(8L, 8L)),
    all(dim(pg_cv_mse) == c(8L, 8L))
  )

  pn_heatmap_path <- file.path(dirname(script_path), sprintf("pn_cv_heatmap_rep%03d.png", rep_i))
  save_cv_heatmap(
    cv_mse = pn_cv_mse,
    x_vals = lam_gamma_vec,
    y_vals = lam_omega_vec,
    best_idx = best_idx_pn,
    out_path = pn_heatmap_path,
    title = sprintf("PN CV Heatmap (Rep %d)", rep_i),
    x_lab = "log10(lambda_gamma)",
    y_lab = "log10(lambda_omega)"
  )

  pg_heatmap_path <- file.path(dirname(script_path), sprintf("pg_cv_heatmap_rep%03d.png", rep_i))
  save_cv_heatmap(
    cv_mse = pg_cv_mse,
    x_vals = lam_gamma_vec,
    y_vals = lam_omega_vec,
    best_idx = best_idx_pg,
    out_path = pg_heatmap_path,
    title = sprintf("PG CV Heatmap (Rep %d)", rep_i),
    x_lab = "log10(lambda_gamma)",
    y_lab = "log10(lambda_omega)"
  )
  
  # Refit & Test
  sc_full <- scale_with(X_tr_raw, Y_tr_raw, rows = seq_len(nrow(X_tr_raw)))
  X_tr <- sc_full$X_tr
  Y_tr <- sc_full$Y_tr
  X_te <- scale(X_te_raw, center = sc_full$centerX, scale = sc_full$scaleX)
  Y_te <- scale(Y_te_raw, center = sc_full$centerY, scale = sc_full$scaleY)
  
  # MRCE Refit
  t0_mrce <- Sys.time()
  mrce_fit <- mrce(
    X = X_tr, Y = Y_tr, lam1 = lam_best_mrce1, lam2 = lam_best_mrce2,
    method = "single", standardize = FALSE, silent = TRUE
  )
  time_mrce <- as.numeric(difftime(Sys.time(), t0_mrce, units = "secs"))
  
  beta_mrce <- mrce_fit$Bhat
  mu_mrce <- as.numeric(mrce_fit$muhat)
  Yhat_mrce <- tcrossprod(matrix(1, nrow = nrow(X_te), ncol = 1), mu_mrce) + X_te %*% beta_mrce
  mse_mrce_test <- mean((Y_te - Yhat_mrce)^2)
  
  # PN Refit
  X_tr_pn <- cbind(Intercept = rep(1, nrow(X_tr)), X_tr)
  X_te_pn <- cbind(Intercept = rep(1, nrow(X_te)), X_te)
  
  t0_pn_mse <- Sys.time()
  pn_fit_mse <- pn_sparse_mvreg(
    Xc = X_tr_pn,
    Y = Y_tr,
    lambda_gamma = lam_best_pn_mse,
    lambda_Omega = lam_best_pn_mse_Omega,
    gamma_init = NULL,
    Omega_init = diag(ncol(Y_tr)),
    max_iter = refit_max_iter,
    tol = refit_tol
  )
  time_pn_mse <- as.numeric(difftime(Sys.time(), t0_pn_mse, units = "secs"))
  
  beta_pn_mse <- pn_fit_mse$gamma %*% solve(pn_fit_mse$Omega)
  mu_pn_mse <- as.numeric(colMeans(Y_tr) - t(beta_pn_mse) %*% colMeans(X_tr_pn))
  Yhat_pn_mse <- tcrossprod(matrix(1, nrow = nrow(X_te_pn), ncol = 1), mu_pn_mse) + X_te_pn %*% beta_pn_mse
  mse_pn_mse_test <- mean((Y_te - Yhat_pn_mse)^2)

  # PG Refit
  t0_pg_mse <- Sys.time()
  pg_fit_mse <- pg_sparse_mvreg(
    Xc = X_tr_pn,
    Y = Y_tr,
    lambda_gamma = lam_best_pg_mse,
    lambda_Omega = lam_best_pg_mse_Omega,
    gamma_init = NULL,
    Omega_init = diag(ncol(Y_tr)),
    max_iter = refit_max_iter,
    tol = refit_tol
  )
  time_pg_mse <- as.numeric(difftime(Sys.time(), t0_pg_mse, units = "secs"))

  beta_pg_mse <- pg_fit_mse$gamma %*% solve(pg_fit_mse$Omega)
  mu_pg_mse <- as.numeric(colMeans(Y_tr) - t(beta_pg_mse) %*% colMeans(X_tr_pn))
  Yhat_pg_mse <- tcrossprod(matrix(1, nrow = nrow(X_te_pn), ncol = 1), mu_pg_mse) + X_te_pn %*% beta_pg_mse
  mse_pg_mse_test <- mean((Y_te - Yhat_pg_mse)^2)

  results_df[rep_i, ] <- c(
    rep_i, mse_mrce_test, mse_pn_mse_test, mse_pg_mse_test,
    time_mrce, time_pn_mse, time_pg_mse
  )
}

t_end <- Sys.time()
total_time <- as.numeric(difftime(t_end, t_start, units = "secs"))

mse_row <- results_df[1, c("mse_mrce", "mse_pn_mse", "mse_pg_mse")]
winner <- names(mse_row)[which.min(as.numeric(mse_row))]

cat("\n=== Final Results (1 Repetition) ===\n")
cat(sprintf("Total Runtime: %.1f seconds\n\n", total_time))

cat("--- Test MSE ---\n")
cat(sprintf("MRCE: %.4f\n", results_df$mse_mrce[1]))
cat(sprintf("PN  : %.4f\n", results_df$mse_pn_mse[1]))
cat(sprintf("PG  : %.4f\n", results_df$mse_pg_mse[1]))
cat(sprintf("Winner by Test MSE: %s\n", winner))

cat("\n--- Fit Time (secs) ---\n")
cat(sprintf("MRCE: %.4f\n", results_df$time_mrce[1]))
cat(sprintf("PN  : %.4f\n", results_df$time_pn_mse[1]))
cat(sprintf("PG  : %.4f\n", results_df$time_pg_mse[1]))

saveRDS(results_df, file.path(dirname(script_path), "compare_results_modified.rds"))
