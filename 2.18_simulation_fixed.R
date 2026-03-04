## =========================================================
##  Shared utilities
## =========================================================

symmetrize <- function(A) 0.5 * (A + t(A))

clip_eig <- function(A, eps = 0) {
  ea <- eigen(symmetrize(A), symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  V %*% diag(d, length(d)) %*% t(V)
}

proj_psd <- function(A) clip_eig(A, eps = 0)

is_pd_by_chol_rank <- function(A, eps = 1e-8) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  R <- tryCatch(
    suppressWarnings(chol(A_sym + diag(eps, n), pivot = TRUE)),
    error = function(e) NULL
  )
  if (is.null(R)) return(FALSE)
  attr(R, "rank") == n
}

soft_thresh <- function(A, tau) sign(A) * pmax(abs(A) - tau, 0)

soft_thresh_offdiag <- function(A, tau) {
  B <- symmetrize(A)
  idx <- matrix(TRUE, nrow(A), ncol(A)); diag(idx) <- FALSE
  B[idx] <- soft_thresh(B[idx], tau)
  B
}

safe_inv_psd <- function(A, eps = 1e-8, jitter_init = 0, jitter_max = 1e-5) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  jitter <- jitter_init
  while (jitter <= jitter_max) {
    chol_try <- tryCatch(
      suppressWarnings(chol(A_sym + diag(jitter, n))),
      error = function(e) NULL
    )
    if (!is.null(chol_try)) {
      invA <- chol2inv(chol_try)
      return(symmetrize(invA))
    }
    jitter <- if (jitter == 0) eps else jitter * 10
  }
  ea <- eigen(A_sym, symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  symmetrize(V %*% (t(V) / d))
}

safe_logdet_spd <- function(A, jitter_init = 0, jitter_max = 1e-5, eps = 1e-10) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  jitter <- jitter_init
  while (jitter <= jitter_max) {
    R <- tryCatch(
      suppressWarnings(chol(A_sym + diag(jitter, n))),
      error = function(e) NULL
    )
    if (!is.null(R)) {
      return(2 * sum(log(diag(R))))
    }
    jitter <- if (jitter == 0) eps else jitter * 10
  }
  ea <- eigen(A_sym, symmetric = TRUE, only.values = TRUE)$values
  sum(log(pmax(ea, eps)))
}

eig_max_symmetric <- function(A) {
  eigvals <- suppressWarnings(eigen(symmetrize(A), symmetric = TRUE, only.values = TRUE)$values)
  max(eigvals)
}

matrix_condition_stats <- function(A, eps = 1e-12) {
  eigvals <- suppressWarnings(eigen(symmetrize(A), symmetric = TRUE, only.values = TRUE)$values)
  if (length(eigvals) == 0 || any(!is.finite(eigvals))) {
    return(list(min_eig = NA_real_, max_eig = NA_real_, kappa = NA_real_))
  }
  min_eig <- min(eigvals)
  max_eig <- max(eigvals)
  eig_clipped <- pmax(eigvals, eps)
  kappa <- max(eig_clipped) / min(eig_clipped)
  list(min_eig = min_eig, max_eig = max_eig, kappa = kappa)
}

matrix_sparsity_stats <- function(gamma, Omega, thresh = 1e-10) {
  nnz_gamma_total <- sum(abs(gamma) > thresh)
  nnz_gamma_pen <- if (nrow(gamma) >= 2) {
    sum(abs(gamma[-1, , drop = FALSE]) > thresh)
  } else {
    nnz_gamma_total
  }
  gamma_total_pen <- if (nrow(gamma) >= 2) (nrow(gamma) - 1) * ncol(gamma) else nrow(gamma) * ncol(gamma)
  off_upper <- upper.tri(Omega, diag = FALSE)
  nnz_omega_offdiag_upper <- sum(abs(Omega[off_upper]) > thresh)
  omega_offdiag_total_upper <- sum(off_upper)
  list(
    nnz_gamma_total = nnz_gamma_total,
    nnz_gamma_pen = nnz_gamma_pen,
    gamma_pen_sparsity = if (gamma_total_pen > 0) nnz_gamma_pen / gamma_total_pen else NA_real_,
    nnz_omega_offdiag_upper = nnz_omega_offdiag_upper,
    omega_offdiag_sparsity = if (omega_offdiag_total_upper > 0) nnz_omega_offdiag_upper / omega_offdiag_total_upper else NA_real_
  )
}

bind_rows_safe <- function(df_list) {
  df_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, df_list)
  if (length(df_list) == 0) return(data.frame())
  row.names(df_list) <- NULL
  do.call(rbind, df_list)
}

## =========================================================
##  Shared simulation setup
## =========================================================

shared_seed <- 20251025
set.seed(shared_seed)
shared_params <- list(
  n = 1000,
  p = 50,
  q = 50,
  ratio_lambda = 1e-4,
  gamma_sparsity = 0.05,
  gamma_mag = c(1, 2),
  omega_kappa = 60,
  omega_sparsity = 0.05
)

generate_sparse_gamma <- function(p_aug, q,
                                  sparsity = 0.05,
                                  nnz_per_col = NULL,
                                  mag_range = c(1, 2),
                                  intercept_zero = TRUE,
                                  signed = TRUE) {
  if (length(mag_range) != 2) stop("mag_range must have length 2")
  lo <- min(mag_range); hi <- max(mag_range)
  gamma <- matrix(0, p_aug, q)
  idx_pool <- if (intercept_zero && p_aug > 1) 2:p_aug else seq_len(p_aug)
  if (length(idx_pool) == 0) return(gamma)
  if (is.null(nnz_per_col)) {
    nnz_per_col <- max(1, round(sparsity * length(idx_pool)))
  } else {
    nnz_per_col <- max(1, min(nnz_per_col, length(idx_pool)))
  }
  for (j in seq_len(q)) {
    nnz <- min(nnz_per_col, length(idx_pool))
    active <- sample(idx_pool, nnz)
    vals <- runif(nnz, lo, hi)
    if (signed) vals <- vals * sample(c(-1, 1), nnz, replace = TRUE)
    gamma[active, j] <- vals
  }
  if (intercept_zero && p_aug >= 1) gamma[1, ] <- 0
  gamma
}

generate_sparse_omega <- function(q,
                                  kappa_target = 50,
                                  offdiag_prob = 0.05,
                                  base_val = 0.5,
                                  ensure_pd = TRUE) {
  if (q == 1) return(matrix(1, 1, 1))
  B <- matrix(0, q, q)
  upper_mask <- matrix(runif(q^2) < offdiag_prob, q, q)
  upper_mask[lower.tri(upper_mask, diag = TRUE)] <- FALSE
  B[upper_mask] <- base_val
  B <- B + t(B)
  eig_B <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(eig_B)
  lambda_max <- max(eig_B)
  if (abs(lambda_min - lambda_max) < 1e-12) {
    delta <- 1
  } else {
    delta <- (lambda_max - kappa_target * lambda_min) / (kappa_target - 1)
  }
  if (!is.finite(delta) || delta <= -lambda_min + 1e-8) {
    delta <- -lambda_min + 1e-4
  }
  Omega <- B + delta * diag(q)
  if (ensure_pd) {
    eig_vals <- suppressWarnings(eigen(symmetrize(Omega), symmetric = TRUE, only.values = TRUE)$values)
    eig_vals <- pmax(eig_vals, 1e-8)
    attr(Omega, "cond_est") <- max(eig_vals) / min(eig_vals)
  }
  Omega
}

p_aug <- shared_params$p + 1
gamma_true <- generate_sparse_gamma(
  p_aug,
  shared_params$q,
  sparsity = shared_params$gamma_sparsity,
  mag_range = shared_params$gamma_mag
)
shared_Beta_true <- gamma_true[-1, , drop = FALSE]
shared_Omega_true <- generate_sparse_omega(
  shared_params$q,
  kappa_target = shared_params$omega_kappa,
  offdiag_prob = shared_params$omega_sparsity
)
shared_Sigma_true <- safe_inv_psd(shared_Omega_true)
shared_truth <- list(
  gamma = gamma_true,
  Omega = shared_Omega_true,
  cond_number = attr(shared_Omega_true, "cond_est")
)

shared_X <- matrix(rnorm(shared_params$n * shared_params$p), shared_params$n, shared_params$p)
suppressPackageStartupMessages(library(MASS))
shared_Y <- shared_X %*% shared_Beta_true + MASS::mvrnorm(
  shared_params$n,
  mu = rep(0, shared_params$q),
  Sigma = shared_Sigma_true
)
shared_Xc <- cbind(1, scale(shared_X, center = TRUE, scale = FALSE))
shared_p_aug <- ncol(shared_Xc)

Sxy_shared <- crossprod(shared_Xc, shared_Y) / shared_params$n
Syy_shared <- crossprod(shared_Y) / shared_params$n
lambda_gamma_0 <- 2 * max(abs(Sxy_shared[-1, , drop = FALSE]))
tmp_off_shared <- Syy_shared - diag(diag(Syy_shared))
lambda_Omega_0 <- 2 * max(abs(tmp_off_shared))
path_len <- 60
lambda_gamma_seq <- exp(seq(log(lambda_gamma_0),
                            log(lambda_gamma_0 * shared_params$ratio_lambda),
                            length.out = path_len))
lambda_Omega_seq <- exp(seq(log(lambda_Omega_0),
                            log(lambda_Omega_0 * shared_params$ratio_lambda),
                            length.out = path_len))
lambda_path <- data.frame(
  step = seq_len(path_len),
  lambda_gamma = lambda_gamma_seq,
  lambda_Omega = lambda_Omega_seq
)

compute_null_corner_init <- function(Sxx, Sxy, Syy, unpenalized = 1L, eps = 1e-8) {
  p_aug <- nrow(Sxx); q <- ncol(Sxy)
  gamma_init <- matrix(0, p_aug, q)
  U <- unpenalized
  Sxx_UU <- Sxx[U, U, drop = FALSE]
  Sxy_U <- Sxy[U, , drop = FALSE]
  Sxx_UU_inv <- solve(Sxx_UU)
  B <- t(Sxy_U) %*% Sxx_UU_inv %*% Sxy_U
  diag_term <- diag(Syy - B)
  diag_term <- pmax(diag_term, eps)
  Omega_diag <- 1 / diag_term
  Omega_init <- diag(Omega_diag, nrow = length(Omega_diag))
  gamma_U <- Sxx_UU_inv %*% Sxy_U %*% Omega_init
  gamma_init[U, ] <- gamma_U
  list(gamma = gamma_init, Omega = Omega_init)
}

tol_env <- suppressWarnings(as.numeric(Sys.getenv("SOLVER_TOL", unset = NA_character_)))
if (is.na(tol_env)) tol_env <- 1e-5
maxit_env <- suppressWarnings(as.integer(Sys.getenv("SOLVER_MAXIT", unset = NA_character_)))
if (is.na(maxit_env) || maxit_env <= 0) maxit_env <- 1000

ctrl <- list(
  lambda = list(gamma = lambda_gamma_seq[1], Omega = lambda_Omega_seq[1]),
  maxit = maxit_env,
  tol = tol_env,
  tol_obj = tol_env,
  tol_gm = tol_env,
  use_relative_obj = TRUE
)

prox_tol_env <- suppressWarnings(as.numeric(Sys.getenv("PROX_TOL", unset = NA_character_)))
if (is.na(prox_tol_env)) prox_tol_env <- 0.1 * tol_env
prox_max_admm_env <- suppressWarnings(as.integer(Sys.getenv("PROX_MAX_ADMM", unset = NA_character_)))
if (is.na(prox_max_admm_env) || prox_max_admm_env <= 0) prox_max_admm_env <- 10000
prox_eta_env <- suppressWarnings(as.numeric(Sys.getenv("PROX_ETA", unset = NA_character_)))
if (is.na(prox_eta_env) || prox_eta_env < 0) prox_eta_env <- 1e-6
prox_ctrl <- list(max_admm = prox_max_admm_env, tol = prox_tol_env, eta = prox_eta_env)

pn_eta_env <- suppressWarnings(as.numeric(Sys.getenv("PN_PROX_ETA", unset = NA_character_)))
if (is.na(pn_eta_env) || pn_eta_env < 0) pn_eta_env <- prox_eta_env
pn_prox_ctrl <- list(max_admm = prox_max_admm_env, tol = prox_tol_env, eta = pn_eta_env)

# PN outer monotone safeguard controls.
pn_ls_max_halving_env <- suppressWarnings(as.integer(Sys.getenv("PN_LS_MAX_HALVING", unset = NA_character_)))
if (is.na(pn_ls_max_halving_env) || pn_ls_max_halving_env < 0) pn_ls_max_halving_env <- 3L
pn_ls_min_alpha_env <- suppressWarnings(as.numeric(Sys.getenv("PN_LS_MIN_ALPHA", unset = NA_character_)))
if (is.na(pn_ls_min_alpha_env) || pn_ls_min_alpha_env <= 0) pn_ls_min_alpha_env <- 1e-8
pn_ls_ctrl <- list(max_halving = pn_ls_max_halving_env, min_alpha = pn_ls_min_alpha_env)

# PN inner ADMM: adaptive rho enabled.
pn_admm_ctrl <- list(
  adaptive = TRUE,
  rho = 1,
  rho_mu = 10,
  rho_tau_inc = 2,
  rho_tau_dec = 2,
  rho_min = 1e-8,
  rho_max = 1e8
)

armijo_c <- 0.01
pn_inner_admm_tol_fixed <- 1e-7
pn_outer_progress_every_env <- suppressWarnings(as.integer(Sys.getenv("PN_OUTER_PROGRESS_EVERY", unset = NA_character_)))
if (is.na(pn_outer_progress_every_env) || pn_outer_progress_every_env <= 0) pn_outer_progress_every_env <- 50L

null_corner_init <- compute_null_corner_init(
  crossprod(shared_Xc) / shared_params$n,
  Sxy_shared,
  Syy_shared,
  unpenalized = 1L
)

shared_data <- list(
  X = shared_Xc,
  Y = shared_Y,
  lambda = ctrl$lambda,
  lambda_path = lambda_path,
  gamma_init = null_corner_init$gamma,
  Omega_init = null_corner_init$Omega,
  params = shared_params,
  truth = shared_truth
)

## =========================================================
##  Objective and derivatives
## =========================================================

.ensure_Oinv <- function(Omega, Oinv = NULL) {
  if (is.null(Oinv)) safe_inv_psd(Omega) else symmetrize(Oinv)
}
.ensure_beta <- function(gamma, Omega, Oinv = NULL, beta = NULL) {
  if (!is.null(beta)) return(beta)
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  gamma %*% Oinv_eff
}

smooth_obj <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL, beta = NULL) {
  logdet <- safe_logdet_spd(Omega)
  term1 <- -logdet
  term3 <- -2 * sum(Sxy * gamma)
  term4 <- sum(Omega * Syy)
  
  if (!is.null(beta)) {
    BSB  <- t(beta) %*% Sxx %*% beta
    term2 <- sum(BSB * Omega)
  } else {
    Oinv_eff <- .ensure_Oinv(Omega, Oinv)
    M <- t(gamma) %*% Sxx %*% gamma
    term2 <- sum(Oinv_eff * M)
  }
  term1 + term2 + term3 + term4
}

loss_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n,
                    lambda_gamma, lambda_Omega, lambda_Omega_diag = 0,
                    Oinv = NULL, beta = NULL) {
  smooth_obj(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv, beta = beta) +
    penalty_obj(gamma, Omega, lambda_gamma, lambda_Omega, lambda_Omega_diag = lambda_Omega_diag)
}

penalty_obj <- function(gamma, Omega, lambda_gamma, lambda_Omega, lambda_Omega_diag = 0) {
  pen_g <- lambda_gamma * sum(abs(gamma[-1, , drop = FALSE]))
  pen_o <- lambda_Omega * (sum(abs(Omega)) - sum(abs(diag(Omega)))) +
    lambda_Omega_diag * sum(abs(diag(Omega)))
  pen_g + pen_o
}

nll_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n, Oinv = NULL, beta = NULL) {
  smooth_obj(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv, beta = beta)
}

grad_gamma <- function(gamma, Omega, Sxx, Sxy, Oinv = NULL, beta = NULL) {
  if (is.null(beta)) {
    Oinv_eff <- .ensure_Oinv(Omega, Oinv)
    2 * (Sxx %*% gamma %*% Oinv_eff - Sxy)
  } else {
    2 * (Sxx %*% beta - Sxy)
  }
}

grad_Omega <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  if (is.null(beta)) {
    M <- t(gamma) %*% Sxx %*% gamma
    Syy - Oinv_eff %*% M %*% Oinv_eff - Oinv_eff
  } else {
    BSB <- t(beta) %*% Sxx %*% beta
    Syy - BSB - Oinv_eff
  }
}

H_blocks <- function(gamma, Omega, Sxx, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  beta_eff <- .ensure_beta(gamma, Omega, Oinv_eff, beta)
  
  Hgg <- 2 * Sxx
  Hgo <- -2 * Sxx %*% beta_eff
  Hog <- t(Hgo)
  Hoo <- Oinv_eff + 2 * t(beta_eff) %*% Sxx %*% beta_eff
  list(Hgg = Hgg, Hgo = Hgo, Hog = Hog, Hoo = Hoo)
}

H_matrix <- function(Hb) {
  rbind(cbind(Hb$Hgg, Hb$Hgo),
        cbind(Hb$Hog, Hb$Hoo))
}

pack_xi <- function(gamma, Omega) rbind(gamma, Omega)
unpack_xi <- function(Xi, p) list(
  gamma = Xi[1:p, , drop = FALSE],
  Omega = symmetrize(Xi[(p + 1):nrow(Xi), , drop = FALSE])
)

local_norm <- function(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv_eff, beta = beta)
  H  <- H_matrix(Hb)
  Xi_d <- pack_xi(d_gamma, d_Omega)
  val <- sum(Xi_d * (H %*% Xi_d %*% Oinv_eff))
  sqrt(max(val, 0))
}

## =========================================================
##  Prox operator for Omega
## =========================================================

prox_psd_offdiag_l1 <- function(V,
                                tau,
                                mu       = 1,
                                rho      = 1,
                                max_admm = 1000,
                                tol      = 1e-5,
                                eig_floor = 0,
                                Omega_init = NULL,
                                A_init = NULL) {
  
  Z <- symmetrize(V)
  offdiag_l1 <- function(M) {
    idx <- matrix(TRUE, nrow(M), ncol(M))
    diag(idx) <- FALSE
    sum(abs(M[idx]))
  }
  dual_feasible_Y <- function(A, mu, tau) {
    Y <- -(1 / mu) * A
    diag(Y) <- 0
    idx <- matrix(TRUE, nrow(Y), ncol(Y))
    diag(idx) <- FALSE
    Y[idx] <- pmin(pmax(Y[idx], -tau), tau)
    symmetrize(Y)
  }
  primal_value <- function(Omega) {
    0.5 * sum((Omega - Z)^2) + tau * offdiag_l1(Omega)
  }
  dual_value <- function(Y) {
    K_star <- clip_eig(Z - Y, eps = eig_floor)
    sum(Y * K_star) + 0.5 * sum((K_star - Z)^2)
  }

  omega_hat <- symmetrize(soft_thresh_offdiag(Z, tau / rho))
  q <- nrow(omega_hat)
  if (is_pd_by_chol_rank(omega_hat - eig_floor * diag(q))) {
    attr(omega_hat, "state") <- list(Omega = omega_hat, A = matrix(0, q, q))
    return(omega_hat)
  }

  omega_hat <- clip_eig(omega_hat, eps = eig_floor)
  attr(omega_hat, "state") <- list(Omega = omega_hat, A = matrix(0, q, q))

  Omega <- if (is.null(Omega_init)) omega_hat else symmetrize(Omega_init)
  A     <- if (is.null(A_init)) matrix(0, nrow(Omega), ncol(Omega)) else symmetrize(A_init)
  
  for (l in seq_len(max_admm)) {
    K <- clip_eig(Omega + mu * A, eps = eig_floor)
    
    T <- (K + mu * (Z - A)) / (1 + mu)
    Omega_new <- symmetrize(
      soft_thresh_offdiag(T, (mu * tau) / ((1 + mu) * rho))
    )
    
    A <- A - (K - Omega_new) / mu
    
    p_val <- primal_value(Omega_new)
    Y_feas <- dual_feasible_Y(A, mu, tau)
    d_val <- dual_value(Y_feas)
    gap <- max(p_val - d_val, 0)
    if (is.finite(gap) && is.finite(p_val)) {
      if (gap / (1 + abs(p_val)) <= tol) {
        K_out <- clip_eig(Omega_new, eps = eig_floor)
        attr(K_out, "state") <- list(Omega = K_out, A = A)
        return(K_out)
      }
    }
    Omega <- Omega_new
  }
  
  K_out <- clip_eig(symmetrize(Omega), eps = eig_floor)
  attr(K_out, "state") <- list(Omega = K_out, A = A)
  K_out
}

## =========================================================
##  Stopping helpers
## =========================================================

rel_obj_change <- function(prev, curr) {
  abs(curr - prev) / (1 + abs(prev))
}

oracle_gap_ok <- function(curr, oracle, gap_tol, use_relative_gap = FALSE) {
  if (is.null(oracle) || is.null(gap_tol) || !is.finite(oracle) || !is.finite(curr)) return(FALSE)
  gap <- curr - oracle
  if (!use_relative_gap) return(gap <= gap_tol)
  gap / max(1, abs(oracle)) <= gap_tol
}

## =========================================================
##  Prox-Gradient with improved stopping
## =========================================================

pg_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 200, tol = 1e-4, L0 = 1,
                            track_loss = TRUE,
                            ctrl = NULL, warm_state = NULL,
                            oracle_value = NULL, gap_tol = NULL, use_relative_gap = FALSE) {
  
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit
    if (!is.null(ctrl$tol)) tol <- ctrl$tol
  }
  
  tol_obj <- if (!is.null(ctrl$tol_obj)) ctrl$tol_obj else tol
  tol_gm  <- if (!is.null(ctrl$tol_gm))  ctrl$tol_gm  else tol
  use_rel_obj <- isTRUE(ctrl$use_relative_obj)
  
  n <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q)          else symmetrize(Omega_init)
  Oinv  <- safe_inv_psd(Omega)
  
  t0 <- Sys.time()
  loss_curr <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n,
                       lambda_gamma, lambda_Omega, Oinv = Oinv)
  
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    loss_hist[1] <- loss_curr
    time_hist[1] <- 0
  }
  
  L_prev <- L0
  converged <- FALSE
  stop_reason <- NA_character_
  backtrack_stalls <- 0L
  pg_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  
  for (m in seq_len(max_iter)) {
    Lm <- L_prev
    loss_prev <- loss_curr
    accepted_step <- FALSE
    bt_iter <- 0L
    
    repeat {
      bt_iter <- bt_iter + 1L
      if (bt_iter > 100L) {
        backtrack_stalls <- backtrack_stalls + 1L
        L_prev <- Lm
        break
      }
      Gg_m <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv)
      Go_m <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv)
      
      Vg_m <- gamma - Gg_m / Lm
      Vo_m <- symmetrize(Omega - Go_m / Lm)
      
      s_gamma_m <- soft_thresh(Vg_m, lambda_gamma / Lm)
      s_gamma_m[1, ] <- Vg_m[1, ]
      
      s_Omega_m <- prox_psd_offdiag_l1(
        Vo_m,
        tau = lambda_Omega / Lm,
        Omega_init = if (is.null(pg_psd_state)) Omega else pg_psd_state$Omega,
        A_init = if (is.null(pg_psd_state)) NULL else pg_psd_state$A,
        max_admm = prox_ctrl$max_admm,
        tol = prox_ctrl$tol,
        eig_floor = prox_ctrl$eta
      )
      pg_psd_state <- attr(s_Omega_m, "state")
      
      d_gamma_m <- s_gamma_m - gamma
      d_Omega_m <- s_Omega_m - Omega
      
      lambda_m <- local_norm(d_gamma_m, d_Omega_m, gamma, Omega, Sxx, Oinv = Oinv)
      beta_m   <- sqrt(Lm) * sqrt(sum(d_gamma_m^2) + sum(d_Omega_m^2))
      
      if (lambda_m^2 / max(beta_m^2, .Machine$double.eps) + lambda_m > 1) {
        Lm <- 2 * Lm
        next
      }
      
      step_norm <- sqrt(sum(d_gamma_m^2) + sum(d_Omega_m^2))
      
      alpha_m <- (beta_m^2) / (lambda_m * (lambda_m + beta_m^2))
      gamma_old <- gamma
      Omega_old <- Omega
      
      gamma_candidate <- gamma_old + alpha_m * d_gamma_m
      Omega_candidate <- symmetrize(Omega_old + alpha_m * d_Omega_m)
      
      Oinv_try <- tryCatch(safe_inv_psd(Omega_candidate), error = function(e) NULL)
      if (is.null(Oinv_try)) {
        backtrack_stalls <- backtrack_stalls + 1L
        L_prev <- Lm
        break
      }
      
      loss_try <- loss_fn(gamma_candidate, Omega_candidate, Sxx, Sxy, Syy, n,
                          lambda_gamma, lambda_Omega, Oinv = Oinv_try)
      if (!is.finite(loss_try)) {
        backtrack_stalls <- backtrack_stalls + 1L
        L_prev <- Lm
        break
      }
      
      gamma <- gamma_candidate
      Omega <- Omega_candidate
      Oinv  <- Oinv_try
      loss_curr <- loss_try
      L_prev <- Lm
      accepted_step <- TRUE
      
      ## Gradient mapping norm uses proximal point at (gamma_old, Omega_old, Lm)
      gm_gamma <- Lm * (gamma_old - s_gamma_m)
      gm_Omega <- Lm * (Omega_old - s_Omega_m)
      gm_norm  <- sqrt(sum(gm_gamma^2) + sum(gm_Omega^2))
      
      ## Oracle gap stopping if requested
      if (oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE
        stop_reason <- "oracle gap reached"
        break
      }
      
      ## Objective stabilization stopping (used only when no oracle gap is supplied)
      if (is.null(gap_tol)) {
        if (use_rel_obj) {
          if (rel_obj_change(loss_prev, loss_curr) <= tol_obj) {
            converged <- TRUE
            stop_reason <- "objective stabilized"
            break
          }
        } else {
          if (abs(loss_curr - loss_prev) <= tol_obj) {
            converged <- TRUE
            stop_reason <- "objective stabilized"
            break
          }
        }
      }
      
      break
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- loss_curr
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    
    if (converged) break
    if (!accepted_step) next
  }
  
  elapsed_time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- list(gamma = gamma, Omega = Omega, iters = m,
              converged = converged,
              stop_reason = ifelse(is.na(stop_reason), "max_iter reached", stop_reason),
              elapsed_time_sec = elapsed_time_sec,
              Sxx = Sxx, Sxy = Sxy, Syy = Syy, n = n,
              backtrack_stalls = backtrack_stalls,
              psd_state = pg_psd_state)
  if (track_loss) {
    out$loss_hist <- loss_hist[1:(m + 1)]
    out$time_hist <- time_hist[1:(m + 1)]
  }
  out
}

## =========================================================
##  Prox-Newton subproblem ADMM
## =========================================================

pn_subproblem_admm <- function(gamma, Omega, Sxx, Sxy, Syy,
                               lambda_gamma, lambda_Omega,
                               Oinv = NULL, beta = NULL,
                               rho = 1, mu = 1,
                               max_admm = 10000, admm_tol = 1e-4,
                               eig_floor = 1e-12,
                               adaptive_rho = TRUE,
                               rho_mu = 10,
                               rho_tau_inc = 2,
                               rho_tau_dec = 2,
                               rho_min = 1e-8,
                               rho_max = 1e8,
                               state = NULL,
                               n = NULL,
                               track_trace = FALSE) {
  p <- nrow(gamma); q <- ncol(gamma)
  
  beta_eff <- .ensure_beta(gamma, Omega, Oinv, beta)
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  
  Gg <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv_eff, beta = beta_eff)
  Go <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv_eff, beta = beta_eff)
  P  <- pack_xi(Gg, Go)
  
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv_eff, beta = beta_eff)
  H  <- symmetrize(H_matrix(Hb))
  
  Xi <- pack_xi(gamma, Omega)
  Gm <- P - H %*% Xi %*% Oinv_eff
  
  eigH <- eigen(symmetrize(H), symmetric = TRUE)
  UH   <- eigH$vectors
  lamH <- pmax(eigH$values, eig_floor)
  
  eigS <- eigen(symmetrize(Oinv_eff), symmetric = TRUE)
  US   <- eigS$vectors
  lamS <- pmax(eigS$values, eig_floor)
  
  Xi_var <- Xi
  Z_var  <- Xi
  Gam    <- matrix(0, p + q, q)
  rho_curr <- max(min(if (is.finite(rho) && rho > 0) rho else 1, rho_max), rho_min)
  if (!is.null(state)) {
    if (!is.null(state$Xi)) Xi_var <- state$Xi
    if (!is.null(state$Z))  Z_var  <- state$Z
    if (!is.null(state$Gam)) Gam   <- state$Gam
    if (!is.null(state$rho) && is.finite(state$rho) && state$rho > 0) {
      rho_curr <- max(min(state$rho, rho_max), rho_min)
    }
  }
  prox_state <- if (!is.null(state)) state$psd_state else NULL
  
  Den <- outer(lamH, lamS, "*") + rho_curr
  admm_converged <- FALSE
  k_last <- 0L
  inner_trace_list <- if (track_trace) vector("list", max_admm) else NULL
  
  for (k in 1:max_admm) {
    k_last <- k
    rho_iter <- rho_curr
    Z_prev <- Z_var
    R  <- rho_iter * Z_var - Gam - Gm
    Rt <- crossprod(UH, R) %*% US
    Y  <- Rt / Den
    Xi_var <- UH %*% Y %*% t(US)
    
    Xi_parts <- unpack_xi(Xi_var + Gam / rho_iter, p)
    Zg <- soft_thresh(Xi_parts$gamma, lambda_gamma / rho_iter)
    Zg[1, ] <- Xi_parts$gamma[1, ]
    
    Zo <- prox_psd_offdiag_l1(
      Xi_parts$Omega,
      tau = lambda_Omega / rho_iter,
      mu = mu,
      Omega_init = if (is.null(prox_state)) Omega else prox_state$Omega,
      A_init = if (is.null(prox_state)) NULL else prox_state$A,
      max_admm = pn_prox_ctrl$max_admm,
      tol = pn_prox_ctrl$tol,
      eig_floor = pn_prox_ctrl$eta
    )
    prox_state <- attr(Zo, "state")
    Z_new <- pack_xi(Zg, Zo)
    
    Gam <- Gam + rho_iter * (Xi_var - Z_new)
    
    r_norm <- sqrt(sum((Xi_var - Z_new)^2))
    s_norm <- rho_iter * sqrt(sum((Z_new - Z_prev)^2))
    Z_var  <- Z_new
    
    scale_pri <- max(sqrt(length(Xi_var)), sqrt(sum(Xi_var^2)), sqrt(sum(Z_var^2)))
    scale_dual <- max(sqrt(length(Gam)), sqrt(sum(Gam^2)))
    eps_pri <- admm_tol * scale_pri
    eps_dual <- admm_tol * scale_dual

    if (track_trace) {
      inner_loss <- NA_real_
      if (!is.null(n) && is.finite(n) && n > 0) {
        Oinv_zo <- tryCatch(safe_inv_psd(Zo), error = function(e) NULL)
        if (!is.null(Oinv_zo)) {
          inner_loss <- loss_fn(Zg, Zo, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega, Oinv = Oinv_zo)
        }
      }
      inner_trace_list[[k]] <- data.frame(
        inner_iter = k,
        training_loss = inner_loss,
        r_norm = r_norm,
        s_norm = s_norm,
        residual_sum = r_norm + s_norm,
        eps_pri = eps_pri,
        eps_dual = eps_dual,
        primal_dual_gap = NA_real_,
        admm_tol_target = admm_tol,
        rho = rho_iter,
        stringsAsFactors = FALSE
      )
    }
    
    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      admm_converged <- TRUE
      break
    }

    # Boyd-style residual balancing for adaptive ADMM penalty.
    if (adaptive_rho) {
      rho_new <- rho_iter
      if (r_norm > rho_mu * s_norm) {
        rho_new <- min(rho_iter * rho_tau_inc, rho_max)
      } else if (s_norm > rho_mu * r_norm) {
        rho_new <- max(rho_iter / rho_tau_dec, rho_min)
      }
      if (abs(rho_new - rho_iter) > .Machine$double.eps) {
        # Keep scaled dual variable u = y/rho invariant when rho changes.
        Gam <- Gam * (rho_new / rho_iter)
        rho_curr <- rho_new
        Den <- outer(lamH, lamS, "*") + rho_curr
      }
    }

  }
  
  out_sub <- unpack_xi(Z_var, p)
  attr(out_sub$Omega, "state") <- list(
    psd_state = prox_state,
    Xi = Xi_var,
    Z = Z_var,
    Gam = Gam,
    rho = rho_curr
  )
  out_sub$inner_iters <- k_last
  out_sub$inner_converged <- admm_converged
  out_sub$inner_trace <- if (track_trace && k_last > 0) {
    bind_rows_safe(inner_trace_list[seq_len(k_last)])
  } else {
    NULL
  }
  out_sub
}

## =========================================================
##  Prox-Newton outer loop with optional oracle stopping
## =========================================================

pn_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 100, tol = 1e-4, track_loss = TRUE,
                            ctrl = NULL, warm_state = NULL,
                            oracle_value = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                            track_diagnostics = FALSE) {
  
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit
    if (!is.null(ctrl$tol)) tol <- ctrl$tol
  }
  
  tol_obj <- if (!is.null(ctrl$tol_obj)) ctrl$tol_obj else tol
  use_rel_obj <- isTRUE(ctrl$use_relative_obj)
  
  n <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q) else symmetrize(Omega_init)
  Oinv <- safe_inv_psd(Omega)
  
  t0 <- Sys.time()
  loss_curr <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n,
                       lambda_gamma, lambda_Omega, Oinv = Oinv)
  
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    loss_hist[1] <- loss_curr
    time_hist[1] <- 0
  }
  
  if (track_diagnostics) {
    outer_trace_list <- vector("list", max_iter + 1)
    inner_trace_list <- vector("list", max_iter)
    cond0 <- matrix_condition_stats(Omega)
    sparse0 <- matrix_sparsity_stats(gamma, Omega)
    gap0 <- if (!is.null(oracle_value) && is.finite(oracle_value)) loss_curr - oracle_value else NA_real_
    rel_gap0 <- if (is.finite(gap0) && !is.null(oracle_value) && is.finite(oracle_value)) {
      gap0 / max(1, abs(oracle_value))
    } else {
      NA_real_
    }
    outer_trace_list[[1]] <- data.frame(
      outer_iter = 0L,
      training_loss = loss_curr,
      delta_loss = NA_real_,
      relative_obj_change = NA_real_,
      oracle_gap_abs = gap0,
      oracle_gap_rel = rel_gap0,
      local_norm = NA_real_,
      step_size = NA_real_,
      line_search_halvings = NA_integer_,
      line_search_success = NA,
      line_search_forced = NA,
      armijo_lhs = NA_real_,
      armijo_rhs = NA_real_,
      inner_iters = 0L,
      inner_converged = NA,
      admm_tol_used = NA_real_,
      direction_norm_fro = NA_real_,
      refinement_triggered = NA,
      refinement_rounds = NA_integer_,
      subproblem_attempts = NA_integer_,
      min_eig = cond0$min_eig,
      max_eig = cond0$max_eig,
      kappa = cond0$kappa,
      nnz_gamma_total = sparse0$nnz_gamma_total,
      nnz_gamma_pen = sparse0$nnz_gamma_pen,
      gamma_pen_sparsity = sparse0$gamma_pen_sparsity,
      nnz_omega_offdiag_upper = sparse0$nnz_omega_offdiag_upper,
      omega_offdiag_sparsity = sparse0$omega_offdiag_sparsity,
      elapsed_sec = 0,
      stringsAsFactors = FALSE
    )
  }
  
  converged <- FALSE
  stop_reason <- NA_character_
  lambda_hist <- numeric(max_iter)
  pn_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  
  for (m in seq_len(max_iter)) {
    loss_prev <- loss_curr
    beta <- gamma %*% Oinv
    
    admm_tol_curr <- pn_inner_admm_tol_fixed
    refinement_triggered <- FALSE
    refinement_rounds <- 0L
    total_inner_iters_m <- 0L
    all_inner_converged_m <- TRUE
    subproblem_attempts <- 1L
    total_ls_halvings <- NA_integer_
    ls_success <- FALSE
    line_search_forced <- FALSE
    alpha <- NA_real_
    ls_halvings <- NA_integer_
    armijo_lhs_used <- NA_real_
    armijo_rhs_used <- NA_real_
    gamma_new <- gamma
    Omega_new <- Omega
    Oinv_new <- Oinv
    loss_new <- loss_prev
    last_alpha_try <- NA_real_
    lam <- NA_real_
    d_norm <- NA_real_
    sub <- NULL
    inner_trace_m_accum <- list()
    inner_trace_m_idx <- 0L
    
    sub <- pn_subproblem_admm(gamma, Omega, Sxx, Sxy, Syy,
                              lambda_gamma, lambda_Omega,
                              Oinv = Oinv, beta = beta,
                              rho = pn_admm_ctrl$rho, mu = 1,
                              max_admm = pn_prox_ctrl$max_admm,
                              admm_tol = admm_tol_curr,
                              adaptive_rho = pn_admm_ctrl$adaptive,
                              rho_mu = pn_admm_ctrl$rho_mu,
                              rho_tau_inc = pn_admm_ctrl$rho_tau_inc,
                              rho_tau_dec = pn_admm_ctrl$rho_tau_dec,
                              rho_min = pn_admm_ctrl$rho_min,
                              rho_max = pn_admm_ctrl$rho_max,
                              state = pn_psd_state,
                              n = n,
                              track_trace = track_diagnostics)
    pn_psd_state <- attr(sub$Omega, "state")
    
    if (!is.null(sub$inner_iters) && is.finite(sub$inner_iters)) {
      total_inner_iters_m <- total_inner_iters_m + sub$inner_iters
    }
    if (!isTRUE(sub$inner_converged)) all_inner_converged_m <- FALSE
    
    d_gamma <- sub$gamma - gamma
    d_Omega <- sub$Omega - Omega
    d_norm <- sqrt(sum(d_gamma^2) + sum(d_Omega^2))
    lam <- local_norm(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = Oinv)
    alpha_init <- (1 + lam)^(-1)
    
        ## ---------------------------------------------------------
    ## Step size (PDF Algorithm 2): strictly self-concordant damped Newton
    ## alpha = 1 / (1 + nu), where nu = local_norm(d).
    ##
    ## The paper does NOT require Armijo backtracking. We only keep a very light
    ## monotone safeguard (halve alpha) to protect against inexact subproblem
    ## solutions / numerical issues. We NEVER "forced accept" an uphill step.
    ## ---------------------------------------------------------

    alpha_try <- alpha_init
    ls_halvings <- 0L
    armijo_lhs_used <- NA_real_
    armijo_rhs_used <- NA_real_
    ls_success <- FALSE
    line_search_forced <- FALSE

    for (ls_try in 0:pn_ls_ctrl$max_halving) {
      if (alpha_try < pn_ls_ctrl$min_alpha) break

      gamma_try <- gamma + alpha_try * d_gamma
      Omega_try <- symmetrize(Omega + alpha_try * d_Omega)

      Oinv_try <- tryCatch(safe_inv_psd(Omega_try), error = function(e) NULL)
      if (is.null(Oinv_try)) {
        alpha_try <- alpha_try * 0.5
        ls_halvings <- ls_halvings + 1L
        next
      }

      loss_try <- loss_fn(gamma_try, Omega_try, Sxx, Sxy, Syy, n,
                          lambda_gamma, lambda_Omega, Oinv = Oinv_try)

      if (is.finite(loss_try) && (loss_try <= loss_prev + 2e-12)) {
        alpha <- alpha_try
        ls_success <- TRUE
        gamma_new <- gamma_try
        Omega_new <- Omega_try
        Oinv_new <- Oinv_try
        loss_new <- loss_try
        break
      }

      alpha_try <- alpha_try * 0.5
      ls_halvings <- ls_halvings + 1L
    }

    total_ls_halvings <- as.integer(ls_halvings)
    last_alpha_try <- alpha_try

    pn_step_failed <- FALSE
    if (!ls_success) {
      ## No descent step found: stop rather than accepting an uphill step.
      alpha <- 0
      gamma_new <- gamma
      Omega_new <- Omega
      Oinv_new <- Oinv
      loss_new <- loss_prev
      stop_reason <- "PN step failed to decrease objective"
      pn_step_failed <- TRUE
    }

if (track_diagnostics && !is.null(sub$inner_trace) && nrow(sub$inner_trace) > 0) {
      inner_trace_attempt <- sub$inner_trace
      inner_trace_attempt$outer_iter <- m
      inner_trace_attempt$sub_attempt <- subproblem_attempts
      inner_trace_attempt$admm_tol_target <- admm_tol_curr
      inner_trace_attempt$direction_norm_fro <- d_norm
      names(inner_trace_attempt)[names(inner_trace_attempt) == "training_loss"] <- "inner_training_loss"
      inner_trace_m_idx <- inner_trace_m_idx + 1L
      inner_trace_m_accum[[inner_trace_m_idx]] <- inner_trace_attempt
    }
    
    lambda_hist[m] <- lam
    
    if (track_diagnostics && inner_trace_m_idx > 0L) {
      inner_trace_list[[m]] <- bind_rows_safe(inner_trace_m_accum)
    }
    
    gamma <- gamma_new
    Omega <- Omega_new
    Oinv <- Oinv_new
    loss_curr <- loss_new
    
    if (!pn_step_failed) {
      if (oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE
        stop_reason <- "oracle gap reached"
      } else if (is.null(gap_tol)) {
        if (use_rel_obj) {
          if (rel_obj_change(loss_prev, loss_curr) <= tol_obj) {
            converged <- TRUE
            stop_reason <- "objective stabilized"
          }
        } else {
          if (abs(loss_curr - loss_prev) <= tol_obj) {
            converged <- TRUE
            stop_reason <- "objective stabilized"
          }
        }
      }
    }

    if (track_diagnostics) {
      gap_abs <- if (!is.null(oracle_value) && is.finite(oracle_value)) loss_curr - oracle_value else NA_real_
      gap_rel <- if (is.finite(gap_abs) && !is.null(oracle_value) && is.finite(oracle_value)) {
        gap_abs / max(1, abs(oracle_value))
      } else {
        NA_real_
      }
      cond_m <- matrix_condition_stats(Omega)
      sparse_m <- matrix_sparsity_stats(gamma, Omega)
      outer_trace_list[[m + 1]] <- data.frame(
        outer_iter = m,
        training_loss = loss_curr,
        delta_loss = loss_curr - loss_prev,
        relative_obj_change = rel_obj_change(loss_prev, loss_curr),
        oracle_gap_abs = gap_abs,
        oracle_gap_rel = gap_rel,
        local_norm = lam,
        step_size = alpha,
        line_search_halvings = as.integer(ls_halvings),
        line_search_success = ls_success,
        line_search_forced = line_search_forced,
        armijo_lhs = armijo_lhs_used,
        armijo_rhs = armijo_rhs_used,
        inner_iters = total_inner_iters_m,
        inner_converged = all_inner_converged_m,
        admm_tol_used = admm_tol_curr,
        direction_norm_fro = d_norm,
        refinement_triggered = refinement_triggered,
        refinement_rounds = refinement_rounds,
        subproblem_attempts = subproblem_attempts,
        min_eig = cond_m$min_eig,
        max_eig = cond_m$max_eig,
        kappa = cond_m$kappa,
        nnz_gamma_total = sparse_m$nnz_gamma_total,
        nnz_gamma_pen = sparse_m$nnz_gamma_pen,
        gamma_pen_sparsity = sparse_m$gamma_pen_sparsity,
        nnz_omega_offdiag_upper = sparse_m$nnz_omega_offdiag_upper,
        omega_offdiag_sparsity = sparse_m$omega_offdiag_sparsity,
        elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
        stringsAsFactors = FALSE
      )
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- loss_curr
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    
    if (m == 1L || (pn_outer_progress_every_env > 0L && (m %% pn_outer_progress_every_env == 0L))) {
      cat(sprintf("[PN outer] iter %d/%d | loss=%.6e | rel_change=%.3e | ls_ok=%s | forced=%s | inner=%d\n",
                  m, max_iter, loss_curr, rel_obj_change(loss_prev, loss_curr),
                  ifelse(ls_success, "TRUE", "FALSE"),
                  ifelse(line_search_forced, "TRUE", "FALSE"),
                  total_inner_iters_m))
      flush.console()
    }
    
    if (pn_step_failed) break
    if (converged) break
  }
  
  elapsed_time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- list(
    gamma = gamma,
    Omega = Omega,
    iters = m,
    converged = converged,
    stop_reason = ifelse(is.na(stop_reason), "max_iter reached", stop_reason),
    elapsed_time_sec = elapsed_time_sec,
    Sxx = Sxx,
    Sxy = Sxy,
    Syy = Syy,
    n = n,
    lambda_hist = lambda_hist[seq_len(m)],
    psd_state = pn_psd_state
  )
  if (track_loss) {
    out$loss_hist <- loss_hist[1:(m + 1)]
    out$time_hist <- time_hist[1:(m + 1)]
  }
  if (track_diagnostics) {
    out$outer_trace <- bind_rows_safe(outer_trace_list[seq_len(m + 1)])
    out$inner_trace <- bind_rows_safe(inner_trace_list[seq_len(m)])
    if (!is.null(out$outer_trace) && nrow(out$outer_trace) > 0) {
      out$outer_inner_map <- out$outer_trace[out$outer_trace$outer_iter > 0, c("outer_iter", "inner_iters", "inner_converged"), drop = FALSE]
    } else {
      out$outer_inner_map <- data.frame()
    }
  }
  out
}

## =========================================================
##  Lambda-path runners with optional oracle stopping
## =========================================================

run_pg_path <- function(data_list, ctrl, ref_path = NULL,
                        warm_starts = NULL,
                        oracle_vec = NULL, gap_tol = NULL, use_relative_gap = FALSE) {
  lambda_path <- data_list$lambda_path
  gamma_init <- data_list$gamma_init
  Omega_init <- data_list$Omega_init
  warm_state <- NULL
  use_shared_warm <- !is.null(warm_starts)
  if (use_shared_warm && length(warm_starts) != nrow(lambda_path)) {
    stop("warm_starts length mismatch")
  }
  results <- vector("list", nrow(lambda_path))
  
  for (i in seq_len(nrow(lambda_path))) {
    cat(sprintf("[PG path] step %d/%d\n", i, nrow(lambda_path)))
    flush.console()
    ctrl_i <- ctrl
    ctrl_i$lambda$gamma <- lambda_path$lambda_gamma[i]
    ctrl_i$lambda$Omega <- lambda_path$lambda_Omega[i]
    if (use_shared_warm) {
      gamma_init_step <- warm_starts[[i]]$gamma
      Omega_init_step <- warm_starts[[i]]$Omega
      warm_state_step <- NULL
    } else {
      gamma_init_step <- gamma_init
      Omega_init_step <- Omega_init
      warm_state_step <- warm_state
      if (i == 1 && !is.null(ref_path) && length(ref_path) >= 1) {
        seed_fit <- ref_path[[1]]
        if (!is.null(seed_fit$gamma)) gamma_init_step <- seed_fit$gamma
        if (!is.null(seed_fit$Omega)) Omega_init_step <- seed_fit$Omega
        if (!is.null(seed_fit$psd_state)) warm_state_step <- list(psd_state = seed_fit$psd_state)
      }
    }
    
    oracle_value <- if (!is.null(oracle_vec)) oracle_vec[i] else NULL
    
    fit <- pg_sparse_mvreg(data_list$X, data_list$Y,
                           lambda_path$lambda_gamma[i],
                           lambda_path$lambda_Omega[i],
                           gamma_init = gamma_init_step,
                           Omega_init = Omega_init_step,
                           max_iter = ctrl$maxit, tol = ctrl$tol, L0 = 1,
                           track_loss = TRUE,
                           ctrl = ctrl_i, warm_state = warm_state_step,
                           oracle_value = oracle_value, gap_tol = gap_tol,
                           use_relative_gap = use_relative_gap)
    
    fit$lambda_gamma <- lambda_path$lambda_gamma[i]
    fit$lambda_Omega <- lambda_path$lambda_Omega[i]
    fit$path_index <- lambda_path$step[i]
    results[[i]] <- fit
    if (!use_shared_warm) {
      gamma_init <- fit$gamma
      Omega_init <- fit$Omega
      warm_state <- list(psd_state = fit$psd_state)
    }
  }
  results
}

run_pn_path <- function(data_list, ctrl,
                        warm_starts = NULL,
                        oracle_vec = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                        track_diagnostics = FALSE) {
  lambda_path <- data_list$lambda_path
  gamma_init <- data_list$gamma_init
  Omega_init <- data_list$Omega_init
  warm_state <- NULL
  use_shared_warm <- !is.null(warm_starts)
  if (use_shared_warm && length(warm_starts) != nrow(lambda_path)) {
    stop("warm_starts length mismatch")
  }
  results <- vector("list", nrow(lambda_path))
  
  for (i in seq_len(nrow(lambda_path))) {
    cat(sprintf("[PN path] step %d/%d\n", i, nrow(lambda_path)))
    flush.console()
    ctrl_i <- ctrl
    ctrl_i$lambda$gamma <- lambda_path$lambda_gamma[i]
    ctrl_i$lambda$Omega <- lambda_path$lambda_Omega[i]
    if (use_shared_warm) {
      gamma_init_step <- warm_starts[[i]]$gamma
      Omega_init_step <- warm_starts[[i]]$Omega
      warm_state_step <- NULL
    } else {
      gamma_init_step <- gamma_init
      Omega_init_step <- Omega_init
      warm_state_step <- warm_state
    }
    oracle_value <- if (!is.null(oracle_vec)) oracle_vec[i] else NULL
    
    fit <- pn_sparse_mvreg(data_list$X, data_list$Y,
                           lambda_path$lambda_gamma[i],
                           lambda_path$lambda_Omega[i],
                           gamma_init = gamma_init_step,
                           Omega_init = Omega_init_step,
                           max_iter = ctrl$maxit, tol = ctrl$tol, track_loss = TRUE,
                           ctrl = ctrl_i, warm_state = warm_state_step,
                           oracle_value = oracle_value, gap_tol = gap_tol,
                           use_relative_gap = use_relative_gap,
                           track_diagnostics = track_diagnostics)
    
    fit$lambda_gamma <- lambda_path$lambda_gamma[i]
    fit$lambda_Omega <- lambda_path$lambda_Omega[i]
    fit$path_index <- lambda_path$step[i]
    results[[i]] <- fit
    if (!use_shared_warm) {
      gamma_init <- fit$gamma
      Omega_init <- fit$Omega
      warm_state <- list(psd_state = fit$psd_state)
    }
  }
  results
}

build_shared_warm_starts <- function(data_list, oracle_path) {
  n <- nrow(data_list$lambda_path)
  if (length(oracle_path) < n) stop("oracle_path length mismatch")
  warm_starts <- vector("list", n)
  warm_starts[[1]] <- list(gamma = data_list$gamma_init, Omega = data_list$Omega_init)
  if (n > 1) {
    for (i in 2:n) {
      warm_starts[[i]] <- list(gamma = oracle_path[[i - 1]]$gamma,
                               Omega = oracle_path[[i - 1]]$Omega)
    }
  }
  warm_starts
}

## =========================================================
##  Oracle path utility for fair comparisons
## =========================================================

compute_oracle_path <- function(data_list, ctrl_oracle) {
  pn_oracle <- run_pn_path(data_list, ctrl_oracle)
  oracle_vals <- vapply(pn_oracle, function(fit) tail(fit$loss_hist, 1), numeric(1))
  list(pn_path = pn_oracle, oracle_vals = oracle_vals)
}

## =========================================================
##  Benchmark wrapper
## =========================================================

benchmark_pg_vs_pn <- function(data_list, ctrl,
                               use_oracle = FALSE,
                               gap_tol = NULL,
                               use_relative_gap = FALSE,
                               ctrl_oracle = NULL,
                               track_pn_diagnostics = FALSE) {
  
  oracle_vals <- NULL
  pn_oracle_path <- NULL
  warm_starts <- NULL
  if (use_oracle) {
    if (is.null(ctrl_oracle)) stop("ctrl_oracle must be provided when use_oracle is TRUE")
    oracle_out <- compute_oracle_path(data_list, ctrl_oracle)
    pn_oracle_path <- oracle_out$pn_path
    oracle_vals <- oracle_out$oracle_vals
    warm_starts <- build_shared_warm_starts(data_list, pn_oracle_path)
  }
  
  pn_path <- run_pn_path(data_list, ctrl, warm_starts = warm_starts,
                         oracle_vec = oracle_vals, gap_tol = gap_tol,
                         use_relative_gap = use_relative_gap,
                         track_diagnostics = track_pn_diagnostics)
  pg_path <- run_pg_path(data_list, ctrl, warm_starts = warm_starts,
                         oracle_vec = oracle_vals, gap_tol = gap_tol,
                         use_relative_gap = use_relative_gap)
  
  pn_fit <- pn_path[[length(pn_path)]]
  pg_fit <- pg_path[[length(pg_path)]]
  
  summarize_path <- function(res_list, method_name) {
    do.call(rbind, lapply(res_list, function(fit) {
      data.frame(method = method_name,
                 path_index = fit$path_index,
                 lambda_gamma = fit$lambda_gamma,
                 lambda_Omega = fit$lambda_Omega,
                 elapsed_sec = fit$elapsed_time_sec,
                 iters = fit$iters,
                 converged = fit$converged,
                 stop_reason = fit$stop_reason,
                 final_loss = tail(fit$loss_hist, 1),
                 stringsAsFactors = FALSE)
    }))
  }
  
  path_summary <- rbind(
    summarize_path(pn_path, "Prox-Newton"),
    summarize_path(pg_path, "Prox-Gradient")
  )
  
  list(pg = pg_fit,
       pn = pn_fit,
       pg_path = pg_path,
       pn_path = pn_path,
       oracle_path = pn_oracle_path,
       oracle_vals = oracle_vals,
       path_summary = path_summary)
}

build_lambda_summary <- function(res, solver_tol) {
  summarize_one <- function(fit, method_name) {
    inner_total <- NA_real_
    if (!is.null(fit$outer_inner_map) && nrow(fit$outer_inner_map) > 0) {
      inner_total <- sum(fit$outer_inner_map$inner_iters, na.rm = TRUE)
    }
    data.frame(
      path_index = fit$path_index,
      lambda_gamma = fit$lambda_gamma,
      lambda_Omega = fit$lambda_Omega,
      method = method_name,
      outer_iters = fit$iters,
      inner_iters_total = inner_total,
      converged = fit$converged,
      stop_reason = fit$stop_reason,
      elapsed_sec = fit$elapsed_time_sec,
      final_loss = tail(fit$loss_hist, 1),
      solver_tol = solver_tol,
      stringsAsFactors = FALSE
    )
  }
  pn_df <- bind_rows_safe(lapply(res$pn_path, summarize_one, method_name = "Prox-Newton"))
  pg_df <- bind_rows_safe(lapply(res$pg_path, summarize_one, method_name = "Prox-Gradient"))
  out <- bind_rows_safe(list(pn_df, pg_df))
  out$inner_iters_total[out$method == "Prox-Gradient"] <- NA_real_
  out
}

build_pn_outer_trace <- function(pn_path, solver_tol) {
  bind_rows_safe(lapply(pn_path, function(fit) {
    if (is.null(fit$outer_trace) || nrow(fit$outer_trace) == 0) return(NULL)
    df <- fit$outer_trace
    df$path_index <- fit$path_index
    df$lambda_gamma <- fit$lambda_gamma
    df$lambda_Omega <- fit$lambda_Omega
    df$solver_tol <- solver_tol
    df
  }))
}

build_pn_inner_trace <- function(pn_path, solver_tol) {
  bind_rows_safe(lapply(pn_path, function(fit) {
    if (is.null(fit$inner_trace) || nrow(fit$inner_trace) == 0) return(NULL)
    df <- fit$inner_trace
    df$path_index <- fit$path_index
    df$lambda_gamma <- fit$lambda_gamma
    df$lambda_Omega <- fit$lambda_Omega
    df$solver_tol <- solver_tol
    df
  }))
}

build_pn_outer_inner_map <- function(pn_path, solver_tol) {
  bind_rows_safe(lapply(pn_path, function(fit) {
    if (is.null(fit$outer_inner_map) || nrow(fit$outer_inner_map) == 0) return(NULL)
    df <- fit$outer_inner_map
    df$path_index <- fit$path_index
    df$lambda_gamma <- fit$lambda_gamma
    df$lambda_Omega <- fit$lambda_Omega
    df$solver_tol <- solver_tol
    df
  }))
}

build_run_config_snapshot <- function(solver_tol, prox_tol, solver_maxit, prox_max_admm,
                                      prox_eta, pn_prox_eta, gap_tol, use_oracle,
                                      use_relative_gap, pn_ls_max_halving,
                                      pn_ls_min_alpha, pn_admm_rho,
                                      pn_inner_admm_tol,
                                      pn_ls_armijo_c,
                                      pn_admm_adaptive) {
  data.frame(
    seed = shared_seed,
    n = shared_params$n,
    p = shared_params$p,
    q = shared_params$q,
    path_len = path_len,
    ratio_lambda = shared_params$ratio_lambda,
    gamma_sparsity = shared_params$gamma_sparsity,
    gamma_mag = paste(shared_params$gamma_mag, collapse = ","),
    gamma_mag_lo = min(shared_params$gamma_mag),
    gamma_mag_hi = max(shared_params$gamma_mag),
    omega_kappa = shared_params$omega_kappa,
    omega_sparsity = shared_params$omega_sparsity,
    solver_tol = solver_tol,
    prox_tol = prox_tol,
    solver_maxit = solver_maxit,
    prox_max_admm = prox_max_admm,
    prox_eta = prox_eta,
    pn_prox_eta = pn_prox_eta,
    gap_tol = gap_tol,
    use_oracle = use_oracle,
    use_relative_gap = use_relative_gap,
    pn_ls_max_halving = pn_ls_max_halving,
    pn_ls_min_alpha = pn_ls_min_alpha,
    pn_admm_rho = pn_admm_rho,
    pn_inner_admm_tol = pn_inner_admm_tol,
    pn_ls_armijo_c = pn_ls_armijo_c,
    pn_admm_adaptive = pn_admm_adaptive,
    stringsAsFactors = FALSE
  )
}

## =========================================================
##  Example run and plots
## =========================================================

suppressPackageStartupMessages(library(ggplot2))

format_tol_tag <- function(tol) {
  tag <- format(tol, scientific = TRUE)
  tag <- gsub("\\+", "", tag)
  gsub("\\.", "p", tag)
}

solver_tols <- c(1e-5)
solver_maxit <- 1000L
prox_max_admm <- 1000L
prox_eta <- 1e-8
pn_prox_eta <- 1e-8

settings_table <- data.frame(
  solver_tol = solver_tols,
  prox_tol = 0.1 * solver_tols,
  solver_maxit = solver_maxit,
  prox_max_admm = prox_max_admm,
  prox_eta = prox_eta,
  pn_prox_eta = pn_prox_eta,
  pn_inner_admm_tol = pn_inner_admm_tol_fixed,
  pn_ls_armijo_c = armijo_c,
  n = shared_params$n,
  p = shared_params$p,
  q = shared_params$q,
  path_len = path_len,
  stringsAsFactors = FALSE
)
print(settings_table)

method_levels <- c("Prox-Newton", "Prox-Gradient")
all_results <- list()
time_tables <- list()
lambda_summary_tables <- list()
pn_outer_trace_tables <- list()
pn_inner_trace_tables <- list()
pn_outer_inner_map_tables <- list()
config_snapshot_tables <- list()

for (solver_tol in solver_tols) {
  cat(sprintf("\n=== Running solver_tol = %.1e (prox_tol = %.1e) ===\n",
              solver_tol, 0.1 * solver_tol))
  flush.console()
  ctrl <- list(
    lambda = list(gamma = lambda_gamma_seq[1], Omega = lambda_Omega_seq[1]),
    maxit = solver_maxit,
    tol = solver_tol,
    tol_obj = solver_tol,
    tol_gm = solver_tol,
    use_relative_obj = TRUE
  )

  prox_ctrl <<- list(max_admm = prox_max_admm, tol = 0.1 * solver_tol, eta = prox_eta)
  pn_prox_ctrl <<- list(max_admm = prox_max_admm, tol = 0.1 * solver_tol, eta = pn_prox_eta)

  # Stop by objective stabilization only (no oracle-gap stopping).
  use_oracle_mode <- FALSE
  use_relative_gap_mode <- FALSE
  gap_tol_mode <- NULL
  ctrl_oracle <- NULL

  res <- benchmark_pg_vs_pn(
    shared_data, ctrl,
    use_oracle = use_oracle_mode,
    gap_tol = gap_tol_mode,
    use_relative_gap = use_relative_gap_mode,
    ctrl_oracle = ctrl_oracle,
    track_pn_diagnostics = TRUE
  )
  res$path_summary$solver_tol <- solver_tol
  all_results[[as.character(solver_tol)]] <- res
  tol_tag <- format_tol_tag(solver_tol)
  
  lambda_summary_df <- build_lambda_summary(res, solver_tol)
  pn_outer_trace_df <- build_pn_outer_trace(res$pn_path, solver_tol)
  pn_inner_trace_df <- build_pn_inner_trace(res$pn_path, solver_tol)
  pn_outer_inner_map_df <- build_pn_outer_inner_map(res$pn_path, solver_tol)
  config_snapshot_df <- build_run_config_snapshot(
    solver_tol = solver_tol,
    prox_tol = 0.1 * solver_tol,
    solver_maxit = solver_maxit,
    prox_max_admm = prox_max_admm,
    prox_eta = prox_eta,
    pn_prox_eta = pn_prox_eta,
    gap_tol = if (is.null(gap_tol_mode)) NA_real_ else gap_tol_mode,
    use_oracle = use_oracle_mode,
    use_relative_gap = use_relative_gap_mode,
    pn_ls_max_halving = pn_ls_ctrl$max_halving,
    pn_ls_min_alpha = pn_ls_ctrl$min_alpha,
    pn_admm_rho = pn_admm_ctrl$rho,
    pn_inner_admm_tol = pn_inner_admm_tol_fixed,
    pn_ls_armijo_c = armijo_c,
    pn_admm_adaptive = pn_admm_ctrl$adaptive
  )
  
  lambda_summary_tables[[as.character(solver_tol)]] <- lambda_summary_df
  pn_outer_trace_tables[[as.character(solver_tol)]] <- pn_outer_trace_df
  pn_inner_trace_tables[[as.character(solver_tol)]] <- pn_inner_trace_df
  pn_outer_inner_map_tables[[as.character(solver_tol)]] <- pn_outer_inner_map_df
  config_snapshot_tables[[as.character(solver_tol)]] <- config_snapshot_df
  
  write.csv(lambda_summary_df, sprintf("lambda_summary_tol_%s.csv", tol_tag), row.names = FALSE)
  write.csv(pn_outer_trace_df, sprintf("pn_outer_trace_tol_%s.csv", tol_tag), row.names = FALSE)
  write.csv(pn_inner_trace_df, sprintf("pn_inner_trace_tol_%s.csv", tol_tag), row.names = FALSE)
  write.csv(pn_outer_inner_map_df, sprintf("pn_outer_inner_map_tol_%s.csv", tol_tag), row.names = FALSE)
  write.csv(config_snapshot_df, sprintf("run_config_snapshot_tol_%s.csv", tol_tag), row.names = FALSE)
  saveRDS(list(
    lambda_summary = lambda_summary_df,
    pn_outer_trace = pn_outer_trace_df,
    pn_inner_trace = pn_inner_trace_df,
    pn_outer_inner_map = pn_outer_inner_map_df,
    run_config_snapshot = config_snapshot_df
  ), file = sprintf("pn_diagnostics_tol_%s.rds", tol_tag))

  loss_plot_df <- res$path_summary
  loss_plot_df$method <- factor(loss_plot_df$method, levels = method_levels)
  loss_plot_df$log_lambda <- log(loss_plot_df$lambda_gamma)

  final_loss_plot <- ggplot(loss_plot_df,
                            aes(x = log_lambda, y = final_loss,
                                colour = method, shape = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2, stroke = 0.9) +
    scale_shape_manual(values = c("Prox-Newton" = 4, "Prox-Gradient" = 16)) +
    labs(x = expression(log(lambda)),
         y = "Training Loss",
         colour = NULL,
         shape = NULL) +
    theme_minimal(base_size = 13)

  iter_plot <- ggplot(loss_plot_df,
                      aes(x = log_lambda, y = iters,
                          colour = method, shape = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2, stroke = 0.9) +
    scale_shape_manual(values = c("Prox-Newton" = 4, "Prox-Gradient" = 16)) +
    labs(x = expression(log(lambda[gamma])),
         y = "number of iterations",
         colour = NULL,
         shape = NULL) +
    theme_minimal(base_size = 13)

  ggsave(paste0("pn_pg_final_loss_vs_lambda_tol_", tol_tag, ".png"),
         final_loss_plot, width = 9, height = 4.5, dpi = 300)
  ggsave(paste0("pn_pg_iters_vs_lambda_tol_", tol_tag, ".png"),
         iter_plot, width = 9, height = 4.5, dpi = 300)

  print(final_loss_plot)
  print(iter_plot)

  total_time_by_method <- aggregate(elapsed_sec ~ method, data = res$path_summary, sum)
  total_time_by_method$solver_tol <- solver_tol
  time_tables[[as.character(solver_tol)]] <- total_time_by_method
  print(total_time_by_method)

  saveRDS(res, file = sprintf("res_tol_%s.rds", tol_tag))
  saveRDS(list(
    pn_path = res$pn_path,
    pg_path = res$pg_path,
    lambda_path = shared_data$lambda_path
  ), file = sprintf("lambda_path_estimators_tol_%s.rds", tol_tag))
}

time_summary <- do.call(rbind, time_tables)
row.names(time_summary) <- NULL
time_summary <- time_summary[, c("solver_tol", "method", "elapsed_sec")]
names(time_summary) <- c("solver_tol", "method", "total_elapsed_sec")
print(time_summary)
write.csv(time_summary, "time_method_summary.csv", row.names = FALSE)

lambda_summary_all <- bind_rows_safe(lambda_summary_tables)
pn_outer_trace_all <- bind_rows_safe(pn_outer_trace_tables)
pn_inner_trace_all <- bind_rows_safe(pn_inner_trace_tables)
pn_outer_inner_map_all <- bind_rows_safe(pn_outer_inner_map_tables)
run_config_snapshot_all <- bind_rows_safe(config_snapshot_tables)

write.csv(lambda_summary_all, "lambda_summary_all.csv", row.names = FALSE)
write.csv(pn_outer_trace_all, "pn_outer_trace_all.csv", row.names = FALSE)
write.csv(pn_inner_trace_all, "pn_inner_trace_all.csv", row.names = FALSE)
write.csv(pn_outer_inner_map_all, "pn_outer_inner_map_all.csv", row.names = FALSE)
write.csv(run_config_snapshot_all, "run_config_snapshot_all.csv", row.names = FALSE)

saveRDS(list(
  lambda_summary_all = lambda_summary_all,
  pn_outer_trace_all = pn_outer_trace_all,
  pn_inner_trace_all = pn_inner_trace_all,
  pn_outer_inner_map_all = pn_outer_inner_map_all,
  run_config_snapshot_all = run_config_snapshot_all
), file = "pn_diagnostics_all.rds")
