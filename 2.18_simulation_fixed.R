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

safe_inv_psd <- function(A, eps = NULL, jitter_init = NULL, jitter_max = NULL) {
  if (is.null(eps)) eps <- GLOBAL_NUMERICS$safe_inv_eps
  if (is.null(jitter_init)) jitter_init <- GLOBAL_NUMERICS$safe_inv_jitter_init
  if (is.null(jitter_max)) jitter_max <- GLOBAL_NUMERICS$safe_inv_jitter_max
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

safe_logdet_spd <- function(A, jitter_init = NULL, jitter_max = NULL, eps = NULL) {
  if (is.null(jitter_init)) jitter_init <- GLOBAL_NUMERICS$safe_logdet_jitter_init
  if (is.null(jitter_max)) jitter_max <- GLOBAL_NUMERICS$safe_logdet_jitter_max
  if (is.null(eps)) eps <- GLOBAL_NUMERICS$safe_logdet_eps
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

matrix_condition_stats <- function(A, eps = NULL) {
  if (is.null(eps)) eps <- GLOBAL_NUMERICS$matrix_condition_eps
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

get_env_numeric <- function(name, default, integer = FALSE) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw) || !nzchar(raw)) return(default)
  parsed <- suppressWarnings(if (integer) as.integer(raw) else as.numeric(raw))
  if (!is.finite(parsed) || is.na(parsed)) return(default)
  parsed
}

flatten_named_list <- function(x, prefix = character()) {
  rows <- list()
  row_idx <- 0L

  recurse <- function(obj, path) {
    if (is.list(obj) && length(obj) > 0) {
      nms <- names(obj)
      if (is.null(nms)) nms <- rep("", length(obj))
      for (i in seq_along(obj)) {
        nm <- nms[[i]]
        if (!nzchar(nm)) nm <- sprintf("[%d]", i)
        recurse(obj[[i]], c(path, nm))
      }
      return(invisible(NULL))
    }
    row_idx <<- row_idx + 1L
    key <- paste(path, collapse = ".")
    value <- if (is.null(obj) || length(obj) == 0) {
      "NULL"
    } else if (length(obj) > 1) {
      paste(as.character(obj), collapse = ",")
    } else {
      as.character(obj)
    }
    rows[[row_idx]] <<- data.frame(
      parameter = key,
      value = value,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  recurse(x, prefix)
  bind_rows_safe(rows)
}

print_sim_config <- function(cfg) {
  cat("=== Active SIM_CONFIG ===\n")
  print(cfg)
  cat("=========================\n")
  flush.console()
}

build_parameter_catalog <- function(cfg) {
  flat <- flatten_named_list(cfg)
  flat$timestamp_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  flat
}

get_default_sim_config <- function() {
  solver_tol <- get_env_numeric("SOLVER_TOL", 1e-7)
  solver_maxit <- get_env_numeric("SOLVER_MAXIT", 5000L, integer = TRUE)
  pg_min_iter <- get_env_numeric("PG_MIN_ITER", 100L, integer = TRUE)
  if (!is.finite(pg_min_iter) || is.na(pg_min_iter)) pg_min_iter <- 100L
  pg_min_iter <- as.integer(pg_min_iter)
  if (pg_min_iter < 0L) pg_min_iter <- 0L
  if (pg_min_iter > solver_maxit) {
    warning(sprintf("PG_MIN_ITER (%d) exceeds SOLVER_MAXIT (%d); clipping to SOLVER_MAXIT.", pg_min_iter, solver_maxit))
    pg_min_iter <- as.integer(solver_maxit)
  }
  prox_max_admm <- get_env_numeric("PROX_MAX_ADMM", 1000L, integer = TRUE)
  prox_tol <- get_env_numeric("PROX_TOL", 1e-8)
  prox_eta <- get_env_numeric("PROX_ETA", 1e-8)
  pn_prox_eta <- get_env_numeric("PN_PROX_ETA", prox_eta)
  pn_tol_admm_base <- get_env_numeric("PN_TOL_ADMM_BASE", 1e-8)
  pn_tol_admm_floor <- get_env_numeric("PN_TOL_ADMM_FLOOR", 1e-10)
  omega_tol_floor <- get_env_numeric("OMEGA_TOL_FLOOR", 1e-10)
  pn_outer_progress_every <- get_env_numeric("PN_OUTER_PROGRESS_EVERY", 50L, integer = TRUE)
  pn_ls_max_halving <- get_env_numeric("PN_LS_MAX_HALVING", 3L, integer = TRUE)
  pn_ls_min_alpha <- get_env_numeric("PN_LS_MIN_ALPHA", 1e-8)
  oracle_fixed_iters <- get_env_numeric("ORACLE_FIXED_ITERS", 100L, integer = TRUE)
  oracle_run_flag <- get_env_numeric("ORACLE_RUN_FOR_EXPORTS", 1L, integer = TRUE)
  oracle_reuse_flag <- get_env_numeric("ORACLE_REUSE_EXISTING", 0L, integer = TRUE)
  oracle_reuse_file <- Sys.getenv("ORACLE_REUSE_FILE", unset = "")
  pg_omega_mode <- Sys.getenv("PG_OMEGA_STOP_MODE", unset = "current")
  pn_omega_mode <- Sys.getenv("PN_OMEGA_STOP_MODE", unset = "duality_gap")
  pn_sub_mode <- Sys.getenv("PN_SUB_ADMM_STOP_MODE", unset = "residual_plus_proxy")

  list(
    simulation = list(
      seed = 20251025L,
      n = 1000L,
      p = 50L,
      q = 50L,
      gamma_sparsity = 0.05,
      gamma_mag = c(1, 2),
      gamma_signed = TRUE,
      intercept_zero = TRUE,
      center_X = TRUE,
      omega_kappa = 60,
      omega_sparsity = 0.05,
      omega_base_val = 0.5
    ),
    lambda_path = list(
      path_len = 60L,
      ratio_lambda = 1e-4
    ),
    outer_pg = list(
      max_iter = solver_maxit,
      min_iter = pg_min_iter,
      tol_PG_sub = solver_tol,
      tol_obj = solver_tol,
      tol_gm = solver_tol,
      use_relative_obj = TRUE,
      L0 = 1,
      stop_modes = c("relative_change", "oracle_change", "l2_step_norm"),
      stop_logic = "any",
      l2_norm_tol = solver_tol
    ),
    outer_pn = list(
      max_iter = solver_maxit,
      tol_PN_sub = solver_tol,
      stop_modes = c("local_norm", "relative_change", "oracle_change"),
      stop_logic = "any",
      line_search = list(
        max_halving = max(pn_ls_max_halving, 0L),
        min_alpha = max(pn_ls_min_alpha, 1e-16),
        descent_eps = 2e-12
      ),
      progress_every = max(1L, pn_outer_progress_every)
    ),
    pn_subproblem_admm = list(
      max_iter = prox_max_admm,
      tol_PN_ADMM_base = pn_tol_admm_base,
      tol_PN_ADMM_floor = pn_tol_admm_floor,
      refine_trigger_halvings = 3L,
      refine_factor = 10,
      max_refinement_rounds = 3L,
      stop_mode = pn_sub_mode,
      adaptive_rho = list(
        enabled = TRUE,
        rho = 1,
        rho_mu = 10,
        rho_tau_inc = 2,
        rho_tau_dec = 2,
        rho_min = 1e-8,
        rho_max = 1e8
      )
    ),
    omega_prox_admm_pg = list(
      max_iter = prox_max_admm,
      tol_Omega_ADMM = prox_tol,
      tol_Omega_ADMM_floor = omega_tol_floor,
      mu = 1,
      stop_mode = pg_omega_mode,
      eig_floor = prox_eta
    ),
    omega_prox_admm_pn = list(
      max_iter = prox_max_admm,
      tol_Omega_ADMM = prox_tol,
      tol_Omega_ADMM_floor = omega_tol_floor,
      mu = 1,
      stop_mode = pn_omega_mode,
      eig_floor = pn_prox_eta
    ),
    oracle_eval = list(
      fixed_iters = max(1L, oracle_fixed_iters),
      run_for_exports = as.logical(oracle_run_flag > 0),
      reuse_existing = as.logical(oracle_reuse_flag > 0),
      reuse_file = oracle_reuse_file,
      nondecrease_eps = 1e-12,
      loss_definition = "best_within_fixed_iters",
      scope = "whole_path"
    ),
    numerics = list(
      safe_inv_eps = 1e-8,
      safe_inv_jitter_init = 0,
      safe_inv_jitter_max = 1e-5,
      safe_logdet_jitter_init = 0,
      safe_logdet_jitter_max = 1e-5,
      safe_logdet_eps = 1e-10,
      matrix_condition_eps = 1e-12,
      pn_subproblem_eig_floor = 1e-12,
      null_corner_eps = 1e-8
    ),
    benchmark = list(
      use_oracle_mode = FALSE,
      use_relative_gap_mode = FALSE,
      gap_tol_mode = NULL
    ),
    reporting = list(
      solver_tols = c(solver_tol),
      method_levels = c("Prox-Newton", "Prox-Gradient"),
      write_parameter_catalog = TRUE
    )
  )
}

GLOBAL_NUMERICS <- list(
  safe_inv_eps = 1e-8,
  safe_inv_jitter_init = 0,
  safe_inv_jitter_max = 1e-5,
  safe_logdet_jitter_init = 0,
  safe_logdet_jitter_max = 1e-5,
  safe_logdet_eps = 1e-10,
  matrix_condition_eps = 1e-12
)

set_numeric_controls <- function(cfg) {
  if (is.null(cfg$numerics)) return(invisible(NULL))
  for (nm in names(cfg$numerics)) {
    GLOBAL_NUMERICS[[nm]] <<- cfg$numerics[[nm]]
  }
  invisible(NULL)
}

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
    eps <- GLOBAL_NUMERICS$matrix_condition_eps
    eig_vals <- suppressWarnings(eigen(symmetrize(Omega), symmetric = TRUE, only.values = TRUE)$values)
    eig_vals <- pmax(eig_vals, eps)
    attr(Omega, "cond_est") <- max(eig_vals) / min(eig_vals)
  }
  Omega
}

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

build_shared_data <- function(cfg) {
  sim <- cfg$simulation
  p_aug <- sim$p + 1L

  set.seed(sim$seed)
  gamma_true <- generate_sparse_gamma(
    p_aug,
    sim$q,
    sparsity = sim$gamma_sparsity,
    mag_range = sim$gamma_mag,
    intercept_zero = isTRUE(sim$intercept_zero),
    signed = isTRUE(sim$gamma_signed)
  )
  Beta_true <- gamma_true[-1, , drop = FALSE]
  Omega_true <- generate_sparse_omega(
    sim$q,
    kappa_target = sim$omega_kappa,
    offdiag_prob = sim$omega_sparsity,
    base_val = sim$omega_base_val
  )
  Sigma_true <- safe_inv_psd(Omega_true)
  truth <- list(
    gamma = gamma_true,
    Omega = Omega_true,
    cond_number = attr(Omega_true, "cond_est")
  )

  X <- matrix(rnorm(sim$n * sim$p), sim$n, sim$p)
  X_use <- if (isTRUE(sim$center_X)) scale(X, center = TRUE, scale = FALSE) else X
  Xc <- cbind(1, X_use)

  suppressPackageStartupMessages(library(MASS))
  Y <- X %*% Beta_true + MASS::mvrnorm(sim$n, mu = rep(0, sim$q), Sigma = Sigma_true)

  Sxy <- crossprod(Xc, Y) / sim$n
  Syy <- crossprod(Y) / sim$n
  lambda_gamma_0 <- 2 * max(abs(Sxy[-1, , drop = FALSE]))
  tmp_off <- Syy - diag(diag(Syy))
  lambda_Omega_0 <- 2 * max(abs(tmp_off))
  path_len <- cfg$lambda_path$path_len
  ratio_lambda <- cfg$lambda_path$ratio_lambda
  lambda_gamma_seq <- exp(seq(log(lambda_gamma_0),
                              log(lambda_gamma_0 * ratio_lambda),
                              length.out = path_len))
  lambda_Omega_seq <- exp(seq(log(lambda_Omega_0),
                              log(lambda_Omega_0 * ratio_lambda),
                              length.out = path_len))
  lambda_path <- data.frame(
    step = seq_len(path_len),
    lambda_gamma = lambda_gamma_seq,
    lambda_Omega = lambda_Omega_seq
  )

  null_corner_init <- compute_null_corner_init(
    crossprod(Xc) / sim$n,
    Sxy,
    Syy,
    unpenalized = 1L,
    eps = cfg$numerics$null_corner_eps
  )

  shared_data <- list(
    X = Xc,
    Y = Y,
    lambda = list(gamma = lambda_gamma_seq[1], Omega = lambda_Omega_seq[1]),
    lambda_path = lambda_path,
    gamma_init = null_corner_init$gamma,
    Omega_init = null_corner_init$Omega,
    params = sim,
    truth = truth
  )

  list(
    shared_data = shared_data,
    params = sim,
    lambda_gamma_seq = lambda_gamma_seq,
    lambda_Omega_seq = lambda_Omega_seq,
    path_len = path_len
  )
}

SIM_CONFIG <- get_default_sim_config()
set_numeric_controls(SIM_CONFIG)
shared_seed <- SIM_CONFIG$simulation$seed
shared_params <- SIM_CONFIG$simulation
sim_setup <- build_shared_data(SIM_CONFIG)
shared_data <- sim_setup$shared_data
lambda_gamma_seq <- sim_setup$lambda_gamma_seq
lambda_Omega_seq <- sim_setup$lambda_Omega_seq
path_len <- sim_setup$path_len

if (isTRUE(SIM_CONFIG$reporting$write_parameter_catalog)) {
  param_catalog <- build_parameter_catalog(SIM_CONFIG)
  write.csv(param_catalog, "simulation_parameter_catalog.csv", row.names = FALSE)
  saveRDS(param_catalog, "simulation_parameter_catalog.rds")
}
print_sim_config(SIM_CONFIG)

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
                                cfg_omega_admm = NULL,
                                stop_mode = NULL,
                                Omega_init = NULL,
                                A_init = NULL) {
  cfg <- if (is.null(cfg_omega_admm)) {
    list(max_iter = 1000L, tol_Omega_ADMM = 1e-6, mu = 1, stop_mode = "duality_gap", eig_floor = 0)
  } else {
    cfg_omega_admm
  }
  max_admm <- as.integer(if (!is.null(cfg$max_iter)) cfg$max_iter else 1000L)
  tol <- if (!is.null(cfg$tol_Omega_ADMM)) cfg$tol_Omega_ADMM else 1e-6
  mu <- if (!is.null(cfg$mu)) cfg$mu else 1
  eig_floor <- if (!is.null(cfg$eig_floor)) cfg$eig_floor else 0
  mode <- if (!is.null(stop_mode)) stop_mode else if (!is.null(cfg$stop_mode)) cfg$stop_mode else "duality_gap"

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

  omega_hat <- symmetrize(soft_thresh_offdiag(Z, tau))
  q <- nrow(omega_hat)
  if (is_pd_by_chol_rank(omega_hat - eig_floor * diag(q))) {
    attr(omega_hat, "state") <- list(
      Omega = omega_hat,
      A = matrix(0, q, q),
      mode = mode,
      converged = TRUE,
      iters = 0L,
      final_metric = 0
    )
    return(omega_hat)
  }

  omega_hat <- clip_eig(omega_hat, eps = eig_floor)
  Omega <- if (is.null(Omega_init)) omega_hat else symmetrize(Omega_init)
  A <- if (is.null(A_init)) matrix(0, nrow(Omega), ncol(Omega)) else symmetrize(A_init)
  prev_Omega <- Omega
  final_metric <- NA_real_
  converged <- FALSE
  iters <- 0L

  for (l in seq_len(max_admm)) {
    iters <- l
    K <- clip_eig(Omega + mu * A, eps = eig_floor)
    T <- (K + mu * (Z - A)) / (1 + mu)
    Omega_new <- symmetrize(soft_thresh_offdiag(T, (mu * tau) / (1 + mu)))
    A <- A - (K - Omega_new) / mu

    rel_change <- sqrt(sum((Omega_new - prev_Omega)^2)) / max(1, sqrt(sum(prev_Omega^2)))
    p_val <- primal_value(Omega_new)
    Y_feas <- dual_feasible_Y(A, mu, tau)
    d_val <- dual_value(Y_feas)
    gap <- max(p_val - d_val, 0)
    gap_rel <- if (is.finite(p_val)) gap / (1 + abs(p_val)) else Inf
    final_metric <- if (identical(mode, "duality_gap")) gap_rel else rel_change

    if (identical(mode, "duality_gap")) {
      if (is.finite(gap_rel) && gap_rel <= tol) {
        converged <- TRUE
        Omega <- Omega_new
        break
      }
    } else {
      if (is.finite(rel_change) && rel_change <= tol) {
        converged <- TRUE
        Omega <- Omega_new
        break
      }
    }

    prev_Omega <- Omega
    Omega <- Omega_new
  }

  K_out <- clip_eig(symmetrize(Omega), eps = eig_floor)
  attr(K_out, "state") <- list(
    Omega = K_out,
    A = A,
    mode = mode,
    converged = converged,
    iters = iters,
    final_metric = final_metric
  )
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

evaluate_pn_stop_modes <- function(loss_prev, loss_curr, local_norm_value,
                                   oracle_value, tol_pn_sub,
                                   stop_modes = c("local_norm", "relative_change", "oracle_change"),
                                   stop_logic = "any") {
  rel_change <- (loss_prev - loss_curr) / max(1, abs(loss_prev))
  oracle_change <- if (!is.null(oracle_value) && is.finite(oracle_value)) {
    (loss_curr - oracle_value) / max(1, abs(oracle_value))
  } else {
    Inf
  }
  checks <- list(
    local_norm = is.finite(local_norm_value) && local_norm_value <= tol_pn_sub,
    relative_change = is.finite(rel_change) && rel_change >= 0 && rel_change <= tol_pn_sub,
    oracle_change = is.finite(oracle_change) && oracle_change <= tol_pn_sub
  )
  active_modes <- unique(stop_modes)
  active_modes <- active_modes[active_modes %in% names(checks)]
  if (length(active_modes) == 0) active_modes <- "relative_change"
  active_vals <- vapply(active_modes, function(m) checks[[m]], logical(1))
  should_stop <- if (identical(stop_logic, "all")) all(active_vals) else any(active_vals)
  hit_modes <- active_modes[active_vals]
  list(
    should_stop = should_stop,
    hit_modes = hit_modes,
    rel_change = rel_change,
    oracle_change = oracle_change
  )
}

evaluate_pg_stop_modes <- function(loss_prev, loss_curr, l2_step_norm,
                                   oracle_value, tol_pg_sub,
                                   stop_modes = c("relative_change", "oracle_change", "l2_step_norm"),
                                   stop_logic = "any",
                                   l2_norm_tol = NULL) {
  rel_change <- (loss_prev - loss_curr) / max(1, abs(loss_prev))
  oracle_change <- if (!is.null(oracle_value) && is.finite(oracle_value)) {
    (loss_curr - oracle_value) / max(1, abs(oracle_value))
  } else {
    Inf
  }
  l2_tol <- if (is.null(l2_norm_tol) || !is.finite(l2_norm_tol)) tol_pg_sub else l2_norm_tol
  checks <- list(
    relative_change = is.finite(rel_change) && rel_change >= 0 && rel_change <= tol_pg_sub,
    oracle_change = is.finite(oracle_change) && oracle_change <= tol_pg_sub,
    l2_step_norm = is.finite(l2_step_norm) && l2_step_norm <= l2_tol
  )
  active_modes <- unique(stop_modes)
  active_modes <- active_modes[active_modes %in% names(checks)]
  if (length(active_modes) == 0) active_modes <- "relative_change"
  active_vals <- vapply(active_modes, function(m) checks[[m]], logical(1))
  should_stop <- if (identical(stop_logic, "all")) all(active_vals) else any(active_vals)
  hit_modes <- active_modes[active_vals]
  list(
    should_stop = should_stop,
    hit_modes = hit_modes,
    rel_change = rel_change,
    oracle_change = oracle_change,
    l2_step_norm = l2_step_norm
  )
}

## =========================================================
##  Prox-Gradient with improved stopping
## =========================================================

pg_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 200, tol = 1e-4, L0 = 1,
                            track_loss = TRUE,
                            ctrl = NULL, warm_state = NULL,
                            oracle_value = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                            cfg = SIM_CONFIG) {
  
  if (!is.null(cfg$outer_pg)) {
    max_iter <- cfg$outer_pg$max_iter
    tol <- cfg$outer_pg$tol_PG_sub
    L0 <- cfg$outer_pg$L0
  }
  min_iter_pg <- if (!is.null(cfg$outer_pg$min_iter)) as.integer(cfg$outer_pg$min_iter) else 100L
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit
    if (!is.null(ctrl$min_iter)) min_iter_pg <- as.integer(ctrl$min_iter)
    if (!is.null(ctrl$tol)) tol <- ctrl$tol
    if (!is.null(ctrl$L0)) L0 <- ctrl$L0
  }
  if (!is.finite(min_iter_pg) || is.na(min_iter_pg)) min_iter_pg <- 100L
  min_iter_pg <- as.integer(min_iter_pg)
  if (min_iter_pg < 0L) min_iter_pg <- 0L
  if (min_iter_pg > max_iter) {
    warning(sprintf("PG min_iter (%d) exceeds max_iter (%d); clipping to max_iter.", min_iter_pg, max_iter))
    min_iter_pg <- as.integer(max_iter)
  }
  
  pg_stop_modes <- if (!is.null(cfg$outer_pg$stop_modes)) cfg$outer_pg$stop_modes else c("relative_change")
  pg_stop_logic <- if (!is.null(cfg$outer_pg$stop_logic)) cfg$outer_pg$stop_logic else "any"
  l2_norm_tol <- if (!is.null(cfg$outer_pg$l2_norm_tol)) cfg$outer_pg$l2_norm_tol else tol
  if (!is.null(ctrl$stop_modes)) pg_stop_modes <- ctrl$stop_modes
  if (!is.null(ctrl$stop_logic)) pg_stop_logic <- ctrl$stop_logic
  if (!is.null(ctrl$l2_norm_tol)) l2_norm_tol <- ctrl$l2_norm_tol
  pg_omega_cfg <- cfg$omega_prox_admm_pg
  
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
        cfg_omega_admm = pg_omega_cfg,
        Omega_init = if (is.null(pg_psd_state)) Omega else pg_psd_state$Omega,
        A_init = if (is.null(pg_psd_state)) NULL else pg_psd_state$A
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
      
      alpha_m <- if (lambda_m <= .Machine$double.eps) {
        1
      } else {
        (beta_m^2) / (lambda_m * (lambda_m + beta_m^2))
      }
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

      l2_step_norm <- sqrt(sum(d_gamma_m^2) + sum(d_Omega_m^2))
      stop_eval <- evaluate_pg_stop_modes(
        loss_prev = loss_prev,
        loss_curr = loss_curr,
        l2_step_norm = l2_step_norm,
        oracle_value = oracle_value,
        tol_pg_sub = tol,
        stop_modes = pg_stop_modes,
        stop_logic = pg_stop_logic,
        l2_norm_tol = l2_norm_tol
      )
      if (m >= min_iter_pg && isTRUE(stop_eval$should_stop)) {
        converged <- TRUE
        stop_reason <- sprintf("PG stop modes reached: %s", paste(stop_eval$hit_modes, collapse = "|"))
        break
      }

      ## Optional extra oracle-gap stop (independent of stop-mode evaluator).
      if (m >= min_iter_pg && oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE
        stop_reason <- "oracle gap reached"
        break
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
                               max_admm = 10000, admm_tol = NULL,
                               eig_floor = 1e-12,
                               adaptive_rho = TRUE,
                               rho_mu = 10,
                               rho_tau_inc = 2,
                               rho_tau_dec = 2,
                               rho_min = 1e-8,
                               rho_max = 1e8,
                               state = NULL,
                               n = NULL,
                               track_trace = FALSE,
                               cfg = SIM_CONFIG) {
  p <- nrow(gamma); q <- ncol(gamma)

  sub_cfg <- cfg$pn_subproblem_admm
  rho_cfg <- sub_cfg$adaptive_rho
  max_admm <- if (!is.null(sub_cfg$max_iter)) sub_cfg$max_iter else max_admm
  if (is.null(admm_tol)) {
    admm_tol <- if (!is.null(sub_cfg$tol_PN_ADMM_base)) sub_cfg$tol_PN_ADMM_base else 1e-4
  }
  eig_floor <- if (!is.null(cfg$numerics$pn_subproblem_eig_floor)) cfg$numerics$pn_subproblem_eig_floor else eig_floor
  adaptive_rho <- isTRUE(rho_cfg$enabled)
  rho <- if (!is.null(rho_cfg$rho)) rho_cfg$rho else rho
  rho_mu <- if (!is.null(rho_cfg$rho_mu)) rho_cfg$rho_mu else rho_mu
  rho_tau_inc <- if (!is.null(rho_cfg$rho_tau_inc)) rho_cfg$rho_tau_inc else rho_tau_inc
  rho_tau_dec <- if (!is.null(rho_cfg$rho_tau_dec)) rho_cfg$rho_tau_dec else rho_tau_dec
  rho_min <- if (!is.null(rho_cfg$rho_min)) rho_cfg$rho_min else rho_min
  rho_max <- if (!is.null(rho_cfg$rho_max)) rho_cfg$rho_max else rho_max
  stop_mode <- if (!is.null(sub_cfg$stop_mode)) sub_cfg$stop_mode else "residual_plus_proxy"
  
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
  prev_inner_loss <- NA_real_
  normalized_gap_proxy_last <- Inf
  
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
      cfg_omega_admm = cfg$omega_prox_admm_pn,
      Omega_init = if (is.null(prox_state)) Omega else prox_state$Omega,
      A_init = if (is.null(prox_state)) NULL else prox_state$A
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

    inner_loss <- NA_real_
    if (!is.null(n) && is.finite(n) && n > 0) {
      Oinv_zo <- tryCatch(safe_inv_psd(Zo), error = function(e) NULL)
      if (!is.null(Oinv_zo)) {
        inner_loss <- loss_fn(Zg, Zo, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega, Oinv = Oinv_zo)
      }
    }
    progress_proxy <- if (is.finite(prev_inner_loss) && is.finite(inner_loss)) {
      abs(prev_inner_loss - inner_loss) / max(1, abs(prev_inner_loss))
    } else {
      NA_real_
    }
    prev_inner_loss <- inner_loss
    normalized_gap_proxy <- max(r_norm / max(eps_pri, .Machine$double.eps),
                                s_norm / max(eps_dual, .Machine$double.eps))
    normalized_gap_proxy_last <- normalized_gap_proxy

    if (track_trace) {
      inner_trace_list[[k]] <- data.frame(
        inner_iter = k,
        training_loss = inner_loss,
        r_norm = r_norm,
        s_norm = s_norm,
        residual_sum = r_norm + s_norm,
        eps_pri = eps_pri,
        eps_dual = eps_dual,
        primal_dual_gap = normalized_gap_proxy,
        progress_proxy = progress_proxy,
        admm_tol_target = admm_tol,
        rho = rho_iter,
        sub_stop_mode = stop_mode,
        omega_stop_mode = if (!is.null(prox_state$mode)) prox_state$mode else NA_character_,
        omega_stop_metric = if (!is.null(prox_state$final_metric)) prox_state$final_metric else NA_real_,
        stringsAsFactors = FALSE
      )
    }

    residual_ok <- (r_norm <= eps_pri && s_norm <= eps_dual)
    proxy_ok <- is.finite(progress_proxy) && progress_proxy <= admm_tol
    if (identical(stop_mode, "residual_plus_proxy")) {
      if (residual_ok && (is.na(progress_proxy) || proxy_ok)) {
        admm_converged <- TRUE
        break
      }
    } else if (residual_ok) {
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
  out_sub$normalized_gap_proxy <- normalized_gap_proxy_last
  out_sub$subproblem_stop_mode <- stop_mode
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
                            track_diagnostics = FALSE,
                            cfg = SIM_CONFIG,
                            force_max_iter = NULL,
                            ignore_stop = FALSE) {

  pn_cfg <- cfg$outer_pn
  pn_sub_cfg <- cfg$pn_subproblem_admm
  line_search_cfg <- pn_cfg$line_search

  max_iter <- pn_cfg$max_iter
  tol <- pn_cfg$tol_PN_sub
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit
    if (!is.null(ctrl$tol)) tol <- ctrl$tol
  }
  if (!is.null(force_max_iter)) max_iter <- as.integer(force_max_iter)

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
    refinement_trace_list <- vector("list", max_iter * (pn_sub_cfg$max_refinement_rounds + 1L))
    refinement_trace_idx <- 0L
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
      stop_mode_hit = NA_character_,
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
  lambda_hist <- rep(NA_real_, max_iter)
  pn_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  m_last <- 0L

  for (m in seq_len(max_iter)) {
    m_last <- m
    loss_prev <- loss_curr
    beta <- gamma %*% Oinv

    admm_tol_curr <- pn_sub_cfg$tol_PN_ADMM_base
    refinement_triggered <- FALSE
    refinement_rounds <- 0L
    total_inner_iters_m <- 0L
    all_inner_converged_m <- TRUE
    subproblem_attempts <- 0L
    total_ls_halvings <- NA_integer_
    ls_success <- FALSE
    line_search_forced <- FALSE
    alpha <- NA_real_
    ls_halvings <- 0L
    armijo_lhs_used <- NA_real_
    armijo_rhs_used <- NA_real_
    gamma_new <- gamma
    Omega_new <- Omega
    Oinv_new <- Oinv
    loss_new <- loss_prev
    lam <- NA_real_
    d_norm <- NA_real_
    pn_step_failed <- FALSE
    stop_mode_hit <- NA_character_
    inner_trace_m_accum <- list()
    inner_trace_m_idx <- 0L

    max_attempts <- pn_sub_cfg$max_refinement_rounds + 1L
    for (attempt in seq_len(max_attempts)) {
      subproblem_attempts <- attempt
      sub <- pn_subproblem_admm(
        gamma, Omega, Sxx, Sxy, Syy,
        lambda_gamma, lambda_Omega,
        Oinv = Oinv, beta = beta,
        admm_tol = admm_tol_curr,
        state = pn_psd_state,
        n = n,
        track_trace = track_diagnostics,
        cfg = cfg
      )
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

      alpha_try <- alpha_init
      ls_halvings <- 0L
      ls_success <- FALSE
      for (ls_try in 0:line_search_cfg$max_halving) {
        if (alpha_try < line_search_cfg$min_alpha) break
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
        if (is.finite(loss_try) && (loss_try <= loss_prev + line_search_cfg$descent_eps)) {
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

      if (track_diagnostics && !is.null(sub$inner_trace) && nrow(sub$inner_trace) > 0) {
        inner_trace_attempt <- sub$inner_trace
        inner_trace_attempt$outer_iter <- m
        inner_trace_attempt$sub_attempt <- attempt
        inner_trace_attempt$admm_tol_target <- admm_tol_curr
        inner_trace_attempt$direction_norm_fro <- d_norm
        names(inner_trace_attempt)[names(inner_trace_attempt) == "training_loss"] <- "inner_training_loss"
        inner_trace_m_idx <- inner_trace_m_idx + 1L
        inner_trace_m_accum[[inner_trace_m_idx]] <- inner_trace_attempt
      }

      if (track_diagnostics) {
        refinement_trace_idx <- refinement_trace_idx + 1L
        refinement_trace_list[[refinement_trace_idx]] <- data.frame(
          outer_iter = m,
          sub_attempt = attempt,
          admm_tol = admm_tol_curr,
          line_search_halvings = ls_halvings,
          line_search_success = ls_success,
          inner_iters = sub$inner_iters,
          inner_converged = sub$inner_converged,
          stringsAsFactors = FALSE
        )
      }

      if (ls_success) break

      next_tol <- max(admm_tol_curr / pn_sub_cfg$refine_factor, pn_sub_cfg$tol_PN_ADMM_floor)
      can_refine <- (ls_halvings >= pn_sub_cfg$refine_trigger_halvings) &&
        (attempt <= pn_sub_cfg$max_refinement_rounds) &&
        (next_tol < admm_tol_curr - .Machine$double.eps)
      if (can_refine) {
        refinement_triggered <- TRUE
        refinement_rounds <- refinement_rounds + 1L
        admm_tol_curr <- next_tol
        next
      }

      pn_step_failed <- TRUE
      if (!ignore_stop) stop_reason <- "PN step failed to decrease objective"
      break
    }

    if (!ls_success) {
      alpha <- 0
      gamma_new <- gamma
      Omega_new <- Omega
      Oinv_new <- Oinv
      loss_new <- loss_prev
      if (ignore_stop) {
        pn_step_failed <- FALSE
      } else {
        pn_step_failed <- TRUE
      }
    }

    lambda_hist[m] <- lam
    if (track_diagnostics && inner_trace_m_idx > 0L) {
      inner_trace_list[[m]] <- bind_rows_safe(inner_trace_m_accum)
    }

    gamma <- gamma_new
    Omega <- Omega_new
    Oinv <- Oinv_new
    loss_curr <- loss_new

    if (!pn_step_failed && !ignore_stop) {
      stop_eval <- evaluate_pn_stop_modes(
        loss_prev = loss_prev,
        loss_curr = loss_curr,
        local_norm_value = lam,
        oracle_value = oracle_value,
        tol_pn_sub = tol,
        stop_modes = pn_cfg$stop_modes,
        stop_logic = pn_cfg$stop_logic
      )
      if (isTRUE(stop_eval$should_stop)) {
        converged <- TRUE
        stop_mode_hit <- paste(stop_eval$hit_modes, collapse = "|")
        stop_reason <- sprintf("PN stop modes reached: %s", stop_mode_hit)
      } else if (oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE
        stop_reason <- "oracle gap reached"
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
        line_search_halvings = as.integer(total_ls_halvings),
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
        stop_mode_hit = stop_mode_hit,
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

    if (m == 1L || (pn_cfg$progress_every > 0L && (m %% pn_cfg$progress_every == 0L))) {
      cat(sprintf("[PN outer] iter %d/%d | loss=%.6e | rel_change=%.3e | ls_ok=%s | inner=%d | admm_tol=%.1e\n",
                  m, max_iter, loss_curr, rel_obj_change(loss_prev, loss_curr),
                  ifelse(ls_success, "TRUE", "FALSE"),
                  total_inner_iters_m, admm_tol_curr))
      flush.console()
    }

    if (!ignore_stop && pn_step_failed) break
    if (!ignore_stop && converged) break
  }

  elapsed_time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.na(stop_reason)) {
    stop_reason <- if (isTRUE(ignore_stop)) "fixed iteration oracle run completed" else "max_iter reached"
  }
  out <- list(
    gamma = gamma,
    Omega = Omega,
    iters = m_last,
    converged = converged,
    stop_reason = stop_reason,
    elapsed_time_sec = elapsed_time_sec,
    Sxx = Sxx,
    Sxy = Sxy,
    Syy = Syy,
    n = n,
    lambda_hist = lambda_hist[seq_len(max(m_last, 1L))],
    psd_state = pn_psd_state
  )
  if (track_loss) {
    out$loss_hist <- loss_hist[1:(m_last + 1)]
    out$time_hist <- time_hist[1:(m_last + 1)]
  }
  if (track_diagnostics) {
    out$outer_trace <- bind_rows_safe(outer_trace_list[seq_len(m_last + 1)])
    out$inner_trace <- bind_rows_safe(inner_trace_list[seq_len(max(m_last, 1L))])
    out$refinement_trace <- bind_rows_safe(refinement_trace_list[seq_len(max(refinement_trace_idx, 1L))])
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
                        oracle_vec = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                        cfg = SIM_CONFIG) {
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
                           use_relative_gap = use_relative_gap,
                           cfg = cfg)
    
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
                        track_diagnostics = FALSE,
                        cfg = SIM_CONFIG,
                        force_max_iter = NULL,
                        ignore_stop = FALSE) {
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
                           track_diagnostics = track_diagnostics,
                           cfg = cfg,
                           force_max_iter = force_max_iter,
                           ignore_stop = ignore_stop)
    
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

run_pn_oracle_fixed_iters <- function(data_list, cfg = SIM_CONFIG, fixed_iters = NULL,
                                      track_diagnostics = FALSE) {
  fixed_k <- if (is.null(fixed_iters)) cfg$oracle_eval$fixed_iters else as.integer(fixed_iters)
  ctrl_oracle <- list(
    lambda = list(gamma = data_list$lambda_path$lambda_gamma[1], Omega = data_list$lambda_path$lambda_Omega[1]),
    maxit = fixed_k,
    tol = cfg$outer_pn$tol_PN_sub
  )
  pn_oracle_path <- run_pn_path(
    data_list,
    ctrl_oracle,
    track_diagnostics = track_diagnostics,
    cfg = cfg,
    force_max_iter = fixed_k,
    ignore_stop = TRUE
  )

  eps_nondec <- cfg$oracle_eval$nondecrease_eps
  oracle_loss_summary <- bind_rows_safe(lapply(pn_oracle_path, function(fit) {
    loss_hist <- fit$loss_hist
    oracle_loss <- min(loss_hist)
    updates <- length(loss_hist) - 1L
    data.frame(
      path_index = fit$path_index,
      lambda_gamma = fit$lambda_gamma,
      lambda_Omega = fit$lambda_Omega,
      fixed_iters = fixed_k,
      updates_observed = updates,
      updates_expected = fixed_k,
      fixed_iter_ok = as.logical(updates == fixed_k),
      oracle_loss_best100 = oracle_loss,
      oracle_loss_final100 = tail(loss_hist, 1),
      stringsAsFactors = FALSE
    )
  }))

  nondecrease_counts <- bind_rows_safe(lapply(pn_oracle_path, function(fit) {
    loss_hist <- fit$loss_hist
    updates <- length(loss_hist) - 1L
    nd_count <- if (updates > 0) {
      sum(loss_hist[-1] >= loss_hist[-length(loss_hist)] - eps_nondec)
    } else {
      0L
    }
    data.frame(
      path_index = fit$path_index,
      lambda_gamma = fit$lambda_gamma,
      lambda_Omega = fit$lambda_Omega,
      fixed_iters = fixed_k,
      updates_observed = updates,
      nondecrease_count = nd_count,
      nondecrease_rate = if (updates > 0) nd_count / updates else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  aggregate_row <- data.frame(
    path_index = NA_integer_,
    lambda_gamma = NA_real_,
    lambda_Omega = NA_real_,
    fixed_iters = fixed_k,
    updates_observed = sum(nondecrease_counts$updates_observed, na.rm = TRUE),
    nondecrease_count = sum(nondecrease_counts$nondecrease_count, na.rm = TRUE),
    nondecrease_rate = sum(nondecrease_counts$nondecrease_count, na.rm = TRUE) /
      max(1, sum(nondecrease_counts$updates_observed, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  list(
    pn_path = pn_oracle_path,
    oracle_vals = oracle_loss_summary$oracle_loss_best100,
    oracle_loss_summary = oracle_loss_summary,
    nondecrease_counts = nondecrease_counts,
    nondecrease_counts_with_total = bind_rows_safe(list(nondecrease_counts, aggregate_row))
  )
}

load_oracle_from_existing <- function(lambda_path, solver_tol, tol_tag, cfg = SIM_CONFIG) {
  oracle_cfg <- cfg$oracle_eval
  loss_candidates <- character(0)
  if (!is.null(oracle_cfg$reuse_file) && nzchar(oracle_cfg$reuse_file)) {
    loss_candidates <- c(loss_candidates, oracle_cfg$reuse_file)
  }
  loss_candidates <- unique(c(
    loss_candidates,
    sprintf("oracle_pn_100iter_loss_summary_tol_%s.csv", tol_tag),
    "oracle_pn_100iter_loss_summary.csv",
    "oracle_pn_100iter_loss_summary_tol_1e-07.csv"
  ))

  read_loss <- function(path) {
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    req <- c("path_index", "oracle_loss_best100")
    if (!all(req %in% names(df))) return(NULL)
    if ("solver_tol" %in% names(df)) {
      idx <- which(is.finite(df$solver_tol) & abs(df$solver_tol - solver_tol) <= .Machine$double.eps^0.5)
      if (length(idx) > 0) df <- df[idx, , drop = FALSE]
    }
    df <- df[!is.na(df$path_index), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- df[order(df$path_index), , drop = FALSE]
    df <- df[!duplicated(df$path_index), , drop = FALSE]
    if (!all(lambda_path$step %in% df$path_index)) return(NULL)
    df
  }

  loss_df <- NULL
  source_loss <- NA_character_
  for (cand in loss_candidates) {
    df_try <- read_loss(cand)
    if (!is.null(df_try)) {
      loss_df <- df_try
      source_loss <- cand
      break
    }
  }
  if (is.null(loss_df)) {
    return(list(
      oracle_vals = NULL,
      oracle_loss_df = data.frame(),
      oracle_nondec_df = data.frame(),
      source_file = NA_character_
    ))
  }

  nondec_candidates <- unique(c(
    if (!is.na(source_loss)) sub("loss_summary", "nondecrease_counts", source_loss) else character(0),
    sprintf("oracle_pn_100iter_nondecrease_counts_tol_%s.csv", tol_tag),
    "oracle_pn_100iter_nondecrease_counts.csv",
    "oracle_pn_100iter_nondecrease_counts_tol_1e-07.csv"
  ))
  nondec_df <- data.frame()
  for (cand in nondec_candidates) {
    if (!file.exists(cand)) next
    nd_try <- tryCatch(read.csv(cand, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(nd_try)) next
    if ("solver_tol" %in% names(nd_try)) {
      idx <- which(is.finite(nd_try$solver_tol) & abs(nd_try$solver_tol - solver_tol) <= .Machine$double.eps^0.5)
      if (length(idx) > 0) nd_try <- nd_try[idx, , drop = FALSE]
    }
    nondec_df <- nd_try
    break
  }

  oracle_vals <- loss_df$oracle_loss_best100[match(lambda_path$step, loss_df$path_index)]
  if (any(!is.finite(oracle_vals))) {
    return(list(
      oracle_vals = NULL,
      oracle_loss_df = data.frame(),
      oracle_nondec_df = data.frame(),
      source_file = NA_character_
    ))
  }
  list(
    oracle_vals = oracle_vals,
    oracle_loss_df = loss_df,
    oracle_nondec_df = nondec_df,
    source_file = source_loss
  )
}

compute_oracle_path <- function(data_list, cfg = SIM_CONFIG, ctrl_oracle = NULL) {
  if (is.null(ctrl_oracle)) {
    return(run_pn_oracle_fixed_iters(data_list, cfg = cfg, fixed_iters = cfg$oracle_eval$fixed_iters))
  }
  pn_oracle <- run_pn_path(data_list, ctrl_oracle, cfg = cfg)
  oracle_vals <- vapply(pn_oracle, function(fit) min(fit$loss_hist), numeric(1))
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
                               track_pn_diagnostics = FALSE,
                               cfg = SIM_CONFIG,
                               oracle_vals_override = NULL,
                               oracle_path_override = NULL,
                               warm_starts_override = NULL) {
  
  oracle_vals <- NULL
  pn_oracle_path <- NULL
  warm_starts <- NULL
  if (!is.null(oracle_vals_override)) {
    oracle_vals <- oracle_vals_override
    pn_oracle_path <- oracle_path_override
    warm_starts <- warm_starts_override
  } else if (use_oracle) {
    oracle_out <- compute_oracle_path(data_list, cfg = cfg, ctrl_oracle = ctrl_oracle)
    pn_oracle_path <- oracle_out$pn_path
    oracle_vals <- oracle_out$oracle_vals
    warm_starts <- build_shared_warm_starts(data_list, pn_oracle_path)
  }
  
  pn_path <- run_pn_path(data_list, ctrl, warm_starts = warm_starts,
                         oracle_vec = oracle_vals, gap_tol = gap_tol,
                         use_relative_gap = use_relative_gap,
                         track_diagnostics = track_pn_diagnostics,
                         cfg = cfg)
  pg_path <- run_pg_path(data_list, ctrl, warm_starts = warm_starts,
                         oracle_vec = oracle_vals, gap_tol = gap_tol,
                         use_relative_gap = use_relative_gap,
                         cfg = cfg)
  
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

build_pn_refinement_trace <- function(pn_path, solver_tol) {
  bind_rows_safe(lapply(pn_path, function(fit) {
    if (is.null(fit$refinement_trace) || nrow(fit$refinement_trace) == 0) return(NULL)
    df <- fit$refinement_trace
    df$path_index <- fit$path_index
    df$lambda_gamma <- fit$lambda_gamma
    df$lambda_Omega <- fit$lambda_Omega
    df$solver_tol <- solver_tol
    df
  }))
}

build_pn_omega_mode_summary <- function(pn_path, solver_tol) {
  bind_rows_safe(lapply(pn_path, function(fit) {
    if (is.null(fit$inner_trace) || nrow(fit$inner_trace) == 0) return(NULL)
    tab <- aggregate(
      cbind(omega_stop_metric, inner_iter) ~ omega_stop_mode + sub_stop_mode,
      data = fit$inner_trace,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    tab$path_index <- fit$path_index
    tab$lambda_gamma <- fit$lambda_gamma
    tab$lambda_Omega <- fit$lambda_Omega
    tab$solver_tol <- solver_tol
    tab
  }))
}

build_pg_omega_mode_summary <- function(pg_path, solver_tol) {
  bind_rows_safe(lapply(pg_path, function(fit) {
    st <- fit$psd_state
    if (is.null(st)) return(NULL)
    data.frame(
      path_index = fit$path_index,
      lambda_gamma = fit$lambda_gamma,
      lambda_Omega = fit$lambda_Omega,
      omega_stop_mode = if (!is.null(st$mode)) st$mode else NA_character_,
      omega_stop_metric = if (!is.null(st$final_metric)) st$final_metric else NA_real_,
      omega_inner_iters = if (!is.null(st$iters)) st$iters else NA_real_,
      omega_inner_converged = if (!is.null(st$converged)) st$converged else NA,
      solver_tol = solver_tol,
      stringsAsFactors = FALSE
    )
  }))
}

build_run_config_snapshot <- function(cfg, solver_tol,
                                      gap_tol, use_oracle, use_relative_gap) {
  sim <- cfg$simulation
  lam <- cfg$lambda_path
  pg <- cfg$outer_pg
  pn <- cfg$outer_pn
  sub <- cfg$pn_subproblem_admm
  pg_omega <- cfg$omega_prox_admm_pg
  pn_omega <- cfg$omega_prox_admm_pn
  num <- cfg$numerics
  data.frame(
    seed = shared_seed,
    n = sim$n,
    p = sim$p,
    q = sim$q,
    path_len = lam$path_len,
    ratio_lambda = lam$ratio_lambda,
    gamma_sparsity = sim$gamma_sparsity,
    gamma_mag = paste(sim$gamma_mag, collapse = ","),
    gamma_mag_lo = min(sim$gamma_mag),
    gamma_mag_hi = max(sim$gamma_mag),
    omega_kappa = sim$omega_kappa,
    omega_sparsity = sim$omega_sparsity,
    solver_tol = solver_tol,
    max_iter_PG_outer = pg$max_iter,
    min_iter_PG_outer = pg$min_iter,
    max_iter_PN_outer = pn$max_iter,
    tol_PG_sub = pg$tol_PG_sub,
    tol_PN_sub = pn$tol_PN_sub,
    tol_PN_ADMM_base = sub$tol_PN_ADMM_base,
    tol_PN_ADMM_floor = sub$tol_PN_ADMM_floor,
    max_iter_PN_sub = sub$max_iter,
    max_iter_Omega_ADMM_PG = pg_omega$max_iter,
    max_iter_Omega_ADMM_PN = pn_omega$max_iter,
    tol_Omega_ADMM_PG = pg_omega$tol_Omega_ADMM,
    tol_Omega_ADMM_PN = pn_omega$tol_Omega_ADMM,
    stop_mode_Omega_ADMM_PG = pg_omega$stop_mode,
    stop_mode_Omega_ADMM_PN = pn_omega$stop_mode,
    stop_mode_PN_sub = sub$stop_mode,
    eig_floor_Omega_PG = pg_omega$eig_floor,
    eig_floor_Omega_PN = pn_omega$eig_floor,
    eig_floor_PN_subproblem = num$pn_subproblem_eig_floor,
    safe_inv_eps = num$safe_inv_eps,
    safe_inv_jitter_init = num$safe_inv_jitter_init,
    safe_inv_jitter_max = num$safe_inv_jitter_max,
    safe_logdet_eps = num$safe_logdet_eps,
    safe_logdet_jitter_init = num$safe_logdet_jitter_init,
    safe_logdet_jitter_max = num$safe_logdet_jitter_max,
    matrix_condition_eps = num$matrix_condition_eps,
    null_corner_eps = num$null_corner_eps,
    gap_tol = gap_tol,
    use_oracle = use_oracle,
    use_relative_gap = use_relative_gap,
    pn_ls_max_halving = pn$line_search$max_halving,
    pn_ls_min_alpha = pn$line_search$min_alpha,
    pn_ls_descent_eps = pn$line_search$descent_eps,
    pn_admm_rho = sub$adaptive_rho$rho,
    pn_admm_adaptive = sub$adaptive_rho$enabled,
    pn_admm_rho_mu = sub$adaptive_rho$rho_mu,
    pn_admm_rho_tau_inc = sub$adaptive_rho$rho_tau_inc,
    pn_admm_rho_tau_dec = sub$adaptive_rho$rho_tau_dec,
    pn_refine_trigger_halvings = sub$refine_trigger_halvings,
    pn_refine_factor = sub$refine_factor,
    pn_refine_max_rounds = sub$max_refinement_rounds,
    oracle_fixed_iters = cfg$oracle_eval$fixed_iters,
    oracle_nondecrease_eps = cfg$oracle_eval$nondecrease_eps,
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

solver_tols <- SIM_CONFIG$reporting$solver_tols
settings_table <- data.frame(
  solver_tol = solver_tols,
  max_iter_PG_outer = SIM_CONFIG$outer_pg$max_iter,
  min_iter_PG_outer = SIM_CONFIG$outer_pg$min_iter,
  max_iter_PN_outer = SIM_CONFIG$outer_pn$max_iter,
  max_iter_PN_sub = SIM_CONFIG$pn_subproblem_admm$max_iter,
  tol_PN_ADMM_base = SIM_CONFIG$pn_subproblem_admm$tol_PN_ADMM_base,
  tol_PN_ADMM_floor = SIM_CONFIG$pn_subproblem_admm$tol_PN_ADMM_floor,
  tol_Omega_ADMM_PG = SIM_CONFIG$omega_prox_admm_pg$tol_Omega_ADMM,
  tol_Omega_ADMM_PN = SIM_CONFIG$omega_prox_admm_pn$tol_Omega_ADMM,
  n = shared_params$n,
  p = shared_params$p,
  q = shared_params$q,
  path_len = path_len,
  stringsAsFactors = FALSE
)
print(settings_table)

stop_mode_cases <- list(
  case_relative_change = list(
    label = "relative_change",
    pn_stop_modes = c("relative_change"),
    pg_stop_modes = c("relative_change")
  ),
  case_oracle_change = list(
    label = "oracle_change",
    pn_stop_modes = c("oracle_change"),
    pg_stop_modes = c("oracle_change")
  ),
  case_norm_mode = list(
    label = "norm_mode",
    pn_stop_modes = c("local_norm"),
    pg_stop_modes = c("l2_step_norm")
  )
)

ensure_data_frame <- function(df) {
  if (is.null(df)) data.frame(stringsAsFactors = FALSE) else df
}

add_case_column <- function(df, case_name) {
  df <- ensure_data_frame(df)
  if (nrow(df) == 0) {
    df$case <- character(0)
  } else {
    df$case <- case_name
  }
  df
}

method_levels <- SIM_CONFIG$reporting$method_levels
all_results <- list()
time_tables <- list()
lambda_summary_tables <- list()
pn_outer_trace_tables <- list()
pn_inner_trace_tables <- list()
pn_outer_inner_map_tables <- list()
pn_refinement_tables <- list()
pn_omega_mode_tables <- list()
pg_omega_mode_tables <- list()
config_snapshot_tables <- list()
oracle_loss_tables <- list()
oracle_nondec_tables <- list()

for (solver_tol in solver_tols) {
  cat(sprintf("\n=== Running solver_tol = %.1e ===\n", solver_tol))
  flush.console()
  tol_tag <- format_tol_tag(solver_tol)

  cfg_tol <- SIM_CONFIG
  cfg_tol$outer_pg$tol_PG_sub <- solver_tol
  cfg_tol$outer_pg$tol_obj <- solver_tol
  cfg_tol$outer_pg$tol_gm <- solver_tol
  cfg_tol$outer_pg$l2_norm_tol <- solver_tol
  cfg_tol$outer_pn$tol_PN_sub <- solver_tol
  set_numeric_controls(cfg_tol)

  oracle_vals_shared <- NULL
  oracle_loss_df <- data.frame()
  oracle_nondec_df <- data.frame()
  if (isTRUE(cfg_tol$oracle_eval$reuse_existing)) {
    oracle_loaded <- load_oracle_from_existing(
      lambda_path = shared_data$lambda_path,
      solver_tol = solver_tol,
      tol_tag = tol_tag,
      cfg = cfg_tol
    )
    oracle_vals_shared <- oracle_loaded$oracle_vals
    oracle_loss_df <- oracle_loaded$oracle_loss_df
    oracle_nondec_df <- oracle_loaded$oracle_nondec_df
    if (is.null(oracle_vals_shared)) {
      stop("Oracle reuse requested, but no compatible oracle file was found.")
    }
    cat(sprintf("Reusing oracle from existing file: %s\n", oracle_loaded$source_file))
    flush.console()
  } else if (isTRUE(cfg_tol$oracle_eval$run_for_exports)) {
    oracle_run <- run_pn_oracle_fixed_iters(
      shared_data,
      cfg = cfg_tol,
      fixed_iters = cfg_tol$oracle_eval$fixed_iters,
      track_diagnostics = FALSE
    )
    oracle_vals_shared <- oracle_run$oracle_vals
    oracle_loss_df <- oracle_run$oracle_loss_summary
    oracle_nondec_df <- oracle_run$nondecrease_counts_with_total
  } else {
    stop("No oracle source available. Enable ORACLE_REUSE_EXISTING=1 or ORACLE_RUN_FOR_EXPORTS=1.")
  }

  if (nrow(oracle_loss_df) > 0) {
    oracle_loss_df$solver_tol <- solver_tol
    oracle_loss_df$oracle_case <- "case_oracle_change"
  }
  if (nrow(oracle_nondec_df) > 0) {
    oracle_nondec_df$solver_tol <- solver_tol
    oracle_nondec_df$oracle_case <- "case_oracle_change"
  }
  oracle_loss_tables[[as.character(solver_tol)]] <- oracle_loss_df
  oracle_nondec_tables[[as.character(solver_tol)]] <- oracle_nondec_df
  write.csv(oracle_loss_df, sprintf("oracle_pn_100iter_loss_summary_tol_%s.csv", tol_tag), row.names = FALSE)
  write.csv(oracle_nondec_df, sprintf("oracle_pn_100iter_nondecrease_counts_tol_%s.csv", tol_tag), row.names = FALSE)

  case_time_rows <- list()
  for (case_name in names(stop_mode_cases)) {
    case_spec <- stop_mode_cases[[case_name]]
    cat(sprintf("  -> case: %s | PN=%s | PG=%s\n",
                case_name,
                paste(case_spec$pn_stop_modes, collapse = "|"),
                paste(case_spec$pg_stop_modes, collapse = "|")))
    flush.console()

    cfg_case <- cfg_tol
    cfg_case$outer_pn$stop_modes <- case_spec$pn_stop_modes
    cfg_case$outer_pn$stop_logic <- "any"
    cfg_case$outer_pg$stop_modes <- case_spec$pg_stop_modes
    cfg_case$outer_pg$stop_logic <- "any"
    cfg_case$outer_pg$l2_norm_tol <- solver_tol
    set_numeric_controls(cfg_case)

    ctrl_case <- list(
      lambda = list(gamma = lambda_gamma_seq[1], Omega = lambda_Omega_seq[1]),
      maxit = cfg_case$outer_pn$max_iter,
      tol = solver_tol,
      L0 = cfg_case$outer_pg$L0,
      stop_modes = cfg_case$outer_pg$stop_modes,
      stop_logic = cfg_case$outer_pg$stop_logic,
      l2_norm_tol = cfg_case$outer_pg$l2_norm_tol
    )

    oracle_vec_case <- if (identical(case_name, "case_oracle_change")) oracle_vals_shared else NULL

    res <- benchmark_pg_vs_pn(
      shared_data, ctrl_case,
      use_oracle = FALSE,
      gap_tol = cfg_case$benchmark$gap_tol_mode,
      use_relative_gap = cfg_case$benchmark$use_relative_gap_mode,
      ctrl_oracle = NULL,
      track_pn_diagnostics = TRUE,
      cfg = cfg_case,
      oracle_vals_override = oracle_vec_case
    )
    res$path_summary$solver_tol <- solver_tol
    res$path_summary$case <- case_name
    all_results[[paste0(tol_tag, "_", case_name)]] <- res

    lambda_summary_df <- build_lambda_summary(res, solver_tol)
    pn_outer_trace_df <- build_pn_outer_trace(res$pn_path, solver_tol)
    pn_inner_trace_df <- build_pn_inner_trace(res$pn_path, solver_tol)
    pn_outer_inner_map_df <- build_pn_outer_inner_map(res$pn_path, solver_tol)
    pn_refinement_df <- build_pn_refinement_trace(res$pn_path, solver_tol)
    pn_omega_mode_df <- build_pn_omega_mode_summary(res$pn_path, solver_tol)
    pg_omega_mode_df <- build_pg_omega_mode_summary(res$pg_path, solver_tol)
    config_snapshot_df <- build_run_config_snapshot(
      cfg = cfg_case,
      solver_tol = solver_tol,
      gap_tol = if (is.null(cfg_case$benchmark$gap_tol_mode)) NA_real_ else cfg_case$benchmark$gap_tol_mode,
      use_oracle = identical(case_name, "case_oracle_change"),
      use_relative_gap = cfg_case$benchmark$use_relative_gap_mode
    )

    lambda_summary_df <- add_case_column(lambda_summary_df, case_name)
    pn_outer_trace_df <- add_case_column(pn_outer_trace_df, case_name)
    pn_inner_trace_df <- add_case_column(pn_inner_trace_df, case_name)
    pn_outer_inner_map_df <- add_case_column(pn_outer_inner_map_df, case_name)
    pn_refinement_df <- add_case_column(pn_refinement_df, case_name)
    pn_omega_mode_df <- add_case_column(pn_omega_mode_df, case_name)
    pg_omega_mode_df <- add_case_column(pg_omega_mode_df, case_name)
    config_snapshot_df <- add_case_column(config_snapshot_df, case_name)

    result_key <- paste0(tol_tag, "_", case_name)
    lambda_summary_tables[[result_key]] <- lambda_summary_df
    pn_outer_trace_tables[[result_key]] <- pn_outer_trace_df
    pn_inner_trace_tables[[result_key]] <- pn_inner_trace_df
    pn_outer_inner_map_tables[[result_key]] <- pn_outer_inner_map_df
    pn_refinement_tables[[result_key]] <- pn_refinement_df
    pn_omega_mode_tables[[result_key]] <- pn_omega_mode_df
    pg_omega_mode_tables[[result_key]] <- pg_omega_mode_df
    config_snapshot_tables[[result_key]] <- config_snapshot_df

    write.csv(lambda_summary_df, sprintf("lambda_summary_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(pn_outer_trace_df, sprintf("pn_outer_trace_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(pn_inner_trace_df, sprintf("pn_inner_trace_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(pn_outer_inner_map_df, sprintf("pn_outer_inner_map_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(pn_refinement_df, sprintf("pn_adaptive_tol_refinement_trace_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(pn_omega_mode_df, sprintf("pn_omega_admm_mode_summary_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(pg_omega_mode_df, sprintf("pg_omega_admm_mode_summary_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    write.csv(config_snapshot_df, sprintf("run_config_snapshot_%s_tol_%s.csv", case_name, tol_tag), row.names = FALSE)
    saveRDS(list(
      lambda_summary = lambda_summary_df,
      pn_outer_trace = pn_outer_trace_df,
      pn_inner_trace = pn_inner_trace_df,
      pn_outer_inner_map = pn_outer_inner_map_df,
      pn_adaptive_tol_refinement_trace = pn_refinement_df,
      pn_omega_admm_mode_summary = pn_omega_mode_df,
      pg_omega_admm_mode_summary = pg_omega_mode_df,
      run_config_snapshot = config_snapshot_df
    ), file = sprintf("pn_diagnostics_%s_tol_%s.rds", case_name, tol_tag))

    loss_plot_df <- res$path_summary
    loss_plot_df$method <- factor(loss_plot_df$method, levels = method_levels)
    loss_plot_df$log_lambda <- log(loss_plot_df$lambda_gamma)

    total_time_by_method <- aggregate(elapsed_sec ~ method, data = res$path_summary, sum)
    total_time_by_method$solver_tol <- solver_tol
    total_time_by_method$case <- case_name
    total_time_by_method$tol_outer <- solver_tol
    total_time_by_method$tol_inner_base <- cfg_case$pn_subproblem_admm$tol_PN_ADMM_base
    total_time_by_method$tol_inner_floor <- cfg_case$pn_subproblem_admm$tol_PN_ADMM_floor
    case_time_rows[[case_name]] <- total_time_by_method
    time_tables[[paste0(tol_tag, "_", case_name)]] <- total_time_by_method

    pn_time <- total_time_by_method$elapsed_sec[total_time_by_method$method == "Prox-Newton"]
    pg_time <- total_time_by_method$elapsed_sec[total_time_by_method$method == "Prox-Gradient"]
    if (length(pn_time) == 0) pn_time <- NA_real_
    if (length(pg_time) == 0) pg_time <- NA_real_
    time_subtitle <- sprintf("Total path time (sec): PN=%.3f, PG=%.3f", pn_time[1], pg_time[1])

    final_loss_plot <- ggplot(loss_plot_df,
                              aes(x = log_lambda, y = final_loss,
                                  colour = method, shape = method, group = method)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.2, stroke = 0.9) +
      scale_shape_manual(values = c("Prox-Newton" = 4, "Prox-Gradient" = 16)) +
      labs(
        title = sprintf("Training loss vs log(lambda): %s", case_name),
        subtitle = time_subtitle,
        x = expression(log(lambda)),
        y = "Training Loss",
        colour = NULL,
        shape = NULL
      ) +
      theme_minimal(base_size = 13)

    iter_plot <- ggplot(loss_plot_df,
                        aes(x = log_lambda, y = iters,
                            colour = method, shape = method, group = method)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.2, stroke = 0.9) +
      scale_shape_manual(values = c("Prox-Newton" = 4, "Prox-Gradient" = 16)) +
      labs(
        title = sprintf("Iterations vs log(lambda): %s", case_name),
        subtitle = time_subtitle,
        x = expression(log(lambda)),
        y = "Number of iterations",
        colour = NULL,
        shape = NULL
      ) +
      theme_minimal(base_size = 13)

    ggsave(
      sprintf("pn_pg_final_loss_vs_loglambda_%s_tol_%s.png", case_name, tol_tag),
      final_loss_plot,
      width = 9,
      height = 4.5,
      dpi = 300
    )
    ggsave(
      sprintf("pn_pg_iters_vs_loglambda_%s_tol_%s.png", case_name, tol_tag),
      iter_plot,
      width = 9,
      height = 4.5,
      dpi = 300
    )

    print(final_loss_plot)
    print(iter_plot)

    saveRDS(res, file = sprintf("res_%s_tol_%s.rds", case_name, tol_tag))
    saveRDS(list(
      pn_path = res$pn_path,
      pg_path = res$pg_path,
      lambda_path = shared_data$lambda_path
    ), file = sprintf("lambda_path_estimators_%s_tol_%s.rds", case_name, tol_tag))
  }

  case_time_df <- bind_rows_safe(case_time_rows)
  if (!is.null(case_time_df) && nrow(case_time_df) > 0) {
    case_time_export <- case_time_df[, c("case", "method", "elapsed_sec", "tol_outer", "tol_inner_base", "tol_inner_floor"), drop = FALSE]
    names(case_time_export) <- c("case", "method", "total_elapsed_sec", "tol_outer", "tol_inner_base", "tol_inner_floor")
    write.csv(
      case_time_export,
      sprintf("solution_path_time_by_stop_mode_tol_%s.csv", tol_tag),
      row.names = FALSE
    )
    print(case_time_export)
  }
}

time_summary <- bind_rows_safe(time_tables)
if (is.null(time_summary)) time_summary <- data.frame()
if (nrow(time_summary) > 0) {
  time_summary <- time_summary[, c("solver_tol", "case", "method", "elapsed_sec", "tol_outer", "tol_inner_base", "tol_inner_floor"), drop = FALSE]
  names(time_summary) <- c("solver_tol", "case", "method", "total_elapsed_sec", "tol_outer", "tol_inner_base", "tol_inner_floor")
}
print(time_summary)
write.csv(time_summary, "time_method_summary.csv", row.names = FALSE)

lambda_summary_all <- bind_rows_safe(lambda_summary_tables)
pn_outer_trace_all <- bind_rows_safe(pn_outer_trace_tables)
pn_inner_trace_all <- bind_rows_safe(pn_inner_trace_tables)
pn_outer_inner_map_all <- bind_rows_safe(pn_outer_inner_map_tables)
pn_adaptive_tol_refinement_trace_all <- bind_rows_safe(pn_refinement_tables)
pn_omega_admm_mode_summary_all <- bind_rows_safe(pn_omega_mode_tables)
pg_omega_admm_mode_summary_all <- bind_rows_safe(pg_omega_mode_tables)
oracle_pn_100iter_loss_summary_all <- bind_rows_safe(oracle_loss_tables)
oracle_pn_100iter_nondecrease_counts_all <- bind_rows_safe(oracle_nondec_tables)
run_config_snapshot_all <- bind_rows_safe(config_snapshot_tables)

write.csv(lambda_summary_all, "lambda_summary_all.csv", row.names = FALSE)
write.csv(pn_outer_trace_all, "pn_outer_trace_all.csv", row.names = FALSE)
write.csv(pn_inner_trace_all, "pn_inner_trace_all.csv", row.names = FALSE)
write.csv(pn_outer_inner_map_all, "pn_outer_inner_map_all.csv", row.names = FALSE)
write.csv(pn_adaptive_tol_refinement_trace_all, "pn_adaptive_tol_refinement_trace.csv", row.names = FALSE)
write.csv(pn_omega_admm_mode_summary_all, "pn_omega_admm_mode_summary.csv", row.names = FALSE)
write.csv(pg_omega_admm_mode_summary_all, "pg_omega_admm_mode_summary.csv", row.names = FALSE)
write.csv(oracle_pn_100iter_loss_summary_all, "oracle_pn_100iter_loss_summary.csv", row.names = FALSE)
write.csv(oracle_pn_100iter_nondecrease_counts_all, "oracle_pn_100iter_nondecrease_counts.csv", row.names = FALSE)
write.csv(run_config_snapshot_all, "run_config_snapshot_all.csv", row.names = FALSE)

saveRDS(list(
  lambda_summary_all = lambda_summary_all,
  pn_outer_trace_all = pn_outer_trace_all,
  pn_inner_trace_all = pn_inner_trace_all,
  pn_outer_inner_map_all = pn_outer_inner_map_all,
  pn_adaptive_tol_refinement_trace_all = pn_adaptive_tol_refinement_trace_all,
  pn_omega_admm_mode_summary_all = pn_omega_admm_mode_summary_all,
  pg_omega_admm_mode_summary_all = pg_omega_admm_mode_summary_all,
  oracle_pn_100iter_loss_summary_all = oracle_pn_100iter_loss_summary_all,
  oracle_pn_100iter_nondecrease_counts_all = oracle_pn_100iter_nondecrease_counts_all,
  run_config_snapshot_all = run_config_snapshot_all
), file = "pn_diagnostics_all.rds")
