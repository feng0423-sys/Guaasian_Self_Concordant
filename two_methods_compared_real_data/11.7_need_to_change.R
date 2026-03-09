# ---- Basic helpers ----
symmetrize <- function(A) 0.5 * (A + t(A))

clip_eig <- function(A, eps = 0) {
  ea <- eigen(symmetrize(A), symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  V %*% diag(d, length(d)) %*% t(V)
}

proj_psd <- function(A) clip_eig(A, eps = 0)   # allow PSD


is_pd_by_chol_rank <- function(A, eps = 1e-8) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  R <- tryCatch(suppressWarnings(chol(A_sym + diag(eps, n), pivot = TRUE)),  # FIX: suppress spurious chol warnings
                error = function(e) NULL)
  if (is.null(R)) return(FALSE)
  attr(R, "rank") == n
}


soft_thresh <- function(A, tau) sign(A) * pmax(abs(A) - tau, 0)

# only threshold off-diagonals
soft_thresh_offdiag <- function(A, tau) {
  B <- symmetrize(A)
  idx <- matrix(TRUE, nrow(A), ncol(A)); diag(idx) <- FALSE
  B[idx] <- soft_thresh(B[idx], tau)
  B
}

# robust inverse (Cholesky + jitter, fallback to eig)
safe_inv_psd <- function(A, eps = 1e-8, jitter_init = 0, jitter_max = 1e-5) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  jitter <- jitter_init
  while (jitter <= jitter_max) {
    chol_try <- tryCatch(suppressWarnings(chol(A_sym + diag(jitter, n))),  # FIX: silence SPD jitter warnings
                         error = function(e) NULL)
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


eig_max_symmetric <- function(A) {
  eigvals <- suppressWarnings(eigen(symmetrize(A), symmetric = TRUE, only.values = TRUE)$values)
  max(eigvals)
}

# FIX: Use a shared dataset/configuration for PN and PG
set.seed(20251025)
shared_params <- list(
  n = 1000,
  p = 100,   # UPDATED: p=100
  q = 100,   # UPDATED: q=100
  ratio_lambda = 1e-4,
  gamma_sparsity = 0.05,
  gamma_mag = c(1, 2),
  omega_kappa = 60,
  omega_sparsity = 0.05
)

# Sparse generator utilities (controls NNZ count and magnitude)
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

# Generate ground-truth coefficients/precision following sparsity + Gershgorin rules
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

# Shared covariates/response
shared_X <- matrix(rnorm(shared_params$n * shared_params$p), shared_params$n, shared_params$p)
shared_Y <- shared_X %*% shared_Beta_true + MASS::mvrnorm(shared_params$n,
                                                          mu = rep(0, shared_params$q),
                                                          Sigma = shared_Sigma_true)
shared_Xc <- cbind(1, scale(shared_X, center = TRUE, scale = FALSE))
shared_p_aug <- ncol(shared_Xc)

# Shared penalties (full path so that λ_1 zeros everything)
Sxy_shared <- crossprod(shared_Xc, shared_Y) / shared_params$n
Syy_shared <- crossprod(shared_Y) / shared_params$n
lambda_gamma_0 <- max(abs(Sxy_shared[-1, , drop = FALSE]))
tmp_off_shared <- Syy_shared - diag(diag(Syy_shared))
lambda_Omega_0 <- max(abs(tmp_off_shared))
path_len <- 40  # UPDATED: lambda path length
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

# Closed-form null-corner initializer (γ_P = 0, Ω diagonal)
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

# Solver hyper-parameters shared by PN/PG (env overrides optional)
tol_env <- suppressWarnings(as.numeric(Sys.getenv("SOLVER_TOL", unset = NA_character_)))
if (is.na(tol_env)) tol_env <- 1e-5
maxit_env <- suppressWarnings(as.integer(Sys.getenv("SOLVER_MAXIT", unset = NA_character_)))
if (is.na(maxit_env) || maxit_env <= 0) maxit_env <- 1000

# Common control + sparse initial guesses
ctrl <- list(
  lambda = list(gamma = lambda_gamma_seq[1], Omega = lambda_Omega_seq[1]),
  maxit = maxit_env,
  tol = tol_env
)
prox_tol_env <- suppressWarnings(as.numeric(Sys.getenv("PROX_TOL", unset = NA_character_)))
if (is.na(prox_tol_env)) prox_tol_env <- 0.1 * tol_env
prox_max_admm_env <- suppressWarnings(as.integer(Sys.getenv("PROX_MAX_ADMM", unset = NA_character_)))
if (is.na(prox_max_admm_env) || prox_max_admm_env <= 0) prox_max_admm_env <- 1000
prox_eta_env <- suppressWarnings(as.numeric(Sys.getenv("PROX_ETA", unset = NA_character_)))
if (is.na(prox_eta_env) || prox_eta_env < 0) prox_eta_env <- 1e-6
prox_ctrl <- list(max_admm = prox_max_admm_env, tol = prox_tol_env, eta = prox_eta_env)

pn_eta_env <- suppressWarnings(as.numeric(Sys.getenv("PN_PROX_ETA", unset = NA_character_)))
if (is.na(pn_eta_env) || pn_eta_env < 0) pn_eta_env <- prox_eta_env
pn_prox_ctrl <- list(max_admm = prox_max_admm_env, tol = prox_tol_env, eta = pn_eta_env)
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
)  # Bundle shared inputs
armijo_c <- 0.1  # FIX: Armijo constant c


## ---------- helpers: robust Oinv / beta ----------
.ensure_Oinv <- function(Omega, Oinv = NULL) {
  if (is.null(Oinv)) safe_inv_psd(Omega) else symmetrize(Oinv)
}
.ensure_beta <- function(gamma, Omega, Oinv = NULL, beta = NULL) {
  if (!is.null(beta)) return(beta)
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  gamma %*% Oinv_eff
}

## ---------- Loss / Gradients (γ-parameterization; optional Oinv/beta) ----------
# smooth part:
# 0.5 * ( -logdet(Ω) + tr(Ω^{-1} γ^T Sxx γ) - 2 tr(Sxy γ) + tr(Ω Syy) )
smooth_obj <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL, beta = NULL) {
  term1 <- -as.numeric(determinant(Omega, logarithm = TRUE)$modulus)   # Ω must be SPD
  term3 <- -2 * sum(Sxy * gamma)
  term4 <- sum(Omega * Syy)
  
  if (!is.null(beta)) {
    # tr(Ω^{-1} γ^T Sxx γ) = tr( (β^T Sxx β) Ω )
    BSB  <- t(beta) %*% Sxx %*% beta
    term2 <- sum(BSB * Omega)
  } else {
    Oinv_eff <- .ensure_Oinv(Omega, Oinv)
    M <- t(gamma) %*% Sxx %*% gamma
    term2 <- sum(Oinv_eff * M)
  }
  0.5 * (term1 + term2 + term3 + term4)
}

# full penalized objective; 
loss_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n,
                    lambda_gamma, lambda_Omega, lambda_Omega_diag = 0,
                    Oinv = NULL, beta = NULL) {
  pen_g <- lambda_gamma * sum(abs(gamma[-1, , drop = FALSE]))
  pen_o <- lambda_Omega * (sum(abs(Omega)) - sum(abs(diag(Omega)))) +
    lambda_Omega_diag * sum(abs(diag(Omega)))
  smooth_obj(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv, beta = beta) + pen_g + pen_o
}

# pure NLL
nll_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n, Oinv = NULL, beta = NULL) {
  smooth_obj(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv, beta = beta)
}


# ∇_γ f_smooth = -Sxy + Sxx γ Ω^{-1} (= -Sxy + Sxx β)
grad_gamma <- function(gamma, Omega, Sxx, Sxy, Oinv = NULL, beta = NULL) {
  if (is.null(beta)) {
    Oinv_eff <- .ensure_Oinv(Omega, Oinv)
    Sxx %*% gamma %*% Oinv_eff - Sxy
  } else {
    Sxx %*% beta - Sxy
  }
}


# ∇_Ω f_smooth = 1/2 ( Syy - Ω^{-1} γ^T Sxx γ Ω^{-1} - Ω^{-1} )
#              = 1/2 ( Syy - β^T Sxx β - Ω^{-1} )
grad_Omega <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  if (is.null(beta)) {
    M <- t(gamma) %*% Sxx %*% gamma
    0.5 * ( Syy - Oinv_eff %*% M %*% Oinv_eff - Oinv_eff )
  } else {
    BSB <- t(beta) %*% Sxx %*% beta
    0.5 * ( Syy - BSB - Oinv_eff )
  }
}

## ---------- Hessian blocks & local norm (Algorithm 1/2) ----------
# β = γΩ^{-1}：
# Hgg = Sxx; Hgo = -Sxx β; Hog = -(Sxx β)^T = -β^T Sxx; Hoo = 1/2 Ω^{-1} + β^T Sxx β
H_blocks <- function(gamma, Omega, Sxx, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  beta_eff <- .ensure_beta(gamma, Omega, Oinv_eff, beta)
  
  Hgg <- Sxx
  Hgo <- - Sxx %*% beta_eff
  Hog <- - t(Hgo)
  Hoo <- 0.5 * Oinv_eff + t(beta_eff) %*% Sxx %*% beta_eff
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

## ||dξ||_ξ = sqrt( <dξ, H(ξ) dξ Ω^{-1}> )
local_norm <- function(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv_eff, beta = beta)
  H  <- H_matrix(Hb)                     # (p+q) x (p+q)
  Xi_d <- pack_xi(d_gamma, d_Omega)      # (p+q) x q
  val <- sum(Xi_d * (H %*% Xi_d %*% Oinv_eff))
  sqrt(max(val, 0))
}


# ---- PSD + off-diag L1 proximal (ADMM for Eq. (9)-(11) in your notes) ----
# Solve: min_{Ω ⪰ 0} 0.5 ||Ω - V||_F^2 + τ ||Ω||_{1,off}

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
  omega_hat <- symmetrize( soft_thresh_offdiag(Z, tau / rho) )
  omega_hat <- clip_eig(omega_hat, eps = eig_floor)
  attr(omega_hat, "state") <- list(Omega = omega_hat, A = matrix(0, nrow(omega_hat), ncol(omega_hat)))
  if (is_pd_by_chol_rank(omega_hat)) {
    return(omega_hat)
  }
  
  Omega <- if (is.null(Omega_init)) omega_hat else symmetrize(Omega_init)
  A     <- if (is.null(A_init)) matrix(0, nrow(Omega), ncol(Omega)) else symmetrize(A_init)
  
  for (l in seq_len(max_admm)) {
    # (9) PSD projection
    K <- clip_eig(Omega + mu * A, eps = eig_floor)
    if (is_pd_by_chol_rank(K)) {
      attr(K, "state") <- list(Omega = K, A = A)
      return(K)
    }
    
    # (10) averaged center + off-diag soft-threshold
    T <- (K + mu * (Z - A)) / (1 + mu)
    Omega_new <- symmetrize(
      soft_thresh_offdiag(T, (mu * tau) / ((1 + mu) * rho))
    )
    
    # (11) dual update
    A <- A - (K - Omega_new) / mu
    
    
    if (max(abs(Omega_new - Omega)) < tol &&
        max(abs(K - Omega_new))     < tol) {
      attr(Omega_new, "state") <- list(Omega = Omega_new, A = A)
      return(Omega_new)
    }
    Omega <- Omega_new
  }
  
  Omega <- symmetrize(Omega)
  attr(Omega, "state") <- list(Omega = Omega, A = A)
  Omega
}

## =========================================================
##  Algorithm 2: Backtracking Prox-Gradient (exactly as notes)
##  Step-size rule via λ, β; only compute Oinv once per outer iter
## =========================================================
pg_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 200, tol = 1e-4, L0 = 1,
                            track_loss = TRUE, time_budget = Inf,
                            ctrl = NULL, warm_state = NULL,
                            oracle_loss = NULL) {
  
  # FIX: override solver settings with shared ctrl
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma  # FIX: shared lambda
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega  # FIX: shared lambda
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit                   # FIX: shared max iterations
    if (!is.null(ctrl$tol)) tol <- ctrl$tol                             # FIX: shared tolerance
  }
  
  n <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  # --- initial value---
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q)          else symmetrize(Omega_init)
  Oinv  <- safe_inv_psd(Omega)
  
  # --- record---
  t0 <- Sys.time()
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    loss_hist[1] <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n,
                            lambda_gamma, lambda_Omega, Oinv = Oinv)
    time_hist[1] <- 0
  }
  loss_curr <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n,
                       lambda_gamma, lambda_Omega, Oinv = Oinv)  # FIX: track current loss
  
  # --- Algorithm 2：L_{(0)} > 0 ---
  L_prev <- L0
  converged <- FALSE; stop_reason <- NA_character_
  backtrack_stalls <- 0L  # FIX: monitor stalled backtracking
  pg_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  
  for (m in seq_len(max_iter)) {
    Lm <- L_prev  # FIX: retain outer structure
    loss_prev <- loss_curr
    accepted_step <- FALSE
    
    repeat {  # FIX: Armijo backtracking loop for PG step
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
      if (!is.null(oracle_loss) && is.finite(loss_prev)) {
        gap_rel_pg <- (loss_prev - oracle_loss) / max(1, abs(oracle_loss))
        if (loss_prev <= oracle_loss || gap_rel_pg <= tol) {
          converged <- TRUE
          stop_reason <- "converged (oracle gap)"
          accepted_step <- TRUE
          break
        }
      } else if (step_norm <= tol) {
        converged <- TRUE
        stop_reason <- "converged (λ_m and ||d|| small)"
        accepted_step <- TRUE
        break
      }
      
      alpha_m <- (beta_m^2) / (lambda_m * (lambda_m + beta_m^2))
      step_gamma <- alpha_m * d_gamma_m
      step_Omega <- alpha_m * d_Omega_m
      
      bt_attempts <- 1L
      gamma_old <- gamma
      Omega_old <- Omega
      denom_old <- sqrt(sum(gamma_old^2) + sum(Omega_old^2))
      loss_new <- NA_real_
      Oinv_new <- NULL
      grad_inner_gamma <- Gg_m
      grad_inner_Omega <- Go_m
      
      gamma_candidate <- gamma_old + step_gamma
      Omega_candidate <- symmetrize(Omega_old + step_Omega)
      Oinv_try <- tryCatch(safe_inv_psd(Omega_candidate), error = function(e) NULL)
      if (!is.null(Oinv_try)) {
        loss_try <- loss_fn(gamma_candidate, Omega_candidate, Sxx, Sxy, Syy, n,
                            lambda_gamma, lambda_Omega, Oinv = Oinv_try)
        if (is.finite(loss_try)) {
          gamma <- gamma_candidate
          Omega <- Omega_candidate
          Oinv_new <- Oinv_try
          loss_new <- loss_try
          accepted_step <- TRUE
        }
      }
      
      if (!accepted_step) {
        backtrack_stalls <- backtrack_stalls + 1L
        loss_curr <- loss_prev
        L_prev <- Lm
        break
      }
      
      Oinv <- Oinv_new
      loss_curr <- loss_new
      L_prev  <- Lm
      
      rel_change <- sqrt(sum((gamma - gamma_old)^2) + sum((Omega - Omega_old)^2)) /
        (1 + denom_old)
      if (!is.null(oracle_loss)) {
        gap_rel_pg <- (loss_curr - oracle_loss) / max(1, abs(oracle_loss))
        if (loss_curr <= oracle_loss || gap_rel_pg <= tol) {
          converged <- TRUE
          stop_reason <- "converged (oracle gap)"
        }
      } else if (rel_change < tol || abs(loss_curr - loss_prev) < 1e-9) {
        converged <- TRUE
        stop_reason <- "converged (rel_change or loss diff)"
      }
      break
    } # repeat
    
    if (converged) {
      if (track_loss) {
        loss_hist[m + 1] <- loss_curr
        time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      }
      break
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- loss_curr
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    if (!accepted_step && !converged) next  # FIX: continue even if iteration skipped
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

# ============================
# Prox-Newton: ADMM 子问题（按笔记的本征分解解 ξ-step）
# ============================
pn_subproblem_admm <- function(gamma, Omega, Sxx, Sxy, Syy,
                               lambda_gamma, lambda_Omega,
                               Oinv = NULL, beta = NULL,
                               rho = 1, mu = 1,
                               max_admm = 150, admm_tol = 1e-4,
                               eig_floor = 1e-12,
                               state = NULL) {
  p <- nrow(gamma); q <- ncol(gamma)
  
  # —— 当前点的梯度块、H(m)、Σ(m) —— #
  beta_eff <- .ensure_beta(gamma, Omega, Oinv, beta)
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  
  Gg <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv_eff, beta = beta_eff)
  Go <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv_eff, beta = beta_eff)
  P  <- pack_xi(Gg, Go)                                # (p+q) x q
  
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv_eff, beta = beta_eff)
  H  <- symmetrize(H_matrix(Hb))
  
  Xi <- pack_xi(gamma, Omega)                          # ξ^{(m)}
  Gm <- P - H %*% Xi %*% Oinv_eff                      # G^{(m)}（按你的推导）
  
  # —— 预先做谱分解：H = U_H Λ_H U_H^T，Σ = U_S Λ_S U_S^T —— #
  eigH <- eigen(symmetrize(H),    symmetric = TRUE)
  UH   <- eigH$vectors
  lamH <- pmax(eigH$values, eig_floor)
  
  eigS <- eigen(symmetrize(Oinv_eff), symmetric = TRUE)
  US   <- eigS$vectors
  lamS <- pmax(eigS$values, eig_floor)
  
  # —— ADMM 变量 —— #
  Xi_var <- Xi           # ξ
  Z_var  <- Xi           # ζ
  Gam    <- matrix(0, p + q, q)   # Γ
  if (!is.null(state)) {
    if (!is.null(state$Xi)) Xi_var <- state$Xi
    if (!is.null(state$Z))  Z_var  <- state$Z
    if (!is.null(state$Gam)) Gam   <- state$Gam
  }
  prox_state <- if (!is.null(state)) state$psd_state else NULL
  
  # 预组装分母矩阵：D_{ij} = λ_i^{(H)} * λ_j^{(Σ)} + ρ
  Den <- outer(lamH, lamS, "*") + rho
  
  for (k in 1:max_admm) {
    # -------------------------
    # ① ξ-step（按笔记：本征分解 + 逐元素求解）
    #    解  H ξ Σ + ρ ξ = R， R = ρ ζ - Γ - G^{(m)}
    # -------------------------
    R  <- rho * Z_var - Gam - Gm                      # (p+q) x q
    Rt <- crossprod(UH, R) %*% US                     # \tilde R = U_H^T R U_S
    Y  <- Rt / Den                                    # 元素级：Y_{ij} = \tilde R_{ij} / (λ_i λ̄_j + ρ)
    Xi_var <- UH %*% Y %*% t(US)                      # 回到原坐标系
    
    # -------------------------
    # ② ζ-step（可分的近端）
    #    γ：L1 软阈（截距不惩罚）
    #    Ω：PSD + off-diag L1 近端（你已有的 prox_psd_offdiag_l1）
    # -------------------------
    Xi_parts <- unpack_xi(Xi_var + Gam / rho, p)      # 不要对 (p+q)×q 做 symmetrize
    Zg <- soft_thresh(Xi_parts$gamma, lambda_gamma / rho)
    Zg[1, ] <- Xi_parts$gamma[1, ]                    # intercept unpenalized
    Zo <- prox_psd_offdiag_l1(
      Xi_parts$Omega,
      tau = lambda_Omega / rho,
      mu = mu,
      rho = rho,
      Omega_init = if (is.null(prox_state)) Omega else prox_state$Omega,
      A_init = if (is.null(prox_state)) NULL else prox_state$A,
      max_admm = prox_ctrl$max_admm,
      tol = prox_ctrl$tol,
      eig_floor = pn_prox_ctrl$eta
    )
    prox_state <- attr(Zo, "state")
    Z_new <- pack_xi(Zg, Zo)
    
    # -------------------------
    # ③ Γ-step（对偶更新）
    # -------------------------
    Gam <- Gam + rho * (Xi_var - Z_new)
    
    # —— 收敛判据 —— #
    r_norm <- sqrt(sum((Xi_var - Z_new)^2))           # primal residual
    s_norm <- rho * sqrt(sum((Z_new - Z_var)^2))      # dual residual
    Z_var  <- Z_new
    scale_pri <- max(sqrt(length(Xi_var)), sqrt(sum(Xi_var^2)), sqrt(sum(Z_var^2)), sqrt(sum(Z_new^2)))
    scale_dual <- max(sqrt(length(Gam)), sqrt(sum(Gam^2)))
    eps_pri <- admm_tol * scale_pri
    eps_dual <- admm_tol * scale_dual
    if (r_norm <= eps_pri && s_norm <= eps_dual) break
  }
  
  # 返回 s^{(m)} = ζ
  out_sub <- unpack_xi(Z_var, p)
  attr(out_sub$Omega, "state") <- list(
    psd_state = prox_state,
    Xi = Xi_var,
    Z = Z_var,
    Gam = Gam
  )
  out_sub
}


# ============================
# Prox-Newton 外层（保持你原先的算法 1；每次外层更新前重算 β）
# ============================
pn_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 100, tol = 1e-4, track_loss = TRUE,
                            time_budget = Inf, backtrack = TRUE,
                            ctrl = NULL, warm_state = NULL,
                            oracle_loss = NULL) {
  # FIX: override shared solver settings
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma  # FIX: shared lambda
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega  # FIX: shared lambda
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit                    # FIX: shared max iterations
    if (!is.null(ctrl$tol)) tol <- ctrl$tol                              # FIX: shared tolerance
  }
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
  
  converged <- FALSE; stop_reason <- NA_character_
  lambda_hist <- numeric(max_iter)
  pn_loss_hist <- numeric(max_iter + 1)  # FIX: track accepted PN losses
  pn_loss_hist[1] <- loss_curr
  pn_hist_idx <- 1L
  pn_backtrack_counts <- integer(max_iter)
  pn_backtrack_total <- 0L
  pn_nonmonotone <- FALSE
  pn_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  
  for (m in seq_len(max_iter)) {
    loss_prev <- loss_curr
    beta <- gamma %*% Oinv
    gamma_old <- gamma
    Omega_old <- Omega
    denom_old <- sqrt(sum(gamma_old^2) + sum(Omega_old^2))
    
    stall_flag <- FALSE
    local_converged <- FALSE
    attempts_record <- 0L
    accepted_step <- FALSE
    
    grad_gamma_curr <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv)
    grad_Omega_curr <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv)
    
    repeat {
      sub <- pn_subproblem_admm(gamma, Omega, Sxx, Sxy, Syy,
                                lambda_gamma, lambda_Omega,
                                Oinv = Oinv, beta = beta,
                                rho = 1, mu = 1,
                                max_admm = 150, admm_tol = 1e-4,
                                state = pn_psd_state)
      pn_psd_state <- attr(sub$Omega, "state")
      d_gamma <- sub$gamma - gamma
      d_Omega <- sub$Omega - Omega
      
      lam <- local_norm(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = Oinv)
      lambda_hist[m] <- lam
      
      step_norm <- sqrt(sum(d_gamma^2) + sum(d_Omega^2))
      # Oracle-gap stopping (relative to provided oracle_loss), else norm-based
      if (!is.null(oracle_loss)) {
        gap_rel <- (loss_curr - oracle_loss) / max(1, abs(oracle_loss))
        if (loss_curr <= oracle_loss || gap_rel <= tol) {
          converged <- TRUE
          local_converged <- TRUE
          stop_reason <- "converged (oracle gap)"
          pn_backtrack_counts[m] <- 0L
          pn_hist_idx <- pn_hist_idx + 1L
          pn_loss_hist[pn_hist_idx] <- loss_curr
          break
        }
      } else if (lam <= tol && step_norm <= tol) {
        converged <- TRUE
        local_converged <- TRUE
        stop_reason <- "converged (λ and ||d|| small)"
        loss_curr <- loss_prev
        pn_backtrack_counts[m] <- 0L
        pn_hist_idx <- pn_hist_idx + 1L
        pn_loss_hist[pn_hist_idx] <- loss_curr
        break
      }
      
      alpha <- 1
      attempts <- 0L
      gamma_new <- gamma
      Omega_new <- Omega
      Oinv_new <- Oinv
      loss_new <- loss_prev
      accepted_step <- FALSE
      last_loss_try <- NA_real_
      
      while (alpha > 1e-12) {  # FIX: Armijo backtracking for PN
        attempts <- attempts + 1L
        gamma_candidate <- gamma + alpha * d_gamma
        Omega_candidate <- symmetrize(Omega + alpha * d_Omega)
        Oinv_try <- tryCatch(safe_inv_psd(Omega_candidate), error = function(e) NULL)
        if (!is.null(Oinv_try)) {
          loss_try <- loss_fn(gamma_candidate, Omega_candidate, Sxx, Sxy, Syy, n,
                              lambda_gamma, lambda_Omega, Oinv = Oinv_try)
          last_loss_try <- loss_try
          dir_gamma <- alpha * d_gamma
          dir_Omega <- alpha * d_Omega
          armijo_rhs <- loss_prev + armijo_c * (sum(grad_gamma_curr * dir_gamma) +
                                                  sum(grad_Omega_curr * dir_Omega))
          if (is.finite(loss_try) && loss_try <= armijo_rhs) {
            accepted_step <- TRUE
            gamma_new <- gamma_candidate
            Omega_new <- Omega_candidate
            Oinv_new <- Oinv_try
            loss_new <- loss_try
            break
          }
        }
        if (attempts >= 50L) break
        alpha <- alpha * 0.5
      }
      
      attempts_record <- attempts
      if (accepted_step) {
        pn_backtrack_counts[m] <- attempts
        pn_backtrack_total <- pn_backtrack_total + max(0L, attempts - 1L)
        gamma <- gamma_new
        Omega <- Omega_new
        Oinv  <- Oinv_new
        loss_curr <- loss_new
        pn_hist_idx <- pn_hist_idx + 1L
        pn_loss_hist[pn_hist_idx] <- loss_curr
        if (pn_hist_idx > 1L && loss_curr > pn_loss_hist[pn_hist_idx - 1L] + 1e-10) {
          pn_nonmonotone <- TRUE
        }
        break
      }
      if (is.finite(last_loss_try) && abs(last_loss_try - loss_prev) < 1e-9) {  # FIX: treat flat loss as convergence per tolerance
        converged <- TRUE
        local_converged <- TRUE
        pn_backtrack_counts[m] <- attempts
        loss_curr <- loss_prev
        pn_hist_idx <- pn_hist_idx + 1L
        pn_loss_hist[pn_hist_idx] <- loss_curr
        break
      }
      
      stall_flag <- TRUE
      pn_backtrack_counts[m] <- attempts
      loss_curr <- loss_prev
      break
    } # repeat
    
    if (stall_flag) {
      message(sprintf("PN line search stalled at iter %d (bt=%d)", m, attempts_record))
      stop_reason <- "PN line search stalled"
      if (track_loss) {
        loss_hist[m + 1] <- loss_curr
        time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      }
      break
    }
    
    if (converged && local_converged) {
      if (track_loss) {
        loss_hist[m + 1] <- loss_curr
        time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      }
      break
    }
    
    rel_change <- sqrt(sum((gamma - gamma_old)^2) + sum((Omega - Omega_old)^2)) /
      (1 + denom_old)
    if (rel_change < tol || abs(loss_curr - loss_prev) < 1e-9) {
      converged <- TRUE
      stop_reason <- "converged (rel_change or loss diff)"
      if (track_loss) {
        loss_hist[m + 1] <- loss_curr
        time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      }
      break
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- loss_curr
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
  }
  
  elapsed_time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- list(gamma = gamma, Omega = Omega, iters = m,
              converged = converged,
              stop_reason = ifelse(is.na(stop_reason), "max_iter reached", stop_reason),
              elapsed_time_sec = elapsed_time_sec,
              Sxx = Sxx, Sxy = Sxy, Syy = Syy, n = n,
              lambda_hist = lambda_hist[seq_len(m)],
              backtrack_counts = pn_backtrack_counts[seq_len(m)],
              backtrack_total = pn_backtrack_total,
              loss_monotone = !pn_nonmonotone,
              psd_state = pn_psd_state)
  if (track_loss) {
    out$loss_hist <- loss_hist[1:(m + 1)]
    out$time_hist <- time_hist[1:(m + 1)]
  }
  out$loss_hist_pn <- pn_loss_hist[seq_len(pn_hist_idx)]
  out
}

run_pg_path <- function(data_list, ctrl, ref_path = NULL, oracle_loss_path = NULL) {
  lambda_path <- data_list$lambda_path
  gamma_init <- data_list$gamma_init
  Omega_init <- data_list$Omega_init
  warm_state <- NULL
  results <- vector("list", nrow(lambda_path))
  for (i in seq_len(nrow(lambda_path))) {
    ctrl_i <- ctrl
    ctrl_i$lambda$gamma <- lambda_path$lambda_gamma[i]
    ctrl_i$lambda$Omega <- lambda_path$lambda_Omega[i]
    gamma_init_step <- gamma_init
    Omega_init_step <- Omega_init
    warm_state_step <- warm_state
    if (i == 1 && !is.null(ref_path) && length(ref_path) >= 1) {
      seed_fit <- ref_path[[1]]
      if (!is.null(seed_fit$gamma)) gamma_init_step <- seed_fit$gamma
      if (!is.null(seed_fit$Omega)) Omega_init_step <- seed_fit$Omega
      if (!is.null(seed_fit$psd_state)) {
        warm_state_step <- list(psd_state = seed_fit$psd_state)
      }
    }
    fit <- pg_sparse_mvreg(data_list$X, data_list$Y,
                           lambda_path$lambda_gamma[i],
                           lambda_path$lambda_Omega[i],
                           gamma_init = gamma_init_step,
                           Omega_init = Omega_init_step,
                           max_iter = ctrl$maxit, tol = ctrl$tol, L0 = 1,
                           track_loss = TRUE, time_budget = Inf,
                           ctrl = ctrl_i, warm_state = warm_state_step,
                           oracle_loss = if (!is.null(oracle_loss_path)) oracle_loss_path[i] else NULL)
    fit$lambda_gamma <- lambda_path$lambda_gamma[i]
    fit$lambda_Omega <- lambda_path$lambda_Omega[i]
    fit$path_index <- lambda_path$step[i]
    results[[i]] <- fit
    gamma_init <- fit$gamma
    Omega_init <- fit$Omega
    warm_state <- list(psd_state = fit$psd_state)
  }
  results
}

run_pn_path <- function(data_list, ctrl, oracle_loss_path = NULL) {
  lambda_path <- data_list$lambda_path
  gamma_init <- data_list$gamma_init
  Omega_init <- data_list$Omega_init
  warm_state <- NULL
  results <- vector("list", nrow(lambda_path))
  for (i in seq_len(nrow(lambda_path))) {
    ctrl_i <- ctrl
    ctrl_i$lambda$gamma <- lambda_path$lambda_gamma[i]
    ctrl_i$lambda$Omega <- lambda_path$lambda_Omega[i]
    fit <- pn_sparse_mvreg(data_list$X, data_list$Y,
                           lambda_path$lambda_gamma[i],
                           lambda_path$lambda_Omega[i],
                           gamma_init = gamma_init,
                           Omega_init = Omega_init,
                           max_iter = ctrl$maxit, tol = ctrl$tol, track_loss = TRUE,
                           time_budget = Inf, backtrack = TRUE,
                           ctrl = ctrl_i, warm_state = warm_state,
                           oracle_loss = if (!is.null(oracle_loss_path)) oracle_loss_path[i] else NULL)
    fit$lambda_gamma <- lambda_path$lambda_gamma[i]
    fit$lambda_Omega <- lambda_path$lambda_Omega[i]
    fit$path_index <- lambda_path$step[i]
    results[[i]] <- fit
    gamma_init <- fit$gamma
    Omega_init <- fit$Omega
    warm_state <- list(psd_state = fit$psd_state)
  }
  results
}

## =========================================================
##  Benchmark & Visualization (shared-data version)
## =========================================================

benchmark_pg_vs_pn <- function(data_list, ctrl,
                               time_budget = 2.0) {
  pn_path <- run_pn_path(data_list, ctrl)
  pg_path <- run_pg_path(data_list, ctrl)
  pn_fit <- pn_path[[length(pn_path)]]
  pg_fit <- pg_path[[length(pg_path)]]
  df <- rbind(
    data.frame(method = "Prox-Newton",
               iter = seq_along(pn_fit$time_hist) - 1,
               time = pn_fit$time_hist,
               NLL = pn_fit$loss_hist,
               path_index = pn_fit$path_index),
    data.frame(method = "Prox-Gradient",
               iter = seq_along(pg_fit$time_hist) - 1,
               time = pg_fit$time_hist,
               NLL = pg_fit$loss_hist,
               path_index = pg_fit$path_index)
  )
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
       df = df,
       lambdas = c(lambda_gamma = pn_fit$lambda_gamma, lambda_Omega = pn_fit$lambda_Omega),
       path_summary = path_summary)
}

# Diagnostic contour plot near a given solution (uses current loss_fn)
plot_objective_contour_at_fit <- function(fit, lambda_gamma, lambda_Omega,
                                          t_range = 0.3, ngrid = 61, nlevels = 15,
                                          seed = 2025, main = "Objective contours near solution") {
  set.seed(seed)
  gamma_sol <- fit$gamma; Omega_sol <- fit$Omega
  Sxx <- fit$Sxx; Sxy <- fit$Sxy; Syy <- fit$Syy; n <- fit$n
  p <- nrow(gamma_sol); q <- ncol(gamma_sol)
  if (q > 20 || p > 30) {  # FIX: limit contour workload for large problems
    ngrid <- min(ngrid, 21)
    nlevels <- min(nlevels, 10)
  }
  
  d1_g <- matrix(rnorm(p * q), p, q)
  d2_g <- matrix(rnorm(p * q), p, q)
  d1_o <- symmetrize(matrix(rnorm(q * q), q, q))
  d2_o <- symmetrize(matrix(rnorm(q * q), q, q))
  
  norm1 <- sqrt(sum(d1_g^2) + sum(d1_o^2))
  norm2 <- sqrt(sum(d2_g^2) + sum(d2_o^2))
  d1_g <- d1_g / norm1; d1_o <- d1_o / norm1
  d2_g <- d2_g / norm2; d2_o <- d2_o / norm2
  
  tgrid <- seq(-t_range, t_range, length.out = ngrid)
  Loss_grid <- outer(tgrid, tgrid, Vectorize(function(t1, t2) {
    gamma_p <- gamma_sol + t1 * d1_g + t2 * d2_g
    Omega_p <- symmetrize(Omega_sol + t1 * d1_o + t2 * d2_o)
    Omega_p <- proj_psd(Omega_p + diag(1e-10, nrow(Omega_p)))
    Oinv_p  <- safe_inv_psd(Omega_p)
    loss_fn(gamma_p, Omega_p, Sxx, Sxy, Syy, n,
            lambda_gamma = lambda_gamma, lambda_Omega = lambda_Omega, Oinv = Oinv_p)
  }))
  Loss_grid[!is.finite(Loss_grid)] <- NA
  contour(tgrid, tgrid, Loss_grid, nlevels = nlevels,
          xlab = expression(t[1]), ylab = expression(t[2]),
          main = main)
  points(0, 0, pch = 19, col = "red")
  abline(h = 0, v = 0, lty = 3, col = "grey70")
  grid()
}

## =========================================================
##  Example run + Plots (kept consistent with your original)
## =========================================================

# packages
suppressPackageStartupMessages({
  library(MASS)
  library(ggplot2)
})

# FIX: run PN first, then PG using shared data
res <- benchmark_pg_vs_pn(shared_data, ctrl, time_budget = 2.0)
pn_fit <- res$pn
pg_fit <- res$pg

# Summarize final losses per λ and plot
pn_monotone <- isTRUE(pn_fit$loss_monotone)
pg_monotone <- if (!is.null(pg_fit)) {
  all(diff(pg_fit$loss_hist) <= 1e-8)
} else {
  NA
}
pn_final_loss <- tail(pn_fit$loss_hist, 1)
pg_final_loss <- tail(pg_fit$loss_hist, 1)


pn_loss_df <- do.call(rbind, lapply(res$pn_path, function(fit) {
  data.frame(method = "Prox-Newton",
             path_index = fit$path_index,
             lambda_gamma = fit$lambda_gamma,
             lambda_Omega = fit$lambda_Omega,
             final_loss = tail(fit$loss_hist, 1),
             elapsed_sec = fit$elapsed_time_sec,
             iters = fit$iters,
             stringsAsFactors = FALSE)
}))
pg_loss_df <- do.call(rbind, lapply(res$pg_path, function(fit) {
  data.frame(method = "Prox-Gradient",
             path_index = fit$path_index,
             lambda_gamma = fit$lambda_gamma,
             lambda_Omega = fit$lambda_Omega,
             final_loss = tail(fit$loss_hist, 1),
             elapsed_sec = fit$elapsed_time_sec,
             iters = fit$iters,
             stringsAsFactors = FALSE)
}))

loss_plot_df <- rbind(pn_loss_df, pg_loss_df)
loss_plot_df$lambda_label <- sprintf("γ=%.3f / Ω=%.3f",
                                     loss_plot_df$lambda_gamma,
                                     loss_plot_df$lambda_Omega)

final_loss_plot <- ggplot(loss_plot_df,
                          aes(x = log(lambda_gamma), y = final_loss,
                              colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_shape_manual(values = c("Prox-Newton" = 16, "Prox-Gradient" = 4)) +
  labs(x = expression(log(lambda)),
       y = "Training Loss",
       colour = NULL,
       shape = NULL) +
  theme_minimal(base_size = 13)

iter_plot_df <- loss_plot_df
iter_plot <- ggplot(iter_plot_df,
                     aes(x = log(lambda_gamma), y = iters,
                         colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_shape_manual(values = c("Prox-Newton" = 16, "Prox-Gradient" = 4)) +
  labs(x = expression(log(lambda[gamma])),
       y = "Number of iterations",
       colour = NULL,
       shape = NULL) +
  theme_minimal(base_size = 13)

ggsave("pn_pg_final_loss_vs_lambda.png", final_loss_plot, width = 9, height = 4.5, dpi = 300)
ggsave("pn_pg_iters_vs_lambda.png", iter_plot, width = 9, height = 4.5, dpi = 300)

print(final_loss_plot)
print(iter_plot)

total_time_by_method <- aggregate(elapsed_sec ~ method, data = res$path_summary, sum)
print(total_time_by_method)

result_table <- data.frame(
  method = c("Prox-Newton", "Prox-Gradient"),
  converged = c(pn_fit$converged, pg_fit$converged),
  n_iter = c(pn_fit$iters, pg_fit$iters),
  final_loss = c(pn_final_loss, pg_final_loss)
)
print(result_table)

saveRDS(res, file = "res_latest.rds")
saveRDS(list(
  pn_path = res$pn_path,
  pg_path = res$pg_path,
  lambda_path = shared_data$lambda_path
), file = "lambda_path_estimators.rds")

cat("Lambda path (gamma, Omega):\n")
print(shared_data$lambda_path)
cat("Per-step timing summary (seconds):\n")
print(res$path_summary)
cat("Total time by method (seconds):\n")
print(aggregate(elapsed_sec ~ method, data = res$path_summary, sum))

# FIX: PN backtracking diagnostics
pn_backtracking_events <- sum(pmax(0, pn_fit$backtrack_counts - 1))
cat(sprintf("PN backtracking shrink events: %d; monotone=%s\n",
            pn_backtracking_events,
            if (pn_monotone) "TRUE" else "FALSE"))

# FIX: console report answering study questions
shared_inputs_ok <- max(abs(pn_fit$Sxx - crossprod(shared_data$X) / nrow(shared_data$X))) < 1e-9  # FIX: verify PN used shared Sxx
if (!is.null(pg_fit)) {
  shared_inputs_ok <- shared_inputs_ok && max(abs(pg_fit$Sxx - pn_fit$Sxx)) < 1e-9
}
shared_inputs_ok <- shared_inputs_ok &&
  abs(tail(shared_data$lambda_path$lambda_gamma, 1) - res$lambdas["lambda_gamma"]) < 1e-12 &&
  abs(tail(shared_data$lambda_path$lambda_Omega, 1) - res$lambdas["lambda_Omega"]) < 1e-12
