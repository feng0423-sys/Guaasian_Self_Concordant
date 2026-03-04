# PG-01：Loss（scaled self-concordant）+ 设计矩阵生成 + 稀疏真值代码

本文件对应你要求的 **第 1 个 md**：  
- self-concordant 的 Gaussian loss（文档 `main_Dec2_revised_gemini.pdf` Eq (2.2)）  
- 设计矩阵 `X` / `Xc`（含 intercept）的生成方式  
- 生成不同 sparsity 的真值（β / γ / Ω）以及常用统计量

> 记号对齐（和 pdf 一致）：  
> - `β ∈ R^{p_aug×q}`：原始回归系数（含 intercept 行）  
> - `Ω ∈ S_{++}^q`：precision matrix  
> - **自然参数**：`γ := β Ω`（同 pdf）  
> - `ξ := (γ, Ω)`  
> - `D_δ := { (γ,Ω): Ω ⪰ δ I_q }`

---

## 1. 数据与经验散度矩阵（Sxx/Sxy/Syy）

给定数据矩阵：

- `Xc ∈ R^{n×p_aug}`（通常 `p_aug = p + 1`，第一列是 intercept 1）
- `Y  ∈ R^{n×q}`

定义经验散度（和代码 `crossprod()/n` 一致）：

\[
\widehat\Sigma_{XX} = \frac{1}{n} X_c^\top X_c \quad (p_{aug}\times p_{aug})
\]
\[
\widehat\Sigma_{XY} = \frac{1}{n} X_c^\top Y \quad (p_{aug}\times q)
\]
\[
\widehat\Sigma_{YY} = \frac{1}{n} Y^\top Y \quad (q\times q)
\]

R 代码（建议封装成函数）：

```r
scr_scatter <- function(Xc, Y) {
  n <- nrow(Xc)
  list(
    n = n,
    Sxx = crossprod(Xc) / n,
    Sxy = crossprod(Xc, Y) / n,
    Syy = crossprod(Y) / n
  )
}
```

---

## 2. Self-concordant 的 scaled loss（Eq (2.2)）

pdf 先定义 reparameterized negative log-likelihood \(\tilde L\)，然后使用缩放后的 loss：

\[
L(\gamma,\Omega)
= -\log\det(\Omega)
+ \mathrm{tr}\big(\Omega^{-1}\gamma^\top \widehat\Sigma_{XX}\gamma\big)
- 2\,\mathrm{tr}\big(\widehat\Sigma_{XY}^\top\gamma\big)
+ \mathrm{tr}\big(\Omega\widehat\Sigma_{YY}\big).
\]

这正是你 `2.18_simulation.R` 里 `smooth_obj()` 采用的形式（注意它是 scaled 版本，不带 1/2）。

R 实现（数值上要安全 logdet / inverse）：

```r
scr_smooth_loss <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL) {
  if (is.null(Oinv)) Oinv <- safe_inv_spd(Omega)
  M <- t(gamma) %*% Sxx %*% gamma  # q×q
  -safe_logdet_spd(Omega) +
    sum(Oinv * M) -
    2 * sum(Sxy * gamma) +
    sum(Omega * Syy)
}
```

---

## 3. penalty（SCR-MRCE）与目标函数

SCR-MRCE 的 composite 目标（pdf Eq (3.1)-(3.2)）：

\[
F(\gamma,\Omega) = L(\gamma,\Omega) + \lambda_\gamma\|\gamma_P\|_1 + \lambda_\Omega \|\Omega\|_{1,\mathrm{off}}
\]

其中：
- \(U\subset\{1,\dots,p_{aug}\}\)：不惩罚的行（通常 \(U=\{1\}\) intercept）
- \(P=\{1,\dots,p_{aug}\}\setminus U\) ：惩罚的行
- \(\|\Omega\|_{1,\mathrm{off}}=\sum_{i\neq j}|\Omega_{ij}|\)

实现：

```r
scr_penalty_obj <- function(gamma, Omega, lambda_gamma, lambda_Omega, U = 1L) {
  p_aug <- nrow(gamma)
  U <- sort(unique(as.integer(U)))
  P <- setdiff(seq_len(p_aug), U)

  pen_gamma <- if (length(P) > 0) sum(abs(gamma[P, , drop = FALSE])) else 0
  pen_Omega <- sum(abs(Omega)) - sum(abs(diag(Omega)))

  lambda_gamma * pen_gamma + lambda_Omega * pen_Omega
}
```

---

## 4. 设计矩阵 X 的生成（simulation 版本）

### 4.1 常用做法：高斯设计 + mean-centering + intercept

```r
X <- matrix(rnorm(n * p), n, p)
X <- scale(X, center = TRUE, scale = FALSE)  # mean-center predictors
Xc <- cbind(1, X)                             # add intercept
```

这样会使：
- intercept 与其他列近似正交：\(\widehat\Sigma_{XX}[P,U]\approx 0\)
- 有利于 λmax 公式简化（见 PG-03）

---

## 5. 稀疏真值：β / Ω 的生成（并得到 γ=βΩ）

### 5.1 生成稀疏 β（列稀疏或元素稀疏）

下面沿用你代码的“每列固定 nnz 个激活行”的思路，但变量命名改成 **beta_true** 更准确：

```r
scr_generate_sparse_beta <- function(p_aug, q,
                                     sparsity = 0.05,
                                     nnz_per_col = NULL,
                                     mag_range = c(1, 2),
                                     intercept_zero = TRUE,
                                     signed = TRUE) {
  beta <- matrix(0, p_aug, q)
  idx_pool <- if (intercept_zero && p_aug > 1) 2:p_aug else seq_len(p_aug)

  if (is.null(nnz_per_col)) {
    nnz_per_col <- max(1, round(sparsity * length(idx_pool)))
  } else {
    nnz_per_col <- max(1, min(nnz_per_col, length(idx_pool)))
  }

  for (j in seq_len(q)) {
    active <- sample(idx_pool, nnz_per_col)
    vals <- runif(nnz_per_col, min(mag_range), max(mag_range))
    if (signed) vals <- vals * sample(c(-1,1), nnz_per_col, replace = TRUE)
    beta[active, j] <- vals
  }

  if (intercept_zero && p_aug >= 1) beta[1, ] <- 0
  beta
}
```

### 5.2 生成稀疏 Ω（SPD + 控制 condition number）

你的 `generate_sparse_omega()` 采用：先随机稀疏 offdiag，再加 diagonal shift 使其 SPD 并控制 κ。沿用即可：

```r
scr_generate_sparse_omega <- function(q,
                                      kappa_target = 50,
                                      offdiag_prob = 0.05,
                                      base_val = 0.5) {
  B <- matrix(0, q, q)
  mask <- matrix(runif(q^2) < offdiag_prob, q, q)
  mask[lower.tri(mask, diag = TRUE)] <- FALSE
  B[mask] <- base_val
  B <- B + t(B)

  eig_B <- eigen(symmetrize(B), symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(eig_B); lambda_max <- max(eig_B)

  delta <- (lambda_max - kappa_target * lambda_min) / (kappa_target - 1)
  if (!is.finite(delta) || delta <= -lambda_min + 1e-8) delta <- -lambda_min + 1e-4

  Omega <- B + delta * diag(q)
  eig_O <- eigen(symmetrize(Omega), symmetric = TRUE, only.values = TRUE)$values
  attr(Omega, "cond_est") <- max(eig_O) / min(eig_O)
  Omega
}
```

### 5.3 用 (β_true, Ω_true) 得到 γ_true

\[
\gamma_{\text{true}} = \beta_{\text{true}}\Omega_{\text{true}}.
\]

```r
gamma_true <- beta_true %*% Omega_true
```

> ⚠️ 你原始 `2.18_simulation.R` 里把 `gamma_true[-1,]` 当作 `Beta_true` 用来生成 Y，  
> 这会让 “真值 γ” 和 “生成数据的 β” 不一致（除非 Ω=I）。  
> 建议以后 simulation 用上面这种一致方式：先生成 β、Ω，再令 γ=βΩ。

---

## 6. 生成 Y：多响应高斯回归

模型：\(Y_i|x_i\sim N_q(\beta^\top x_i,\Sigma)\)，其中 \(\Sigma=\Omega^{-1}\)。

```r
Sigma_true <- safe_inv_spd(Omega_true)
E <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = Sigma_true)
Y <- Xc %*% beta_true + E
```

---

## 7. 常用 sparsity / condition diagnostics

```r
matrix_sparsity_stats(gamma_hat, Omega_hat)
matrix_condition_stats(Omega_hat)
```

---

## 8. 给 Codex 的实现提示

- **不要混淆 β 与 γ**：优化变量是 (γ,Ω)，输出 β 用 `β = γ Ω^{-1}`。
- 经验散度用 `crossprod()/n`，维度严格对齐。
- 任何地方出现 `Omega` 都要 `symmetrize()`，并保证 `Omega ⪰ δ I`（δ 是 safeguard）。
