# PG-03：λγ,max / λΩ,max（Lemma 3.1）+ null-corner init + regularization path

本文件对应你要求的 **第 3 个 md**：  
- 给定 `Xc,Y`（或 `Sxx,Sxy,Syy`）计算 `lambda_gamma_max` 与 `lambda_Omega_max`  
- null-corner（全零角点）解析初始化（γ(0), Ω(0)）  
- 构造二维 regularization path（log spacing + warm-start）

> 参考：`main_Dec2_revised_gemini.pdf` Section 3.4, Lemma 3.1.

---

## 1. 为什么要算 λmax？

SCR-MRCE 的目标：

\[
F(\gamma,\Omega)=L(\gamma,\Omega)+\lambda_\gamma\|\gamma_P\|_1+\lambda_\Omega\|\Omega\|_{1,\mathrm{off}}.
\]

当 \(\lambda_\gamma,\lambda_\Omega\) 足够大时，最优解会落在一个“null corner”：
- 惩罚块完全为 0：\(\gamma_P=0\)
- \(\Omega\) 的 off-diagonal 也被压成 0（即 Ω 是对角）
- 不惩罚的 γU（如 intercept 行）仍然有解析解

我们想把路径的起点设在这个“严格为零”的角点，**这样整个路径计算既稳定又可 warm-start**。

---

## 2. 记号：U/P 行划分

- `U ⊂ {1,…,p_aug}`：不惩罚行（通常 intercept：`U={1}`）
- `P := {1,…,p_aug} \ U`：惩罚行（predictors）

---

## 3. Lemma 3.1：λγ,max 与 λΩ,max 的公式

令：

\[
B_0 := S_{XY}[U,\cdot]^\top\, S_{XX}[U,U]^{-1}\, S_{XY}[U,\cdot] \in \mathbb{R}^{q\times q}
\]

\[
B_\gamma := S_{XY}[P,\cdot] - S_{XX}[P,U]\, S_{XX}[U,U]^{-1}\, S_{XY}[U,\cdot] \in \mathbb{R}^{|P|\times q}.
\]

那么阈值为：

\[
\lambda_{\gamma,\max} = 2 \|B_\gamma\|_{\max}.
\]

\[
\lambda_{\Omega,\max} = \|(S_{YY} - B_0)_{\mathrm{off}}\|_{\max}.
\]

> 注意：这里用的是 **scaled loss L**，所以 \(\lambda_{\gamma,\max}\) 前面是 2。  
> 你旧代码里 `lambda_Omega_0 <- 2*max(abs(offdiag(Syy)))` 与 Lemma 3.1 不一致（应当无 2，并且要用 \(S_{YY}-B_0\)）。

---

## 4. null-corner 解析初始化

当 \(\lambda_\gamma \ge \lambda_{\gamma,\max}\) 且 \(\lambda_\Omega \ge \lambda_{\Omega,\max}\)，全局唯一最优解是：

\[
\Omega^{(0)} = \mathrm{diag}(S_{YY}-B_0)^{-1},
\]
\[
\gamma_U^{(0)} = S_{XX}[U,U]^{-1} S_{XY}[U,\cdot]\;\Omega^{(0)},
\]
\[
\gamma_P^{(0)} = 0.
\]

---

## 5. 代码：从 (Sxx,Sxy,Syy) 计算 λmax

建议封装成：

```r
scr_lambda_max <- function(Sxx, Sxy, Syy, U = 1L) {
  U <- sort(unique(as.integer(U)))
  P <- setdiff(seq_len(nrow(Sxx)), U)

  Sxx_UU <- Sxx[U, U, drop=FALSE]
  Sxy_U  <- Sxy[U, , drop=FALSE]
  Sxx_UU_inv <- solve(Sxx_UU)

  B0 <- t(Sxy_U) %*% Sxx_UU_inv %*% Sxy_U

  if (length(P) > 0) {
    Bgamma <- Sxy[P, , drop=FALSE] - Sxx[P, U, drop=FALSE] %*% Sxx_UU_inv %*% Sxy_U
    lambda_gamma_max <- 2 * max(abs(Bgamma))
  } else {
    lambda_gamma_max <- 0
  }

  R <- Syy - B0
  R_off <- R; diag(R_off) <- 0
  lambda_Omega_max <- max(abs(R_off))

  list(lambda_gamma_max=lambda_gamma_max, lambda_Omega_max=lambda_Omega_max,
       B0=B0, U=U, P=P)
}
```

---

## 6. 代码：null-corner init

```r
scr_null_corner_init <- function(Sxx, Sxy, Syy, U = 1L, eps = 1e-10) {
  lm <- scr_lambda_max(Sxx, Sxy, Syy, U = U)
  B0 <- lm$B0
  R  <- Syy - B0

  d <- pmax(diag(R), eps)
  Omega0 <- diag(1/d, nrow=ncol(Syy))

  gamma0 <- matrix(0, nrow(Sxx), ncol(Sxy))
  gamma0[lm$U, ] <- solve(Sxx[lm$U,lm$U,drop=FALSE]) %*% Sxy[lm$U,,drop=FALSE] %*% Omega0

  list(gamma=gamma0, Omega=Omega0)
}
```

---

## 7. 简化情形：predictors mean-centered

如果你像 simulation 一样做：

```r
X <- scale(X, center=TRUE, scale=FALSE)
Xc <- cbind(1, X)
```

则经验上 \(S_{XX}[P,U]\approx 0\)，因此：

- \(B_\gamma \approx S_{XY}[P,\cdot]\)
- \(\lambda_{\gamma,\max} \approx 2 \max|S_{XY}[P,\cdot]|\)

---

## 8. Regularization path（二维）

建议对 `(lambda_gamma, lambda_Omega)` 同时做 log-spaced：

```r
path_len <- 60
ratio <- 1e-4

lambda_gamma_seq <- exp(seq(log(lambda_gamma_max),
                            log(lambda_gamma_max * ratio),
                            length.out = path_len))
lambda_Omega_seq <- exp(seq(log(lambda_Omega_max),
                            log(lambda_Omega_max * ratio),
                            length.out = path_len))
```

然后：
- 第 1 个点用 null-corner init（解析解）
- 后续点 warm-start：把上一个点的 `(gamma,Omega)` 作为下一个点的初值；ADMM dual state 也 warm-start（对 PN 尤其重要）。

---

## 9. 给 Codex 的实现提示

- λmax 公式必须与 **scaled loss** 的梯度一致（否则 path 的起点不是真正的 “all zeros corner”）。
- Ω 的 λmax 依赖 \(S_{YY}-B_0\)，不是只看 \(S_{YY}\)。
- 若 `U` 不止 intercept 一行，公式仍成立（直接把 U 当成多个行索引）。
