# PG-02：Gradient / Hessian / Local norm（self-concordant 几何）

本文件对应你要求的 **第 2 个 md**：  
- 梯度公式（pdf Eq (2.3)）  
- Hessian block matrix（pdf Eq (2.4)）  
- local norm / Newton decrement（pdf Eq (2.7)，以及 PN 用到的 ν）  
- R 实现时的 stack 方式与维度检查

---

## 1. 维度与 stack 约定

- `gamma ∈ R^{p_aug×q}`
- `Omega ∈ S_{++}^q`
- `Oinv := Omega^{-1} ∈ S_{++}^q`
- `Sxx ∈ R^{p_aug×p_aug}`, `Sxy ∈ R^{p_aug×q}`, `Syy ∈ R^{q×q}`

pdf 把 \(\xi=(\gamma,\Omega)\) 也写成 stacked matrix：

\[
\xi \;\equiv\;
\begin{pmatrix}
\gamma\\
\Omega
\end{pmatrix}\in\mathbb{R}^{(p_{aug}+q)\times q}.
\]

在 R 里建议用：

```r
pack_xi   <- function(gamma, Omega) rbind(gamma, Omega)
unpack_xi <- function(Xi, p_aug) list(
  gamma = Xi[1:p_aug, , drop=FALSE],
  Omega = Xi[(p_aug+1):nrow(Xi), , drop=FALSE]
)
```

---

## 2. 梯度（pdf Eq (2.3)）

scaled loss：

\[
L(\gamma,\Omega)= -\log\det(\Omega)
+ \mathrm{tr}(\Omega^{-1}\gamma^\top S_{xx}\gamma)
-2\,\mathrm{tr}(S_{xy}^\top\gamma)
+ \mathrm{tr}(\Omega S_{yy})
\]

对应梯度：

\[
\nabla_\gamma L(\xi)
= -2 S_{xy} + 2 S_{xx}\gamma\Omega^{-1}
\]

\[
\nabla_\Omega L(\xi)
= S_{yy} - \Omega^{-1} - \Omega^{-1}\gamma^\top S_{xx}\gamma \Omega^{-1}.
\]

R 实现：

```r
grad_gamma <- function(gamma, Omega, Sxx, Sxy, Oinv = NULL) {
  if (is.null(Oinv)) Oinv <- safe_inv_spd(Omega)
  2 * (Sxx %*% gamma %*% Oinv - Sxy)
}

grad_Omega <- function(gamma, Omega, Sxx, Syy, Oinv = NULL) {
  if (is.null(Oinv)) Oinv <- safe_inv_spd(Omega)
  M <- t(gamma) %*% Sxx %*% gamma
  Syy - Oinv - Oinv %*% M %*% Oinv
}
```

> ✅ 这和你 `2.18_simulation.R` 的 `grad_gamma/grad_Omega` 一致（它们实现的是 scaled loss 的梯度）。

---

## 3. Hessian block matrix H(ξ)（pdf Eq (2.4)）

pdf 给出 Hessian quadratic form：

\[
\nabla^2 L[\xi](d\xi,d\xi)
=\mathrm{tr}\big(d\xi^\top H(\xi)\, d\xi \,\Omega^{-1}\big),
\]

其中 block matrix \(H(\xi)\in\mathbb{R}^{(p_{aug}+q)\times(p_{aug}+q)}\) 是：

\[
H(\xi)=
\begin{pmatrix}
2S_{xx} & -2S_{xx}\gamma\Omega^{-1}\\
-2(\,S_{xx}\gamma\Omega^{-1}\,)^\top &
\Omega^{-1}+2\Omega^{-1}\gamma^\top S_{xx}\gamma\Omega^{-1}
\end{pmatrix}.
\]

R 里建议按 block 组装：

```r
H_blocks <- function(gamma, Omega, Sxx, Oinv = NULL) {
  if (is.null(Oinv)) Oinv <- safe_inv_spd(Omega)
  beta_eff <- gamma %*% Oinv  # β = γ Ω^{-1}

  Hgg <- 2 * Sxx
  Hgo <- -2 * Sxx %*% beta_eff
  Hog <- t(Hgo)
  Hoo <- Oinv + 2 * t(beta_eff) %*% Sxx %*% beta_eff

  list(Hgg=Hgg, Hgo=Hgo, Hog=Hog, Hoo=Hoo)
}

H_matrix <- function(Hb) {
  rbind(cbind(Hb$Hgg, Hb$Hgo),
        cbind(Hb$Hog, Hb$Hoo))
}
```

### 3.1 关键实现警告：`Hog` 只能取 `t(Hgo)`（不要多一个负号）

在你的 `2.18_simulation.R` 里曾经出现过：

```r
Hgo <- -2 * Sxx %*% (gamma %*% Oinv)
Hog <- - t(Hgo)      # ❌ bug
```

这是一个**符号错误**：因为 `Hgo` 本身已经带了负号（来自公式 \(H_{\gamma\Omega}=-2S_{xx}\gamma\Omega^{-1}\)），再加一个负号会把下左块的符号翻转。随后如果你再做：

```r
H <- symmetrize(H_matrix(Hb))
```

则交叉块会被平均成：

\[
\tfrac12\left(H_{\gamma\Omega} + H_{\Omega\gamma}^\top\right)
=\tfrac12\left(Hgo + (-Hgo)\right)=0,
\]

导致 Hessian 退化成**分块对角阵**，直接丢失 \(\gamma\) 与 \(\Omega\) 的耦合二阶信息。

后果：
- **PN**：Newton 方向会明显跑偏（因为二次模型少了交叉二阶导）；
- **PG/PN**：local norm（Newton decrement）里的交叉项被抵消，\(\nu\) 被系统性低估，显式 damping 步长会偏大。

✅ 正确写法只有一种：

```r
Hog <- t(Hgo)        # ✅ correct
```

> 说明：`symmetrize(H)` 可以保留（用于消除数值误差），但它不能替代你在组装 block 时保持 `Hog=t(Hgo)`。


---

## 4. local norm 与 Newton decrement（pdf Eq (2.7)）

对增量 \(d\xi=(d\gamma,d\Omega)\)：

\[
\|d\xi\|_\xi :=
\sqrt{\nabla^2 L[\xi](d\xi,d\xi)}
=
\sqrt{\langle d\xi,\; H(\xi)d\xi\Omega^{-1}\rangle}.
\]

R 实现（stacked 版本）：

```r
local_norm <- function(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = NULL) {
  if (is.null(Oinv)) Oinv <- safe_inv_spd(Omega)
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv)
  H  <- H_matrix(Hb)

  dXi <- rbind(d_gamma, d_Omega)          # (p_aug+q)×q
  val <- sum(dXi * (H %*% dXi %*% Oinv))  # Frobenius inner product
  sqrt(max(val, 0))
}
```

---

## 5. Codex 实现细节（常见坑）

1. **对称性**：任何 `Omega` 更新后都 `symmetrize(Omega)`；H 也建议 `symmetrize(H)`.
2. **Oinv 只用 SPD inverse**：用 `safe_inv_spd()`（Cholesky + jitter 或 eigen fallback）。
3. **维度确认**：  
   - `Sxx %*% gamma` 是 `p_aug×q`  
   - `gamma %*% Oinv` 是 `p_aug×q`  
   - `t(gamma)%*%Sxx%*%gamma` 是 `q×q`
4. **local norm 是半范数**：高维时 `Sxx` 可能奇异，因此 \(\|\cdot\|_\xi\) 可能是 semi-norm；实现中要允许 `val≈0`.
