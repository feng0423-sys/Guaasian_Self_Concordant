# PG-05：prox 子问题的闭式解（γ）+ Ω-prox 的 inner ADMM（Appendix B）

本文件对应你要求的 **第 5 个 md**：  
- 内层 prox：γ 的 closed form  
- Ω 的 prox（off-diagonal ℓ1 + PSD 约束）用 inner ADMM 怎么算  
- 这套 Ω-prox 既用于 PG，也用于 PN 的外层 ADMM 的 z-update

> 参考：`main_Dec2_revised_gemini.pdf` Appendix B（Ω-block proximal operator）  
> 以及 Algorithm 1 的 prox 子问题（Eq (3.3)）。

---

## 1. PG 的 prox 子问题为什么可分块？

PG 子问题（Eq 3.3）：

\[
s^{(m)}=\arg\min_{\xi\in D_\delta}
\langle\nabla L(\xi^{(m)}), \xi-\xi^{(m)}\rangle
+\frac{L(m)}{2}\|\xi-\xi^{(m)}\|_F^2
+\Phi(\xi)
\]

其中：
\[
\Phi(\gamma,\Omega)=\lambda_\gamma\|\gamma_P\|_1+\lambda_\Omega\|\Omega\|_{1,\mathrm{off}}
\]

由于 Frobenius 二次项是分块的，并且 penalty 也对 (γ,Ω) 分离，所以 prox 子问题分解为：

- γ-block：标准 soft-threshold（除 U 行不惩罚）
- Ω-block：带 PSD 约束的 off-diagonal soft-threshold（需 ADMM）

---

## 2. γ-block prox：闭式解（soft-threshold）

设

\[
V_\gamma = \gamma^{(m)} - \frac{1}{L(m)}\nabla_\gamma L(\xi^{(m)}).
\]

则

\[
s_\gamma
=\arg\min_\gamma \frac{1}{2}\|\gamma-V_\gamma\|_F^2 + \frac{\lambda_\gamma}{L(m)}\|\gamma_P\|_1
\]

闭式解就是对 penalized rows \(P\) 做 soft-threshold：

\[
[s_\gamma]_{P,\cdot} = \mathrm{SoftThresh}\!\left([V_\gamma]_{P,\cdot},\; \lambda_\gamma/L(m)\right),
\qquad
[s_\gamma]_{U,\cdot} = [V_\gamma]_{U,\cdot}.
\]

R 代码：

```r
scr_prox_gamma_l1 <- function(Vg, tau, U = 1L) {
  p_aug <- nrow(Vg)
  U <- sort(unique(as.integer(U)))
  P <- setdiff(seq_len(p_aug), U)

  G <- Vg
  if (length(P) > 0) G[P, ] <- soft_thresh(Vg[P,,drop=FALSE], tau)
  if (length(U) > 0) G[U, ] <- Vg[U,,drop=FALSE]
  G
}
```

---

## 3. Ω-block prox：带 PSD 约束的 off-diagonal ℓ1

PG 里：

\[
V_\Omega = \Omega^{(m)} - \frac{1}{L(m)}\nabla_\Omega L(\xi^{(m)}),
\]

Ω-prox 要解：

\[
\min_{\Omega \succeq \delta I}\;
\frac{1}{2}\|\Omega - V_\Omega\|_F^2
+\frac{\lambda_\Omega}{L(m)}\|\Omega\|_{1,\mathrm{off}}.
\]

这是 Appendix B 的主问题（把 \(\tau\lambda_\Omega\) 记为 `tau_lambda`）：

\[
\min_{\Omega\succeq \delta I}\;
\frac{1}{2}\|\Omega - V\|_F^2
+ \tau_\lambda \|\Omega\|_{1,\mathrm{off}}.
\]

---

## 4. Appendix B：inner ADMM 分裂与更新公式

引入 splitting variable \(K\)：

\[
\min_{\Omega,K}\;
\frac{1}{2}\|\Omega - V\|_F^2
+\tau_\lambda\|\Omega\|_{1,\mathrm{off}}
+\iota_{\{K\succeq\delta I\}}(K)
\quad s.t.\quad \Omega=K.
\]

scaled dual \(A\)，罚参数 \(\mu>0\)：

### (i) Ω-update（闭式）

\[
\Omega^{(\ell+1)}
=\arg\min_{\Omega}
\frac{1}{2}\|\Omega - V\|_F^2
+\tau_\lambda\|\Omega\|_{1,\mathrm{off}}
+\frac{\mu}{2}\|\Omega - K^{(\ell)} + A^{(\ell)}\|_F^2.
\]

闭式解（利用 soft-threshold 的齐次性）：

\[
\Omega^{(\ell+1)}
= \frac{1}{1+\mu}\,
\mathrm{SoftThresh}_{off}\!\Big(V+\mu(K^{(\ell)}-A^{(\ell)}),\;\tau_\lambda\Big).
\]

> 等价写法：  
> \(\Omega^{(\ell+1)} = \mathrm{SoftThresh}_{off}((V+\mu(K-A))/(1+\mu),\;\tau_\lambda/(1+\mu))\).

### (ii) K-update（shifted PSD projection）

\[
K^{(\ell+1)} = \Pi_{\Omega\succeq\delta I}(\Omega^{(\ell+1)} + A^{(\ell)}).
\]

实现是 eigenvalue clipping：

如果 \(W = U\mathrm{diag}(\lambda)U^\top\)，则

\[
\Pi_{\Omega\succeq\delta I}(W)=U\mathrm{diag}(\max\{\lambda,\delta\})U^\top.
\]

### (iii) Dual update

\[
A^{(\ell+1)} = A^{(\ell)} + \Omega^{(\ell+1)} - K^{(\ell+1)}.
\]

---

## 5. stopping rule：primal/dual residual 或 duality gap（Appendix B）

Appendix B 也给了 duality gap：

- primal：
\[
p(\Omega)=\frac{1}{2}\|\Omega - V\|_F^2 + \tau_\lambda\|\Omega\|_{1,\mathrm{off}}.
\]

- dual feasible set：
\[
D_{\tau_\lambda}=\{Y\in S^q:\;Y_{ii}=0,\;|Y_{ij}|\le \tau_\lambda\;\forall i\ne j\}.
\]

- dual objective（闭式）：
\[
K^\star(Y)=\Pi_{S_+^q}(V-Y),\quad
d(Y)=\langle Y,K^\star(Y)\rangle+\frac{1}{2}\|K^\star(Y)-V\|_F^2.
\]

- gap：
\[
gap(\Omega,Y)=p(\Omega)-d(Y)\ge 0,
\quad
\frac{gap}{1+|p|}\le\varepsilon_{gap} \Rightarrow \text{stop}.
\]

---

## 6. R 代码骨架（建议直接照 Appendix B 实现）

```r
scr_prox_omega_offdiag_l1_psd <- function(V, tau_lambda, delta=1e-8,
                                         mu=1, maxit=2000, gap_tol=1e-7,
                                         warm_state=NULL) {
  # warm_state: list(Omega,K,A)
  # returns: list(Omega, state, iters, converged)
}
```

> ✅ 关键点：**tau_lambda 已经是 penalty 系数**（比如 λΩ/L 或 λΩ/ρ），Ω-prox 内部不要再额外除以 ρ。  
> 你旧代码在 PN 的 z-update 里把 `tau=lambda_Omega/rho` 同时又传了 `rho=rho`，导致阈值变成 λΩ/ρ²（明显过小）。

---

## 7. 给 Codex 的实现提示

1. soft-threshold 仅对 off-diagonal：diag 不惩罚、保持不变。  
2. PSD 约束用 eigen clipping，`delta` 就是 safeguard。  
3. inner ADMM 强烈建议 warm-start（把上一轮的 (Omega,K,A) 接着用），能显著省迭代。
