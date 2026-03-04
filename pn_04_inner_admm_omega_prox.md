# PN-04：Ω-prox 的 inner ADMM（Appendix B）—— PN/PG 共用

本文件对应你要求的：「inner ADMM 怎么解」部分，具体是 **Ω-block proximal operator**。

> 参考：`main_Dec2_revised_gemini.pdf` Appendix B（尤其 Eq (B.1)-(B.7)）。

---

## 1. Ω-prox 出现在哪里？

- **PG**：每次 prox step 都要解一次 Ω-prox（参数 `tau_lambda = lambda_Omega / L(m)`）
- **PN 子问题 outer ADMM**：每次 ζ-update 的 Ω-block 又要解一次 Ω-prox  
  （参数 `tau_lambda = lambda_Omega / rho`）

因此 Ω-prox 的实现如果缩放写错，会同时污染 PG/PN。

---

## 2. Ω-prox 的数学形式（带 safeguard δ）

给定对称矩阵 \(V\)，解：

\[
\min_{\Omega\succeq \delta I}\;
\frac12\|\Omega - V\|_F^2
+ \tau_\lambda \|\Omega\|_{1,\mathrm{off}}.
\]

其中：
- \(\tau_\lambda\) 是 penalty 系数（例如 λΩ/L 或 λΩ/ρ）
- \(\|\Omega\|_{1,\mathrm{off}}=\sum_{i\ne j}|\Omega_{ij}|\)
- \(\delta>0\) 用于保证 SPD（数值稳定）

---

## 3. Appendix B 的 inner ADMM（带 splitting K）

分裂：

\[
\min_{\Omega,K}\;
\frac12\|\Omega - V\|_F^2
+\tau_\lambda \|\Omega\|_{1,\mathrm{off}}
+\iota_{\{K\succeq\delta I\}}(K)
\quad \text{s.t.}\quad \Omega=K.
\]

scaled dual \(A\)，罚参数 \(\mu\)：

### 3.1 Ω-update（闭式）

\[
\Omega^{(\ell+1)}
= \frac{1}{1+\mu}\,
\mathrm{SoftThresh}_{off}\!\left(
V+\mu(K^{(\ell)}-A^{(\ell)}),
\tau_\lambda
\right).
\]

### 3.2 K-update（shifted PSD projection）

\[
K^{(\ell+1)}=\Pi_{\Omega\succeq\delta I}\big(\Omega^{(\ell+1)}+A^{(\ell)}\big).
\]

实现：对称化 + eigenvalue clipping below δ。

### 3.3 Dual update

\[
A^{(\ell+1)} = A^{(\ell)} + \Omega^{(\ell+1)} - K^{(\ell+1)}.
\]

---

## 4. stopping：dual residual / primal residual 或 duality gap

### 4.1 residual（工程最常用）

- \(r^{(\ell)}=\Omega^{(\ell)}-K^{(\ell)}\)
- \(s^{(\ell)}=\mu(K^{(\ell)}-K^{(\ell-1)})\)

stop if \(\|r\|_F\le\varepsilon\) 且 \(\|s\|_F\le\varepsilon\)。

### 4.2 duality gap（pdf 推荐）

dual feasible set：

\[
D_{\tau_\lambda}=\{Y\in S^q:\;Y_{ii}=0,\;|Y_{ij}|\le\tau_\lambda\}.
\]

dual closed form：

\[
K^\star(Y)=\Pi_{\Omega\succeq\delta I}(V-Y),
\qquad
d(Y)=\langle Y,K^\star(Y)\rangle+\frac12\|K^\star(Y)-V\|_F^2.
\]

primal：

\[
p(\Omega)=\frac12\|\Omega-V\|_F^2+\tau_\lambda\|\Omega\|_{1,off}.
\]

gap = p - d ≥ 0，若 gap/(1+|p|) ≤ eps 则停。

---

## 5. PN 里最常见的缩放错误（会导致 loss 不降）

在 PN 的 outer ADMM ζ-update 中，Ω-block prox 应该是：

\[
\min_{\Omega\succeq\delta I}\;
\frac12\|\Omega-V\|^2 + \frac{\lambda_\Omega}{\rho}\|\Omega\|_{1,off}.
\]

也就是说 **tau_lambda = lambda_Omega / rho**。

工程上实现有两种等价写法：

- 写法 A：prox 函数固定解 `0.5||Ω-V||^2 + tau_lambda ||Ω||_1`  
  → 调用时 `tau_lambda=lambda_Omega/rho`

- 写法 B：prox 函数解 ` (rho/2)||Ω-V||^2 + lambda_Omega||Ω||_1`  
  → 调用时传 rho，同时 prox 内部要把阈值设成 `lambda_Omega/rho`

你旧代码在 PN 里属于 “写法 A + 又额外除一次 rho”，结果阈值变 λΩ/ρ²，这会让 Ω 变得过密、并可能破坏 PN 方向性质，从而出现 loss 不下降。

---

## 6. 与本项目的模块化实现对照

我已经把 Appendix B 的 inner ADMM 生成在：

- `scr_code/scr_pn_inner_omega_admm.R`

核心函数：

```r
scr_prox_omega_offdiag_l1_psd(V, tau_lambda,
                             delta=..., mu=..., maxit=..., gap_tol=...,
                             warm_state=...)
```

默认用 duality gap（也支持 residual stopping）。

---

## 7. Codex 实现提示

1) Ω-prox 是整个算法最耗时部分：必须 warm-start（state = {Omega,K,A}）。  
2) 任何时候都 `symmetrize()`，并投影到 \(\Omega\succeq\delta I\)。  
3) 若 q 较大，eigen decomposition 是瓶颈；可以考虑只在需要时做或用更快的近似（后续再优化）。
