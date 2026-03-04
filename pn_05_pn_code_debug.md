# PN-05：对 `2.18_simulation.R` 里 PN 实现的对照检查（为什么 loss 不降？）

你说的现象：「理论上 loss 是下降的，但 simulation 里 loss 没下降」——  
从你 `2.18_simulation.R` 的实现看，确实存在 **至少三个非常关键的实现偏差**，足以解释 loss 不单调。

本文件把 pdf（Algorithm 2/3 + Appendix B）逐条对照你的代码，并指出需要修的地方。

---

## A. pdf 的“理论保证”依赖哪些条件？

在 `main_Dec2_revised_gemini.pdf` 的 Prox-Newton（Algorithm 2）里，常见的单调下降结论通常需要：

1) 子问题 \(s^{(m)}\)（Eq 3.6）**解得足够准确**（或满足一定 inexact 条件）  
2) 用 damped 步长
\[
\alpha^{(m)}=\frac{1}{1+\nu^{(m)}}
\]
（Eq 3.7），其中 \(\nu^{(m)}\) 是 Newton decrement  
3) Ω 始终在 \(D_\delta\) 内（Ω ⪰ δI）

如果子问题没按正确的 penalty 缩放求解，或者 Hessian / local norm 的实现本身有符号错误（导致 ν 计算失真），再加上 line-search/acceptance 允许 loss 增加，单调性就没了。

---

## B. 关键 bug 0：Hessian 交叉块符号写反，`symmetrize()` 会把交叉二阶导“完美清零”

这点是**毁灭性**的：它会把 Hessian 的非对角块（\(\gamma\) 与 \(\Omega\) 的交叉二阶导）直接抹掉，导致 PN 的二次模型变成错误的分块对角近似。

### B.1 你代码里的写法（确实有 bug）

在 `2.18_simulation.R` 的 `H_blocks()`（约第 400 行）：

```r
Hgo <- -2 * Sxx %*% beta_eff
Hog <- - t(Hgo)      # ❌ bug
```

但按 pdf Eq (2.4)，交叉块必须满足 Hessian 对称：

- \(H_{\gamma\Omega} = -2\,\hat\Sigma_{XX}\,\gamma\,\Omega^{-1}\)（本身就带负号）
- \(H_{\Omega\gamma} = (H_{\gamma\Omega})^\top\)

所以正确实现只能是：

```r
Hog <- t(Hgo)        # ✅ correct
```

### B.2 为什么会“清零”（你现在看到的 loss 不降，很大一部分来自这里）

你在 PN 子问题里还会做：

```r
Hm <- symmetrize(H_matrix(Hb))
```

当 `Hog = -t(Hgo)` 时，

\[
H_{12}^{sym}=\tfrac12\left(Hgo + t(Hog)\right)
=\tfrac12\left(Hgo + (-Hgo)\right)=0,
\]

导致 Hessian **退化成分块对角阵**，交叉二阶导被完美抵消。

### B.3 后果

- **PN**：二次模型缺少 \(\gamma\)-\(\Omega\) 耦合项，Newton 方向会跑偏，damped 步长 \(\alpha=1/(1+\nu)\) 的下降保证也会随之失效；  
- **PG/PN**：local norm（Newton decrement）中的交叉项被抵消，\(\nu\) 被系统性低估 → 步长偏大；  
- 叠加你后面 Armijo 的符号 bug + forced accept，loss 极容易上升。

✅ 修复后，你再去看 PN 的 loss 曲线，通常会立刻“正常很多”。

---

## C. 关键 bug 1：PN 子问题 outer ADMM 的 Ω-prox 缩放错了（λΩ/ρ 被写成 λΩ/ρ²）

### C.1 pdf 需要的 Ω-block prox（Algorithm 3 的 ζ-update）

在 outer ADMM（Algorithm 3）里，ζ-update 的 Ω-block 是：

\[
\Omega^{k+1}
=\arg\min_{\Omega\succeq\delta I}
\frac12\|\Omega - V\|_F^2 + \frac{\lambda_\Omega}{\rho}\|\Omega\|_{1,off}
\]

也就是说，Ω-prox 里的 penalty 系数必须是：

\[
\tau_\lambda = \lambda_\Omega/\rho.
\]

### C.2 你代码的写法

在 `pn_subproblem_admm()` 里（约第 770 行附近）：

```r
Zo <- prox_psd_offdiag_l1(
  Xi_parts$Omega,
  tau = lambda_Omega / rho_iter,
  mu = mu,
  rho = rho_iter,
  ...
)
```

而你的 `prox_psd_offdiag_l1()` 内部又用到了 `tau/rho`（阈值）  
→ 实际阈值变成 **(lambda_Omega / rho_iter) / rho_iter = lambda_Omega / rho_iter^2**。

这会导致：
- Ω 的 offdiag shrink 太弱（Ω 更密）
- PN 子问题并没有在正确的目标下求解
- 方向 d 不再对应 pdf 的 Newton step，自然就可能 loss 上升

✅ 修复（两种等价方案，选一种就行）：

- 方案 1（推荐，最不容易错）  
  把 Ω-prox 写成固定解：
  \[
  \min \frac12\|\Omega-V\|^2+\tau_\lambda\|\Omega\|_{1,off}
  \]
  则调用时只传 `tau_lambda=lambda_Omega/rho`，并且 **不要再传 rho**。

- 方案 2  
  若你想 Ω-prox 解：
  \[
  \min \frac{\rho}{2}\|\Omega-V\|^2+\lambda_\Omega\|\Omega\|_{1,off}
  \]
  则调用传 `(lambda_Omega, rho)`，prox 内阈值用 `lambda_Omega/rho`，并且 primal/dual gap 也要相应缩放。

我在新的模块化实现里采用了方案 1：见 `scr_code/scr_pn_solver.R` 与 `scr_code/scr_pn_inner_omega_admm.R`。

---

## D. 关键 bug 2：你的 Armijo 条件把梯度项符号写反了，导致“允许 loss 增加”

你在 `pn_sparse_mvreg()` 的 line-search 里写了：

```r
armijo_rhs_try <- armijo_c * (pen_try - pen_prev - alpha_try * grad_dir_ip)
armijo_lhs_try <- loss_try - loss_prev
accept if armijo_lhs_try <= armijo_rhs_try
```

其中 `grad_dir_ip = <∇smooth, d>`。

标准 Armijo（smooth 情况）是：

\[
F(x+\alpha d)\le F(x)+c\alpha \langle\nabla f(x),d\rangle
\]

对 composite（f+g）常用改写是：

\[
F(x+\alpha d)\le F(x)+c\left(\alpha \langle\nabla f(x),d\rangle + g(x+\alpha d)-g(x)\right)
\]

也就是说右侧应该是：

```r
armijo_rhs = c * ( alpha*grad_dir_ip + (pen_try - pen_prev) )
```

而你写成了 `- alpha*grad_dir_ip`。  
由于 descent 方向通常 `grad_dir_ip < 0`，所以 `-alpha*grad_dir_ip > 0`，右侧变成正的，结果：

- 即使 `loss_try - loss_prev > 0`（loss 上升），也可能满足不等式  
- 于是 line-search 会“通过”，loss 就不单调了

✅ 修复：把 rhs 改成 `+ alpha_try * grad_dir_ip`（不是减号）。

---

## E. 关键策略问题：line-search 失败时你“强行接受 best α”，即使 best 仍然增大

你代码在失败时：

```r
if (!line_search_ok) {
  # forced accept best loss among tried alphas (even if > loss_prev)
  gamma <- best_gamma_try
  Omega <- best_Omega_try
}
```

这会直接破坏 monotone。  
如果你希望“理论一致：loss 必须下降”，建议改成：

- 若所有尝试的 α 都无法使 loss 下降：**不更新**，并返回/停止/或重新解子问题（更高精度）

我在 `scr_code/scr_pn_solver.R` 里做的是：
- 默认 `monotone=TRUE`
- 若 backtracking 到 `min_alpha` 仍不降，则 stop_reason = `"no_descent"`，保持当前点不动

---

## F. 其他可能导致 PN 不稳的点（次要，但值得注意）

1) **子问题精度不够**：  
`pn_inner_admm_tol_fixed=1e-7` 但 `max_admm=100`，可能经常达不到 tol。  
建议：  
- 用相对 residual tol；或  
- 外层迭代早期 tol 宽、后期逐渐变严；或  
- 提高 max_admm / warm-start dual state

2) **λmax 的 Ω 部分用错**：会影响路径起点（见 PG-03）。

---

## G. 我提供的修正版（按 pdf）在哪里？

- `scr_code/scr_pn_solver.R`：Algorithm 2 + Algorithm 3（外层 ADMM）  
- `scr_code/scr_pn_inner_omega_admm.R`：Appendix B Ω-prox inner ADMM  
- 关键修复点：
  - Hessian 交叉块：强制 `Hog = t(Hgo)`，并且 `symmetrize(H)` 仅用于数值对称化（避免交叉项被清零）
  - Ω-prox 缩放：严格用 `tau_lambda=lambda_Omega/rho`（不再出现 /rho²）
  - 步长：按 pdf 用显式 damping `alpha=1/(1+nu)`；不再用 Armijo（避免符号/重复计 penalty 的坑）
  - 默认 monotone safeguard：若 loss 不降则继续 halve α；若仍不降则停止并不更新（不 forced accept）

---
