# PG-06：对 `2.18_simulation.R` 里 PG 实现的 code review（潜在错误/改进点）

本文件用于回答你要求的「还要 check 一下这个 PG 算法有啥错误吗」以及「补充 PG 的东西」。

结论先说：  
你 `pg_sparse_mvreg()` 的 **主框架**（梯度、prox、ν/β、L backtracking、α(m)）和 pdf 的 Algorithm 1 基本一致；  
但有一些 **容易造成数值异常或路径起点不准** 的点建议修。

---

## 1. 与 pdf 完全对齐的点（✅）

1) 使用的是 **scaled loss L**（对应 pdf Eq (2.2)），并且 `grad_gamma/grad_Omega` 的实现与 pdf Eq (2.3) 一致。  
但要注意：你代码里的 Hessian 交叉块写法有符号 bug（`Hog <- -t(Hgo)`），会让 `local_norm` 不再对应 pdf Eq (2.7)（见 2.1）。

2) Backtracking 判别：
\[
\nu^2/\beta^2+\nu>1 \Rightarrow L\leftarrow 2L
\]
你的实现一致。

3) 步长：
\[
\alpha=\frac{\beta^2}{\nu(\nu+\beta^2)}
\]
你的实现一致（但缺少 ν=0 的保护，见下）。

---

## 2. 可能的“真 bug”/不稳点（⚠️）

### 2.1 Hessian 交叉块符号 bug：`Hog <- -t(Hgo)` 会把交叉二阶导清零

你在 `H_blocks()` 里写了：

```r
Hgo <- -2 * Sxx %*% beta_eff
Hog <- - t(Hgo)     # ❌ bug
```

但按 pdf Eq (2.4) 必须满足 Hessian 对称：**`Hog = t(Hgo)`**。  
因为 `Hgo` 本身已经带负号（\(-2S_{xx}\gamma\Omega^{-1}\)），再多一个负号会把交叉块整体翻转。

更糟的是，你在 PN 子问题里还会做：

```r
H <- symmetrize(H_matrix(Hb))
```

当 `Hog=-t(Hgo)` 时，`symmetrize()` 会把交叉块平均成：

\[
\tfrac12\left(Hgo + t(Hog)\right)=\tfrac12\left(Hgo + (-Hgo)\right)=0,
\]

于是 Hessian 退化成**分块对角阵**，完全丢失 \(\gamma\)-\(\Omega\) 的耦合二阶导。

后果：
- **PG**：`local_norm()` 的交叉项被抵消 → \(\nu\) 被系统性低估 → \(\alpha(m)\) 偏大、backtracking 判别偏松；  
- **PN**：更严重，Newton 方向跑偏，loss 很容易不降甚至上升。

✅ 修复（必须改）：

```r
Hog <- t(Hgo)        # ✅ correct
```

并且 `symmetrize(H)` 只作为数值对称化手段，不能用来“修复”块结构错误。

### 2.2 α(m) 在 ν=0 时会 NaN

你代码直接：

```r
alpha_m <- (beta_m^2) / (lambda_m * (lambda_m + beta_m^2))
```

若 `lambda_m == 0`（也就是 local_norm 近似 0），会出现 0/0。  
建议改成（我在新的 `scr_pg_solver.R` 里已处理）：

```r
if (lambda_m <= .Machine$double.eps) alpha <- 1 else alpha <- ...
```

### 2.3 λΩ,max 的公式（path 起点）与你 pdf Lemma 3.1 不一致

你用了：

```r
lambda_Omega_0 <- 2 * max(abs(offdiag(Syy)))
```

pdf Lemma 3.1 是：

\[
\lambda_{\Omega,\max}=\|(S_{YY}-B_0)_{\mathrm{off}}\|_{\max}
\]
（无 2，并且要减去 \(B_0=S_{XY}[U]^T S_{XX}[U,U]^{-1} S_{XY}[U]\)）

这会导致：
- path 的“第一点”不一定真的是 **offdiag(Ω)=0** 的角点
- warm-start 时可能出现不必要的振荡

建议按 PG-03 修改。

### 2.4 PG 的 stopping 没有使用 gm_norm（只是算了）

你算了 gradient mapping norm：

```r
gm_gamma <- Lm * (gamma_old - s_gamma_m)
gm_Omega <- Lm * (Omega_old - s_Omega_m)
gm_norm  <- ...
```

但没有用 `gm_norm <= tol_gm` 作为收敛判断。  
如果 oracle gap 不用，你只用 objective stabilized，可能会：
- 停得太早（objective变化小但离最优还远）
- 或停得太晚（objective变化小但 gm_norm 已经很小）

建议：把 `gm_norm <= tol_gm` 纳入 stopping（特别是没有 oracle_value 时）。

---

## 3. “不是 bug，但建议改”的工程点（🛠️）

### 3.1 backtracking loop 里重复计算梯度

`Gg_m/Go_m` 在 backtracking 的 repeat 中反复算，但当前点没变，可以提到外面节省时间。

### 3.2 全局变量 prox_ctrl

你在 `pg_sparse_mvreg()` 里直接用全局 `prox_ctrl`。  
建议把 prox 控制项作为函数参数传入（避免不同实验互相污染）。

### 3.3 Omega 维度/对称性

你每次 `Vo_m <- symmetrize(Omega - Go/L)` 很好；  
建议更新 `Omega_candidate` 后再 `symmetrize()` 并做一次 `proj_psd_shifted`（你已经做 `safe_inv_psd` 检查了，但投影更稳）。

---

## 4. Simulation 里一个“命名/真值一致性”问题（会影响评估）

你 simulation 里：

```r
gamma_true <- generate_sparse_gamma(p_aug, q, ...)
Beta_true  <- gamma_true[-1, ]
Y <- X %*% Beta_true + noise(Sigma_true)
```

这里 `gamma_true` 实际上更像 `beta_true_aug`（回归系数），而不是 reparam 的 γ=βΩ。  
这不会影响 PG/PN 优化是否收敛，但会影响你后续对 “gamma 真值” 的比较。

建议按 PG-01 的方式统一：
- 先生成 `beta_true` 和 `Omega_true`
- 再令 `gamma_true = beta_true %*% Omega_true`
- 用 `Y = Xc %*% beta_true + noise`

---

## 5. 我在新的模块化 R 文件里做了哪些修正？

见 `scr_code/scr_pg_solver.R`：
- 修复 Hessian 交叉块符号：`Hog <- t(Hgo)`（避免交叉项被 `symmetrize()` 抵消，导致 local norm 失真）  
- ν≈0 时 α=1 的保护  
- Ω 投影到 `Omega ⪰ delta I` 的数值稳定保护  
- 把 prox / stop / settings 都集中管理（避免全局变量）

同时 λmax 和 null-corner init 见 `scr_code/scr_settings.R`（按 pdf Lemma 3.1）。

---
