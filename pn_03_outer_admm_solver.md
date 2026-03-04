# PN-03：PN 子问题的 outer ADMM（Algorithm 3）+ ξ-update Sylvester 方程怎么解

本文件对应你要求的：「PN subproblem outer ADMM 怎么解」。

> 参考：`main_Dec2_revised_gemini.pdf` Section 3.3（Eq (3.8)-(3.12), Algorithm 3）。

---

## 1. outer ADMM 的 splitting

PN 子问题（Eq 3.6）：

\[
\min_{\xi\in D_\delta}\;
\frac12\mathrm{tr}(\xi^\top H \xi \Sigma) + \mathrm{tr}(\xi^\top G) + \Phi(\xi).
\]

做变量分裂：引入 \(\zeta\)（pdf 用 ζ）：

\[
\min_{\xi,\zeta\in D_\delta}
\frac12\mathrm{tr}(\xi^\top H \xi \Sigma) + \mathrm{tr}(\xi^\top G) + \Phi(\zeta)
\quad s.t.\quad \xi=\zeta.
\]

---

## 2. scaled ADMM 迭代（Algorithm 3）

令 scaled dual 为 \(U\)，罚参数为 \(\rho>0\)：

1) **ξ-update**：

\[
\xi^{k+1}=
\arg\min_\xi\;
\frac12\mathrm{tr}(\xi^\top H \xi \Sigma)
+\mathrm{tr}(\xi^\top G)
+\frac{\rho}{2}\|\xi-\zeta^{k}+U^{k}\|_F^2.
\]

2) **ζ-update**（prox）：

\[
\zeta^{k+1}=
\arg\min_{\zeta\in D_\delta}\;
\Phi(\zeta)+\frac{\rho}{2}\|\xi^{k+1}-\zeta+U^k\|_F^2.
\]

这一步是 blockwise prox：
- γ：soft-threshold（阈值 λγ/ρ，U 行不惩罚）
- Ω：PSD+offdiag ℓ1 prox（阈值 λΩ/ρ，需要 inner ADMM，见 PN-04）

3) **dual update**：

\[
U^{k+1}=U^k+\xi^{k+1}-\zeta^{k+1}.
\]

---

## 3. ξ-update 推导到 Sylvester（pdf Eq 3.12）

对 ξ-update 求一阶条件（矩阵微积分）得到：

\[
H\,\xi\,\Sigma + \rho \xi = \rho(\zeta^k - U^k) - G.
\]

这就是 Sylvester 型方程（左边是 “两边夹” 结构）。

---

## 4. Sylvester 的高效解法：特征分解 / Kronecker trick

令：

- \(H = U_H \mathrm{diag}(\lambda_H) U_H^\top\)，维度 \((p_{aug}+q)\times(p_{aug}+q)\)
- \(\Sigma = U_\Sigma \mathrm{diag}(\lambda_\Sigma) U_\Sigma^\top\)，维度 \(q\times q\)

设 RHS：

\[
R := \rho(\zeta^k - U^k) - G.
\]

将方程左右分别乘以 \(U_H^\top\) 与 \(U_\Sigma\)，得到元素级别的解：

\[
Y_{ij}=\frac{(U_H^\top R U_\Sigma)_{ij}}{\lambda_{H,i}\lambda_{\Sigma,j}+\rho},
\]

然后：

\[
\xi^{k+1}=U_H\,Y\,U_\Sigma^\top.
\]

R 伪代码：

```r
eigH <- eigen(symmetrize(H), symmetric=TRUE)
UH <- eigH$vectors
lamH <- pmax(eigH$values, eig_floor)   # clip for stability

eigS <- eigen(symmetrize(Sigma), symmetric=TRUE)
US <- eigS$vectors
lamS <- pmax(eigS$values, 1e-12)

RHS <- rho * (Z - U) - G
Rt  <- crossprod(UH, RHS) %*% US

Den <- outer(lamH, lamS, "*") + rho
Y   <- Rt / Den
Xi_new <- UH %*% Y %*% t(US)
```

> ⚠️ 这里的 `symmetrize(H)` 只能用于**数值对称化**（消除浮点误差）。  
> 前提是你组装 Hessian block 时已经满足 `Hog = t(Hgo)`。如果误写成 `Hog = -t(Hgo)`，`symmetrize()` 会把交叉块平均成 0，导致 Sylvester 方程对应的二次模型被改坏（H 退化为分块对角）。


---

## 5. 收敛判据与 ρ 自适应

标准 ADMM residual：

- primal residual：\(r^k=\xi^k-\zeta^k\)
- dual residual：\(s^k=\rho(\zeta^k-\zeta^{k-1})\)

停止条件：\(\|r^k\|_F\le\varepsilon,\;\|s^k\|_F\le\varepsilon\)

自适应 ρ（常用 heuristic）：

- 若 \(\|r\|>\mu\|s\|\)：ρ ← ρ * τ_inc
- 若 \(\|s\|>\mu\|r\|\)：ρ ← ρ / τ_dec

---

## 6. 给 Codex 的实现提示（重要）

1. ξ-update 要用 **Sigma = Omega^{-1}**（pdf 里叫 Σ(m)），不要误用 Omega。  
2. ζ-update 的 Ω prox 参数必须是 **λΩ/ρ**（见 PN-04 的 bug 分析）。  
3. H 的最小特征值可能为 0（Sxx 奇异），对 `lamH` 做 clipping 能避免除 0。  
4. Warm-start（上一 outer Newton iteration 的 ADMM state）会明显加速。
