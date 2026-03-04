# PN-02：PN 子问题（Eq (3.6)）+ H(m), Σ(m), G(m) 怎么构造

本文件对应你要求的：「PN subproblem 在解决什么、outer ADMM 要解的子问题是什么」。

> 参考：`main_Dec2_revised_gemini.pdf` Section 3.2（Eq (3.5)-(3.7)）。

---

## 1. 从二阶近似到 PN 子问题

在 \(\xi^{(m)}\) 附近：

\[
L(\xi)\approx L(\xi^{(m)}) +
\langle\nabla L(\xi^{(m)}),\xi-\xi^{(m)}\rangle
+\frac12\nabla^2L[\xi^{(m)}](\xi-\xi^{(m)},\xi-\xi^{(m)}).
\]

把常数项丢掉，PN 的“prox Newton step”是解：

\[
s^{(m)}=\arg\min_{\xi\in D_\delta}\;
\frac12\nabla^2L[\xi^{(m)}](\xi,\xi)
+\langle G^{(m)},\xi\rangle
+\Phi(\xi)
\]
其中 \(G^{(m)}\) 是把展开后关于 \(\xi\) 的线性项整理出来得到的矩阵（见下面）。

pdf 把 Hessian quadratic form写成：

\[
\nabla^2L[\xi^{(m)}](d\xi,d\xi)
=
\mathrm{tr}\big(d\xi^\top H^{(m)} d\xi \Sigma^{(m)}\big)
\]
其中 \(\Sigma^{(m)}=(\Omega^{(m)})^{-1}\)。

---

## 2. PN 子问题的最终形式（pdf Eq (3.6)）

把变量 \(\xi\) 视为 stacked matrix（维度 \((p_{aug}+q)\times q\)）后，pdf Eq (3.6) 写为：

\[
s^{(m)}=
\arg\min_{\xi\in D_\delta}\;
\frac12\mathrm{tr}\left(\xi^\top H^{(m)}\xi\Sigma^{(m)}\right)
+\mathrm{tr}\left(\xi^\top G^{(m)}\right)
+\Phi(\xi)
\]

其中：

- \(H^{(m)} := H(\xi^{(m)})\)（见 PG-02 的 Hessian block matrix）
- \(\Sigma^{(m)} := (\Omega^{(m)})^{-1}\)
- \(G^{(m)} := \nabla L(\xi^{(m)}) - H^{(m)}\xi^{(m)}\Sigma^{(m)}\)

> 你代码里的对应项：  
> `Hm <- H_matrix(H_blocks(...))`  
> `Xi <- pack_xi(gamma,Omega)`  
> `Gm <- Grad_stack - Hm %*% Xi %*% Oinv`

---

## 3. 计算 H(m)

给定 \(S_{xx},\gamma^{(m)},\Omega^{(m)}\)，令 \(\Sigma^{(m)}=\Omega^{-1}\)，则

\[
H^{(m)}=
\begin{pmatrix}
2S_{xx} & -2S_{xx}\gamma^{(m)}\Sigma^{(m)}\\
-2(S_{xx}\gamma^{(m)}\Sigma^{(m)})^\top &
\Sigma^{(m)}+2\Sigma^{(m)}(\gamma^{(m)})^\top S_{xx}\gamma^{(m)}\Sigma^{(m)}
\end{pmatrix}
\]

R 实现（同 PG-02）：

```r
Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv)
H  <- symmetrize(H_matrix(Hb))
```

> ⚠️ **重要**：`H_blocks()` 必须保证交叉块满足 `Hog = t(Hgo)`（不要再额外加负号）。  
> 如果误写成 `Hog = -t(Hgo)`，再配合 `symmetrize(H_matrix(Hb))` 会把交叉块平均成 0（完美抵消），使 Hessian 退化为分块对角阵，从而破坏 PN 的二阶耦合信息与下降性质。


---

## 4. 计算 G(m)

先算梯度：

```r
Gg <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv)
Go <- grad_Omega(gamma, Omega, Sxx, Syy, Oinv = Oinv)
Grad_stack <- pack_xi(Gg, Go)
```

再算：

```r
Xi <- pack_xi(gamma, Omega)
Gm <- Grad_stack - H %*% Xi %*% Oinv
```

这就对应 pdf 的：

\[
G^{(m)}=\nabla L(\xi^{(m)})-H^{(m)}\xi^{(m)}\Sigma^{(m)}.
\]

---

## 5. PN 方向与 Newton decrement

解出 \(s^{(m)}\) 后：

\[
d^{(m)} = s^{(m)} - \xi^{(m)}
\]

Newton decrement（用于 damped step）：

\[
\nu^{(m)} = \|d^{(m)}\|_{\xi^{(m)}}.
\]

---

## 6. 给 Codex 的实现提示

- PN 子问题是凸的（在 \(D_\delta\) 内），但规模大：必须用 ADMM / 结构化求解。
- `Gm` 的公式非常容易写错：是 **Grad_stack - H Xi Σ**，不是别的组合。
- stack 的 `Xi` 维度是 `(p_aug+q)×q`，确保矩阵乘法合法：
  - `H` 是 `(p_aug+q)×(p_aug+q)`
  - `Xi` 是 `(p_aug+q)×q`
  - `Oinv` 是 `q×q`
