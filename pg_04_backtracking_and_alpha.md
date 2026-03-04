# PG-04：L(m) backtracking + α(m)（Eq (3.4)）+ 停止准则（Algorithm 1）

本文件对应你要求的 **第 4 个 md**：  
- α(m) 怎么算  
- backtracking 用到的参数：`L(m)`、`ν(m)`、`β(m)` 以及判别条件  
- 这些量在 R 代码里的计算位置/形式

> 参考：`main_Dec2_revised_gemini.pdf` Section 3.1, Algorithm 1, Eq (3.3)-(3.4).

---

## 1. PG 每一步在做什么？

在当前迭代点 \(\xi^{(m)}=(\gamma^{(m)},\Omega^{(m)})\)：

1) 先用一个 quadratic majorizer 构造 prox 子问题：

\[
s^{(m)}
=\arg\min_{\xi\in D_\delta}
\langle \nabla L(\xi^{(m)}),\xi-\xi^{(m)}\rangle
+\frac{L(m)}{2}\|\xi-\xi^{(m)}\|_F^2
+\Phi(\xi).
\]
（pdf Eq (3.3)）

2) 得到方向：

\[
d^{(m)} = s^{(m)} - \xi^{(m)}.
\]

3) 计算 self-concordant 几何下的两个标量：

- Newton-like decrement（local norm）：
\[
\nu(m) := \|d^{(m)}\|_{\xi^{(m)}}.
\]

- Frobenius norm 缩放：
\[
\beta(m) := \sqrt{L(m)}\|d^{(m)}\|_F.
\]

4) backtracking：如果条件不满足就把 \(L(m)\) 翻倍。

5) 一旦满足条件，用 closed-form 步长更新：

\[
\alpha(m)=\frac{\beta(m)^2}{\nu(m)\left(\nu(m)+\beta(m)^2\right)}.
\]
（pdf Eq (3.4)）

更新：
\[
\xi^{(m+1)}=\xi^{(m)}+\alpha(m)d^{(m)}.
\]

---

## 2. backtracking 判别条件

Algorithm 1 的核心判别（line 6）是：

\[
\frac{\nu(m)^2}{\beta(m)^2}+\nu(m) \le 1.
\]

若不满足，则：
\[
L(m)\leftarrow 2L(m),
\]
并重新解 prox 子问题（因为 prox 子问题依赖 \(L(m)\)）。

---

## 3. R 里怎么计算 ν(m)、β(m)、α(m)

设已经得到了 `s_gamma`, `s_Omega`，方向：

```r
d_gamma <- s_gamma - gamma
d_Omega <- s_Omega - Omega
```

Frobenius norm：

```r
dnormF <- sqrt(sum(d_gamma^2) + sum(d_Omega^2))
beta   <- sqrt(Lm) * dnormF
```

local norm（调用 PG-02 定义的 local_norm）：

```r
nu <- local_norm(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = Oinv)
```

> ⚠️ `nu` 的正确性高度依赖 Hessian 的交叉块实现：必须满足 `Hog = t(Hgo)`。  
> 如果误写成 `Hog <- -t(Hgo)`，再配合 `symmetrize(H)` 会把交叉块平均成 0，导致 local norm \(\nu\) 被系统性低估，进而让 \(\alpha(m)\) 偏大、backtracking 判别条件失真。


判别：

```r
if (nu^2 / max(beta^2, .Machine$double.eps) + nu > 1) {
  Lm <- 2 * Lm
  next
}
```

α(m)：

```r
if (nu <= .Machine$double.eps) {
  alpha <- 1
} else {
  alpha <- beta^2 / (nu * (nu + beta^2))
}
```

---

## 4. 停止准则（Algorithm 1 line 12）

pdf 给的停止是：

\[
\|d^{(m)}\|_F \le \varepsilon.
\]

R 实现一般直接用 `dnormF <= tol`，或加上 objective stabilization：

```r
if (dnormF <= tol) stop
if (abs(obj_new - obj_old) / (1+abs(obj_old)) <= tol_obj) stop
```

---

## 5. Codex 实现提示（常见数值坑）

1. **ν 很小**：直接设 `alpha=1`（否则 0/0）。
2. **beta^2 很小**：分母 `max(beta^2, eps)` 防止爆炸。
3. **Omega 的可行域**：prox(Ω) 已保证 Ω ⪰ δI；更新是 convex combination，理论上保持可行；但建议每步后 `Omega <- proj_psd_shifted(Omega, delta)` 防漂移。
4. **L(m) 初值**：`L0` 太小会导致频繁 backtracking；太大则步长小。一般设 1 或基于谱范数估计。
