# PN-01：二阶算法在做什么？为什么要用 Damped Prox-Newton？

本文件用于你要求的第二部分：「解释二阶算法都在干什么、解决什么问题、damped PN 在解决什么」。

> 参考：`main_Dec2_revised_gemini.pdf` Section 3.2（Algorithm 2）。

---

## 1. 目标问题（复合 self-concordant）

我们要解：

\[
\min_{\xi=(\gamma,\Omega)\in D_\delta}\;
F(\xi)=L(\xi)+\Phi(\xi)
\]

其中：
- \(L(\xi)\) 是 scaled Gaussian loss（Eq (2.2)），**标准 self-concordant**（Section 2.3）
- \(\Phi(\xi)=\lambda_\gamma\|\gamma_P\|_1+\lambda_\Omega\|\Omega\|_{1,\mathrm{off}}\) 是凸但不可微
- \(D_\delta=\{(\gamma,\Omega):\Omega\succeq \delta I\}\) 是可行域（防止 Ω 奇异）

---

## 2. 为什么一阶 PG 可能慢？为什么 PN 有意义？

### 2.1 关键：loss 的梯度不满足全局 Lipschitz

对 logistic/least-squares 类型的 smooth loss，\(\nabla L\) 常常有 Lipschitz 常数，从而标准 PG 的线性收敛更容易分析。  
但本问题里有：

- \(-\log\det(\Omega)\)
- \(\Omega^{-1}\gamma^\top S_{xx}\gamma\) 里的 \(\Omega^{-1}\)

当 Ω 的最小特征值变小，Hessian 会爆炸；因此 **全局 Lipschitz 常数不存在**，用 “Lipschitz-based PG” 会很保守甚至不成立。

### 2.2 self-concordant 框架：用局部几何替代 Lipschitz 常数

self-concordant 理论允许用局部 Hessian 定义的几何量（local norm / Newton decrement）来控制下降，不需要全局 Lipschitz。

这就是你 pdf 里 PG/PN 都基于：
- local norm \(\|\cdot\|_{\xi}\)
- Newton decrement \(\nu\)

---

## 3. Prox-Newton 的核心思想

在当前点 \(\xi^{(m)}\)，对 smooth part 做二阶近似：

\[
L(\xi)\approx
L(\xi^{(m)})
+\langle\nabla L(\xi^{(m)}),\xi-\xi^{(m)}\rangle
+\frac12 \nabla^2L[\xi^{(m)}](\xi-\xi^{(m)},\xi-\xi^{(m)}).
\]

然后加上非光滑惩罚 \(\Phi(\xi)\)，得到 PN 子问题（见 PN-02）。

直觉上：
- PG 用的是 \(\frac{L}{2}\|\xi-\xi^{(m)}\|_F^2\) 这种 **欧氏** 二次项
- PN 用的是 \(\frac12 \nabla^2L[\xi^{(m)}](\cdot,\cdot)\) 这种 **自适应的 Hessian 二次项**
  → 更贴合局部曲率，通常更快（尤其在靠近最优解时）

---

## 4. 为什么需要 damped（阻尼）？

在经典 Newton 方法里，如果直接用 full step（α=1），当离最优点远时会：
- 可能跑出可行域（这里 Ω 可能不再 SPD）
- 可能导致目标上升（尤其 composite 情况）

self-concordant 理论提供一个 **无需 Lipschitz** 的全局步长：

\[
\alpha^{(m)} = \frac{1}{1+\nu^{(m)}}
\]

⚠️ **重要实现约束（按 pdf Algorithm 2）**：上面这个显式步长 \(\alpha=1/(1+\nu)\) 本身就是 self-concordant 理论给的“下降保证”。  
一般不需要、也不建议再额外套 Armijo line-search。  

如果你为了工程稳健性想加 *monotone safeguard*，建议只做最简单的回溯：  
- 先用 \(\alpha=1/(1+\nu)\)；  
- 若 `F(gamma+αd_gamma, Omega+αd_Omega)` 没下降，就令 `α <- α/2` 再试；  
- 若降不下来就 **不更新并停止/提高子问题精度**（不要 “forced accept”）。  

尤其不要把 Armijo 的梯度项写成 `-α * <∇f, d>`（负负得正会错误接受上升步）。


其中 \(\nu^{(m)}\) 是 Newton decrement（本质是方向在 local norm 下的长度）。

这个 damped 步长保证：
- 在 self-concordant 框架下，**全局收敛**
- 并且通常能保证目标下降（如果 PN 子问题解得足够精确）

---

## 5. PN 与两层 ADMM 的关系（你代码里最容易乱的点）

PN 里有两个层次的“难点”：

1) **PN 子问题（Eq (3.6)）** 本身是一个带 ℓ1 + PSD 约束的凸问题（规模大），所以用 **外层 ADMM**（Algorithm 3）来解。

2) 外层 ADMM 的 z-update 里会出现 Ω 的 proximal mapping：

\[
\min_{\Omega\succeq\delta I}\;\frac12\|\Omega-V\|_F^2+\tau_\lambda\|\Omega\|_{1,off}
\]

这个 Ω-prox 又需要 **inner ADMM**（Appendix B）来解。

所以整体是：
- outer loop：Prox-Newton（Algorithm 2）
- middle loop：PN subproblem outer ADMM（Algorithm 3）
- inner loop：Ω-prox inner ADMM（Appendix B）

下一份文件 PN-02/03/04 会分别拆解。

---
