# Notebook
# D–D* 的 OBE 势：标准化数学式与同位旋因子

本说明整理并标准化 D–D* 体系中常用的一重玻色子交换（OBE）势的坐标空间表达、带单极型（monopole）形状因子后的正则化核，以及三体 D D D*（总同位旋 Itot=1/2）情形下用 Clebsch–Gordan / 6j 重耦合得到的两体同位旋因子。

## 0. 记号与单位

- 取自然单位 $\hbar=c=1$。半经典常数换算：$hc \equiv 197.3269804\ \mathrm{MeV\cdot fm}$。
- 距离：$r$（MeV$^{-1}$）。若用实验常用的 $R$（fm），则 $r = R/hc$。
- 质量与截断：$m$ 表示介子质量，$\Lambda$ 为形状因子（单极）截断，均以 MeV 计。
- 便于数值实现的无量纲变量：
  $$X \equiv m r,\qquad L \equiv \Lambda/m.$$

- 两体同位旋算子（对等双重态）本征值关系：
  $$\tau_1\cdot\tau_2 = 2\,T(T+1) - 3,\quad T=0\Rightarrow -3,\ T=1\Rightarrow +1.$$

## 1. 带单极形状因子的正则化 Yukawa 核

采用动量空间形状因子
$$F(q)=\frac{\Lambda^2-m^2}{\Lambda^2+q^2}$$
傅里叶变换到坐标空间的标准结果（常用规范）为：

- 维度齐次的核（直接用于势能）：
  $$
  \begin{aligned}
  Y(m,\Lambda,r)
  &= \frac{e^{-mr}-e^{-\Lambda r}}{4\pi r}
     -\frac{\Lambda^2-m^2}{8\pi\Lambda}\,e^{-\Lambda r},\\[4pt]
  T(m,\Lambda,r)
  &= \frac{e^{-mr}}{4\pi r}\!\left(1+\frac{3}{mr}+\frac{3}{(mr)^2}\right)
   - \frac{e^{-\Lambda r}}{4\pi r}\!\left(1+\frac{3}{\Lambda r}+\frac{3}{(\Lambda r)^2}\right)\\
  &\quad - \frac{\Lambda^2-m^2}{8\pi\Lambda}\,e^{-\Lambda r}\!\left(1+\frac{1}{\Lambda r}\right).
  \end{aligned}
  $$

- 便于编码的无量纲写法（把总体尺度提到耦合系数上；用上面定义的 $X=mr,\, L=\Lambda/m$）：
  $$
  \begin{aligned}
  \widehat Y(X,L)
  &= \frac{e^{-X}}{4\pi X}
   - \frac{L\,e^{-L X}}{4\pi X}
   - \frac{L^2-1}{8\pi L}\,e^{-L X},\\[4pt]
  \widehat T(X,L)
  &= \frac{e^{-X}}{4\pi X}\!\left(1+\frac{3}{X}+\frac{3}{X^2}\right)
   - \frac{e^{-L X}}{4\pi X}\!\left(1+\frac{3}{L X}+\frac{3}{(L X)^2}\right)\\
  &\quad - \frac{L^2-1}{8\pi L}\,e^{-L X}\!\left(1+\frac{1}{L X}\right).
  \end{aligned}
  $$
  与上式的量纲关系为
  $$Y(m,\Lambda,r)=m\,\widehat Y(X,L),\qquad T(m,\Lambda,r)=m\,\widehat T(X,L).$$

注：
- 上述形式对应“单极形状因子 + 均匀规范”的常见坐标空间表达。它在 $r\to 0$ 时对非解析的 1/r 奇性进行一致正则化，避免了裸 Yukawa 的发散。
- 若你在代码中以“核函数×质量”的方式写势（例如用 $\widehat Y$ 再乘以 $m$），务必保留第二项中的系数 $L$；少了这一因子会显著加深短程势井。

## 2. D–D* 的 OBE 势分解（S 波对角中心项）

将势写成各介子交换的和：
$$V(r)=V_\pi(r)+V_\sigma(r)+V_\rho(r)+V_\omega(r)+\cdots$$

- 标量 $\sigma$ 交换（通常吸引）：
  $$V_\sigma(r) = -g_S^2\,Y(m_\sigma,\Lambda_\sigma,r).$$

- 等矢量 $\rho$ 交换（随同位旋变号）：
  $$V_\rho(r) = \big(\tau_1\cdot\tau_2\big)\,g_\rho^2\,Y(m_\rho,\Lambda_\rho,r).$$

- 同味矢量 $\omega$ 交换（D–D* 情形一般取同号耦合，不带 $\tau$）：
  $$V_\omega(r) = +g_\omega^2\,Y(m_\omega,\Lambda_\omega,r).$$

- 伪标量 $\pi$ 交换：
  - 对于 P–V（$D$–$D^*$）的“对角 S 波中心势”通常为零；
  - 其物理作用主要通过张量核 $T$ 引起的 S–D 耦合与非对角（$PV\leftrightarrow VP$）道间耦合。若需要张量项：
    $$V_\pi^{(\text{tensor})}(r) = C_\pi\ T(m_\pi,\Lambda_\pi,r)\ S_{12},\qquad C_\pi\propto \frac{g^2}{f_\pi^2}\ (\tau_1\cdot\tau_2),$$
    其中 $S_{12}$ 为张量算子。纯 S 波对角元时 $S_{12}$ 的期望值为 0。

把上式换成无量纲核实现（配合前节定义）：
$$
\begin{aligned}
V_\sigma(r) &= -g_S^2\,m_\sigma\,\widehat Y(X_\sigma,L_\sigma),\\
V_\rho(r)   &= (\tau_1\cdot\tau_2)\,g_\rho^2\,m_\rho\,\widehat Y(X_\rho,L_\rho),\\
V_\omega(r) &= +g_\omega^2\,m_\omega\,\widehat Y(X_\omega,L_\omega),
\end{aligned}
\quad
\begin{aligned}
X_\phi &= m_\phi r,\\
L_\phi &= \Lambda_\phi/m_\phi.
\end{aligned}
$$

参数说明：
- $g_S,\ g_\rho,\ g_\omega$ 为重强子–介子耦合（推荐采用重强子手征拉氏量给出的 $g,\ \beta,\ \lambda,\ g_V$ 等组合；切勿直接用核子–介子 $NN$-OBE 的数值）。
- D–D* 与 D–$\bar D^*$ 在矢量交换上的相对号可能不同；本式针对 D–D* 的“同号耦合”常见取法。

## 3. 三体 D D D*（Itot=1/2）中两体同位旋因子（CG/6j）

考虑三粒子同位旋均为 $1/2$：$D(1)=1/2,\ D(2)=1/2,\ D^*(3)=1/2$。给定总同位旋 $I_{\mathrm{tot}}$，对某一对（如 $1$–$3$ 的 D–D* 对）算子的期望值
$$\langle \tau_1\cdot\tau_3\rangle = \sum_{T_{13}=0,1} P(T_{13}\mid T_{12},I_{\mathrm{tot}})\,\big[2T_{13}(T_{13}+1)-3\big],$$
其中
$$P(T_{13}\mid T_{12},I_{\mathrm{tot}})
= (2T_{12}+1)(2T_{13}+1)\,
\begin{Bmatrix}
\frac12 & \frac12 & T_{12}\\
\frac12 & I_{\mathrm{tot}} & T_{13}
\end{Bmatrix}^2,$$
为 Wigner 6j 给出的重耦合概率；$T_{12}$ 是先耦合的 $(D,D)$ 对的同位旋。

常用情形与结果（亦满足总和恒等式 $\sum_{i<j}\tau_i\cdot\tau_j = 2I_{\mathrm{tot}}(I_{\mathrm{tot}}+1)-\tfrac{9}{2}$）：

- $I_{\mathrm{tot}}=\frac12$，若两颗 D 处于 S 波（交换对称）则 $T_{DD}=1$，有
  $$
  \begin{aligned}
  \langle\tau_{DD}\rangle &= +1,\\
  P(T_{13}=0)&=\tfrac34,\quad P(T_{13}=1)=\tfrac14,\\
  \Rightarrow\ \langle\tau_{D D^*}\rangle &= (-3)\times\tfrac34 + (+1)\times\tfrac14 = -2.
  \end{aligned}
  $$

- $I_{\mathrm{tot}}=\frac12$，若两颗 D 处于相对奇波（交换反对称）则 $T_{DD}=0$，有
  $$
  \langle\tau_{DD}\rangle = -3,\qquad
  \langle\tau_{D D^*}\rangle = (-3)\times\tfrac14 + (+1)\times\tfrac34 = 0.
  $$

- $I_{\mathrm{tot}}=\frac32$（仅供参考）：必有 $T_{DD}=1$ 且任意两体同位旋 $T=1$，故
  $$\langle\tau_{DD}\rangle=\langle\tau_{D D^*}\rangle=+1.$$

因此，你的目标（三体 $I_{\mathrm{tot}}=\tfrac12$，两 D 为 S 波）中用于 D–D* 两体势的同位旋因子应取
$$\boxed{\ \langle\tau_{D D^*}\rangle = -2\ }.$$
若在三体框架中显式区分 $T_{DD}=0/1$ 两通道，则应保留完整的 2×2 同位旋算子结构；若已选定 $T_{DD}=1$，则可将其作为 c-number 插入两体势。

## 4. 汇总：可直接用于实现的标准式

- 无量纲核（建议在代码中使用）：
  $$
  \widehat Y(X,L)
  = \frac{e^{-X}}{4\pi X}
   - \frac{L\,e^{-L X}}{4\pi X}
   - \frac{L^2-1}{8\pi L}\,e^{-L X},
  $$
  $$
  \widehat T(X,L)
  = \frac{e^{-X}}{4\pi X}\!\left(1+\frac{3}{X}+\frac{3}{X^2}\right)
   - \frac{e^{-L X}}{4\pi X}\!\left(1+\frac{3}{L X}+\frac{3}{(L X)^2}\right)
   - \frac{L^2-1}{8\pi L}\,e^{-L X}\!\left(1+\frac{1}{L X}\right).
  $$

- D–D* 的 S 波对角中心势（不含 $\pi$ 的对角中心项）：
  $$
  V_{DD^*}(r)=
   -g_S^2\,m_\sigma\,\widehat Y(X_\sigma,L_\sigma)
   +(\tau_1\cdot\tau_2)\,g_\rho^2\,m_\rho\,\widehat Y(X_\rho,L_\rho)
   +g_\omega^2\,m_\omega\,\widehat Y(X_\omega,L_\omega),
  $$
  其中 $\tau_1\cdot\tau_2$ 在你的三体态（$I_{\mathrm{tot}}=\tfrac12$, 两 D S 波）下取 $-2$。

- 若加入 $\pi$ 的张量作用（S–D 耦合或道间耦合）：
  $$
  V_\pi^{(\text{tensor})}(r)=
  C_\pi\ \widehat T(X_\pi,L_\pi)\ m_\pi\ S_{12},\qquad
  C_\pi\propto \frac{g^2}{f_\pi^2}\ (\tau_1\cdot\tau_2).
  $$

以上各式在符号、系数与维度上彼此一致，便于直接替换到你的代码里。若需要，我可以把这些标准式对应到你现有函数（修正 $\widehat Y$ 第二项缺失的 $L$ 因子、移除不必要的 $\pi$ 接触项、设置 $\tau=-2$）并给出一份完整的 Julia 实现。
