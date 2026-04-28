## 1. 运动学

考虑一块厚板 $\Omega$，厚度为 $h$，中面位移为 $u_\alpha$（$\alpha=1,2$），横向挠度为 $w$，法线转角为 $\phi_\alpha$。则任意点的位移可表示为：

$$
\bar{u}_{\alpha}(x_1,x_2,x_3) = u_{\alpha}(x_1,x_2) - x_3\phi_{\alpha},\qquad \bar{u}_3 = w
$$

相应的应变分量为：

$$
\epsilon_{\alpha\beta} = \varepsilon_{\alpha\beta} + x_3\kappa_{\alpha\beta}
$$

其中中面应变和曲率定义为：

$$
\varepsilon_{\alpha\beta} = \frac{1}{2}(u_{\alpha,\beta}+u_{\beta,\alpha}),\qquad
\kappa_{\alpha\beta} = -\frac{1}{2}(\phi_{\alpha,\beta}+\phi_{\beta,\alpha})
$$

横向剪切应变为：

$$
\gamma_{\alpha} = 2\epsilon_{\alpha3} = w_{,\alpha} - \phi_{\alpha}
$$

## 2. 本构关系

对于各向同性线弹性材料（平面应力状态，$\sigma_{33}=0$），应力-应变关系为：

$$
\sigma_{\alpha\beta} = D_{\alpha\beta\gamma\eta}\,\epsilon_{\gamma\eta},\qquad
\tau_{3\alpha} = G\,\gamma_{\alpha}
$$

其中弹性张量 $D_{\alpha\beta\gamma\eta}$ 和剪切模量 $G$ 分别为：

$$
D_{\alpha\beta\gamma\eta} = \frac{E}{1-\nu^2}\Big(\nu\delta_{\alpha\beta}\delta_{\gamma\eta} + \frac{1-\nu}{2}(\delta_{\alpha\gamma}\delta_{\beta\eta}+\delta_{\alpha\eta}\delta_{\beta\gamma})\Big)
$$

$$
G = \frac{E}{2(1+\nu)}
$$

## 3. 合力（内力素）

沿厚度积分得到膜力、弯矩和剪力：

$$
{N}_{\alpha\beta} = \int_{-h/2}^{h/2}\sigma_{\alpha\beta}dx_3 = hD_{\alpha\beta\gamma\eta}\,\varepsilon_{\gamma\eta}
$$

$$
M_{\alpha\beta} = \int_{-h/2}^{h/2}x_3\sigma_{\alpha\beta}dx_3 = \frac{h^3}{12}D_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta}
$$

$$
Q_{\alpha} = k\int_{-h/2}^{h/2}\tau_{3\alpha}dx_3 = kGh\,\gamma_{\alpha}
$$

其中 $k$ 为剪切修正系数（通常取 $5/6$）。

## 4. 平衡方程与边界条件（强形式）

考虑初始面内力 $P_{\alpha\beta}$ 的几何效应（屈曲分析），总势能的一阶变分为零给出：

$$
\delta\Pi = \int_{\Omega}\Big(M_{\alpha\beta}\,\delta\kappa_{\alpha\beta}+Q_{\alpha}\,\delta\gamma_{\alpha}-P_{\alpha\beta}\,w_{,\alpha}\delta w_{,\beta}\Big)d\Omega = 0
$$

对上述变分方程进行分部积分（忽略边界项推导细节），得到域内平衡方程：

$$
\begin{aligned}
& M_{\alpha\beta,\beta} - Q_{\alpha} = 0 \quad \text{在 } \Omega \text{ 内}, \\
& Q_{\alpha,\alpha} - (P_{\alpha\beta} w_{,\beta})_{,\alpha} = 0 \quad \text{在 } \Omega \text{ 内}
\end{aligned}
$$


## 5. Galerkin 弱形式（详细推导）

### 5.1 从强形式出发

将强形式中的每个平衡方程乘以对应的虚位移（测试函数），并在全域 $\Omega$ 上积分：

$$
\int_{\Omega} (M_{\alpha\beta,\beta} - Q_{\alpha})\,\delta\phi_\alpha \,d\Omega = 0,\qquad
\int_{\Omega} \bigl( Q_{\alpha,\alpha} - (P_{\alpha\beta} w_{,\beta})_{,\alpha} \bigr)\,\delta w \,d\Omega = 0
$$

### 5.2 分部积分（转移导数）

对第一个方程中的 $M_{\alpha\beta,\beta}$ 项使用散度定理：

$$
\int_{\Omega} M_{\alpha\beta,\beta}\,\delta\phi_\alpha \,d\Omega = \int_{\Gamma} M_{\alpha\beta} n_{\beta}\,\delta\phi_\alpha \,d\Gamma - \int_{\Omega} M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} \,d\Omega
$$

于是第一个方程变为：

$$
\int_{\Gamma} M_{\alpha\beta} n_{\beta}\,\delta\phi_\alpha \,d\Gamma - \int_{\Omega} M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} \,d\Omega - \int_{\Omega} Q_{\alpha}\,\delta\phi_\alpha \,d\Omega = 0 \tag{1}
$$

对第二个方程中的两项分别分部积分：

$$
\int_{\Omega} Q_{\alpha,\alpha}\,\delta w \,d\Omega = \int_{\Gamma} Q_{\alpha} n_{\alpha}\,\delta w \,d\Gamma - \int_{\Omega} Q_{\alpha}\,\delta w_{,\alpha} \,d\Omega
$$

$$
-\int_{\Omega} (P_{\alpha\beta} w_{,\beta})_{,\alpha}\,\delta w \,d\Omega = -\int_{\Gamma} P_{\alpha\beta} w_{,\beta} n_{\alpha}\,\delta w \,d\Gamma + \int_{\Omega} P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \,d\Omega
$$

将这两部分相加，得到第二个方程的展开式：

$$
\int_{\Gamma} Q_{\alpha} n_{\alpha}\,\delta w \,d\Gamma - \int_{\Omega} Q_{\alpha}\,\delta w_{,\alpha} \,d\Omega - \int_{\Gamma} P_{\alpha\beta} w_{,\beta} n_{\alpha}\,\delta w \,d\Gamma + \int_{\Omega} P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \,d\Omega = 0 \tag{2}
$$

### 5.3 整理域积分与边界积分

将方程 (1) 和 (2) 相加，把所有域积分项移到等式左边，边界积分项移到右边，得到：

$$
\begin{aligned}
& \int_{\Omega} \Bigl( - M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} - Q_{\alpha}\,\delta\phi_\alpha - Q_{\alpha}\,\delta w_{,\alpha} + P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \Bigr) d\Omega \\
= & -\int_{\Gamma} M_{\alpha\beta} n_{\beta}\,\delta\phi_\alpha \,d\Gamma - \int_{\Gamma} Q_{\alpha} n_{\alpha}\,\delta w \,d\Gamma + \int_{\Gamma} P_{\alpha\beta} w_{,\beta} n_{\alpha}\,\delta w \,d\Gamma
\end{aligned}
$$

两边乘以 $-1$ 并重新整理，使内力虚功项为正：

$$
\int_{\Omega} \Bigl( M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} + Q_{\alpha}\,\delta\phi_\alpha + Q_{\alpha}\,\delta w_{,\alpha} - P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \Bigr) d\Omega = \int_{\Gamma} M_{\alpha\beta} n_{\beta}\,\delta\phi_\alpha \,d\Gamma + \int_{\Gamma} Q_{\alpha} n_{\alpha}\,\delta w \,d\Gamma - \int_{\Gamma} P_{\alpha\beta} w_{,\beta} n_{\alpha}\,\delta w \,d\Gamma
$$

利用自然边界条件 $M_{\alpha\beta}n_{\beta} = \bar{M}_\alpha$（在 $\Gamma_M$ 上）和 $(Q_{\alpha} - P_{\alpha\beta} w_{,\beta})n_{\alpha} = \bar{Q}$（在 $\Gamma_Q$ 上），边界积分可合并为：

$$
\int_{\Gamma_M} \bar{M}_\alpha\,\delta\phi_\alpha \,d\Gamma + \int_{\Gamma_Q} \bar{Q}\,\delta w \,d\Gamma
$$

因此弱形式可写为：

$$
\boxed{
\int_{\Omega} \Bigl( M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} + Q_{\alpha}\,\delta\phi_\alpha + Q_{\alpha}\,\delta w_{,\alpha} - P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \Bigr) d\Omega = \int_{\Gamma_M} \bar{M}_\alpha\,\delta\phi_\alpha \,d\Gamma + \int_{\Gamma_Q} \bar{Q}\,\delta w \,d\Gamma
} \tag{3}
$$

$$
\int_{\Omega} \Bigl( M_{\alpha\beta}\,\delta\kappa_{\alpha\beta} + Q_{\alpha}\,\delta\gamma_{\alpha} - P_{\alpha\beta} w_{,\alpha}\,\delta w_{,\beta} \Bigr) d\Omega = \int_{\Gamma_M} \bar{M}_\alpha\,\delta\phi_\alpha \,d\Gamma + \int_{\Gamma_Q} \bar{Q}\,\delta w \,d\Gamma
$$

$$
\boxed{\int_{\Omega}\Bigl(M_{\alpha\beta}\,\delta\kappa_{\alpha\beta} + Q_{\alpha}\,\delta\gamma_{\alpha} - P_{\alpha\beta}\, w_{,\alpha}\,\delta w_{,\beta} \Bigr) d\Omega = 0}
$$

记弯曲刚度贡献的虚功为 $a_b$，剪切刚度贡献的虚功为 $a_s$，几何刚度贡献的虚功为 $a_g$，则上式可简记为：

$$
a_b + a_s - a_g = 0
$$

## 6. 有限元离散

采用标准等参单元，形函数为 $N_K(x)$，节点自由度为 $w_K, \phi_{\alpha K}$：

$$
w^h(x) = \sum_{K} N_K(x) w_K,\qquad
\phi_{\alpha}^h (x)= \sum_{R} N_R(x) \phi_{\alpha R}
$$

相应的离散曲率、剪切应变和挠度梯度为：

$$
\kappa_{\alpha\beta R} = -\frac{1}{2}\bigl(N_{R,\beta}\phi_{\alpha R}+N_{R,\alpha}\phi_{\beta R}\bigr),\qquad
\gamma_{\alpha K} = N_{K,\alpha} w_K - N_K \phi_{\alpha K},\qquad
w_{,\alpha K} = N_{K,\alpha} w_K
$$

### 6.1 弯曲刚度部分

$$
a_b = \int_{\Omega} \frac{h^3}{12} D_{\alpha\beta\gamma\eta} \, \kappa_{\alpha\beta}^h \, \delta\kappa_{\gamma\eta}^h \, d\Omega
$$

将 $\kappa_{\alpha\beta}^h$ 的表达式代入，得到离散形式：

$$
a_b = \sum_{R,S} \delta\phi_{iR}
\left( \int_{\Omega} \frac{h^3}{48} D_{\alpha\beta\gamma\eta}
\bigl( N_{R,\beta}\delta_{\alpha i} + N_{R,\alpha}\delta_{\beta i} \bigr)
\bigl( N_{S,\eta}\delta_{\gamma j} + N_{S,\gamma}\delta_{\eta j} \bigr)
d\Omega \right) \phi_{jS}
$$

定义弯曲刚度子矩阵：

$$
K^b_{RS} = \frac{h^3}{48} \int_{\Omega} D_{\alpha\beta\gamma\eta}
\bigl( N_{R,\beta}\delta_{\alpha i} + N_{R,\alpha}\delta_{\beta i} \bigr)
\bigl( N_{S,\eta}\delta_{\gamma j} + N_{S,\gamma}\delta_{\eta j} \bigr)
d\Omega
$$

因此

$$
a^b = \sum_{K,L} \begin{bmatrix} \delta\phi_{1K} & \delta\phi_{2K} \end{bmatrix} \mathbf{K}^b_{KL} \begin{bmatrix} \phi_{1L} \\ \phi_{2L} \end{bmatrix}
$$

### 6.2 剪切刚度部分

$$
a_s = \int_{\Omega} k G h \, \gamma_{\alpha}^h \, \delta\gamma_{\alpha}^h \, d\Omega
$$

展开后得到：

$$
\begin{aligned}
a_s = \sum_{K,L} \Bigg[ &\delta w_K \left( \int_{\Omega} k G h \, N_{K,\alpha} N_{L,\alpha} \, d\Omega \right) w_L \\
+ &\delta w_K \left( \int_{\Omega} -k G h \, N_{K,\alpha} N_L \, d\Omega \right) \phi_{\alpha L} \\
+ &\delta\phi_{\alpha K} \left( \int_{\Omega} -k G h \, N_K N_{L,\alpha} \, d\Omega \right) w_L \\
+ &\delta\phi_{\alpha K} \left( \int_{\Omega} k G h \, N_K N_L \, \delta_{\alpha\beta} \, d\Omega \right) \phi_{\beta L} \Bigg]
\end{aligned}
$$

剪切刚度子矩阵为：

$$
\mathbf{K}^s_{KL} = kGh \int_{\Omega}
\begin{bmatrix}
N_{K,\alpha}N_{L,\alpha} & -N_{K,\alpha}N_L \\
-N_K N_{L,\beta} & N_K N_L\,\delta_{\alpha\beta}
\end{bmatrix} d\Omega
$$

于是

$$
a_s = \sum_{K,L}
\begin{bmatrix}
\delta w_K & \delta\phi_{\alpha K}
\end{bmatrix}
\mathbf{K}^s_{KL}
\begin{bmatrix}
w_L \\
\phi_{\beta L}
\end{bmatrix}
$$

### 6.3 几何刚度部分

$$
a_g = P_{\alpha\beta}\int_{\Omega} w_{,\alpha}\,\delta w_{,\beta}\,d\Omega
$$

离散后：

$$
a_g = \sum_{K,L} \delta w_K \left(P_{\alpha\beta}\int_{\Omega} \, N_{K,\alpha} N_{L,\beta} \, d\Omega \right) w_L
$$

几何刚度子矩阵：

$$
\mathbf{K}_{KL}^g = P_{\alpha\beta}\int_{\Omega} \, N_{K,\alpha} N_{L,\beta} \, d\Omega
$$

### 6.4 总体离散方程

将上述刚度矩阵组装，得到屈曲特征值问题的有限元方程：

$$
\left(
\begin{bmatrix}
\mathbf{K}_{ww}^{s} & \mathbf{K}_{w\phi}^{s} \\
\mathbf{K}_{\phi w}^{s} & \mathbf{K}^{b}+\mathbf{K}_{\phi\phi}^{s}
\end{bmatrix}
-
\begin{bmatrix}
\mathbf{K}_{ww}^{g} & 0 \\
0 & 0
\end{bmatrix}
\right)
\begin{bmatrix}
\mathbf{d}_w \\
\mathbf{d}_\phi
\end{bmatrix} = \mathbf{0}
$$

该方程用于求解临界载荷参数（隐含在 $P_{\alpha\beta}$ 中）。

---

以上即为完整的中文翻译与第 5 节的详细推导替换。
