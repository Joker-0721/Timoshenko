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

沿厚度积分得到弯矩和剪力：

$$
M_{\alpha\beta} = \int_{-h/2}^{h/2}x_3\sigma_{\alpha\beta}dx_3 = \frac{h^3}{12}D_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta}
$$

$$
Q_{\alpha} = k\int_{-h/2}^{h/2}\tau_{3\alpha}dx_3 = kGh\,\gamma_{\alpha}
$$

其中 $k$ 为剪切修正系数（通常取 $5/6$）。

## 4. 弱形式

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

$$

\int_{\Omega} \Bigl(M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} + Q_{\alpha}\,\delta\phi_\alpha + Q_{\alpha}\,\delta w_{,\alpha} - P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \Bigr) d\Omega = 0
\tag{3}
$$

$$
\int_{\Omega} \Bigl( M_{\alpha\beta}\,\delta\kappa_{\alpha\beta} + Q_{\alpha}\,\delta\gamma_{\alpha} - P_{\alpha\beta} w_{,\alpha}\,\delta w_{,\beta} \Bigr) d\Omega = 0
$$


记弯曲刚度贡献的虚功为 $a_b$，剪切刚度贡献的虚功为 $a_s$，几何刚度贡献的虚功为 $a_g$，则上式可简记为：

$$
a_b + a_s - a_g = 0
$$

## 5. 有限元离散



单元有 $n$ 个节点，形函数 $N_K(x,y)$。位移插值：
$$
w^h = \sum_K N_K w_K,\quad
\phi_x^h = \sum_K N_K \phi_{xK},\quad
\phi_y^h = \sum_K N_K \phi_{yK}.
$$
虚位移同理。导数：
$$
w^h_{,x} = \sum_K N_{K,x} w_K,\quad
w^h_{,y} = \sum_K N_{K,y} w_K,\quad
\phi^h_{x,x} = \sum_K N_{K,x} \phi_{xK},\quad
\phi^h_{x,y} = \sum_K N_{K,y} \phi_{xK},\quad
\phi^h_{y,x} = \sum_K N_{K,x} \phi_{yK},\quad
\phi^h_{y,y} = \sum_K N_{K,y} \phi_{yK}.
$$
数值积分：$\int_{\Omega^e} f d\Omega \approx \sum_g f(\xi_g,\eta_g) w_g$。

自由度编号：节点 $K$ 有 $w_K$（编号 $I$），$\phi_{xK}$（编号 $2I-1$），$\phi_{yK}$（编号 $2I$）。

---

### 5.1弯曲刚度部分

$$
a_b = \int_{\Omega} \frac{h^3}{12} D_{\alpha\beta\gamma\eta} \,\kappa_{\alpha\beta} \,\delta\kappa_{\gamma\eta} \, d\Omega,
$$
$$
\kappa_{\alpha\beta} = -\frac12(\phi_{\alpha,\beta}+\phi_{\beta,\alpha}),\quad
\delta\kappa_{\gamma\eta} = -\frac12(\delta\phi_{\gamma,\eta}+\delta\phi_{\eta,\gamma}).
$$
材料张量（各向同性）：
$$
D_{1111}=D_{2222}=\frac{Eh^3}{12(1-\nu^2)},\quad
D_{1122}=D_{2211}=\frac{E\nu h^3}{12(1-\nu^2)},\quad
D_{1212}=D_{1221}=D_{2112}=D_{2121}=\frac{Eh^3}{24(1+\nu)}.
$$

将 
$\phi_{\alpha,\beta}=\sum_R N_{R,\beta}\phi_{\alpha R}$ 和 $\delta\phi_{\gamma,\eta}=\sum_S N_{S,\eta}\delta\phi_{\gamma S}$ 代入，得：
$$
a_b = \sum_{R,S} \delta\phi_{iR}\, \bigl(\mathbf{K}_{\phi\phi,RS}^b\bigr)_{ij}\, \phi_{jS},
$$
其中
$$
\begin{aligned}
\bigl(\mathbf{K}_{\phi\phi,RS}^b\bigr)_{11} &=
\frac{h^3}{12} \int_{\Omega} \Bigl( D_{1111} N_{R,x} N_{S,x} + D_{1212} N_{R,y} N_{S,y} \Bigr) d\Omega, \\[4pt]
\bigl(\mathbf{K}_{\phi\phi,RS}^b\bigr)_{12} &=
\frac{h^3}{12} \int_{\Omega} \Bigl( D_{1122} N_{R,x} N_{S,y} + D_{1212} N_{R,y} N_{S,x} \Bigr) d\Omega, \\[4pt]
\bigl(\mathbf{K}_{\phi\phi,RS}^b\bigr)_{21} &=
\frac{h^3}{12} \int_{\Omega} \Bigl( D_{1122} N_{R,y} N_{S,x} + D_{1212} N_{R,x} N_{S,y} \Bigr) d\Omega, \\[4pt]
\bigl(\mathbf{K}_{\phi\phi,RS}^b\bigr)_{22} &=
\frac{h^3}{12} \int_{\Omega} \Bigl( D_{1212} N_{R,x} N_{S,x} + D_{1111} N_{R,y} N_{S,y} \Bigr) d\Omega.
\end{aligned}
$$
该子矩阵嵌入整体刚度时，对应自由度 $(2R-1,2R)$ 与 $(2S-1,2S)$。

---

### 5.2 剪切刚度部分

剪切双线性形式：
$$
a_s= \int_{\Omega} k G h \, \gamma_\alpha \, \delta\gamma_\alpha \, d\Omega,\qquad
\gamma_\alpha = w_{,\alpha} - \phi_\alpha.
$$
其中 $k=5/6$，$G = E/(2(1+\nu))$，$kGh = \frac{5}{6}\cdot\frac{Eh}{2(1+\nu)}$。

展开为三个独立双线性形式：

### 5.2.1 挠度-挠度
$$
a_{s}^{ww}(w,\delta w) = \int_{\Omega} k G h \, w_{,\alpha} \,\delta w_{,\alpha} \, d\Omega.
$$
离散得子矩阵：
$$
K_{ww,KL}^s = kGh \int_{\Omega} \bigl( N_{K,x}N_{L,x} + N_{K,y}N_{L,y} \bigr) d\Omega.
$$


### 5.2.2 转角-转角
$$
a_{s}^{\phi\phi}(\phi,\delta\phi) = \int_{\Omega} k G h \, \phi_\alpha \,\delta\phi_\alpha \, d\Omega.
$$
离散得子矩阵：
$$
\mathbf{K}_{\phi\phi,KL}^s = kGh \int_{\Omega} N_K N_L \, d\Omega \; \mathbf{I}_2.
$$

$\mathbf{I}_2$ 是 $2\times 2$ 单位矩阵，表示 $\phi_x$ 和 $\phi_y$ 的耦合关系。

### 5.2.3 耦合

#### 第一项：
$$
a_{s}^{w\delta\phi}(w,\delta\phi) = -kGh \int_{\Omega} w_{,\alpha} \,\delta\phi_\alpha \, d\Omega.
$$
离散后给出子矩阵 $\mathbf{K}_{w\phi}^s$（从 $w$ 到 $\delta\phi$）：
$$
\bigl(\mathbf{K}_{w\phi,KL}^s\bigr)_{1} = -kGh \int_{\Omega} N_{K,x} N_L \, d\Omega,\qquad
\bigl(\mathbf{K}_{w\phi,KL}^s\bigr)_{2} = -kGh \int_{\Omega} N_{K,y} N_L \, d\Omega.
$$

#### 第二项：
$$
a_{s}^{\phi\delta w}(\phi,\delta w) = -kGh \int_{\Omega} \phi_\alpha \,\delta w_{,\alpha} \, d\Omega.
$$
离散后给出子矩阵 $\mathbf{K}_{\phi w}^s$（从 $\phi$ 到 $\delta w$）：
$$
\bigl(\mathbf{K}_{\phi w,KL}^s\bigr)_{1} = -kGh \int_{\Omega} N_K N_{L,x} \, d\Omega,\qquad
\bigl(\mathbf{K}_{\phi w,KL}^s\bigr)_{2} = -kGh \int_{\Omega} N_K N_{L,y} \, d\Omega.
$$


---

### 5.3几何刚度部分

双线性形式 $a_g(w,\delta w)$：
$$
a_g(w,\delta w) =  P_{\alpha\beta} \int_{\Omega}  w_{,\alpha} \,\delta w_{,\beta} \, d\Omega.
$$
离散后子矩阵：
$$
K_{ww,KL}^g = P_{\alpha\beta} \int_{\Omega}  N_{K,\alpha} N_{L,\beta} \, d\Omega.
$$

---

### 5.4 总体离散方程

将所有双线性形式相加，并令总虚功为零，得到：
$$
a_b(\phi,\delta\phi) + a_s^{ww}(w,\delta w) + a_s^{\phi\phi}(\phi,\delta\phi) + a_s^{w\phi}(w,\phi,\delta w,\delta\phi) - a_g(w,\delta w) = 0.
$$

按自由度分组，写成矩阵形式：
则矩阵形式为：
$$
\begin{bmatrix}
K_{ww}^{s} & K_{w\phi1}^{s} & K_{w\phi2}^{s} \\
K_{\phi1 w}^{s} & K_{\phi1\phi1}^{b}+K_{\phi1\phi1}^{s} & K_{\phi1\phi2}^{b}+K_{\phi1\phi2}^{s} \\
K_{\phi2 w}^{s} & K_{\phi2\phi1}^{b}+K_{\phi2\phi1}^{s} & K_{\phi2\phi2}^{b}+K_{\phi2\phi2}^{s}
\end{bmatrix}
\begin{bmatrix}
\mathbf{d}_w \\ \mathbf{d}_{\phi1} \\ \mathbf{d}_{\phi2}
\end{bmatrix}
-
\begin{bmatrix}
K_{ww}^{g} & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
\mathbf{d}_w \\ \mathbf{d}_{\phi1} \\ \mathbf{d}_{\phi2}
\end{bmatrix}
= \mathbf{0}.
$$

令整体线性刚度矩阵 $K$ 和几何刚度矩阵 $K_g$ 为：
$$
K =
\begin{bmatrix}
K_{ww}^{s} & K_{w\phi}^{s} \\
K_{\phi w}^{s} & K_{\phi\phi}^{b} + K_{\phi\phi}^{s}
\end{bmatrix},\qquad
K_g =
\begin{bmatrix}
K_{ww}^{g} & 0 \\
0 & 0
\end{bmatrix}.
$$
则特征方程写为：
$$
[\,K - K_g\,]\,\mathbf{d} = \mathbf{0}.
$$


---
