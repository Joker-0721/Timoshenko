# Mindlin 板线性屈曲有限元分析


## 一、基本假设

板中面为 $x_3 = 0$，厚度 $h$ 常数。  
中面上任一点 $(x_1,x_2)$ 有 5 个独立位移场：
- 中面内位移：$u_\alpha(x_1,x_2)$，$\alpha=1,2$
- 横向挠度：$w(x_1,x_2)$
- 转角：$\phi_\alpha(x_1,x_2)$

板内任意一点 $(x_1,x_2,x_3)$ 的位移：
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

---

## 二、本构关系

由于板很薄，忽略厚度方向正应力：$\sigma_{33}=0$。  
各向同性线弹性材料的平面应力本构：
$$
\sigma_{\alpha\beta} = C_{\alpha\beta\gamma\eta}\; \varepsilon_{\gamma\eta},
$$
其中
$$
C_{\alpha\beta\gamma\eta} = \frac{E}{1-\nu^2}\Bigl[ \nu\,\delta_{\alpha\beta}\delta_{\gamma\eta} + \frac{1-\nu}{2}(\delta_{\alpha\gamma}\delta_{\beta\eta}+\delta_{\alpha\eta}\delta_{\beta\gamma}) \Bigr].
$$

横向剪切应力：
$$
\sigma_{\alpha3} = G \,\gamma_\alpha,\qquad G = \frac{E}{2(1+\nu)}.
$$

---

## 三、弱形式

对于处于平衡的弹性体，虚功方程：
$$
\delta W_{\text{int}} = \delta W_{\text{ext}}.
$$

三维内虚功：
$$
\delta W_{\text{int}} = \int_V \sigma_{ij}\,\delta\varepsilon_{ij}\,dV
= \int_V \bigl( \sigma_{\alpha\beta}\,\delta\varepsilon_{\alpha\beta} + \sigma_{\alpha3}\,\delta\gamma_\alpha \bigr) dV.
$$

将体积分沿厚度方向积分，并代入应变分解 $\delta\varepsilon_{\alpha\beta} = \delta\varepsilon_{\alpha\beta}^{(0)} + x_3\delta\kappa_{\alpha\beta}$，得：
$$
\delta W_{\text{int}} = \int_{\Omega} \Bigl(\underbrace{\int_{-h/2}^{h/2}\sigma_{\alpha\beta}x_3dx_3}_{M_{\alpha\beta}} \delta\kappa_{\alpha\beta}
+ \underbrace{\int_{-h/2}^{h/2}\sigma_{\alpha3}dx_3}_{Q_\alpha} \delta\gamma_\alpha \Bigr) d\Omega.
$$

将本构关系 $\sigma_{\alpha\beta}=C_{\alpha\beta\gamma\eta}(\varepsilon_{\gamma\eta}^{(0)}+x_3\kappa_{\gamma\eta})$ 代入积分，$k=\frac56$：
$$
\begin{aligned}
M_{\alpha\beta} &= \frac{h^3}{12}C_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta} \equiv D_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta},\\
Q_{\alpha} &= k G h \,\gamma_\alpha \equiv S_{\alpha\beta}\,\gamma_\beta,\quad S_{\alpha\beta}=kGh\delta_{\alpha\beta}
\end{aligned}
$$

代入虚功表达式，得到 **仅含中面应变和曲率的双线性形式**：
$$
\delta W_{\text{int}} = \int_{\Omega} \Bigl(\delta\kappa_{\alpha\beta} D_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta}
+ \delta\gamma_\alpha S_{\alpha\beta}\,\gamma_\beta \Bigr) d\Omega.
$$

此即为有限元离散的弱形式基础。

---

## 四、几何刚度

### 4.1 实际载荷与载荷因子
实际面内载荷为 $P = \lambda P_{\text{ref}}$，则实际应力为：
$$
\sigma_{\alpha\beta}^0 = \lambda \,\sigma_{\alpha\beta}^{\text{ref}}.
$$
$\lambda$ 是一个无量纲的**载荷因子**.

### 4.2 几何刚度的来源（非线性应变）
当板已经存在面内应力 $\sigma_{\alpha\beta}^0$ 时，在横向挠曲过程中，非线性应变$\frac12 w_{,\alpha} w_{,\beta}$ 会使初始应力做虚功。此虚功的一阶变分为：
$$
\delta W_{\text{nl}} = \int_{\Omega} \sigma_{\alpha\beta}^0 \; w_{,\alpha}\,\delta w_{,\beta} \, d\Omega.
$$
将 $\sigma_{\alpha\beta}^0 = \lambda \sigma_{\alpha\beta}^{\text{ref}}$ 代入，得：
$$
\delta W_{\text{nl}} = \lambda \int_{\Omega} \sigma_{\alpha\beta}^{\text{ref}} \; w_{,\alpha}\,\delta w_{,\beta} \, d\Omega.
$$

### 4.3 屈曲时的总虚功平衡
在线性屈曲分析中，忽略外载荷的横向分量（即 $\delta W_{\text{ext}}=0$），并假设屈曲发生时横向位移为任意小的扰动。总内虚功包含线性部分（弯曲+剪切）和几何非线性部分：
$$
\delta W_{\text{int,lin}} + \delta W_{\text{nl}} = 0.
$$
即：
$$
\int_{\Omega} \Bigl(\delta\kappa_{\alpha\beta} D_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta}
+ \delta\gamma_\alpha S_{\alpha\beta}\,\gamma_\beta \Bigr) d\Omega
+ \lambda \int_{\Omega} \sigma_{\alpha\beta}^{\text{ref}} \; w_{,\alpha}\,\delta w_{,\beta} \, d\Omega = 0.
$$

---

## 五、有限元離散（雙線性形式）



$$
w_{,\alpha}=\sum_I N_{I,\alpha}w_I,\qquad
\delta w_{,\beta}=\sum_J N_{J,\beta}\delta w_J
$$

$$
\phi_{\alpha,\beta}=\sum_K N_{K,\beta}\phi_{\alpha K},\qquad
\delta\phi_{\alpha,\beta}=\sum_L N_{L,\beta}\delta\phi_{\alpha L}
$$

$$
\begin{aligned}
0
={}&
\int_{\Omega}
\Bigl(
\delta\kappa_{\alpha\beta}D_{\alpha\beta\gamma\eta}\kappa_{\gamma\eta}
+\delta\gamma_\alpha S_{\alpha\beta}\gamma_\beta
\Bigr)d\Omega
+\lambda\int_{\Omega}
\sigma_{\alpha\beta}^{\mathrm{ref}}
w_{,\alpha}\delta w_{,\beta}\,d\Omega
\end{aligned}
$$

$$
\boxed{
K_b(\delta\phi,\phi)
+K_s\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
+\lambda K_g(\delta w,w)
=0
}
$$

$$
\begin{aligned}
K_b(\delta\phi,\phi)
&=
D_{\alpha\beta\gamma\eta}\int_{\Omega}
\delta\kappa_{\alpha\beta}

\kappa_{\gamma\eta}\,d\Omega
\\
&=
\frac14
D_{\alpha\beta\gamma\eta}\sum_{IJ}
\int_{\Omega}
\left(
N_{J,\beta}\delta\phi_{\alpha J}
+N_{J,\alpha}\delta\phi_{\beta J}
\right)

\left(
N_{I,\eta}\phi_{\gamma I}
+N_{I,\gamma}\phi_{\eta I}
\right)d\Omega
\end{aligned}
$$

$$
\begin{aligned}
K_s\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
&=
S_{\alpha\beta}\int_{\Omega}
\delta\gamma_\alpha

\gamma_\beta\,d\Omega
\\
&=
K_s^{ww}(\delta w,w)
+K_s^{\psi\psi}(\delta\phi,\phi)
+K_s^{w\psi}\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
\end{aligned}
$$

$$
\begin{aligned}
K_s^{ww}(\delta w,w)
&=
S_{\alpha\beta}\sum_{IJ}
\int_{\Omega}
N_{J,\alpha}\delta w_J\,

N_{I,\beta}w_I\,d\Omega
\end{aligned}
$$

$$
\begin{aligned}
K_s^{\psi\psi}(\delta\phi,\phi)
&=
S_{\alpha\beta}\sum_{IJ}
\int_{\Omega}
N_J\delta\phi_{\alpha J}\,
N_I\phi_{\beta I}\,d\Omega
\end{aligned}
$$

$$
\begin{aligned}
K_s^{w\psi}\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
&=
-S_{\alpha\beta}\sum_{IJ}
\int_{\Omega}
N_{J,\alpha}\delta w_J\,

N_I\phi_{\beta I}\,d\Omega
\\
&\quad
-S_{\alpha\beta}\sum_{IJ}
\int_{\Omega}
N_J\delta\phi_{\alpha J}\,
N_{I,\beta}w_I\,d\Omega
\end{aligned}
$$

$$
\begin{aligned}
K_g(\delta w,w)
&=
\int_{\Omega}
\sigma_{\alpha\beta}^{\mathrm{ref}}
w_{,\alpha}\delta w_{,\beta}\,d\Omega
\\
&=
\sum_{IJ}
\int_{\Omega}
\sigma_{\alpha\beta}^{\mathrm{ref}}
N_{I,\alpha}N_{J,\beta}
w_I\delta w_J\,d\Omega
\end{aligned}
$$

### 總體離散方程

$$
\boxed{

\left[
K_b
+K_s^{ww}
+K_s^{\psi\psi}
+K_s^{w\psi}
+\lambda K_g
\right]
=0
}
$$


---
