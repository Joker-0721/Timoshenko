# Mindlin 板线性屈曲有限元分析


## 一、基本假设

板中面为 $x_3 = 0$，厚度 $h$ 常数。  
变量：
- 中面内位移：$u_\alpha(x_1,x_2)$，$\alpha=1,2$
- 横向挠度：$w(x_1,x_2)$
- 转角：$\phi_\alpha(x_1,x_2)$

位移场：
$$
\bar{u}_{\alpha}(x_1,x_2,x_3) = u_{\alpha}(x_1,x_2) - x_3\phi_{\alpha},\qquad \bar{u}_3 = w
$$

应变：

$$
\epsilon_{\alpha\beta} = \varepsilon_{\alpha\beta} + x_3\kappa_{\alpha\beta}
$$

中面应变与曲率：

$$
\varepsilon_{\alpha\beta} = \frac{1}{2}(u_{\alpha,\beta}+u_{\beta,\alpha}),\qquad
\kappa_{\alpha\beta} = -\frac{1}{2}(\phi_{\alpha,\beta}+\phi_{\beta,\alpha})
$$

横向剪切应变：

$$
\gamma_{\alpha} = 2\epsilon_{\alpha3} = w_{,\alpha} - \phi_{\alpha}
$$

---

## 二、本构关系

平面应力：$\sigma_{33}=0$。
$$
\sigma_{\alpha\beta} = D_{\alpha\beta\gamma\eta}\; \epsilon_{\gamma\eta},
$$
材料张量：
$$
D_{\alpha\beta\gamma\eta} = \frac{E}{1-\nu^2}\Bigl[ \nu\,\delta_{\alpha\beta}\delta_{\gamma\eta} + \frac{1-\nu}{2}(\delta_{\alpha\gamma}\delta_{\beta\eta}+\delta_{\alpha\eta}\delta_{\beta\gamma}) \Bigr].
$$

横向剪切：
$$
\sigma_{\alpha3} = G \,\gamma_\alpha,\qquad G = \frac{E}{2(1+\nu)}.
$$

---

## 三、弱形式

虚功：
$$
\delta W_{\text{int}} = \delta W_{\text{ext}}.
$$

三维内虚功：
$$
\delta W_{\text{int}} = \int_V \sigma_{ij}\,\delta\varepsilon_{ij}\,dV
= \int_V \bigl( \sigma_{\alpha\beta}\,\delta\varepsilon_{\alpha\beta} + \sigma_{\alpha3}\,\delta\gamma_\alpha \bigr) dV.
$$


代入 $\delta\varepsilon_{\alpha\beta}
=\delta\varepsilon_{\alpha\beta}^{(0)}+x_3\,\delta\kappa_{\alpha\beta}$：
$$
\begin{aligned}
\delta W_{\text{int}}
&=
\int_V
\sigma_{\alpha\beta}
\left(
\delta\varepsilon_{\alpha\beta}^{(0)}
+x_3\,\delta\kappa_{\alpha\beta}
\right)dV
+
\int_V
\sigma_{\alpha3}\,\delta\gamma_\alpha\,dV
\\
&=
\int_{\Omega} \left[\int_{-h/2}^{h/2} \sigma_{\alpha\beta}\,dx_3\right] \delta\varepsilon_{\alpha\beta}^{(0)} \, d\Omega
+
\int_{\Omega} \left[\int_{-h/2}^{h/2} \sigma_{\alpha\beta}\,x_3\,dx_3\right] \delta\kappa_{\alpha\beta} \, d\Omega
+
\int_{\Omega} \left[\int_{-h/2}^{h/2} \sigma_{\alpha3}\,dx_3\right] \delta\gamma_\alpha \, d\Omega .
\end{aligned}
$$

截面力：
- 膜力：$N_{\alpha\beta} = \int_{-h/2}^{h/2} \sigma_{\alpha\beta}\,dx_3$
- 弯矩：$M_{\alpha\beta} = \int_{-h/2}^{h/2} \sigma_{\alpha\beta}\,x_3\,dx_3$
- 剪力：$Q_\alpha = \int_{-h/2}^{h/2} \sigma_{\alpha3}\,dx_3$


#### 二維虚功形式

忽略膜力虚功项：
$$
\delta W_{\text{int}} = \int_{\Omega} \Bigl(\underbrace{\int_{-h/2}^{h/2}\sigma_{\alpha\beta}x_3dx_3}_{M_{\alpha\beta}} \delta\kappa_{\alpha\beta}
+ \underbrace{\int_{-h/2}^{h/2}\sigma_{\alpha3}dx_3}_{Q_\alpha} \delta\gamma_\alpha \Bigr) d\Omega.
$$

厚度积分，$k_s=\frac56$：
$$
\begin{aligned}
M_{\alpha\beta}
&= D^b_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta},
\qquad
D^b_{\alpha\beta\gamma\eta}
=\frac{h^3}{12}D_{\alpha\beta\gamma\eta},
\\
Q_{\alpha}
&= S_{\alpha\beta}\,\gamma_\beta,
\qquad
S_{\alpha\beta}=k_sGh\,\delta_{\alpha\beta}.
\end{aligned}
$$

线性双线性式：
$$
\delta W_{\text{int}} = \int_{\Omega} \Bigl(\delta\kappa_{\alpha\beta} D^b_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta}
+ \delta\gamma_\alpha S_{\alpha\beta}\,\gamma_\beta \Bigr) d\Omega.
$$

---

## 四、几何刚度

### 4.1 实际载荷与载荷因子
压缩为正：
$$
N_{\alpha\beta}^{0} = \lambda \,N_{\alpha\beta}^{\mathrm{ref}}.
$$
$\lambda$ 为载荷因子。

$N_{\alpha\beta}^{\mathrm{ref}}$ 为膜力合力；若给三维应力，先做厚度积分。

### 4.2 几何刚度的来源（非线性应变）
几何虚功：
$$
\delta W_{g}
=
-\int_{\Omega}
N_{\alpha\beta}^0\,
w_{,\alpha}\,\delta w_{,\beta}
\,d\Omega.
$$
代入 $N_{\alpha\beta}^0 = \lambda N_{\alpha\beta}^{\mathrm{ref}}$：
$$
\delta W_{g}
=
-\lambda
\int_{\Omega}
N_{\alpha\beta}^{\mathrm{ref}}\,
w_{,\alpha}\,\delta w_{,\beta}
\,d\Omega.
$$

### 4.3 屈曲时的总虚功平衡
总虚功：
$$
\delta W_{\text{int,lin}} + \delta W_g = 0.
$$
$$
\int_{\Omega} \Bigl(\delta\kappa_{\alpha\beta} D^b_{\alpha\beta\gamma\eta}\,\kappa_{\gamma\eta}
+ \delta\gamma_\alpha S_{\alpha\beta}\,\gamma_\beta \Bigr) d\Omega
- \lambda \int_{\Omega} N_{\alpha\beta}^{\mathrm{ref}} \; w_{,\alpha}\,\delta w_{,\beta} \, d\Omega = 0.
$$

---

## 五、有限元離散（雙線性形式）

以節點基函數 $N_I(x_1,x_2)$ 離散。

### 5.1 位移場離散化

位移與轉角：

$$
w_{,\alpha}=\sum_I N_{I,\alpha}w_I,\qquad
\delta w_{,\beta}=\sum_J N_{J,\beta}\delta w_J
$$

$$
\phi_{\alpha,\beta}=\sum_K N_{K,\beta}\phi_{\alpha K},\qquad
\delta\phi_{\alpha,\beta}=\sum_L N_{L,\beta}\delta\phi_{\alpha L}
$$

### 5.2 應變與曲率的離散表達

曲率與剪切應變：

$$
\kappa_{\alpha\beta} = -\frac{1}{2}(\phi_{\alpha,\beta}+\phi_{\beta,\alpha}) = -\frac{1}{2}\sum_I(N_{I,\beta}\phi_{\alpha I} + N_{I,\alpha}\phi_{\beta I})
$$

$$
\gamma_\alpha = w_{,\alpha} - \phi_\alpha = \sum_I N_{I,\alpha}w_I - \sum_I N_I\phi_{\alpha I}
$$

節點自由度：
$$
\mathbf d_I=
\begin{bmatrix}
w_I\\
\phi_{1I}\\
\phi_{2I}
\end{bmatrix}.
$$

### 5.3 總體虛功平衡（離散形式）

$$
\begin{aligned}
&
\int_{\Omega}
\Bigl(
\delta\kappa_{\alpha\beta}D^b_{\alpha\beta\gamma\eta}\kappa_{\gamma\eta}
+\delta\gamma_\alpha S_{\alpha\beta}\gamma_\beta
\Bigr)d\Omega
-\lambda\int_{\Omega}
N_{\alpha\beta}^{\mathrm{ref}}
w_{,\alpha}\delta w_{,\beta}\,d\Omega =0
\end{aligned}
$$


$$
\boxed{
a_b(\delta\phi,\phi)
+a_s\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
-\lambda K_g(\delta w,w)
=0
}
$$


**彎曲雙線性形式 $a_b$**：
$$
\begin{aligned}
a_b(\delta\phi,\phi)
&:=
\int_{\Omega}
\delta\kappa_{\alpha\beta}
D^b_{\alpha\beta\gamma\eta}
\kappa_{\gamma\eta}\,d\Omega
\\
&=
\sum_{IJ}
\delta\boldsymbol\phi_J^{T}
\mathbf K^{b}_{JI}
\boldsymbol\phi_I
\end{aligned}
$$


$$
\begin{aligned}
a_b(\delta\phi,\phi)
&=
\int_{\Omega}
\left[-\frac12\left(\delta\phi_{\alpha,\beta}+\delta\phi_{\beta,\alpha}\right)\right]
D^b_{\alpha\beta\gamma\eta}
\left[-\frac12\left(\phi_{\gamma,\eta}+\phi_{\eta,\gamma}\right)\right]d\Omega
\\
&=
\frac14\int_{\Omega}
\left(\delta\phi_{\alpha,\beta}+\delta\phi_{\beta,\alpha}\right)
D^b_{\alpha\beta\gamma\eta}
\left(\phi_{\gamma,\eta}+\phi_{\eta,\gamma}\right)d\Omega
\\
&=
\frac14\sum_{I,J}\int_{\Omega}
\left(N_{J,\beta}\delta\phi_{\alpha J}+N_{J,\alpha}\delta\phi_{\beta J}\right)
D^b_{\alpha\beta\gamma\eta}
\left(N_{I,\eta}\phi_{\gamma I}+N_{I,\gamma}\phi_{\eta I}\right)d\Omega \\
&=
\sum_{I,J}\delta\boldsymbol\phi_J^{T}
\mathbf K^b_{JI}
\boldsymbol\phi_I
\end{aligned}
$$

矩陣形式：
$$
\mathbf K_{JI}^{b}
=
\int_{\Omega}
\left(\mathbf B_J^b\right)^T
\mathbf D^b
\mathbf B_I^b\,d\Omega
$$

$\mathbf B_I^b$：
$$
\left(\mathbf B_I^b\boldsymbol\phi_I\right)_{\alpha\beta}
=
-\frac12
\left(
N_{I,\beta}\phi_{\alpha I}
+N_{I,\alpha}\phi_{\beta I}
\right).
$$

**剪切雙線性形式 $a_s$**：
$$
\begin{aligned}
a_s\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
&:=
\int_{\Omega}
\delta\gamma_\alpha
S_{\alpha\beta}
\gamma_\beta\,d\Omega
\\
&=
\sum_{IJ}
\delta\mathbf d_J^{T}
\mathbf K^{s}_{JI}
\mathbf d_I
\end{aligned}
$$

$$
\begin{aligned}
a_s\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
&=
\int_{\Omega}
\left(\delta w_{,\alpha}-\delta\phi_{\alpha}\right)
S_{\alpha\beta}
\left(w_{,\beta}-\phi_{\beta}\right)\,d\Omega
\\
&=
\int_{\Omega}
\left[\left(\sum_J N_{J,\alpha}\delta w_J\right)-\left(\sum_J N_J\delta\phi_{\alpha J}\right)\right]
S_{\alpha\beta}
\left[\left(\sum_I N_{I,\beta}w_I\right)-\left(\sum_I N_I\phi_{\beta I}\right)\right]d\Omega
\\
&=
\sum_{I,J}\int_{\Omega}
\left(N_{J,\alpha}\delta w_J-N_J\delta\phi_{\alpha J}\right)
S_{\alpha\beta}
\left(N_{I,\beta}w_I-N_I\phi_{\beta I}\right)d\Omega
\\
&=:
a_s^{ww}(\delta w,w)
+a_s^{\phi\phi}(\delta\phi,\phi)
+a_s^{w\phi}\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
\\

a_s^{ww}(\delta w,w)
&=
\sum_{I,J}S_{\alpha\beta}\,\int_{\Omega}
N_{J,\alpha}\,N_{I,\beta}
\,\delta w_J w_I\,d\Omega
\\
a_s^{\phi\phi}(\delta\phi,\phi)
&=
\sum_{I,J}\int_{\Omega}
N_J\,S_{\alpha\beta}\,N_I
\,\delta\phi_{\alpha J}\phi_{\beta I}\,d\Omega
\\
a_s^{w\phi}\bigl((\delta w,\delta\phi),(w,\phi)\bigr)
&=
-\sum_{I,J}\int_{\Omega}
N_{J,\alpha}\,S_{\alpha\beta}\,N_I
\,\delta w_J\phi_{\beta I}\,d\Omega
\\
&\quad
-\sum_{I,J}\int_{\Omega}
N_J\,S_{\alpha\beta}\,N_{I,\beta}
\,\delta\phi_{\alpha J}w_I\,d\Omega
\\
&=
\sum_{I,J}
\delta\mathbf d_J^{T}
\mathbf K^{s}_{JI}
\mathbf d_I
\end{aligned}
$$


$$
\mathbf K_{JI}^{s}
=
\int_{\Omega}
\left(\mathbf B_J^s\right)^T
\mathbf S
\mathbf B_I^s\,d\Omega
$$

$\mathbf B_I^s$：
$$
\mathbf B_I^s=
\begin{bmatrix}
N_{I,1} & -N_I & 0\\
N_{I,2} & 0 & -N_I
\end{bmatrix},
\qquad
\mathbf S=k_sGh
\begin{bmatrix}
1 & 0\\
0 & 1
\end{bmatrix}.
$$

**幾何刚度 $K_g$**：

$$
\begin{aligned}
K_g(\delta w,w)
&=
\int_{\Omega}
N_{\alpha\beta}^{\mathrm{ref}}
w_{,\alpha}\delta w_{,\beta}\,d\Omega
\\
&=
\sum_{IJ}
\int_{\Omega}
N_{\alpha\beta}^{\mathrm{ref}}
N_{I,\alpha}N_{J,\beta}
w_I\delta w_J\,d\Omega
\end{aligned}
$$

#### 5.4 幾何剛度（向量形式）

向量形式：
$$
K_g(\delta w,w)=\sum_{IJ}\delta\mathbf w_J^T\,\mathbf K_{JI}^g\,\mathbf w_I,
\qquad
\mathbf K_{JI}^{g}=
\int_{\Omega}
\left(\mathbf B_J^g\right)^T
\mathbf N^{\mathrm{ref}}
\mathbf B_I^g\,d\Omega
$$

$$
\mathbf B_I^g=
\begin{bmatrix}
N_{I,1}\\
N_{I,2}
\end{bmatrix},
\qquad
\mathbf N^{\mathrm{ref}}=
\begin{bmatrix}
N_{11}^{\mathrm{ref}} & N_{12}^{\mathrm{ref}}\\
N_{12}^{\mathrm{ref}} & N_{22}^{\mathrm{ref}}
\end{bmatrix}
$$

### 總體離散方程

屈曲特徵值問題：

$$
\boxed{

\left[
\mathbf K^b
+\mathbf K^s
-\lambda\mathbf K^g
\right]\mathbf d=\mathbf 0
}
$$

等價形式：
$$
\boxed{
\left(\mathbf K^b+\mathbf K^s\right)\mathbf d
=
\lambda\mathbf K^g\mathbf d
}
$$
