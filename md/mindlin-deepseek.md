
## 1 位移场与几何关系

### 1.1 位移场（Timoshenko/Engesser 板理论）

设板中面位于 $z=0$，厚度为 $h$。一阶剪切变形理论下板内任意点位移：



$u_1(x_1,x_2,z) = u(x_1,x_2) + z\,\psi_1(x_1,x_2)\\ $
$u_2(x_1,x_2,z) = v(x_1,x_2) + z\,\psi_2(x_1,x_2)\\$
$u_3(x_1,x_2,z) = w(x_1,x_2)$



其中 $u, v, w$ 为中面位移，$\psi_1, \psi_2$ 为中面法线独立转角。

### 1.2 应变-位移关系

线性应变张量 $\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$，分离面内与横向：

**面内应变：**

$$
\varepsilon_{\alpha\beta} = \varepsilon_{\alpha\beta}^{(0)} + z\,\kappa_{\alpha\beta}
\quad (\alpha,\beta = 1,2)
$$
膜应变与曲率定义为：

$$
\varepsilon_{\alpha\beta}^{(0)} = \frac{1}{2}(u_{\alpha,\beta} + u_{\beta,\alpha}),\qquad
\kappa_{\alpha\beta} = \frac{1}{2}(\psi_{\alpha,\beta} + \psi_{\beta,\alpha})
$$
**横向剪切应变：**

$$
\gamma_{\alpha 3} = u_{\alpha,3} + u_{3,\alpha} = \psi_\alpha + w_{,\alpha}
$$
---

## 2 本构关系

各向同性线弹性材料，平面应力状态：

**面内应力：**

$$
\sigma_{\alpha\beta} = \frac{E}{1-\nu^2}\Big[(1-\nu)\varepsilon_{\alpha\beta} + \nu\,\varepsilon_{\gamma\gamma}\,\delta_{\alpha\beta}\Big]
$$
**横向剪切应力：**

$$
\tau_{\alpha 3} = k_s G\,\gamma_{\alpha 3}
$$
其中 $k_s$ 为剪切修正系数，$G = \dfrac{E}{2(1+\nu)}$。

---

## 3 势能泛函

系统总势能由应变能与外力势能组成。

### 3.1 应变能

应变能密度沿厚度积分，得单位中面面积应变能：

$$
\tilde{U} = \frac{1}{2}\Big[ N_{\alpha\beta}\varepsilon_{\alpha\beta}^{(0)} + M_{\alpha\beta}\kappa_{\alpha\beta} + Q_\alpha \gamma_{\alpha 3} \Big]
$$
内力定义：

$$
N_{\alpha\beta} = \int_{-h/2}^{h/2} \sigma_{\alpha\beta}\,dz,\quad
M_{\alpha\beta} = \int_{-h/2}^{h/2} \sigma_{\alpha\beta}\,z\,dz,\quad
Q_\alpha = \int_{-h/2}^{h/2} \tau_{\alpha 3}\,dz
$$
代入本构关系得内力-广义应变关系：

$$
\begin{aligned}
N_{\alpha\beta} &= C_{\alpha\beta\gamma\delta}^{(m)}\,\varepsilon_{\gamma\delta}^{(0)} \\
M_{\alpha\beta} &= C_{\alpha\beta\gamma\delta}^{(b)}\,\kappa_{\gamma\delta} \\
Q_\alpha &= C_{\alpha\beta}^{(s)}\,\gamma_{\beta 3}
\end{aligned}
$$
刚度系数：

$$
C_{\alpha\beta\gamma\delta}^{(m)} = h\,D_{\alpha\beta\gamma\delta},\quad
C_{\alpha\beta\gamma\delta}^{(b)} = \frac{h^3}{12}\,D_{\alpha\beta\gamma\delta},\quad
C_{\alpha\beta}^{(s)} = k_s G h\,\delta_{\alpha\beta}
$$
$$
D_{\alpha\beta\gamma\delta} = \frac{E}{1-\nu^2}\Big[ \nu\,\delta_{\alpha\beta}\delta_{\gamma\delta} + \frac{1-\nu}{2}(\delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma}) \Big]
$$
总应变能：

$$
U = \frac{1}{2}\int_\Omega \Big[ C_{\alpha\beta\gamma\delta}^{(m)}\varepsilon_{\alpha\beta}^{(0)}\varepsilon_{\gamma\delta}^{(0)} + C_{\alpha\beta\gamma\delta}^{(b)}\kappa_{\alpha\beta}\kappa_{\gamma\delta} + C_{\alpha\beta}^{(s)}\gamma_{\alpha 3}\gamma_{\beta 3} \Big] d\Omega
$$
### 3.2 外力势能（屈曲荷载）

设板受面内压力 $N_{\alpha\beta}^0$（受压为正），von Kármán 非线性应变贡献外力功：

$$
W_{\text{ext}} = \frac{1}{2}\int_\Omega N_{\alpha\beta}^0\, w_{,\alpha}w_{,\beta}\, d\Omega
$$
### 3.3 总势能泛函

$$
\Pi = U - W_{\text{ext}} = \frac{1}{2}\int_\Omega \Big[ C_{\alpha\beta\gamma\delta}^{(m)}\varepsilon_{\alpha\beta}^{(0)}\varepsilon_{\gamma\delta}^{(0)} + C_{\alpha\beta\gamma\delta}^{(b)}\kappa_{\alpha\beta}\kappa_{\gamma\delta} + C_{\alpha\beta}^{(s)}\gamma_{\alpha 3}\gamma_{\beta 3} - N_{\alpha\beta}^0\, w_{,\alpha}w_{,\beta} \Big] d\Omega
$$
---

## 4 强形式（控制微分方程）的推导

由最小势能原理 $\delta\Pi = 0$，取一阶变分：

$$
\begin{aligned}
\delta\Pi = \int_\Omega \Big[ 
& N_{\alpha\beta}\,\delta\varepsilon_{\alpha\beta}^{(0)} 
+ M_{\alpha\beta}\,\delta\kappa_{\alpha\beta} 
+ Q_\alpha\,(\delta\psi_\alpha + \delta w_{,\alpha}) 
- N_{\alpha\beta}^0\, w_{,\alpha}\delta w_{,\beta}
\Big] d\Omega = 0
\end{aligned}
$$
将变分写成偏导形式并利用对称性进行**分部积分**，将微分从虚位移转移至内力。

**对 $\delta u_\alpha$ 项：**

$$
\int_\Omega N_{\alpha\beta}\,\delta u_{\alpha,\beta}\,d\Omega = \int_\Gamma N_{\alpha\beta}n_\beta\,\delta u_\alpha\,d\Gamma - \int_\Omega N_{\alpha\beta,\beta}\,\delta u_\alpha\,d\Omega
$$
**对 $\delta \psi_\alpha$ 项：**

$$
\int_\Omega \Big[ M_{\alpha\beta}\,\delta\psi_{\alpha,\beta} + Q_\alpha\,\delta\psi_\alpha \Big] d\Omega = \int_\Gamma M_{\alpha\beta}n_\beta\,\delta\psi_\alpha\,d\Gamma - \int_\Omega (M_{\alpha\beta,\beta} - Q_\alpha)\,\delta\psi_\alpha\,d\Omega
$$
**对 $\delta w$ 项：**

$$
\int_\Omega \Big[ Q_\alpha\,\delta w_{,\alpha} - N_{\alpha\beta}^0\, w_{,\alpha}\delta w_{,\beta} \Big] d\Omega = \int_\Gamma \big( Q_\alpha - N_{\alpha\beta}^0 w_{,\beta} \big) n_\alpha\,\delta w\,d\Gamma - \int_\Omega \big( Q_{\alpha,\alpha} - (N_{\alpha\beta}^0 w_{,\beta})_{,\alpha} \big)\,\delta w\,d\Omega
$$
由 $\delta u_\alpha, \delta w, \delta \psi_\alpha$ 在域内的任意性，得**强形式平衡方程**：

$$
\boxed{
\begin{aligned}
& N_{\alpha\beta,\beta} = 0 \quad \text{在 } \Omega \text{ 内} \\
& M_{\alpha\beta,\beta} - Q_\alpha = 0 \quad \text{在 } \Omega \text{ 内} \\
& Q_{\alpha,\alpha} - (N_{\alpha\beta}^0 w_{,\beta})_{,\alpha} = 0 \quad \text{在 } \Omega \text{ 内}
\end{aligned}
}
$$
自然边界条件：
- 力边界 $\Gamma_N$ 上：$N_{\alpha\beta}n_\beta = \bar{N}_\alpha$，$M_{\alpha\beta}n_\beta = \bar{M}_\alpha$，$(Q_\alpha - N_{\alpha\beta}^0 w_{,\beta})n_\alpha = \bar{Q}$。

---

## 5 弱形式（虚功方程）及其双线性形式

将第 4 节变分积分式**保留原始形式，不作分部积分**，即得弱形式：

$$
\boxed{
\int_\Omega \Big[ 
N_{\alpha\beta}\,\delta\varepsilon_{\alpha\beta}^{(0)} 
+ M_{\alpha\beta}\,\delta\kappa_{\alpha\beta} 
+ Q_\alpha\,(\delta\psi_\alpha + \delta w_{,\alpha}) 
- N_{\alpha\beta}^0\, w_{,\alpha}\delta w_{,\beta}
\Big] d\Omega = 0
}
$$
### 5.1 双线性形式分解

弱形式可视为由四个独立的双线性泛函之和组成。记位移场变分为 $\delta\mathbf{u} = (\delta u, \delta v, \delta w, \delta \psi_1, \delta \psi_2)$，真实位移场为 $\mathbf{u}$。各双线性项如下：

**（1）膜应变双线性形式 $a_m(\mathbf{u}, \delta\mathbf{u})$：**

$$
a_m(\mathbf{u}, \delta\mathbf{u}) = \int_\Omega C_{\alpha\beta\gamma\delta}^{(m)} \varepsilon_{\alpha\beta}^{(0)}(\mathbf{u})\, \delta\varepsilon_{\gamma\delta}^{(0)}(\delta\mathbf{u}) \, d\Omega
$$
**（2）弯曲双线性形式 $a_b(\mathbf{u}, \delta\mathbf{u})$：**

$$
a_b(\mathbf{u}, \delta\mathbf{u}) = \int_\Omega C_{\alpha\beta\gamma\delta}^{(b)} \kappa_{\alpha\beta}(\mathbf{u})\, \delta\kappa_{\gamma\delta}(\delta\mathbf{u}) \, d\Omega
$$
**（3）剪切双线性形式 $a_s(\mathbf{u}, \delta\mathbf{u})$：**

$$
a_s(\mathbf{u}, \delta\mathbf{u}) = \int_\Omega C_{\alpha\beta}^{(s)} \gamma_{\alpha 3}(\mathbf{u})\, \delta\gamma_{\beta 3}(\delta\mathbf{u}) \, d\Omega
$$
**（4）几何刚度双线性形式 $a_g(\mathbf{u}, \delta\mathbf{u})$（屈曲项）：**

$$
a_g(\mathbf{u}, \delta\mathbf{u}) = \int_\Omega N_{\alpha\beta}^0 \, w_{,\alpha}(\mathbf{u})\, \delta w_{,\beta}(\delta\mathbf{u}) \, d\Omega
$$
弱形式即：

$$
a_m(\mathbf{u}, \delta\mathbf{u}) + a_b(\mathbf{u}, \delta\mathbf{u}) + a_s(\mathbf{u}, \delta\mathbf{u}) - \lambda\, a_g(\mathbf{u}, \delta\mathbf{u}) = 0
$$
其中 $\lambda$ 为荷载因子，$N_{\alpha\beta}^0 = \lambda \bar{N}_{\alpha\beta}^0$。

---

## 6 有限元离散：代入位移插值

将板域 $\Omega$ 离散为 $n_e$ 个单元。单元内位移场用形函数插值：

$$
\mathbf{u}^e(\mathbf{x}) = \mathbf{N}^e(\mathbf{x})\,\mathbf{d}^e
$$
其中 $\mathbf{u}^e = [u,\,v,\,w,\,\psi_1,\,\psi_2]^T$，$\mathbf{d}^e$ 为节点自由度向量。形函数矩阵按节点 $i=1,\dots,n_{\text{node}}$ 组装：

$$
\mathbf{N}^e =
\begin{bmatrix}
N_i & 0 & 0 & 0 & 0 \\
0 & N_i & 0 & 0 & 0 \\
0 & 0 & N_i & 0 & 0 \\
0 & 0 & 0 & N_i & 0 \\
0 & 0 & 0 & 0 & N_i
\end{bmatrix}_{5 \times 5n_{\text{node}}}
$$
### 6.1 广义应变-节点位移矩阵

由几何关系与插值公式，定义各应变-位移矩阵：

**膜应变向量** $\boldsymbol{\varepsilon}^{(0)} = [\varepsilon_{11}^{(0)},\ \varepsilon_{22}^{(0)},\ 2\varepsilon_{12}^{(0)}]^T$：

$$
\boldsymbol{\varepsilon}^{(0)} = \mathbf{B}_m^e\,\mathbf{d}^e, \quad
\mathbf{B}_m^e = 
\begin{bmatrix}
N_{i,1} & 0 & 0 & 0 & 0 \\
0 & N_{i,2} & 0 & 0 & 0 \\
N_{i,2} & N_{i,1} & 0 & 0 & 0
\end{bmatrix}
$$
**曲率向量** $\boldsymbol{\kappa} = [\kappa_{11},\ \kappa_{22},\ 2\kappa_{12}]^T$：

$$
\boldsymbol{\kappa} = \mathbf{B}_b^e\,\mathbf{d}^e, \quad
\mathbf{B}_b^e = 
\begin{bmatrix}
0 & 0 & 0 & N_{i,1} & 0 \\
0 & 0 & 0 & 0 & N_{i,2} \\
0 & 0 & 0 & N_{i,2} & N_{i,1}
\end{bmatrix}
$$
**剪切应变向量** $\boldsymbol{\gamma} = [\gamma_{13},\ \gamma_{23}]^T$：

$$
\boldsymbol{\gamma} = \mathbf{B}_s^e\,\mathbf{d}^e, \quad
\mathbf{B}_s^e = 
\begin{bmatrix}
0 & 0 & N_{i,1} & N_i & 0 \\
0 & 0 & N_{i,2} & 0 & N_i
\end{bmatrix}
$$
**挠度梯度向量** $\nabla w = [w_{,1},\ w_{,2}]^T$：

$$
\nabla w = \mathbf{B}_g^e\,\mathbf{d}^e, \quad
\mathbf{B}_g^e = 
\begin{bmatrix}
0 & 0 & N_{i,1} & 0 & 0 \\
0 & 0 & N_{i,2} & 0 & 0
\end{bmatrix}
$$
### 6.2 单元刚度矩阵的显式表达

将上述离散化代入各双线性形式，得单元矩阵。

**膜刚度矩阵：**

$$
\mathbf{K}_m^e = \int_{\Omega_e} (\mathbf{B}_m^e)^T \mathbf{D}_m \mathbf{B}_m^e \, d\Omega
$$
其中 $\mathbf{D}_m = h\,\mathbf{D}_0$，平面应力弹性矩阵 $\mathbf{D}_0$ 为：

$$
\mathbf{D}_0 = \frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{bmatrix}
$$
**弯曲刚度矩阵：**

$$
\mathbf{K}_b^e = \int_{\Omega_e} (\mathbf{B}_b^e)^T \mathbf{D}_b \mathbf{B}_b^e \, d\Omega,\quad
\mathbf{D}_b = \frac{h^3}{12}\mathbf{D}_0
$$
**剪切刚度矩阵：**

$$
\mathbf{K}_s^e = \int_{\Omega_e} (\mathbf{B}_s^e)^T \mathbf{D}_s \mathbf{B}_s^e \, d\Omega,\quad
\mathbf{D}_s = k_s G h
\begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}
$$
**几何刚度矩阵：**

$$
\mathbf{K}_g^e = \int_{\Omega_e} (\mathbf{B}_g^e)^T \mathbf{N}^0 \mathbf{B}_g^e \, d\Omega
$$
其中初始面内应力矩阵（以受压为正）：

$$
\mathbf{N}^0 = 
\begin{bmatrix}
N_{11}^0 & N_{12}^0 \\
N_{21}^0 & N_{22}^0
\end{bmatrix}
$$
### 6.3 全局组装与特征值方程

单元刚度矩阵为：

$$
\mathbf{K}^e = \mathbf{K}_m^e + \mathbf{K}_b^e + \mathbf{K}_s^e
$$
组装全局矩阵：

$$
\mathbf{K} = \operatorname*{\mathbf{A}}_{e=1}^{n_e} \mathbf{K}^e,\qquad
\mathbf{K}_G = \operatorname*{\mathbf{A}}_{e=1}^{n_e} \mathbf{K}_g^e
$$
总势能离散形式：

$$
\Pi(\mathbf{d}) = \frac{1}{2}\mathbf{d}^T\mathbf{K}\mathbf{d} - \frac{\lambda}{2}\mathbf{d}^T\mathbf{K}_G\mathbf{d}
$$
由 $\delta\Pi = 0$ 得屈曲特征值方程：

$$
\boxed{
(\mathbf{K} - \lambda \mathbf{K}_G)\,\mathbf{d} = \mathbf{0}
}
$$
临界屈曲荷载因子为最小特征值 $\lambda_{\text{cr}}$。

