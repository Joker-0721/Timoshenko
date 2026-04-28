# Mindlin-Reissner 薄板屈曲之變分弱形式與有限元離散化

## 1. 基本運動學假設

$$
e_{ij}^{(m)} = e_{ij} + \frac{1}{2} u_{k,i} u_{k,j} - \left( \frac{1}{2} - \frac{m}{4} \right) e_{ki} e_{kj}
$$

$$
m = 2 \quad \text{(Engesser 假設)}
$$

## 2. 位移場與運動學

$$
u_1(x_1, x_2, x_3) = u^0_1(x_1, x_2) + x_3 \theta_1(x_1, x_2)
$$

$$
u_2(x_1, x_2, x_3) = u^0_2(x_1, x_2) + x_3 \theta_2(x_1, x_2)
$$

$$
u_3(x_1, x_2, x_3) = w(x_1, x_2)
$$

$$
\mathbf{u} = [w, \theta_1, \theta_2]^T
$$

$$
\kappa_{\alpha\beta} = \frac{1}{2} (\theta_{\alpha,\beta} + \theta_{\beta,\alpha})
$$

$$
\gamma_{\alpha 3} = \theta_\alpha + w_{,\alpha}
$$

## 3. 本構關係

$$
M_{\alpha\beta} = D_{\alpha\beta\gamma\delta} \kappa_{\gamma\delta}
$$

$$
Q_\alpha = H_{\alpha\gamma} \gamma_{\gamma 3}
$$

$$
D_{\alpha\beta\gamma\delta} = \frac{h^3}{12} C_{\alpha\beta\gamma\delta}
$$

$$
H_{\alpha\gamma} = \kappa_s h C_{\alpha 3\gamma 3}
$$

## 4. 強形式 (控制微分方程)

$$
M_{\alpha\beta,\beta} - Q_\alpha = 0
$$

$$
Q_{\alpha,\alpha} + (N^0_{\alpha\beta} w_{,\beta})_{,\alpha} + q = 0 \quad (N^0_{\alpha\beta,\alpha} = 0)
$$

$$
D_{\alpha\beta\gamma\delta} \theta_{\gamma,\delta\beta} - H_{\alpha\gamma} (\theta_\gamma + w_{,\gamma}) = 0
$$

$$
H_{\alpha\gamma} (\theta_{\gamma,\alpha} + w_{,\gamma\alpha}) + N^0_{\alpha\beta} w_{,\alpha\beta} + q = 0
$$

## 5. 弱形式 (變分形式)

$$
\int_A (M_{\alpha\beta,\beta} - Q_\alpha) \delta \theta_\alpha dA + \int_A (Q_{\alpha,\alpha} + N^0_{\alpha\beta} w_{,\alpha\beta} + q) \delta w dA = 0
$$

$$
\int_A M_{\alpha\beta,\beta} \delta \theta_\alpha dA = \oint_{\partial A} M_{\alpha\beta} n_\beta \delta \theta_\alpha ds - \int_A M_{\alpha\beta} \delta \theta_{\alpha,\beta} dA
$$

$$
\int_A Q_{\alpha,\alpha} \delta w dA = \oint_{\partial A} Q_\alpha n_\alpha \delta w ds - \int_A Q_\alpha \delta w_{,\alpha} dA
$$

$$
\int_A N^0_{\alpha\beta} w_{,\alpha\beta} \delta w dA = \oint_{\partial A} N^0_{\alpha\beta} w_{,\beta} n_\alpha \delta w ds - \int_A N^0_{\alpha\beta} w_{,\beta} \delta w_{,\alpha} dA
$$

$$
-\int_A \left( M_{\alpha\beta} \delta \theta_{\alpha,\beta} + Q_\alpha (\delta \theta_\alpha + \delta w_{,\alpha}) + N^0_{\alpha\beta} w_{,\beta} \delta w_{,\alpha} \right) dA + \int_A q \delta w dA + \oint_{\partial A} (...) = 0
$$

$$
\int_A \left( \delta \theta_{\alpha,\beta} D_{\alpha\beta\gamma\delta} \theta_{\gamma,\delta} \right) dA + \int_A (\delta \theta_\alpha + \delta w_{,\alpha}) H_{\alpha\gamma} (\theta_\gamma + w_{,\gamma}) dA + \int_A \delta w_{,\alpha} N^0_{\alpha\beta} w_{,\beta} dA - \int_A \delta w q dA = 0
$$

## 6. 有限元離散化

$$
w(x_1, x_2) = \sum_{I=1}^{n} N_I(x_1, x_2) w_I
$$

$$
\theta_1(x_1, x_2) = \sum_{I=1}^{n} N_I(x_1, x_2) \theta_{1I}
$$

$$
\theta_2(x_1, x_2) = \sum_{I=1}^{n} N_I(x_1, x_2) \theta_{2I}
$$

$$
\mathbf{d}_I = \begin{bmatrix} w_I \\ \theta_{1I} \\ \theta_{2I} \end{bmatrix}
$$

$$
\mathbf{J} = \begin{bmatrix} J_{11} & J_{12} \\ J_{21} & J_{22} \end{bmatrix} = \begin{bmatrix} \sum_{I=1}^n \frac{\partial N_I}{\partial \xi} x_{1I} & \sum_{I=1}^n \frac{\partial N_I}{\partial \xi} x_{2I} \\ \sum_{I=1}^n \frac{\partial N_I}{\partial \eta} x_{1I} & \sum_{I=1}^n \frac{\partial N_I}{\partial \eta} x_{2I} \end{bmatrix}
$$

$$
\begin{bmatrix} N_{I,1} \\ N_{I,2} \end{bmatrix} = \begin{bmatrix} \frac{\partial N_I}{\partial x_1} \\ \frac{\partial N_I}{\partial x_2} \end{bmatrix} = \mathbf{J}^{-1} \begin{bmatrix} \frac{\partial N_I}{\partial \xi} \\ \frac{\partial N_I}{\partial \eta} \end{bmatrix}
$$

## 7. 顯式矩陣推導 ($\mathbf{K}^e \mathbf{d}^e = \mathbf{F}^e$)

$$
\mathbf{K}^e = \mathbf{K}_b^e + \mathbf{K}_s^e + \mathbf{K}_g^e
$$

### 彎曲剛度矩陣 ($\mathbf{K}_b^e$)

$$
\boldsymbol{\kappa} = [ \kappa_{11}, \kappa_{22}, 2\kappa_{12} ]^T
$$

$$
\mathbf{B}_{bI} = \begin{bmatrix} 0 & N_{I,1} & 0 \\ 0 & 0 & N_{I,2} \\ 0 & N_{I,2} & N_{I,1} \end{bmatrix}
$$

$$
\mathbf{D}_b = \begin{bmatrix} D_{11} & D_{12} & 0 \\ D_{12} & D_{22} & 0 \\ 0 & 0 & D_{66} \end{bmatrix}
$$

$$
\mathbf{K}_{b,IJ} = \int_{A^e} \mathbf{B}_{bI}^T \mathbf{D}_b \mathbf{B}_{bJ} dA
$$

$$
\mathbf{K}_{b,IJ} = \int_{A^e} \begin{bmatrix} 0 & 0 & 0 \\ 0 & D_{11} N_{I,1} N_{J,1} + D_{66} N_{I,2} N_{J,2} & D_{12} N_{I,1} N_{J,2} + D_{66} N_{I,2} N_{J,1} \\ 0 & D_{12} N_{I,2} N_{J,1} + D_{66} N_{I,1} N_{J,2} & D_{22} N_{I,2} N_{J,2} + D_{66} N_{I,1} N_{J,1} \end{bmatrix} dA
$$

### 橫向剪切剛度矩陣 ($\mathbf{K}_s^e$)

$$
\boldsymbol{\gamma} = [ \gamma_{13}, \gamma_{23} ]^T
$$

$$
\mathbf{B}_{sI} = \begin{bmatrix} N_{I,1} & N_I & 0 \\ N_{I,2} & 0 & N_I \end{bmatrix}
$$

$$
\mathbf{D}_s = \begin{bmatrix} H_{44} & 0 \\ 0 & H_{55} \end{bmatrix}
$$

$$
\mathbf{K}_{s,IJ} = \int_{A^e} \mathbf{B}_{sI}^T \mathbf{D}_s \mathbf{B}_{sJ} dA
$$

$$
\mathbf{K}_{s,IJ} = \int_{A^e} \begin{bmatrix} H_{44} N_{I,1} N_{J,1} + H_{55} N_{I,2} N_{J,2} & H_{44} N_{I,1} N_J & H_{55} N_{I,2} N_J \\ H_{44} N_I N_{J,1} & H_{44} N_I N_J & 0 \\ H_{55} N_I N_{J,2} & 0 & H_{55} N_I N_J \end{bmatrix} dA
$$

### 幾何剛度矩陣 ($\mathbf{K}_g^e$)

$$
\mathbf{N}^0 = \begin{bmatrix} N^0_{11} & N^0_{12} \\ N^0_{21} & N^0_{22} \end{bmatrix}
$$

$$
\mathbf{g} = [ w_{,1}, w_{,2} ]^T
$$

$$
\mathbf{B}_{gI} = \begin{bmatrix} N_{I,1} & 0 & 0 \\ N_{I,2} & 0 & 0 \end{bmatrix}
$$

$$
\mathbf{K}_{g,IJ} = \int_{A^e} \mathbf{B}_{gI}^T \mathbf{N}^0 \mathbf{B}_{gJ} dA
$$

$$
\mathbf{K}_{g,IJ} = \int_{A^e} \begin{bmatrix} N_{I,1} N^0_{11} N_{J,1} + N_{I,1} N^0_{12} N_{J,2} + N_{I,2} N^0_{21} N_{J,1} + N_{I,2} N^0_{22} N_{J,2} & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix} dA
$$

### 等效節點力向量 ($\mathbf{F}^e$)

$$
\mathbf{F}_I = \int_{A^e} \begin{bmatrix} N_I q \\ 0 \\ 0 \end{bmatrix} dA
$$

## 8. 全局組裝與特徵值問題

$$
\mathbf{K}_E = \sum_{e=1}^{E_{total}} \left( \mathbf{K}_b^{(e)} + \mathbf{K}_s^{(e)} \right)
$$

$$
\mathbf{K}_G = \sum_{e=1}^{E_{total}} \mathbf{K}_g^{(e)}
$$

$$
\mathbf{F}_{ext} = \sum_{e=1}^{E_{total}} \mathbf{F}^{(e)}
$$

$$
(\mathbf{K}_E + \mathbf{K}_G) \mathbf{D} = \mathbf{F}_{ext}
$$

$$
\mathbf{F}_{ext} = \mathbf{0}
$$

$$
(\mathbf{K}_E + \lambda \mathbf{K}_G) \mathbf{D}_{mode} = \mathbf{0}
$$

$$
\det | \mathbf{K}_E + \lambda \mathbf{K}_G | = 0
$$

---

Source: https://gemini.google.com/app/c351366a5ed02a12
Exported at: 2026-04-07T02:12:37.745Z