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

\int_{\Omega} \Bigl( M_{\alpha\beta}\,\delta\phi_{\alpha,\beta} + Q_{\alpha}\,\delta\phi_\alpha + Q_{\alpha}\,\delta w_{,\alpha} - P_{\alpha\beta} w_{,\beta}\,\delta w_{,\alpha} \Bigr) d\Omega = 0
 \tag{3}
$$

$$
\int_{\Omega} \Bigl( M_{\alpha\beta}\,\delta\kappa_{\alpha\beta} + Q_{\alpha}\,\delta\gamma_{\alpha} - P_{\alpha\beta} w_{,\alpha}\,\delta w_{,\beta} \Bigr) d\Omega = 0
$$


记弯曲刚度贡献的虚功为 $a_b$，剪切刚度贡献的虚功为 $a_s$，几何刚度贡献的虚功为 $a_g$，则上式可简记为：

$$
a_b + a_s - a_g = 0
$$

## 6. 有限元离散

采用标准等参单元，形函数为 $N_K(x)$。记节点总自由度向量及其对应的虚位移向量为：

$$
\mathbf{d}_K =
\begin{bmatrix}
w_K \\
\phi_{1K} \\
\phi_{2K}
\end{bmatrix},
\qquad
\delta \mathbf{d}_K =
\begin{bmatrix}
\delta w_K \\
\delta\phi_{1K} \\
\delta\phi_{2K}
\end{bmatrix}
$$

相应的插值为：

$$
w^h(x) = \sum_{K} N_K(x) w_K,\qquad
\phi_{\alpha}^h (x)= \sum_{K} N_K(x) \phi_{\alpha K}
$$

相应的离散曲率、剪切应变和挠度梯度为：

$$
\kappa_{\alpha\beta}^h = -\frac{1}{2}\sum_{K}\bigl(N_{K,\beta}\phi_{\alpha K}+N_{K,\alpha}\phi_{\beta K}\bigr),\qquad
\gamma_{\alpha}^h = \sum_{K}\bigl(N_{K,\alpha} w_K - N_K \phi_{\alpha K}\bigr),\qquad
w_{,\alpha}^h = \sum_{K} N_{K,\alpha} w_K
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

定义弯曲转角子矩阵 $\mathbf{K}_{\phi\phi,R S}^b$ 为：

$$
\bigl(\mathbf{K}_{\phi\phi,R S}^b\bigr)_{ij}
= \frac{h^3}{48} \int_{\Omega} D_{\alpha\beta\gamma\eta}
\bigl( N_{R,\beta}\delta_{\alpha i} + N_{R,\alpha}\delta_{\beta i} \bigr)
\bigl( N_{S,\eta}\delta_{\gamma j} + N_{S,\gamma}\delta_{\eta j} \bigr)
d\Omega
$$

将其嵌入节点总自由度后，弯曲刚度块矩阵可写为：

$$
\mathbf{K}_{RS}^b =
\begin{bmatrix}
0 & \mathbf{0}_{1\times 2} \\
\mathbf{0}_{2\times 1} & \mathbf{K}_{\phi\phi,R S}^b
\end{bmatrix}
$$

因此

$$
a_b = \sum_{R,S} \delta\mathbf{d}_R^{T}\,\mathbf{K}_{RS}^b\,\mathbf{d}_S
$$

### 6.2 剪切刚度部分

$$
a_s = \int_{\Omega} k G h \, \gamma_{\alpha}^h \, \delta\gamma_{\alpha}^h \, d\Omega
$$

展开后可按 $ww$、$w\phi$ 和 $\phi\phi$ 三个独立区块整理为：

$$
K_{ww,KL}^s = kGh \int_{\Omega} N_{K,\alpha} N_{L,\alpha} \, d\Omega
$$

$$
\mathbf{K}_{w\phi,KL}^s
= -kGh \int_{\Omega}
\begin{bmatrix}
N_{K,1}N_L & N_{K,2}N_L
\end{bmatrix}
d\Omega
$$

$$
\mathbf{K}_{\phi\phi,KL}^s
= kGh \int_{\Omega} N_K N_L\,\mathbf{I}_2 \, d\Omega
$$

相应地，

$$
\mathbf{K}_{\phi w,KL}^s = \left(\mathbf{K}_{w\phi,LK}^s\right)^T
$$

因此节点 $(K,L)$ 的剪切刚度块矩阵为：

$$
\mathbf{K}_{KL}^s =
\begin{bmatrix}
K_{ww,KL}^s & \mathbf{K}_{w\phi,KL}^s \\
\mathbf{K}_{\phi w,KL}^s & \mathbf{K}_{\phi\phi,KL}^s
\end{bmatrix}
$$

从而

$$
a_s = \sum_{K,L} \delta\mathbf{d}_K^{T}\,\mathbf{K}_{KL}^s\,\mathbf{d}_L
$$

装配后的整体剪切刚度矩阵满足

$$
\mathbf{K}_{\phi w}^{s} = \left(\mathbf{K}_{w\phi}^{s}\right)^T
$$

### 6.3 几何刚度部分

$$
a_g = P_{\alpha\beta}\int_{\Omega} w_{,\alpha}\,\delta w_{,\beta}\,d\Omega
$$

离散后：

$$
a_g = \sum_{K,L} \delta w_K \left(P_{\alpha\beta}\int_{\Omega} N_{K,\alpha} N_{L,\beta} \, d\Omega \right) w_L
$$

定义几何刚度的 $ww$ 子块为：

$$
K_{ww,KL}^g = P_{\alpha\beta}\int_{\Omega} N_{K,\alpha} N_{L,\beta} \, d\Omega
$$

将其嵌入节点总自由度后，几何刚度块矩阵为：

$$
\mathbf{K}_{KL}^g =
\begin{bmatrix}
K_{ww,KL}^g & \mathbf{0}_{1\times 2} \\
\mathbf{0}_{2\times 1} & \mathbf{0}_{2\times 2}
\end{bmatrix}
$$

因此

$$
a_g = \sum_{K,L} \delta\mathbf{d}_K^{T}\,\mathbf{K}_{KL}^g\,\mathbf{d}_L
$$

### 6.4 总体离散方程

将上述节点块矩阵组装，并记整体自由度向量为

$$
\mathbf{d} =
\begin{bmatrix}
\mathbf{d}_w \\
\mathbf{d}_\phi
\end{bmatrix}
$$

则屈曲特征值问题的有限元方程为：

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
