# `thick_plate.jl` 公式對照說明

本文對照 [`ApproxOperator.jl/src/operation/thick_plate.jl`](</d:/ApproxOperator.jl/src/operation/thick_plate.jl>)，把程式中的 operator 回寫成 Mindlin-Reissner 厚板的弱式、離散式與矩陣 block。重點不是重做一本板殼教材，而是讓讀者能從公式直接查到對應的 Julia operator，反過來也能從 operator 看出它在總矩陣中的位置。

## 1. 模型、未知量與自由度順序

厚板的主未知量採用

$$
w(x,y), \qquad \phi_1(x,y), \qquad \phi_2(x,y),
$$

其中 $w$ 是橫向位移，$\phi_1,\phi_2$ 是截面轉角。

程式在 combined form `∫κMγQdΩ` 中使用每節點三個自由度，順序為

$$
\mathbf d_I = \begin{bmatrix} w_I & \phi_{1I} & \phi_{2I} \end{bmatrix}^{\mathsf T}.
$$

若採 block form，則把整體矩陣拆成

$$
\mathbf d =
\begin{bmatrix}
\mathbf d_w \\
\mathbf d_\phi
\end{bmatrix},
\qquad
\mathbf d_\phi =
\begin{bmatrix}
\mathbf d_{\phi_1} \\
\mathbf d_{\phi_2}
\end{bmatrix}.
$$

對於 mixed / auxiliary operator，還會再引入獨立剪力與彎矩變數

$$
\mathbf Q = \begin{bmatrix} Q_1 & Q_2 \end{bmatrix}^{\mathsf T},
\qquad
\mathbf M = \begin{bmatrix} M_{11} & M_{22} & M_{12} \end{bmatrix}^{\mathsf T}.
$$

## 2. 運動學與材料常數

Mindlin 厚板的曲率與剪應變可寫成

$$
\kappa_{11} = -\phi_{1,1}, \qquad
\kappa_{22} = -\phi_{2,2}, \qquad
\kappa_{12} = -\frac{1}{2}(\phi_{1,2}+\phi_{2,1}),
$$

$$
\gamma_1 = w_{,1} - \phi_1, \qquad
\gamma_2 = w_{,2} - \phi_2.
$$

在等向性厚板中，程式使用的彎曲與剪切剛度常數為

$$
D_{1111} = \frac{Eh^3}{12(1-\nu^2)}, \qquad
D_{1122} = \frac{\nu Eh^3}{12(1-\nu^2)}, \qquad
D_{1212} = \frac{Eh^3}{24(1+\nu)},
$$

$$
D_s = \frac{5}{6}\frac{Eh}{2(1+\nu)}.
$$

這正是 `thick_plate.jl` 內的

- `E` 對應 Young's modulus $E$
- `谓` 對應 Poisson ratio $\nu$
- `h` 對應板厚 $h$
- `D...` 對應彎曲剛度分量
- `D` 或 `D刷` 對應剪切剛度 $D_s$

## 3. 弱式總覽

若只保留位移與轉角主變數，Mindlin 厚板弱式可寫成

$$
a_b(\phi,\delta\phi) + a_s(w,\phi;\delta w,\delta\phi) = l(\delta w,\delta\phi),
$$

其中

$$
a_b = \int_\Omega M_{\alpha\beta}\,\delta\kappa_{\alpha\beta}\,d\Omega,
\qquad
a_s = \int_\Omega Q_\alpha\,\delta\gamma_\alpha\,d\Omega.
$$

用 constitutive relations

$$
M_{\alpha\beta} = D_{\alpha\beta\gamma\eta}\kappa_{\gamma\eta},
\qquad
Q_\alpha = D_s \gamma_\alpha,
$$

可得

$$
a_b = \int_\Omega D_{\alpha\beta\gamma\eta}\kappa_{\alpha\beta}\,\delta\kappa_{\gamma\eta}\,d\Omega,
$$

$$
a_s = \int_\Omega D_s (w_{,\alpha}-\phi_\alpha)(\delta w_{,\alpha}-\delta\phi_\alpha)\,d\Omega.
$$

域內外力與邊界項則對應

$$
l_\Omega = \int_\Omega q\,\delta w\,d\Omega + \int_\Omega m_\alpha\,\delta\phi_\alpha\,d\Omega,
$$

$$
l_\Gamma = \int_{\Gamma_V} \bar V\,\delta w\,d\Gamma + \int_{\Gamma_M}\bar M_\alpha\,\delta\phi_\alpha\,d\Gamma.
$$

若用 penalty 施加 essential boundary conditions，則額外加入

$$
\int_{\Gamma_w}\alpha (w-g)\,\delta w\,d\Gamma,
\qquad
\int_{\Gamma_\phi}\alpha (\phi-\bar\phi)\cdot\delta\phi\,d\Gamma.
$$

## 4. 離散記號

令

$$
w^h = \sum_J N_J w_J,
\qquad
\phi_1^h = \sum_J N_J \phi_{1J},
\qquad
\phi_2^h = \sum_J N_J \phi_{2J},
$$

且記

$$
B_x = N_{,x}, \qquad B_y = N_{,y}.
$$

程式中對應為

- `N = 尉[:𝝭]` 對應形狀函數 $N_J$
- `B_x = 尉[:∂𝝭∂x]` 對應 $N_{J,x}$
- `B_y = 尉[:∂𝝭∂y]` 對應 $N_{J,y}$
- `𝑤 = 尉.𝑤` 對應積分權重與 Jacobian，以下統一記成 $d\Omega$

## 5. 主要 operator 對照

下列 operator 是 `thick_plate.jl` 中最核心的主系統項。

### 5.1 `∫κMγQdΩ`

對應公式是把彎曲與剪切同時展開後，直接一次組成 $[w,\phi_1,\phi_2]$ 的總剛度矩陣：

$$
a = \int_\Omega D_s (w_{,1}-\phi_1)(\delta w_{,1}-\delta\phi_1)\,d\Omega
+ \int_\Omega D_s (w_{,2}-\phi_2)(\delta w_{,2}-\delta\phi_2)\,d\Omega
+ \int_\Omega D_{\alpha\beta\gamma\eta}\kappa_{\alpha\beta}\,\delta\kappa_{\gamma\eta}\,d\Omega.
$$

離散後可寫成

$$
\mathbf K =
\begin{bmatrix}
\mathbf K_{ww}^s & \mathbf K_{w\phi_1}^s & \mathbf K_{w\phi_2}^s \\
\mathbf K_{\phi_1 w}^s & \mathbf K_{\phi_1\phi_1}^{b+s} & \mathbf K_{\phi_1\phi_2}^{b} \\
\mathbf K_{\phi_2 w}^s & \mathbf K_{\phi_2\phi_1}^{b} & \mathbf K_{\phi_2\phi_2}^{b+s}
\end{bmatrix}.
$$

其 block 內容正是程式中的九個子項：

$$
K_{ww}^s = \int_\Omega D_s(N_{I,x}N_{J,x}+N_{I,y}N_{J,y})\,d\Omega,
$$

$$
K_{w\phi_1}^s = -\int_\Omega D_s N_{J}N_{I,x}\,d\Omega,
\qquad
K_{w\phi_2}^s = -\int_\Omega D_s N_{J}N_{I,y}\,d\Omega,
$$

$$
K_{\phi_1 w}^s = -\int_\Omega D_s N_{I}N_{J,x}\,d\Omega,
\qquad
K_{\phi_2 w}^s = -\int_\Omega D_s N_{I}N_{J,y}\,d\Omega,
$$

$$
K_{\phi_1\phi_1}^{b+s} =
\int_\Omega \left(D_{1111}N_{I,x}N_{J,x}+D_{1212}N_{I,y}N_{J,y}+D_sN_IN_J\right)d\Omega,
$$

$$
K_{\phi_1\phi_2}^{b} =
\int_\Omega \left(D_{1122}N_{I,x}N_{J,y}+D_{1212}N_{I,y}N_{J,x}\right)d\Omega,
$$

$$
K_{\phi_2\phi_2}^{b+s} =
\int_\Omega \left(D_{1111}N_{I,y}N_{J,y}+D_{1212}N_{I,x}N_{J,x}+D_sN_IN_J\right)d\Omega.
$$

矩陣位置：

- row/col `3I-2, 3J-2` 是 $w$-$w$
- row/col `3I-1, 3J-1` 是 $\phi_1$-$\phi_1$
- row/col `3I, 3J` 是 $\phi_2$-$\phi_2$
- 其餘四個 off-diagonal block 是 $w$-$\phi$ 與 $\phi_1$-$\phi_2$ 耦合

### 5.2 `∫κκdΩ`

對應純彎曲能量

$$
a_b = \int_\Omega D_{\alpha\beta\gamma\eta}\kappa_{\alpha\beta}\,\delta\kappa_{\gamma\eta}\,d\Omega.
$$

離散後形成 rotation-rotation block：

$$
\mathbf K_{\phi\phi}^b =
\int_\Omega
\begin{bmatrix}
D_{1111}B_x^{\mathsf T}B_x + D_{1212}B_y^{\mathsf T}B_y &
D_{1122}B_x^{\mathsf T}B_y + D_{1212}B_y^{\mathsf T}B_x \\
D_{1122}B_y^{\mathsf T}B_x + D_{1212}B_x^{\mathsf T}B_y &
D_{1111}B_y^{\mathsf T}B_y + D_{1212}B_x^{\mathsf T}B_x
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- 僅佔據 $\phi$ block，即 `2I-1,2J-1` 到 `2I,2J`

### 5.3 `∫wwdΩ`

對應剪切能中的 $w$-$w$ 部分：

$$
a_{ww}^s = \int_\Omega D_s \nabla w \cdot \nabla \delta w\,d\Omega.
$$

離散式為

$$
K_{ww} = \int_\Omega D_s(N_{I,x}N_{J,x}+N_{I,y}N_{J,y})\,d\Omega.
$$

矩陣位置：

- scalar $w$ block，索引為 `k[I,J]`

### 5.4 `∫φφdΩ`

對應剪切能中的 rotation-rotation 部分：

$$
a_{\phi\phi}^s = \int_\Omega D_s\,\phi_\alpha\,\delta\phi_\alpha\,d\Omega.
$$

離散後

$$
\mathbf K_{\phi\phi}^s
= \int_\Omega D_s
\begin{bmatrix}
N_I N_J & 0 \\
0 & N_I N_J
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- `2I-1,2J-1` 對應 $\phi_1$-$\phi_1$
- `2I,2J` 對應 $\phi_2$-$\phi_2$

### 5.5 `∫φwdΩ`

對應剪切能中的 rotation-displacement 耦合項

$$
a_{\phi w}^s = -\int_\Omega D_s\,\phi_\alpha\,\delta w_{,\alpha}\,d\Omega
\quad\text{或其轉置形式}\quad
-\int_\Omega D_s\,\delta\phi_\alpha\,w_{,\alpha}\,d\Omega.
$$

程式使用的 block 寫成

$$
\mathbf K_{\phi w}
= -\int_\Omega D_s
\begin{bmatrix}
N_I N_{J,x} \\
N_I N_{J,y}
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- `2I-1,J` 對應 $\phi_1$-$w$
- `2I,J` 對應 $\phi_2$-$w$

### 5.6 `∫wqdΩ`

對應板面分佈載重：

$$
l_w = \int_\Omega q\,\delta w\,d\Omega.
$$

離散後

$$
f_I^{(w)} = \int_\Omega N_I q\,d\Omega.
$$

矩陣位置：

- 加到 $w$ 的右手邊向量 `f[I]`

### 5.7 `∫φmdΩ`

對應分佈力矩：

$$
l_\phi = \int_\Omega m_\alpha\,\delta\phi_\alpha\,d\Omega.
$$

離散後

$$
\mathbf f_I^{(\phi)} =
\int_\Omega
\begin{bmatrix}
N_I m_1 \\
N_I m_2
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- `f[2I-1]` 對應 $\phi_1$
- `f[2I]` 對應 $\phi_2$

### 5.8 `∫wVdΓ`

對應自然邊界剪力：

$$
l_{\Gamma_V} = \int_{\Gamma_V} \bar V\,\delta w\,d\Gamma.
$$

離散後

$$
f_I^{(V)} = \int_{\Gamma_V} N_I \bar V\,d\Gamma.
$$

矩陣位置：

- 邊界上的 $w$ 右手邊 `f[I]`

### 5.9 `∫φMdΓ`

對應自然邊界彎矩：

$$
l_{\Gamma_M} = \int_{\Gamma_M}\bar M_\alpha\,\delta\phi_\alpha\,d\Gamma.
$$

離散後

$$
\mathbf f_I^{(M)} =
\int_{\Gamma_M}
\begin{bmatrix}
-N_I \bar M_1 \\
-N_I \bar M_2
\end{bmatrix}
d\Gamma.
$$

矩陣位置：

- `f[2I-1]`, `f[2I]`

### 5.10 `∫αwwdΓ`

對應位移的 penalty 邊界條件：

$$
a_{\Gamma_w} = \int_{\Gamma_w}\alpha\, w\,\delta w\,d\Gamma,
\qquad
l_{\Gamma_w} = \int_{\Gamma_w}\alpha\, g\,\delta w\,d\Gamma.
$$

離散後

$$
K_{IJ}^{(\alpha w)} = \int_{\Gamma_w}\alpha N_I N_J\,d\Gamma,
\qquad
f_I^{(\alpha w)} = \int_{\Gamma_w}\alpha N_I g\,d\Gamma.
$$

矩陣位置：

- $w$-$w$ block 與 $w$ 的右手邊

### 5.11 `∫αφφdΓ`

對應轉角的 penalty 邊界條件：

$$
a_{\Gamma_\phi} = \int_{\Gamma_\phi}\alpha\, \delta\phi^{\mathsf T}\mathbf n_\phi \phi\,d\Gamma,
\qquad
l_{\Gamma_\phi} = \int_{\Gamma_\phi}\alpha\, \delta\phi^{\mathsf T}\mathbf n_\phi \bar\phi\,d\Gamma.
$$

這裡 $\mathbf n_\phi$ 在程式中由邊界局部法向分量組成，對應 `n₁n₁`, `n₁n₂`, `n₂n₂` 一類係數。

離散後

$$
\mathbf K_{\Gamma_\phi} =
\int_{\Gamma_\phi}
\alpha N_I
\begin{bmatrix}
n_{11} & n_{12} \\
n_{12} & n_{22}
\end{bmatrix}
N_J\,d\Gamma.
$$

矩陣位置：

- 只作用在 $\phi$-$\phi$ block 與其右手邊

## 6. mixed / auxiliary operator 對照

這一組 operator 對應的是把剪力 $\mathbf Q$ 與彎矩 $\mathbf M$ 視為獨立未知量時的混合形式。它們不是主位移法必需的項，但在 mixed formulation、邊界條件弱施加或後處理中會出現。

### 6.1 `∫Q∇wdΩ`

對應獨立剪力與位移梯度之間的耦合：

$$
\int_\Omega \delta Q_\alpha\, w_{,\alpha}\,d\Omega.
$$

離散後形成矩形 block

$$
\mathbf K_{Qw}
= \int_\Omega
\begin{bmatrix}
N_I N_{J,x} \\
N_I N_{J,y}
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- rows 為 $\mathbf Q$ 的兩個分量
- cols 為 $w$

### 6.2 `∫∇QwdΩ`

這是與 `∫Q∇wdΩ` 成對的共軛項，符號上可理解為

$$
-\int_\Omega Q_\alpha\,\delta w_{,\alpha}\,d\Omega.
$$

離散後是帶負號的 coupling block：

$$
\mathbf K_{\nabla Q,w}
= -\int_\Omega
\begin{bmatrix}
N_{I,x}N_J \\
N_{I,y}N_J
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- 仍是向量場對純量場的矩形 block
- 單空間與雙空間版本都存在

### 6.3 `∫QφdΩ`

對應獨立剪力與轉角的耦合：

$$
-\int_\Omega \delta Q_\alpha\,\phi_\alpha\,d\Omega.
$$

離散後

$$
\mathbf K_{Q\phi}
= -\int_\Omega
\begin{bmatrix}
N_I N_J & 0 \\
0 & N_I N_J
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- rows 為 $\mathbf Q$
- cols 為 $\phi$

### 6.4 `∫QQdΩ`

對應剪力場自己的 compliance 項：

$$
-\int_\Omega \frac{1}{D_s}\,\delta Q_\alpha Q_\alpha\,d\Omega.
$$

離散後

$$
\mathbf K_{QQ}
= -\int_\Omega \frac{1}{D_s}
\begin{bmatrix}
N_I N_J & 0 \\
0 & N_I N_J
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- 只作用在 $\mathbf Q$-$\mathbf Q$ block

### 6.5 `∫QwdΓ`

對應剪力邊界通量與邊界位移的耦合，可理解成

$$
\int_{\Gamma_Q} \delta Q_\alpha n_\alpha\, w\,d\Gamma
\quad\text{或}\quad
-\int_{\Gamma_Q} Q_\alpha n_\alpha\, \delta w\,d\Gamma.
$$

因此程式提供三種版本：

- 只組矩陣
- 組矩陣加右手邊
- 兩個空間的 cross-space 版本

矩陣位置：

- rows 為 $\mathbf Q$
- cols 為 $w$

### 6.6 `∫MMdΩ`

對應彎矩獨立場的 compliance 項：

$$
-\int_\Omega \delta\mathbf M^{\mathsf T}\mathbf C_b\,\mathbf M\,d\Omega,
$$

其中

$$
\mathbf C_b =
\begin{bmatrix}
\dfrac{12}{Eh^3} & -\dfrac{12\nu}{Eh^3} & 0 \\
-\dfrac{12\nu}{Eh^3} & \dfrac{12}{Eh^3} & 0 \\
0 & 0 & \dfrac{24(1+\nu)}{Eh^3}
\end{bmatrix}.
$$

矩陣位置：

- `3I-2,3J-2` 對應 $M_{11}$-$M_{11}$
- `3I-1,3J-1` 對應 $M_{22}$-$M_{22}$
- `3I,3J` 對應 $M_{12}$-$M_{12}$
- 另外存在 $M_{11}$-$M_{22}$ 的耦合項

### 6.7 `∫MφdΓ`

對應彎矩與邊界轉角之間的弱式邊界項。它可視為

$$
-\int_{\Gamma_M}\delta \mathbf M^{\mathsf T}\mathbf T_n \phi\,d\Gamma
$$

以及帶 prescribed rotation 的一致右手邊版本

$$
\int_{\Gamma_M}\delta \mathbf M^{\mathsf T}\mathbf T_n \bar\phi\,d\Gamma.
$$

這裡 $\mathbf T_n$ 是由法向分量組成的邊界轉換矩陣。

矩陣位置：

- rows 為三分量彎矩場 $\mathbf M$
- cols 為兩分量轉角場 $\phi$

### 6.8 `∫∇MφdΩ`

對應彎矩場與轉角梯度的域內耦合：

$$
\int_\Omega \delta\mathbf M^{\mathsf T}\mathbf B_\phi \phi\,d\Omega.
$$

以程式的 index 排列來看，其 block 內容可理解為

$$
\mathbf K_{\nabla M,\phi}
= \int_\Omega
\begin{bmatrix}
B_x^{\mathsf T}N & 0 \\
0 & B_y^{\mathsf T}N \\
B_y^{\mathsf T}N & B_x^{\mathsf T}N
\end{bmatrix}
d\Omega.
$$

矩陣位置：

- rows 為 $M_{11},M_{22},M_{12}$
- cols 為 $\phi_1,\phi_2$

## 7. 誤差範數

### 7.1 `L₂Q`

對應剪力場的相對 $L_2$ 誤差：

$$
\|e_Q\|_{L_2,\text{rel}}
=
\left(
\frac{\int_\Omega (\mathbf Q_h-\mathbf Q)^T D_s^{-1}(\mathbf Q_h-\mathbf Q)\,d\Omega}
{\int_\Omega \mathbf Q^T D_s^{-1}\mathbf Q\,d\Omega}
\right)^{1/2}.
$$

### 7.2 `L₂φ`

對應轉角場的相對 $L_2$ 誤差：

$$
\|e_\phi\|_{L_2,\text{rel}}
=
\left(
\frac{\int_\Omega [(\phi_{1h}-\phi_1)^2 + (\phi_{2h}-\phi_2)^2]\,d\Omega}
{\int_\Omega (\phi_1^2+\phi_2^2)\,d\Omega}
\right)^{1/2}.
$$

### 7.3 `L₂`

對應位移場 $w$ 的相對 $L_2$ 誤差：

$$
\|e_w\|_{L_2,\text{rel}}
=
\left(
\frac{\int_\Omega (w_h-w)^2\,d\Omega}
{\int_\Omega w^2\,d\Omega}
\right)^{1/2}.
$$

## 8. 快速查表

| Operator | 物理意義 | 主要 block |
| --- | --- | --- |
| `∫κMγQdΩ` | 彎曲加剪切的 combined stiffness | `w-φ` 全系統 |
| `∫κκdΩ` | 純彎曲剛度 | `φ-φ` |
| `∫wwdΩ` | 剪切能中的 `w-w` 項 | `w-w` |
| `∫φφdΩ` | 剪切能中的 `φ-φ` 項 | `φ-φ` |
| `∫φwdΩ` | 剪切耦合項 | `φ-w` |
| `∫wqdΩ` | 分佈載重 | `f_w` |
| `∫φmdΩ` | 分佈力矩 | `f_φ` |
| `∫wVdΓ` | 自然剪力邊界 | `f_w` |
| `∫φMdΓ` | 自然彎矩邊界 | `f_φ` |
| `∫αwwdΓ` | 位移 penalty 邊界 | `w-w`, `f_w` |
| `∫αφφdΓ` | 轉角 penalty 邊界 | `φ-φ`, `f_φ` |
| `∫Q∇wdΩ` | 獨立剪力與位移梯度耦合 | `Q-w` |
| `∫∇QwdΩ` | `Q-w` 共軛耦合項 | `Q-w` |
| `∫QφdΩ` | 獨立剪力與轉角耦合 | `Q-φ` |
| `∫QQdΩ` | 剪力 compliance | `Q-Q` |
| `∫QwdΓ` | 剪力通量邊界耦合 | `Q-w`, `f_Q` |
| `∫MMdΩ` | 彎矩 compliance | `M-M` |
| `∫MφdΓ` | 彎矩與邊界轉角耦合 | `M-φ`, `f_M` |
| `∫∇MφdΩ` | 彎矩與轉角梯度耦合 | `M-φ` |
| `L₂Q` | 剪力相對誤差 | 後處理 |
| `L₂φ` | 轉角相對誤差 | 後處理 |
| `L₂` | 位移相對誤差 | 後處理 |

## 9. 總結

`thick_plate.jl` 可以分成三層來讀：

1. 主位移法的 Mindlin 厚板 operator：`∫κMγQdΩ`, `∫κκdΩ`, `∫wwdΩ`, `∫φφdΩ`, `∫φwdΩ`, `∫wqdΩ`, `∫φmdΩ`, `∫wVdΓ`, `∫φMdΓ`, `∫αwwdΓ`, `∫αφφdΓ`
2. mixed formulation 的獨立剪力 / 彎矩 operator：`∫Q∇wdΩ`, `∫∇QwdΩ`, `∫QφdΩ`, `∫QQdΩ`, `∫QwdΓ`, `∫MMdΩ`, `∫MφdΓ`, `∫∇MφdΩ`
3. 後處理與驗證：`L₂Q`, `L₂φ`, `L₂`

因此若你的目的只是組標準厚板主系統，第一層就足夠；若你要理解 `thick_plate.jl` 為什麼同時有 $Q$ 與 $M$ 這類 operator，就需要把它視為保留 mixed 弱式結構的一份操作庫。
