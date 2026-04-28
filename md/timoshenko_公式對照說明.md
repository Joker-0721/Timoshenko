# `timoshenko.jl` 公式對照說明

本文直接跟隨 [`test.md`](</d:/Joker/Timoshenko/test.md>) 的章節順序，把 2D Mindlin 厚板公式改寫成 1D Timoshenko 梁版本，並對照 [`ApproxOperator.jl/src/operation/timoshenko.jl`](</d:/ApproxOperator.jl/src/operation/timoshenko.jl>) 裡的 operator。

如果把 `test.md` 的主線記成

$$
\text{Kinematics} \rightarrow \text{Constitutive relations} \rightarrow \text{Stress resultants}
\rightarrow \text{Strong form} \rightarrow \text{Weak form} \rightarrow \text{Finite element discretization},
$$

那麼這份文件做的事，就是把板的記號

$$
(w,\phi_1,\phi_2,\gamma_\alpha,\kappa_{\alpha\beta},Q_\alpha,M_{\alpha\beta})
$$

縮減為梁的記號

$$
(w,\phi,\gamma,\kappa,Q,M).
$$

## 1 Kinematics

Timoshenko 梁的位移場可視為厚板位移場在 1D 情況下的簡化：

$$
\bar u_1(x,z) = -z\,\phi(x),\qquad \bar u_2 = 0,\qquad \bar u_3(x,z) = w(x).
$$

對應的應變可以拆成彎曲與剪切兩部分：

$$
\epsilon_{11} = -z\,\kappa,
\qquad
\gamma = 2\epsilon_{13} = w_{,x} - \phi,
$$

其中

$$
\kappa = \phi_{,x}.
$$

和 `test.md` 對照時，可以把

$$
\gamma_\alpha = w_{,\alpha}-\phi_\alpha
$$

看成退化成單一方向的

$$
\gamma = w_{,x}-\phi,
$$

而板的曲率張量 $\kappa_{\alpha\beta}$ 則退化成梁的單一曲率 $\kappa$。

## 2 Constitutive relations

梁的彎矩與剪力本構為

$$
M = EI\,\kappa,
\qquad
Q = kGA\,\gamma.
$$

其中

$$
EI = E I,
\qquad
kGA = \kappa_s G A,
\qquad
G = \frac{E}{2(1+\nu)}.
$$

程式中的對應參數為

- `EI` 對應 $EI$
- `kGA` 對應 $kGA$

## 3 Stress resultants

和 `test.md` 的板理論一樣，梁也可以先透過截面積分定義應力合力與應力矩：

$$
M = \int_A z\,\sigma_{11}\,dA = EI\,\kappa,
$$

$$
Q = k\int_A \tau_{13}\,dA = kGA\,\gamma.
$$

在 1D 梁裡，不再另外保留面內膜力 $N_{\alpha\beta}$，只保留彎矩 $M$ 與剪力 $Q$ 兩個截面力。

## 4 Equilibrium equations and boundary conditions (strong form)

若只考慮靜態彎曲，總位能可寫成

$$
\delta \Pi
= \int_0^L \left(M\,\delta\kappa + Q\,\delta\gamma\right)\,dx
- \int_0^L q\,\delta w\,dx
- \int_0^L m\,\delta\phi\,dx
= 0.
$$

若進一步考慮軸力 $P$ 造成的幾何剛度，則可加上

$$
- \int_0^L P\,w_{,x}\,\delta w_{,x}\,dx
$$

形成 Engesser 幾何剛度形式。

對應的強式平衡方程可寫成

$$
M_{,x} + Q = 0 \quad \text{in } (0,L),
$$

$$
Q_{,x} + q = 0 \quad \text{in } (0,L),
$$

若考慮幾何剛度，第二式會額外帶入 $P$ 對位移梯度的貢獻。

自然邊界條件對應為

$$
M = \bar M \quad \text{on } \Gamma_M,
\qquad
Q = \bar V \quad \text{on } \Gamma_V.
$$

## 5 Galerkin weak form

對強式分別乘上測試函數 $\delta\phi$ 與 $\delta w$ 並積分，可得

$$
\int_0^L (M_{,x}+Q)\,\delta\phi\,dx = 0,
\qquad
\int_0^L (Q_{,x}+q)\,\delta w\,dx = 0.
$$

分部積分後，可寫成

$$
\int_0^L M\,\delta\phi_{,x}\,dx
+ \int_0^L Q\,\delta\phi\,dx
= \int_{\Gamma_M}\bar M\,\delta\phi\,d\Gamma,
$$

$$
\int_0^L Q\,\delta w_{,x}\,dx
= \int_0^L q\,\delta w\,dx
+ \int_{\Gamma_V}\bar V\,\delta w\,d\Gamma.
$$

合併後得到

$$
\boxed{
\int_0^L \left(M\,\delta\kappa + Q\,\delta\gamma\right)\,dx
= \int_0^L q\,\delta w\,dx
+ \int_0^L m\,\delta\phi\,dx
+ \int_{\Gamma_V}\bar V\,\delta w\,d\Gamma
- \int_{\Gamma_M}\bar M\,\delta\phi\,d\Gamma
}.
$$

上面這一串

$$
\int_0^L q\,\delta w\,dx
+ \int_0^L m\,\delta\phi\,dx
+ \int_{\Gamma_V}\bar V\,\delta w\,d\Gamma
- \int_{\Gamma_M}\bar M\,\delta\phi\,d\Gamma
$$

不是剛度矩陣 $\mathbf K$ 的一部分，而是一般外力泛函 $l(\delta w,\delta\phi)$。  
它之所以會出現，是因為這裡寫的是最一般的 weak form：允許同時有

- 分佈橫向載重 $q$
- 分佈力矩 $m$
- 剪力自然邊界 $\bar V$
- 彎矩自然邊界 $\bar M$

如果你的目標只是推導純剛度矩陣，或是做沒有外力的特徵值問題，這一整串就會被移到右手邊，或直接設成零。

因此梁版本的 weak form 可以用和 `test.md` 相同的記號寫成

$$
a_b + a_s - a_g = l,
$$

其中

$$
a_b = \int_0^L M\,\delta\kappa\,dx,
\qquad
a_s = \int_0^L Q\,\delta\gamma\,dx.
$$

對 Engesser 幾何剛度：

$$
a_g^{(E)} = P\int_0^L w_{,x}\,\delta w_{,x}\,dx.
$$

## 6 Finite element discretization

令

$$
w^h(x) = \sum_K N_K(x)w_K,
\qquad
\phi^h(x) = \sum_R N_R(x)\phi_R.
$$

則

$$
\kappa_R = N_{R,x}\phi_R,
\qquad
\gamma_K = N_{K,x}w_K - N_K \phi_K,
\qquad
w_{,xK} = N_{K,x} w_K.
$$

程式中的離散記號對應為

- `N = ξ[:𝝭]` 對應 $N_I$
- `B = ξ[:∂𝝭∂x]` 對應 $N_{I,x}$
- `𝑤 = ξ.𝑤` 對應 $dx$

---

彎曲項：

$$
a_b = \int_0^L EI\,\kappa^h\,\delta\kappa^h\,dx
$$

$$
a_b = \sum_{R,S}\delta\phi_R
\left(
\int_0^L EI\,N_{R,x}N_{S,x}\,dx
\right)\phi_S
$$

$$
K_{\phi\phi,RS}^{b} = \int_0^L EI\,N_{R,x}N_{S,x}\,dx
$$

$$
a^b = \sum_{R,S}\delta\phi_R\,K_{\phi\phi,RS}^{b}\,\phi_S
$$

---

剪切項：

$$
a_s = \int_0^L kGA\,\gamma^h\,\delta\gamma^h\,dx
$$

$$
\begin{aligned}
a_s = \sum_{K,L}\Bigg[
&\delta w_K \left(\int_0^L kGA\,N_{K,x}N_{L,x}\,dx\right) w_L \\
+&\delta w_K \left(\int_0^L -kGA\,N_{K,x}N_L\,dx\right)\phi_L \\
+&\delta\phi_K \left(\int_0^L -kGA\,N_K N_{L,x}\,dx\right) w_L \\
+&\delta\phi_K \left(\int_0^L kGA\,N_K N_L\,dx\right)\phi_L
\Bigg]
\end{aligned}
$$

$$
\mathbf K_{KL}^{s}
=
\int_0^L
\begin{bmatrix}
kGA\,N_{K,x}N_{L,x} & -kGA\,N_{K,x}N_L \\
-kGA\,N_K N_{L,x} & kGA\,N_K N_L
\end{bmatrix} dx
$$

$$
a_s = \sum_{K,L}
\begin{bmatrix}
\delta w_K & \delta\phi_K
\end{bmatrix}
\mathbf K_{KL}^{s}
\begin{bmatrix}
w_L \\
\phi_L
\end{bmatrix}
$$

---

Engesser 幾何剛度：

$$
a_g^{(E)} = P\int_0^L w_{,x}\,\delta w_{,x}\,dx
$$

$$
a_g^{(E)} = \sum_{K,L}\delta w_K
\left(
P\int_0^L N_{K,x}N_{L,x}\,dx
\right)w_L
$$

$$
K_{ww,KL}^{g,E} = \int_0^L N_{K,x}N_{L,x}\,dx
$$

因此總矩陣可寫為

$$
\left(
\begin{bmatrix}
\mathbf K_{ww}^{s} & \mathbf K_{w\phi}^{s} \\
\mathbf K_{\phi w}^{s} & \mathbf K^{b}+\mathbf K_{\phi\phi}^{s}
\end{bmatrix}
-
P
\begin{bmatrix}
\mathbf K_{ww}^{g,E} & \mathbf 0 \\
\mathbf 0 & \mathbf 0
\end{bmatrix}
\right)
\begin{bmatrix}
\mathbf d_w \\
\mathbf d_\phi
\end{bmatrix}
= \mathbf 0
$$

## 7 `timoshenko.jl` operator 對照

### 7.1 Material operators

#### `∫κEIκds`

$$
\int_0^L EI\,\phi_{,x}\,\delta\phi_{,x}\,dx
$$

combined ordering 下只作用在 $\phi$-$\phi$ block，也就是 `k[2I,2J]`。

#### `∫γkGAγds`

$$
\int_0^L kGA(w_{,x}-\phi)(\delta w_{,x}-\delta\phi)\,dx
$$

對應 combined form 的四個 block：

- `2I-1,2J-1` 對應 $K_{ww}^s$
- `2I-1,2J` 對應 $K_{w\phi}^s$
- `2I,2J-1` 對應 $K_{\phi w}^s$
- `2I,2J` 對應 $K_{\phi\phi}^s$

#### `∫vqds`

$$
\int_0^L q\,\delta w\,dx
$$

只加到 combined ordering 下的 $w$ 自由度 `f[2I-1]`。

#### `∫κκdΩ`

$$
K_{\phi\phi}^{b} = \int_0^L EI\,N_{I,x}N_{J,x}\,dx
$$

#### `∫φφdΩ`

$$
K_{\phi\phi}^{s} = \int_0^L kGA\,N_I N_J\,dx
$$

#### `∫wwdΩ`

$$
K_{ww}^{s} = \int_0^L kGA\,N_{I,x}N_{J,x}\,dx
$$

#### `∫φwdΩ`

$$
K_{\phi w}^{s} = -\int_0^L kGA\,N_I N_{J,x}\,dx
$$

實務上 `$K_{w\phi}^{s} = (K_{\phi w}^{s})^{\mathsf T}$`，所以現有 block 例子多半以 `kᵠʷ'` 取得上三角耦合塊。

### 7.2 Load and boundary operators

#### `∫wqdΩ`

$$
f_I^{(w)} = \int_0^L N_I q\,dx
$$

#### `∫φmdΩ`

$$
f_I^{(\phi)} = \int_0^L N_I m\,dx
$$

#### `∫αwwdΓ`

$$
K_{IJ}^{(\alpha w)} = \int_{\Gamma_w}\alpha N_I N_J\,d\Gamma,
\qquad
f_I^{(\alpha w)} = \int_{\Gamma_w}\alpha N_I g\,d\Gamma
$$

#### `∫αφφdΓ`

$$
K_{IJ}^{(\alpha\phi)} = \int_{\Gamma_\phi}\alpha N_I N_J\,d\Gamma,
\qquad
f_I^{(\alpha\phi)} = \int_{\Gamma_\phi}\alpha N_I \bar\phi\,d\Gamma
$$

#### `∫wVdΓ`

$$
f_I^{(V)} = \int_{\Gamma_V}N_I \bar V\,d\Gamma
$$

#### `∫φMdΓ`

$$
f_I^{(M)} = -\int_{\Gamma_M}N_I \bar M\,d\Gamma
$$

### 7.3 Geometric stiffness operators

#### `∫wwGdΩ`

$$
K_{ww}^{g,E} = \int_0^L N_{I,x}N_{J,x}\,dx
$$

對應 Engesser 形式中的 $\mathbf K_{ww}^{g}$。

## 8 方程式組裝觀點

若只做靜態彎曲，block system 可寫成

$$
\begin{bmatrix}
\mathbf K_{\phi\phi}^{b}+\mathbf K_{\phi\phi}^{s} & \mathbf K_{\phi w}^{s} \\
(\mathbf K_{\phi w}^{s})^{\mathsf T} & \mathbf K_{ww}^{s}
\end{bmatrix}
\begin{bmatrix}
\mathbf d_\phi \\
\mathbf d_w
\end{bmatrix}
=
\begin{bmatrix}
\mathbf f_\phi \\
\mathbf f_w
\end{bmatrix}.
$$

這正是 [`src/patch_test.jl`](</d:/Joker/Timoshenko/src/patch_test.jl>) 和 [`src/square.jl`](</d:/Joker/Timoshenko/src/square.jl>) 內使用的

$$
[k_{\phi\phi}\; k_{\phi w};\; k_{\phi w}^{\mathsf T}\; k_{ww}]
$$

組裝方式。

若做 Engesser 屈曲，則可寫成

$$
\left(
\begin{bmatrix}
\mathbf K_{\phi\phi}^{b}+\mathbf K_{\phi\phi}^{s} & \mathbf K_{\phi w}^{s} \\
(\mathbf K_{\phi w}^{s})^{\mathsf T} & \mathbf K_{ww}^{s}
\end{bmatrix}
-
P
\begin{bmatrix}
\mathbf 0 & \mathbf 0 \\
\mathbf 0 & \mathbf K_{ww}^{g,E}
\end{bmatrix}
\right)\mathbf d = \mathbf 0.
$$

這和 [`src/beam_Engesser.jl`](</d:/Joker/Timoshenko/src/beam_Engesser.jl>) 的寫法一致。

## 9 誤差範數與解析解

### 9.1 `L₂`

$$
\|e_w\|_{L_2,\mathrm{rel}}
=
\left(
\frac{\int_0^L (w_h-w)^2\,dx}
{\int_0^L w^2\,dx}
\right)^{1/2}
$$

### 9.2 `L₂φ`

$$
\|e_\phi\|_{L_2,\mathrm{rel}}
=
\left(
\frac{\int_0^L (\phi_h-\phi)^2\,dx}
{\int_0^L \phi^2\,dx}
\right)^{1/2}
$$

### 9.3 `w_exact_ss`, `w_exact_cf`

簡支梁：

$$
w_{\mathrm{exact,ss}}(x)
= \frac{q\,x(L^3-2Lx^2+x^3)}{24EI}
+ \frac{q\,x(L-x)}{2kGA}
$$

懸臂梁：

$$
w_{\mathrm{exact,cf}}(x)
= \frac{q\,x^2(6L^2-4Lx+x^2)}{24EI}
+ \frac{q\,x(2L-x)}{2kGA}
$$

### 9.4 `φ_exact_ss`, `φ_exact_cf`

簡支梁：

$$
\phi_{\mathrm{exact,ss}}(x)
= \frac{q\,(L^3-6Lx^2+4x^3)}{24EI}
$$

懸臂梁：

$$
\phi_{\mathrm{exact,cf}}(x)
= \frac{q\,x(3L^2-3Lx+x^2)}{6EI}
$$

## 10 快速查表

| Operator | 公式 | 位置 |
| --- | --- | --- |
| `∫κEIκds` | $\int EI\,\phi_{,x}\,\delta\phi_{,x}\,dx$ | combined $\phi$-$\phi$ |
| `∫γkGAγds` | $\int kGA(w_{,x}-\phi)(\delta w_{,x}-\delta\phi)\,dx$ | combined 全矩陣 |
| `∫vqds` | $\int q\,\delta w\,dx$ | combined `f_w` |
| `∫κκdΩ` | $\int EI\,N_{I,x}N_{J,x}\,dx$ | block $\phi$-$\phi$ |
| `∫φφdΩ` | $\int kGA\,N_I N_J\,dx$ | block $\phi$-$\phi$ |
| `∫wwdΩ` | $\int kGA\,N_{I,x}N_{J,x}\,dx$ | block $w$-$w$ |
| `∫φwdΩ` | $-\int kGA\,N_I N_{J,x}\,dx$ | block $\phi$-$w$ |
| `∫wqdΩ` | $\int N_I q\,dx$ | $f_w$ |
| `∫φmdΩ` | $\int N_I m\,dx$ | $f_\phi$ |
| `∫αwwdΓ` | $\int \alpha N_I N_J\,d\Gamma$, $\int \alpha N_I g\,d\Gamma$ | 邊界 $w$ |
| `∫αφφdΓ` | $\int \alpha N_I N_J\,d\Gamma$, $\int \alpha N_I \bar\phi\,d\Gamma$ | 邊界 $\phi$ |
| `∫wVdΓ` | $\int N_I \bar V\,d\Gamma$ | 自然剪力邊界 |
| `∫φMdΓ` | $-\int N_I \bar M\,d\Gamma$ | 自然彎矩邊界 |
| `∫wwGdΩ` | $\int N_{I,x}N_{J,x}\,dx$ | Engesser 幾何剛度 |
| `L₂` | 位移相對 $L_2$ 誤差 | 後處理 |
| `L₂φ` | 轉角相對 $L_2$ 誤差 | 後處理 |
| `w_exact_ss`, `w_exact_cf` | 解析位移場 | 驗證 |
| `φ_exact_ss`, `φ_exact_cf` | 解析轉角場 | 驗證 |

## 11 總結

若把 `test.md` 當成厚板版本的模板，那麼 `timoshenko.jl` 的對應版本就是：

1. 先用
   $$
   \gamma = w_{,x}-\phi,\qquad \kappa=\phi_{,x}
   $$
   把 2D Mindlin kinematics 縮成 1D 梁。
2. 再用
   $$
   M=EI\kappa,\qquad Q=kGA\gamma
   $$
   定義截面力。
3. 最後把有限元素離散寫成
   $$
   \mathbf K =
   \begin{bmatrix}
   \mathbf K_{\phi\phi}^{b}+\mathbf K_{\phi\phi}^{s} & \mathbf K_{\phi w}^{s} \\
   (\mathbf K_{\phi w}^{s})^{\mathsf T} & \mathbf K_{ww}^{s}
   \end{bmatrix},
   $$
   並在屈曲問題中扣掉 Engesser 幾何剛度矩陣。

也就是說，現在的 `timoshenko.jl` 已經不只是「靜態彎曲 operator 集合」，而是能和 `test.md` 的主線真正對上的一套梁版本 operator。
