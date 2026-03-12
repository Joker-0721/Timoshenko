### 運動方程 (Kinematic Equations)

* 曲率 (Curvature)
  
  $$
  \kappa = \theta_{,s} \quad \text{--- (1)}
  $$
* 剪切 (Shear)
  
  $$
  \gamma = w_{,s} - \theta + \frac{u}{R} \quad \text{--- (2)}
  $$

> 註：\$w\$ 徑向位移 (radial displacement)

* 軸向 (Axial)
  $$
  \varepsilon = u_{,s} - \frac{w}{R} \quad \text{--- (3)}
  $$

> 註：\$u\$ 切向位移 (tangential displacement)

---

### 本構關係 (Constitutive Relations)

$$
M_b = EI\kappa \quad \text{--- (4)}
$$

$$
V = GAk\gamma \quad \text{--- (5)}
$$

> 註：$G$ 剪切模量
> 
> 註：$k$ 剪切修正係數 (Shear correction factor)

$$
N = EA\varepsilon \quad \text{--- (6)}
$$

---

### 平衡方程 (Equilibrium Equations)

$$
M_{b,s} + V = 0
$$

$$
\Rightarrow EI\kappa_{,s} + GAk\gamma = 0 \quad \text{--- (7)}
$$

$$
\gamma = -\frac{EI}{GAk}\kappa_{,s} \quad \text{--- (9)}
$$

$$
V_{,s} + \frac{N}{R} = 0 \quad \text{--- (8)}
$$

$$
\Rightarrow EI\kappa_{,ss} + \frac{EA}{R}\varepsilon = 0 
$$

$$
\varepsilon = \frac{IR}{A}\kappa_{,ss} \quad \text{--- (10)}
$$

#### 整理與回顧 (Recap)

$$
\kappa = \theta_{,s} \quad \text{--- (11)}
$$

$$
\kappa_{,s} = -\frac{GAk}{EI}\gamma \quad \text{--- (12)}
$$

$$
\kappa_{,ss} = \frac{A}{IR}\varepsilon \quad \text{--- (13)}
$$

曲率是轉角的導數：$\kappa = \frac{d\theta}{ds}$。因此，轉角是曲率的積分：

$$\theta = \int \kappa \, ds + C_\theta  \quad \text{--- (14)}$$

### 消除 $u$ 和 $\theta$ (Eliminating $u$ and $\theta$)

**(12) 對 $s$ 求導 (Differentiate Eq. 12 w.r.t s):**

$$
\kappa_{,ss} = -\frac{GAk}{EI}\gamma_{,s}
$$

$$
= -\frac{GAk}{EI} \left( w_{,ss} - \theta_{,s} + \frac{u_{,s}}{R} \right)
$$

$$
= -\frac{GAk}{EI} \left( w_{,ss} - \kappa + \frac{u_{,s}}{R} \right)
$$

(13) 代入 (12)

(註：此處推導是利用 (3) 式 $\varepsilon = u_{,s} - \frac{w}{R}$ 結合 (13) 來替換 $u_{,s}$)

$$
u_{,s} = \frac{IR}{A}\kappa_{,ss} + \frac{w}{R} \quad \text{--- (13)}
$$

$$
\Rightarrow \kappa_{,ss} = -\frac{GAk}{EI} \left( w_{,ss} - \kappa + \frac{1}{R} \left( \frac{IR}{A}\kappa_{,ss} + \frac{w}{R} \right) \right)
$$

$$
\kappa_{,ss} = -\frac{GAk}{EI} \left( w_{,ss} - \kappa + \frac{I}{A}\kappa_{,ss} + \frac{w}{R^2} \right)
$$

**整理後得到 (Rearranging terms):**

$$
w_{,ss} + \frac{w}{R^2} = \kappa - \left[ \frac{EI}{GAk} + \frac{I}{A} \right] \kappa_{,ss} \quad \text{--- (15)}
$$

公式 (15) 是一個二階常係數非齊次微分方程。數學上，這類方程的解由兩部分組成：齊次解 (Homogeneous solution, $w^h$) 加上 特解 (Particular solution, $w^p$) 15。

$$w = w^h + w^p \quad \text{--- (16)}$$

令公式 (15) 右邊為 0

$$\frac{d^2w}{ds^2} + \frac{1}{R^2}w = 0$$
假設$$w(s) = e^{rs}$$

一階導數：$w' = r e^{rs}$
二階導數：$w'' = r^2 e^{rs}$

$$(r^2 e^{rs}) + \frac{1}{R^2} (e^{rs}) = 0$$

提取公因式 $e^{rs}$（因為指數函數 $e^{rs}$ 永遠不為 0，我們可以把它消掉）：$$e^{rs} \left( r^2 + \frac{1}{R^2} \right) = 0$$剩下的括號內部分必須為 0，這就是所謂的特徵方程式：$$r^2 + \frac{1}{R^2} = 0$$

移項：$$r^2 = -\frac{1}{R^2}$$兩邊開根號：$$r = \pm \sqrt{-\frac{1}{R^2}}$$處理負號開根號（引入虛數單位 $i$，其中 $i^2 = -1$）：
$$r = \pm \sqrt{-1} \cdot \sqrt{\frac{1}{R^2}}$$
$$r = \pm i \frac{1}{R}$$

當特徵根為共軛複數 $\pm bi$ 時（這裡是 $\pm i \frac{1}{R}$），通解可以寫成複數指數形式：$$w(s) = C_1 e^{i \frac{s}{R}} + C_2 e^{-i \frac{s}{R}}$$利用歐拉公式 (Euler's Formula) $e^{ix} = \cos x + i \sin x$，這兩個複數指數項可以組合成實數的三角函數形式：

$$w^h = C_{w1} \cos\left(\frac{s}{R}\right) + C_{w2} \sin\left(\frac{s}{R}\right)\quad \text{--- (17)}$$

$$
w^p = w^p(\kappa) \quad \text{--- (18)}
$$

$u$ 的表達式：
$$u = R(\theta - w_{,s} - \frac{EI}{GAK} \kappa_{,s}) \quad \text{--- (19)}$$

## 曲率形函數 (Shape functions for $k$)

$$\kappa = \mathbf{H}_\kappa \mathbf{V} \quad \text{--- (20a)}$$

$$\mathbf{V} = [\kappa_1, \kappa_2, \kappa_3]^T \quad \text{--- (20b)}$$

$$ \mathbf{H}_\kappa = [h_1^\kappa, h_2^\kappa, h_3^\kappa] \quad \text{--- (20c)}$$

通用的拉格朗日插值公式:推導$h_1^\kappa$（對應節點 1）在 $\xi=0$ 時為 1；在 $\xi=1$ 和 $\xi=0.5$ 時為 0。
$$\begin{array}{l}
h_1^\kappa = \frac{(\xi - 1)(\xi - 0.5)}{0.5} \\
= 2 \cdot (\xi - 1)(\xi - 0.5) \\
= (\xi - 1)(2\xi - 1) \\
= (1 - \xi)(1 - 2\xi) \\
= 1 - 2\xi - \xi + 2\xi^2 \\
= \mathbf{1 - 3\xi + 2\xi^2} \\
\end{array} $$

通用的拉格朗日插值公式:$h_2^\kappa$（對應節點 2）在 $\xi=1$ 時為 1；在 $\xi=0$ 和 $\xi=0.5$ 時為 0。
$$\begin{array}{l}
h_2^\kappa = \frac{(\xi - 0)(\xi - 0.5)}{0.5} \\
= 2 \cdot \xi(\xi - 0.5) \\
= \xi(2\xi - 1) \\
= -\xi + 2\xi^2 \\
\end{array} $$

通用的拉格朗日插值公式:$h_3^\kappa$（對應節點 3）在 $\xi=0.5$ 時為 1；在 $\xi=0$ 和 $\xi=1$ 時為 0。
$$\begin{array}{l}
h_3^\kappa = \frac{(\xi - 0)(\xi - 1)}{-0.25} \\
= -4 \cdot \xi(\xi - 1) \\
= 4\xi - 4\xi^2 \\
\end{array} $$

$\xi = s/L$

$$
\begin{cases}
h_1^\kappa = 1 - 3\frac{s}{L} + 2\left(\frac{s}{L}\right)^2 \quad \text{--- (20d)}\\
h_2^\kappa = -\frac{s}{L} + 2\left(\frac{s}{L}\right)^2 \quad \text{--- (20e)}\\
h_3^\kappa = 4\frac{s}{L} - 4\left(\frac{s}{L}\right)^2 \quad \text{--- (20f)}
\end{cases}
$$

## 截面轉角的推導

根據公式 (14):

$$\theta = \int \kappa \, ds + C_\theta $$

將公式 (20a) 的 $\kappa = \mathbf{H}_\kappa \mathbf{V}$ 代入積分式：

$$\theta = \int (\mathbf{H}_\kappa \mathbf{V}) \, ds + C_\theta$$

由於節點值 $\mathbf{V}$ 是常數向量，可提出積分外：

$$\theta = \left( \int \mathbf{H}_\kappa \, ds \right) \mathbf{V} + C_\theta$$

定義新的形狀函數向量 $\mathbf{H}_\theta = \int \mathbf{H}_\kappa \, ds$，

即得到
$$\theta = \mathbf{H}_\theta \mathbf{V} + C_\theta \quad \text{--- (21a)}$$

轉角形狀函數 $\mathbf{H}_\theta$ 的積分細節，對 $h^\kappa$ 的各項對 $s$ 進行積分。注意變數變換：$d(s/L) = \frac{1}{L} ds \Rightarrow ds = L \cdot d(s/L)$。$h_1^\theta$ 的推導:

$$\int h_1^\kappa \, ds = \int \left( 1 - 3\frac{s}{L} + 2\left(\frac{s}{L}\right)^2 \right) \, ds$$

積分結果：

$$= s - \frac{3}{2L}s^2 + \frac{2}{3L^2}s^3$$

提取 $L$ 作為公因數，並重寫為 $s/L$ 的形式：

$$= L \left\{ \frac{s}{L} - \frac{3}{2}\left(\frac{s}{L}\right)^2 + \frac{2}{3}\left(\frac{s}{L}\right)^3 \right\} \quad \text{--- (21b)}$$

$h_2^\theta$ 的推導:

$$\int h_2^\kappa \, ds = \int \left( -\frac{s}{L} + 2\left(\frac{s}{L}\right)^2 \right) \, ds$$

積分結果：

$$= -\frac{1}{2L}s^2 + \frac{2}{3L^2}s^3$$

提取 $L$：

$$= L \left\{ -\frac{1}{2}\left(\frac{s}{L}\right)^2 + \frac{2}{3}\left(\frac{s}{L}\right)^3 \right\} \quad \text{--- (21c)}$$

$h_3^\theta$ 的推導:

$$\int h_3^\kappa \, ds = \int \left( 4\frac{s}{L} - 4\left(\frac{s}{L}\right)^2 \right) \, ds$$

積分結果：

$$= \frac{4}{2L}s^2 - \frac{4}{3L^2}s^3 = \frac{2}{L}s^2 - \frac{4}{3L^2}s^3$$

提取 $L$：

$$= L \left\{ 2\left(\frac{s}{L}\right)^2 - \frac{4}{3}\left(\frac{s}{L}\right)^3 \right\} \quad \text{--- (21d)}$$

$$\begin{cases}
h_1^\theta = L \left\{ \frac{s}{L} - \frac{3}{2}\left(\frac{s}{L}\right)^2 + \frac{2}{3}\left(\frac{s}{L}\right)^3 \right\} \quad \text{--- (21b)} \\
h_2^\theta = L \left\{ -\frac{1}{2}\left(\frac{s}{L}\right)^2 + \frac{2}{3}\left(\frac{s}{L}\right)^3 \right\} \quad \text{--- (21c)} \\
h_3^\theta = L \left\{ 2\left(\frac{s}{L}\right)^2 - \frac{4}{3}\left(\frac{s}{L}\right)^3 \right\} \quad \text{--- (21d)}
\end{cases}$$

## 徑向位移 $w^p$ 的推導

根據公式 (15)和 (16)：
$$w_{,ss} + \frac{1}{R^2} w = \kappa - (\frac{EI}{GAk} + \frac{I}{A}) \kappa_{,ss} \quad (15)$$
$$w = w^h + w^p \quad (16)$$
已知齊次解 $w^h$ 如公式 (17) 所示，現在需要求解特解 $w^p$。

令 $w$ 的特解分量為 $h_1^{wp}$，對應的曲率形狀函數為 $h_1^\kappa$。方程式變成：$$(h_1^{wp})_{,ss} + \frac{1}{R^2}h_1^{wp} = h_1^\kappa - (\text{材料常數}) \cdot (h_1^\kappa)_{,ss}$$

設計形式：為了消掉左邊的 $\frac{1}{R^2}$，我們假設 $h_1^{wp} = R^2 (\dots)$。為了平衡常數項，我們假設它跟 $h_1^\kappa$ 差一個常數 $\lambda$。

所以假設：
$$h_1^{wp} = R^2 (h_1^\kappa - \lambda)$$

根據推導：$h_1^{wp} = R^2 (h_1^\kappa - \lambda)$ 
 
$$代入\begin{cases}
h_1^\kappa = 1 - 3\xi + 2\xi^2 \\
h_2^\kappa = -\xi + 2\xi^2 \\
h_3^\kappa = 4\xi - 4\xi^2
\end{cases}$$

$$\begin{cases}
h_1^{wp} = R^2 \{ 1 - \lambda - 3\xi + 2\xi^2 \} \quad \text{--- (22a)} \\
h_2^{wp} = R^2 \{ -\lambda - \xi + 2\xi^2 \} \quad \text{--- (22b)} \\
h_3^{wp} = R^2 \{ -\lambda + 4\xi - 4\xi^2 \} \quad \text{--- (22c)}
\end{cases}
$$

代入驗證：將假設代入微分方程左邊：
$$\text{左邊} = (R^2 h_1^\kappa)_{,ss} + \frac{1}{R^2} R^2 (h_1^\kappa - \lambda) = R^2(h_1^\kappa)_{,ss} + h_1^\kappa - \lambda$$

讓它等於微分方程右邊：
$$R^2(h_1^\kappa)_{,ss} + h_1^\kappa - \lambda = h_1^\kappa - \left( \frac{EI}{GAk} + \frac{I}{A} \right) (h_1^\kappa)_{,ss}$$

解出 $\lambda$：消去 $h_1^\kappa$，整理項：$$\lambda = R^2(h_1^\kappa)_{,ss} + \left( \frac{EI}{GAk} + \frac{I}{A} \right) (h_1^\kappa)_{,ss}$$

提取公因數 $(h_1^\kappa)_{,ss}$：$$\lambda = \left( R^2 + \frac{EI}{GAk} + \frac{I}{A} \right) (h_1^\kappa)_{,ss} \quad \text{--- (22d)}$$

$h_1^\kappa = 1 - 3\frac{s}{L} + 2(\frac{s}{L})^2$，其二階導數對 $s$ 微分兩次為：$$(h_1^\kappa)_{,ss} = \frac{4}{L^2}$$代入 $\lambda$：$$\lambda = \frac{4}{L^2} \left( R^2 + \frac{EI}{GAk} + \frac{I}{A} \right)$$




步1: $\theta = \int_0^s \kappa \, ds + C_\theta = H_\theta V + C_\theta \quad (21)$

步2: (21) 代入 (19)

$$u = R((H_\theta V + C_\theta) - w_{,s} - \frac{EI}{GAK} \kappa_{,s})$$

步3: 提取形函數
$$u = H_u V + R C_\theta + \underbrace{\text{齊次項}}_{\color{blue}{ C_{u1} \sin\frac{s}{R} - C_{u2} \cos\frac{s}{R}}}$$
$$\rightarrow u = H_u V + R C_\theta + C_{u1} \sin\frac{s}{R} - C_{u2} \cos\frac{s}{R} \quad (23a)$$

步4:
$h_1^\kappa = 1 - 3\eta + 2\eta^2 \longrightarrow 4$

$h_2^\kappa = -\eta + 2\eta^2 \longrightarrow 4$

$h_3^\kappa = 4\eta - 4\eta^2 \longrightarrow -8$

$-\beta \kappa_{,s} = -\beta \frac{4}{L^2} (\kappa_1 + \kappa_2 - 2\kappa_3)$

$\qquad \quad = -\beta \frac{4}{L^2} \kappa_1 - \beta \frac{4}{L^2} \kappa_1 + 2\beta \frac{4}{L^2} \kappa_3$

(註：手寫筆記展開式的第二項下標看似 $\kappa_1$，但根據上一行應為 $\kappa_2$)

RHS 有常數項，特解常數為

$w^p = \text{常數源} \cdot R^2$

$\rightarrow \kappa_1, \kappa_2 : -\beta \frac{4}{L^2} R^2 \quad \kappa_3 : 2\beta \frac{4}{L^2} R^2$

$\Rightarrow \mu = \frac{4}{L^2} (\frac{EI}{GAK} + R^2) \quad (23f)$

$$\begin{cases} 
w = H_{wp} V + C_{w1} \cos\frac{s}{R} + C_{w2} \sin\frac{s}{R} \\ u = H_u V + R C_\theta + C_{u1} \sin\frac{s}{R} - C_{u2} \cos\frac{s}{R} \\ \theta = H_\theta V + C_\theta
\end{cases}$$

分別代入端點值 : $s=0$ 和 $s=L$

at $s=0$ 時 $\cos(\frac{s}{R}) = 1$ , $\sin(\frac{s}{R}) = 0$

at $s=L$ 時 $\cos(\frac{L}{R})$ , $\sin(\frac{L}{R})$

$$\begin{cases} 
w(0) = H_{wp}(0) V + C_{w1} \cdot 1 + C_{w2} \cdot 0 & (a) \\ u(0) = H_u(0) V + R C_\theta + C_{u1} \cdot 0 - C_{u2} \cdot 1 & (b) \\ \theta(0) = H_\theta(0) V + C_\theta & (c) \\ w(L) = H_{wp}(L) V + C_{w1} \cos\frac{L}{R} + C_{w2} \sin\frac{L}{R} & (d) \\ u(L) = H_u(L) V + R C_\theta + C_{u1} \sin\frac{L}{R} - C_{u2} \cos\frac{L}{R} & (E) \\ \theta(L) = H_\theta(L) V + C_\theta & (F)
 \end{cases}$$

$w_2 - w_1 \cos\phi + u_1 \sin\phi - \theta_1 R \sin\phi$

$= H_{wp|L} V + \cancel{C_{w1} \cos\phi} + \cancel{C_{w2} \sin\phi} - (H_{wp|0} V + \cancel{C_{w1}}) \cos\phi + (H_{u|0} V + \cancel{RC_\theta} - \cancel{C_{u2}}) \sin\phi - (H_{\theta|0} V + \cancel{C_\theta}) R \sin\phi$

$= (H_{wp|L} - H_{wp|0} \cos\frac{L}{R} + H_{u|0} \sin\frac{L}{R} - H_{\theta|0} R \sin\frac{L}{R}) V \quad (25a) \quad (\color{blue}{\text{設 } C_{u2} \approx C_{w2}}) $

$u_2 - w_1 \sin\phi - u_1 \cos\phi - \theta_1 R (1 - \cos\phi)$

$= [H_{u|L} - H_{wp|0} \sin\phi - H_{u|0} \cos\phi - H_{\theta|0} R (1 - \cos\phi)] V \quad (25b)$

$\theta_L - \theta_0 = (H_{\theta|L} - H_{\theta|0}) V + (C_\theta - C_\theta)$

$\qquad \quad = H_{\theta|L} V \qquad \qquad \qquad (25c)$

$T_\kappa V = T_u U$

$U = [w_1, u_1, \theta_1, w_2, u_2, \theta_2]^T$

$$T_\kappa = \begin{bmatrix}
 H_{wp|L} - H_{wp|0} \cos\frac{L}{R} + H_{u|0} \sin\frac{L}{R} \\ H_{u|L} - H_{wp|0} \sin\frac{L}{R} - H_{u|0} \cos\frac{L}{R} \\ H_{\theta|L} 
 \end{bmatrix}$$

$$T_u = \begin{bmatrix} 
-\cos\frac{L}{R} & \sin\frac{L}{R} & -R\sin(\frac{L}{R}) & 1 & 0 & 0 \\ -\sin\frac{L}{R} & -\cos\frac{L}{R} & -R(1-\cos\phi) & 0 & 1 & 0 \\ 0 & 0 & -1 & 0 & 0 & 1 
\end{bmatrix}$$

$V = T U$

$T = T_\kappa^{-1} T_u$

$\kappa(s) = H_\kappa U \theta(s) = H_\theta V + C_\theta$

代入 $V = TU$

$\kappa(s) = H_\kappa T U = S_\kappa U$

$\theta(s) = (H_\theta T + H_{\theta|0}) V = S_\theta U \quad \color{blue}{\text{因 }\theta \text{ 元素內有常數偏移 (剛體方包轉)，需 } H_{\theta|0} = [001000]}$

$w(s) = (H_{wp} - H_{wp|0} \cos\frac{s}{R} - H_{u|0} \sin\frac{s}{R}) T + H_{w0} = S_w U$

$u(s) = (H_u - H_{wp|0} \sin\frac{s}{R} - H_{u|0} \cos\frac{s}{R}) T + H_{u0} = S_u U$

總位能公式
$$\pi = \underbrace{\frac{1}{2} \int_0^L EI \kappa^2 ds}_{\text{彎曲}} + \underbrace{\frac{1}{2} \int_0^L GAK \gamma^2 ds}_{\text{剪切}} + \underbrace{\frac{1}{2} \int_0^L EA \varepsilon^2 ds}_{\text{軸向力}} - \underbrace{\int_0^L (m\theta + p_w w + p_u u) ds}_{\text{外荷載}}$$

用端點位移 $U$ 表示

$(28a-d)$

$\kappa(s) = S_{\kappa(s)} U$

$\theta(s) = S_{\theta(s)} U$

$w(s) = S_{w(s)} U$

$u(s) = S_{u(s)} U$

同時應變

曲率 $\kappa = S_\kappa U$

剪切 $\gamma = B_\gamma (\frac{dS_\kappa}{ds}) U \quad (\gamma \sim \kappa_{,s}) \rightarrow B_\gamma = \frac{d}{ds} H_\kappa$

膜 $\varepsilon = Q_\alpha (\frac{d^2S_\kappa}{ds^2}) U \quad (\varepsilon \sim \kappa_{,ss}) \rightarrow Q_\alpha = \frac{d^2}{ds^2} H_\kappa$

最小位能原理應用 ($\delta \pi = 0$)

內能變分:

$\delta (\frac{1}{2} \int E I \kappa^2 ds) = \int E I \kappa \, \delta \kappa \, ds = \int EI (S_\kappa U)^T \delta (S_\kappa U) ds = \delta U^T (EI \int S_\kappa^T S_\kappa ds) U$

$\delta (\frac{1}{2} \int G A K \gamma^2 ds) = \int G A K \gamma \, \delta \gamma \, ds = \int G A K (B_\gamma U)^T \delta (B_\gamma U) ds = \delta U^T (GAK \int B_\gamma^T B_\gamma) U$

$\delta (\frac{1}{2} \int E A \varepsilon^2 ds) = \int E A \varepsilon \, \delta \varepsilon \, ds = \int EA (Q_\varepsilon U)^T \delta (Q_\varepsilon U) \, ds = \delta U^T (EA \int Q_\varepsilon^T Q_\varepsilon)U $

$$\delta U_{int} = \delta U^T [ EI \int_0^L S_\kappa^T S_\kappa ds + GAK \int_0^L B_\gamma^T B_\gamma ds + EA \int_0^L Q_\alpha^T Q_\alpha ds ] U$$

$$= \delta V^T [ EI \int_0^L H_\kappa^T H_\kappa ds + GAK \int_0^L B_\kappa^T B_\kappa ds + EA \int_0^L Q_\kappa^T Q_\kappa ds ] V$$

因 $V=TU, \delta V = T \delta U$ 代入得 $\delta U^T T^T [\dots\dots] T U$
$\rightarrow K = T^T [\dots\dots] T$

外變分 $-\delta W_{ext} = - (\int m_b \delta \kappa ds + \int p_w \delta w ds + \int p_u \delta u ds)$

$$\delta (\int m_b \kappa ds) = \int m_b \delta \kappa ds = m_b \int \delta (S_\kappa U)^T ds = \delta U^T (m_b \int S_\kappa^T ds)$$

$$\delta (\int p_w w ds) = \int p_w \delta w ds = p_w \int \delta (S_w U)^T ds = \delta U^T (p_w \int S_w^T ds)$$

$$\delta (\int n_u u ds) = \int n_u \delta u ds = n_u \int \delta (S_u U)^T ds = \delta U^T (n_u \int S_u^T ds)$$

$$\rightarrow \delta U^T (\int S_\kappa^T m_b ds + \int S_w^T p_w ds + \int S_u^T n_u ds)$$

$$\delta \pi = \delta U_{int} - \delta W_{ext} = 0$$

$$\delta U^T K U = \delta U^T R$$

$$\Rightarrow K U = R$$

$$K = T^T [ EI \int_0^L H_\kappa^T H_\kappa ds + \alpha \int_0^L B_\gamma^T B_\gamma ds + \beta \int_0^L Q_\alpha^T Q_\alpha ds ] T$$

$$R = \int_0^L S_\theta^T m_b ds + \int_0^L S_w^T p_w ds + \int_0^L S_u^T n_u ds$$

$$\alpha = \frac{EI}{GAK} \quad \beta = \frac{EI}{EA}$$

