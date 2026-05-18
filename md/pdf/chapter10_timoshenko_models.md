# 第10章 Timoshenko 三種模型整理

## 第10章快速導覽

| 章節 | 主題 |
|------|------|
| 10.2 | 靜態分析建模：二節點 Timoshenko 樑插值、元素剛度、靜態位移解析式 |
| 10.3 | 自由振動：動能分解、元素剛度與質量矩陣、廣義特徵值問題、無因次頻率 |
| 10.4 | 屈曲分析：屈曲特徵值問題、幾何剛度矩陣、含有效長度的臨界荷載解析式 |

---

## 模型一：靜態彎曲 (10.2)

### 10.2 前置基本公式 (10.1-10.8)

以下為第10章在元素離散前的連續體層級核心公式：

1. 位移場 (10.1)

$$
u = y\theta_z, \quad w = w_0
$$

2. 正應變 (10.2)

$$
\varepsilon_x = \frac{\partial u}{\partial x} = y\frac{\partial \theta_z}{\partial x}
$$

3. 剪應變 (10.3)

$$
\gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial w}{\partial x} = \theta_z + \frac{\partial w}{\partial x}
$$

4. 應變能分解 (10.4)

$$
U = \frac{1}{2}\int_V \sigma_x\varepsilon_x \, dV + \frac{1}{2}\int_V \tau_{xy}\gamma_{xy} \, dV
$$

5. 本構關係 (10.5)-(10.6)

$$
\sigma_x = E\varepsilon_x, \quad \tau_{xy} = kG\gamma_{xy}
$$

6. 剪切模數 (10.7)

$$
G = \frac{E}{2(1+\nu)}
$$

7. 積分到樑軸後的能量式 (10.8)

$$
U = \frac{1}{2}\int_0^L EI_z\left(\frac{d\theta_z}{dx}\right)^2dx
  + \frac{1}{2}\int_0^L kGA\left(\frac{dw}{dx}+\theta_z\right)^2dx
$$

註：剪應變項符號與轉角正向定義有關，常見等價寫法也可寫為 $\left(dw/dx-\theta\right)$。

### 二節點插值 (對應 10.10-10.11)

自然座標 $\xi \in [-1, 1]$：

$$
N_1(\xi) = \frac{1}{2}(1 - \xi), \quad
N_2(\xi) = \frac{1}{2}(1 + \xi)
$$

場變量插值（二節點元素）：

$$
w = \mathbf{N} \cdot \mathbf{w}_e, \quad
\theta_z = \mathbf{N} \cdot \boldsymbol{\theta}_e
$$

### 元素剛度 (10.12 形式)

元素剛度由彎曲項與剪切項組成：

$$
\mathbf{K}_e = \int_{-1}^{1} \left[ \frac{E I_z}{a^2} \left( \frac{d\mathbf{N}}{d\xi} \right)^T \left( \frac{d\mathbf{N}}{d\xi} \right) \right] a \, d\xi
+ \int_{-1}^{1} \left[ kGA \left( \frac{1}{a}\frac{d\mathbf{N}}{d\xi} + \mathbf{N} \right)^T \left( \frac{1}{a}\frac{d\mathbf{N}}{d\xi} + \mathbf{N} \right) \right] a \, d\xi
$$

其中 $a = L_e / 2$，為 $\xi \to x$ 的 Jacobian 尺度。

> **選擇性積分策略**
> - 彎曲項：2-point Gauss
> - 剪切項：1-point Gauss  
>   （用來減輕 shear locking）

### 靜態解析位移參考 (10.13-10.15)

**簡支薄樑最大位移：**

$$
w_{\text{max}} = \frac{5PL^4}{384EI} \quad \text{(對應 10.13)}
$$

**簡支 Timoshenko 樑位移：**

$$
w(x) = \frac{PL^4}{24D}\left( -\frac{x}{L} + 2\left(\frac{x}{L}\right)^3 - \left(\frac{x}{L}\right)^4 \right)
+ \frac{PL^2}{2S}\left( \frac{x}{L} - \left(\frac{x}{L}\right)^3 \right) \quad \text{(對應 10.14)}
$$

**懸臂 Timoshenko 樑位移：**

$$
w(x) = \frac{PL^4}{24D}\left( 6\left(\frac{x}{L}\right)^2 - 4\left(\frac{x}{L}\right)^3 + \left(\frac{x}{L}\right)^4 \right)
+ \frac{PL^2}{2S}\left( 2\frac{x}{L} - \left(\frac{x}{L}\right)^2 \right) \quad \text{(對應 10.15)}
$$

### 變量定義

| 符號 | 意義 |
|------|------|
| $E$ | 楊氏模數 |
| $G$ | 剪切模數 |
| $I$ 或 $I_z$ | 截面二次矩 |
| $A$ | 截面面積 |
| $k$ | 剪切修正係數 |
| $D$ | 彎曲剛度，$D = \dfrac{Eh^3}{12(1-\nu^2)}$ |
| $S$ | 剪切剛度，$S = kGA$ |
| $P$ | 均布載重參數 |
| $L$ | 樑長 |

---

## 模型二：自由振動 (10.3)

### 動能分解 (10.16)

總動能包含平移與轉動兩部分：

$$
T = \frac{1}{2} \int_0^L \rho A \dot{w}^2 \, dx
+ \frac{1}{2} \int_0^L \rho I_z \dot{\theta}_z^2 \, dx
$$

### 元素剛度與質量矩陣 (10.17-10.18)

**剛度矩陣**（與靜態結構形式一致）：

$$
\mathbf{K}_e = \int_{-1}^{1} \left[ \frac{E I_z}{a^2} \left( \frac{d\mathbf{N}}{d\xi} \right)^T \left( \frac{d\mathbf{N}}{d\xi} \right)
+ kGA \left( \frac{1}{a}\frac{d\mathbf{N}}{d\xi} + \mathbf{N} \right)^T \left( \frac{1}{a}\frac{d\mathbf{N}}{d\xi} + \mathbf{N} \right) \right] a \, d\xi
$$

**質量矩陣**：

$$
\mathbf{M}_e = \int_{-1}^{1} \rho A (\mathbf{N}^T \mathbf{N}) \, a \, d\xi
+ \int_{-1}^{1} \rho I_z (\mathbf{N}^T \mathbf{N}) \, a \, d\xi
$$

### 廣義特徵值問題與無因次頻率 (10.19)

組裝並施加邊界條件後：

$$
\left[ \mathbf{K} - \omega^2 \mathbf{M} \right] \mathbf{X} = \mathbf{0}
$$

第10章使用的**無因次頻率參數**：

$$
\bar{\omega} = \omega L^2 \sqrt{\frac{\rho A}{E I_z}}
$$

### 變量定義

| 符號 | 意義 |
|------|------|
| $\rho$ | 質量密度 |
| $\omega$ | 圓頻率 |
| $\mathbf{X}$ | 振型向量 |
| $\mathbf{K}, \mathbf{M}$ | 全域剛度與質量矩陣 |
| $\bar{\omega}$ | 無因次頻率 |

---

## 模型三：屈曲分析 (10.4)

### 屈曲特徵值問題 (10.20)

$$
\left[ \mathbf{K} - \lambda \mathbf{K}_g \right] \mathbf{X} = \mathbf{0}
$$

其中 $\lambda$ 為臨界荷載，$\mathbf{X}$ 為屈曲模態。

### 幾何剛度矩陣 (10.21)

$$
\mathbf{K}_g = \int_0^L P \left( \frac{d\mathbf{N}_w}{dx} \right)^T \left( \frac{d\mathbf{N}_w}{dx} \right) \, dx
$$

### 臨界荷載解析參考 (10.22)

$$
P_{cr} = \frac{\pi^2 E I}{L_{eff}^2 \left( 1 + \frac{\pi^2 E I}{kGA L_{eff}^2} \right)}
$$

**有效長度**：

| 邊界條件 | $L_{eff}$ |
|----------|-----------|
| pinned-pinned | $L$ |
| fixed-fixed | $L/2$ |

---

## 重建說明

> 由於 OCR 文字品質有限，以上公式已整理為可讀的標準工程記號，同時保持第10章原意。

**公式編號對應表：**

| 模型類型 | 公式編號 |
|----------|----------|
| 靜態 | 10.12 - 10.15 |
| 振動 | 10.16 - 10.19 |
| 屈曲 | 10.20 - 10.22 |
