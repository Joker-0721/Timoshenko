# 鐵摩辛柯梁無鎖定（Locking-Free）公式入門

## 鐵摩辛柯梁無鎖定（Locking-Free）公式入門

本文整理與鐵摩辛柯梁（Timoshenko beam）屈曲問題之「無鎖定」有限元素（locking-free FEM）公式相關的基本概念。

## 鐵摩辛柯梁理論

鐵摩辛柯梁理論是對歐拉–伯努力梁（Euler–Bernoulli beam）理論的延伸，能同時考慮剪力變形與截面轉動所造成的彎曲效應。
它對厚梁、或高頻振動等情況特別重要。

其運動學假設為：
截面在變形後仍保持平面，但不一定仍與梁的縱向軸線保持垂直。
控制方程涉及兩個相互獨立的未知量：橫向位移 $w(x)$ 與截面轉角 $\phi(x)$。

## 剪力鎖定（Shear Locking）現象

剪力鎖定是有限元素分析中的一種數值人為效應（numerical artifact）。當對具有顯著剪力變形的問題（例如鐵摩辛柯梁）使用標準低階元素時，常會出現此現象。

**成因：** 當梁非常細長（薄）時，元素難以正確表示「純彎曲」狀態；剪應變被人為地限制為非零，導致結構回應出現「鎖住」且過度剛硬（overly stiff）。

**影響：** 隨著網格細化，計算得到的位移會趨近於零，而非正確解；對薄梁尤其嚴重。

## 無鎖定（Locking-Free）公式

為了克服剪力鎖定，已發展出多種進階公式。
其目標是讓元素對厚梁與薄梁皆能表現正確。常見方法包含：

* **降階積分（Reduced Integration）：** 對元素剛度矩陣中的剪力項採用較低階的積分規則。有時會引發其他問題，例如虛假零能量模態（spurious zero-energy modes）。
* **混合式公式（Mixed Formulations）：** 將位移場與應力／應變場視為彼此獨立的未知量；Hellinger–Reissner 原理常作為此類方法的基礎。
* **假設應變／應力法（Assumed Strain/Stress Methods）：** 在元素內對應變或應力場引入特定假設，使其滿足必要的一致性條件（consistency conditions）。

## 在屈曲分析中的應用

屈曲分析的核心在於求得結構失穩時的臨界荷重。
對鐵摩辛柯梁而言，這通常會導致一個特徵值問題（eigenvalue problem）。
剪力鎖定會顯著高估臨界屈曲荷重，因此在細長構件中（屈曲為主要控制因素）採用無鎖定公式對精確預測尤為關鍵。

### 屈曲分析的能量泛函

屈曲分析可由「位能駐值原理」（principle of stationary potential energy）出發。
梁在軸向壓縮力 $P$ 作用下，其總位能 $\Pi$ 由應變能 $U$ 與外力位能 $V$ 所構成：

$$
\Pi=U-V
$$

### 應變能

應變能 $U$ 可分為彎曲應變能 $U_{b}$ 與剪力應變能 $U_{s}$ 之和。

1.  **彎曲應變能（$U_{b}$）：** 由彎曲造成的儲能。
    $$
    U_{b}=\frac{1}{2}\int_{0}^{L}EI\phi_{,x}^{2}dx
    $$
    其中 $E$ 為楊氏模數、$I$ 為截面二次矩，而 $\phi$ 為截面轉角。

2.  **剪力應變能（$U_s$）：** 由剪力變形造成的儲能。
    $$
    U_{s}=\frac{1}{2}\int_{0}^{L}kAG(w_{,x}-\phi)^{2}dx
    $$
    其中 $k$ 為剪力修正係數、$A$ 為截面面積、$G$ 為剪力模數，而 $w$ 為橫向位移。


### 外力位能

軸向壓縮力 $P$ 的位能（$V$）可視為該力因梁的橫向變形所做的功：
$$
V=\frac{1}{2}\int_{0}^{L}Pw_{,x}^{2}dx
$$

### 總位能泛函

將上述各項合併，鐵摩辛柯梁屈曲的總位能泛函可寫成：
$$
\Pi=\frac{1}{2}\int_{0}^{L}(EI\phi_{,x}^{2}+kAG(w_{,x}-\phi)^{2}-Pw_{,x}^{2})dx
$$

臨界屈曲荷重 $P_{cr}$ 是使得 $\Pi$ 的二次變分不再為正定（positive definite）的最小 $P$，因此會導出一個特徵值問題。

## 變分原理與控制方程

平衡狀態可由位能駐值原理求得：對任何允許的虛位移（virtual displacement），總位能的一次變分必須為零。

$$
\delta\Pi=\delta\frac{1}{2}\int_{0}^{L}(EI\phi_{,x}^{2}+kAG(w_{,x}-\phi)^{2}-Pw_{,x}^{2})dx=0
$$

對獨立未知量 $w$ 與 $\phi$ 取變分可得：
$$
\int_{0}^{L}(\delta\phi_{,x}EI\phi_{,x}+(\delta w_{,x}-\delta\phi)kAG(w_{,x}-\phi)-\delta w_{,x}Pw_{,x})dx=0
$$

將各項展開並依獨立變分 $\delta w$ 與 $\delta\phi$ 分組後，弱式可拆成兩個互相耦合的方程，對應到矩陣的區塊結構：

**對 $w$ 的弱式：**
$$
\int_{0}^{L}(\delta w_{,x}kAGw_{,x}-\delta w_{,x}kAG\phi-\delta w_{,x}Pw_{,x})dx=0
$$

**對 $\phi$ 的弱式：**
$$
\int_{0}^{L}(\delta\phi_{,x}EI\phi_{,x}-\delta\phi kAGw_{,x}+\delta\phi kAG\phi)dx=0
$$

上述分離後的弱式分別對應剛度矩陣的組成：含 $\delta w_{,x}w_{,x}$ 的項貢獻到 $K_{ww}$；含 $\delta\phi_{,x}\phi_{,x}$ 與 $\delta\phi\,\phi$ 的項貢獻到 $K_{\phi\phi}$；而交叉項則貢獻到耦合矩陣 $K_{w\phi}$ 與 $K_{\phi w}$。對這些弱式做分部積分（integration by parts）即可得到控制微分方程與鐵摩辛柯梁屈曲問題的自然邊界條件。

## 有限元素近似

為了以數值方法解屈曲問題，可將梁域離散成有限元素，並以形狀函數近似場變數 $w$ 與 $\phi$：

$$
w\approx w^{h}=\sum_{I=1}^{n_{w}}N_{I}(x)w_{I}
$$
$$
\phi\approx\phi^{h}=\sum_{J=1}^{n_{\phi}}N_{J}(x)\phi_{J}
$$

其中 $N_{I}(x)$ 與 $N_{J}(x)$ 為形狀函數，$w_{I}$ 與 $\phi_{J}$ 為節點值，$n_{w}$ 與 $n_{\phi}$ 分別代表 $w$ 與 $\phi$ 的總自由度數。

其導數近似為：
$$
w_{,x}\approx w_{,x}^{h}=\sum_{I=1}^{n_{w}}N_{I,x}(x)w_{I}
$$
$$
\phi_{,x}\approx\phi_{,x}^{h}=\sum_{J=1}^{n_{\phi}}N_{J,x}(x)\phi_{J}
$$

將上述近似代入變分式，即可得到離散後的特徵值問題：

$$
\begin{bmatrix}K_{ww}&K_{w\phi}\\ K_{\phi w}&K_{\phi\phi}\end{bmatrix}\begin{bmatrix}d_{w}\\ d_{\phi}\end{bmatrix}-P\begin{bmatrix}K_{,ww}^{G}&0\\ 0&0\end{bmatrix}\begin{bmatrix}d_{w}\\ d_{\phi}\end{bmatrix}=\begin{bmatrix}0\\ 0\end{bmatrix}
$$

其中各元素剛度子矩陣定義為：


$$
K_{ww,IJ}=\int_{0}^{L}(kAGN_{I,x}N_{J,x}-PN_{I,x}N_{J,x})dx
$$
$$
K_{\phi\phi,IJ}=\int_{0}^{L}(EIN_{I,x}N_{J,x}+kAGN_{I}N_{J})dx
$$
$$
K_{w\phi,IJ}=K_{\phi w,JI}=-\int_{0}^{L}kAGN_{I,x}N_{J}dx
$$
$$
K_{,ww,IJ}^{G}=\int_{0}^{L}N_{I,x}N_{J,x}dx
$$

此處 $N_{I}$ 與 $N_{J}$ 為形狀函數，$N_{I,x}$ 與 $N_{J,x}$ 為其對 $x$ 的導數，$I,J$ 為控制點（或節點）索引。
全域剛度矩陣由各元素貢獻組裝而成。