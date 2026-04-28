# `mindlin_plate` 測試流程與標準寫法導讀

## 1. 這個 repo 是什麼性質

`mindlin_plate` 不是一個完整的 Julia package。它沒有 `Project.toml`，所以不是走「`using 某個本地套件`」的封裝方式，而是比較典型的研究腳本集合。

它的執行前提是外部環境已經可以直接使用：

- `ApproxOperator`
- `BenchmarkExample`
- `Gmsh`
- 部分腳本還會用到 `TimerOutputs`、`LinearAlgebra`、`WriteVTK`

所以這個 repo 的定位，不是「提供新 API」，而是「依靠既有 API 寫出一組可直接跑的 Mindlin plate mixed 測試腳本」。

這點很重要，因為它代表你閱讀時不要先找主模組或 package 入口，而要直接把它當成：

1. 生成 `.msh` 的腳本
2. 讀 `.msh` 並求解的腳本

來看。

---

## 2. 這個 repo 內最重要的角色分工

這個 repo 的主線可以分成三層：

### 第一層：生成 mesh

- `generateCircle.jl`
- `generatePatchtest.jl`

這兩個檔案本身不做求解，只是呼叫 `BenchmarkExample` 去批量產生 `.msh`。

### 第二層：真正做 mixed 測試

- `circular.jl`
- `patch_test.jl`
- `square.jl`

這三個才是真正做 Mindlin plate 混合形式測試的核心腳本。

### 第三層：輔助分析

- `dofs.jl`

這個檔案是在做自由度數量的推算與比較，不屬於主求解流程。

所以如果你是新手，要先分清楚：

- `generate*.jl` 不是 solver
- `circular / patch_test / square` 才是 solver

---

## 3. 前處理是怎樣接上的

這個 repo 生成 mesh 的能力，實際上不是自己寫幾何，而是借用你本機的 `D:\BenchmarkExample.jl`。

### `generateCircle.jl`

它的寫法很簡單，核心就是：

```julia
import BenchmarkExample: Circular

for n in 3:15
    Circular.generateMsh("msh/circular_tri3_$n.msh", transfinite=n, order=1, quad=false)
end
```

這段代表：

- 直接用 `BenchmarkExample.Circular.generateMsh(...)`
- 批量產生一系列圓板三角形一階網格
- 檔名規則是 `msh/circular_tri3_n.msh`

這裡亦留有註解掉的 `tri6` 版本，表示它原本就考慮過二階元素，只是目前主線先用 `tri3`。

### `generatePatchtest.jl`

這個檔會批量產生方板/patch test 的多種元素族：

- `patchtest_tri3_n.msh`
- `patchtest_tri6_n.msh`
- `patchtest_quad4_n.msh`
- `patchtest_quad8_n.msh`

核心寫法是：

```julia
import BenchmarkExample: PatchTest

PatchTest.generateMsh(..., order=1, quad=false)
PatchTest.generateMsh(..., order=2, quad=false)
PatchTest.generateMsh(..., order=1, quad=true)
PatchTest.generateMsh(..., order=2, quad=true)
```

也就是說，這個 repo 的 mesh 前處理標準手法是：

1. 幾何與 physical groups 放在 `BenchmarkExample`
2. 這邊只負責決定元素型別、階數、以及檔名批量生成規則

### 為什麼它依賴 `D:\BenchmarkExample.jl`

因為 `.msh` 的真正來源不是 `mindlin_plate` 自己寫的，而是：

- `BenchmarkExample.PatchTest.generateMsh`
- `BenchmarkExample.Circular.generateMsh`

這表示幾何定義、physical group 命名、Gmsh 寫檔，都是在 `BenchmarkExample` 裏完成。  
`mindlin_plate` 這邊做的是「取用已生成的 mesh 進行數值實驗」。

---

## 4. 主求解腳本的共同骨架

三個主腳本雖然案例不同，但寫法骨架幾乎一致。這正是你課題內的標準寫作方法。

### Step 1：載入 `ApproxOperator` 及 MindlinPlate operator

它們都先這樣寫：

```julia
using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ...
```

這個寫法非常標準，表示：

- `ApproxOperator` 提供整體框架
- `GmshImport` 負責讀 `.msh`
- `MindlinPlate` 提供弱式 operator

### Step 2：先定材料與解析量

這些腳本不是先建矩陣，而是先定義：

- 材料常數 `E, ν, h`
- 彎曲剛度 `Dᵇ`
- 剪切剛度 `Dˢ`
- 製造解或參考解：
  - `w`
  - `φ₁`, `φ₂`
  - 以及導數、彎矩 `M`、剪力 `Q`
  - 域內荷載 `q`
  - 轉角方程右端 `m₁`, `m₂`

這套寫法反映一個很明確的觀念：

- 先把「理論場」定義好
- 再把它交給 `prescribe!`

也就是說，數學模型先於矩陣組裝。

### Step 3：用 Gmsh 開 `.msh`

三個腳本都會做：

```julia
gmsh.initialize()
gmsh.open("msh/xxx.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
```

這一步的資料流是：

`.msh -> physical groups / nodes -> ApproxOperator 內部資料結構`

### Step 4：手動建立 block matrix

這個 repo 的 mixed 方法不是單一剛度矩陣，而是拆成：

- `kʷʷ`
- `kᵠᵠ`
- `kᵠʷ`
- `fʷ`
- `fᵠ`

這就是它的核心特徵。

因為 Mindlin plate 在這裡不是單場問題，而是：

- 一個位移場 `w`
- 一個轉角場 `φ = (φ₁, φ₂)`

所以矩陣分塊寫成：

```julia
[kᵠᵠ  kᵠʷ
 kᵠʷ' kʷʷ]
```

這種 block 寫法就是這份 repo 最有代表性的標準風格。

### Step 5：先做 domain 組裝

它們都會先抓域元素：

```julia
elements = getElements(nodes, entities["Ω"])
prescribe!(elements, ...)
set∇𝝭!(elements)
```

然後再定 operator：

```julia
𝑎ʷʷ = ∫wwdΩ => elements
𝑎ᵠʷ = ∫φwdΩ => elements
𝑎ᵠᵠ = [∫φφdΩ => elements, ∫κκdΩ => elements]
𝑓ʷ = ∫wqdΩ => elements
𝑓ᵠ = ∫φmdΩ => elements
```

然後顯式組裝：

```julia
𝑎ʷʷ(kʷʷ)
𝑎ᵠʷ(kᵠʷ)
𝑎ᵠᵠ(kᵠᵠ)
𝑓ʷ(fʷ)
𝑓ᵠ(fᵠ)
```

這裡很值得注意：

- operator 永遠寫成 `∫... => elements`
- 不會把一切包進一個黑箱 solver
- 組裝是顯式、可分塊、可檢查的

這就是你課題組這套寫作方法的最大特徵之一。

### Step 6：再做 boundary 組裝

邊界不直接和 domain 混在一起，而是另外抓各自的 boundary elements：

- `getElements(nodes, entities["Γ¹"])`
- `getElements(nodes, entities["Γ²"])`
- `getElements(nodes, entities["Γᵉ"])`
- `getElements(nodes, entities["Γˡ"], normal=true)`

接著：

- 邊界元素只做 `set𝝭!`
- 因為邊界項不需要域內導數那一套 shape function gradient

這點非常關鍵：

- 域內組裝：`set∇𝝭!`
- 邊界與誤差：`set𝝭!`

這就是它們在實作上的固定分工。

### Step 7：用 penalty 和自然邊界項組裝

常見 operator 包括：

- `∫αwwdΓ`
- `∫αφφdΓ`
- `∫wVdΓ`
- `∫φMdΓ`

這表示它的邊界處理採混合做法：

- 本質邊界條件：多用 penalty 強制
- 自然邊界條件：用邊界荷載 operator 額外加到右端

### Step 8：直接解 block system

最後是：

```julia
d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]
```

這一步完全保持矩陣結構清楚，不做多餘抽象化。

### Step 9：把結果回填到節點

它們不會只把 `d` 留成一條長向量，而是：

```julia
push!(nodes, :d=>..., :d₁=>..., :d₂=>...)
```

也就是：

- `:d` 對應板位移
- `:d₁`, `:d₂` 對應兩個轉角自由度

這樣做的好處是：

- 誤差計算可以直接用 nodes 內的欄位
- 後處理也沿用同一份資料流

### Step 10：再另外做 L2 驗證

最後都會重新抓一份誤差評估元素：

```julia
elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u=>w, :φ₁=>φ₁, :φ₂=>φ₂, ...)
set𝝭!(elements)
L₂_w = L₂(elements)
L₂_φ = L₂φ(elements)
```

這表示：

- 求解用的 mesh / 積分規則是一回事
- 誤差評估另走一條獨立流程

這也是很典型的研究腳本寫法。

---

## 5. 三個主腳本各自做什麼

## 5.1 `patch_test.jl`

這是最標準、最像模板的版本。

它的特色是：

- 讀的是 `msh/patchtest.msh`
- 使用非常簡單的製造解：
  - `w = 1 + x + y`
  - `φ₁ = 1 + x + y`
  - `φ₂ = 1 + x + y`
- 四條邊 `Γ¹ ~ Γ⁴` 全部用 penalty 強制
- 沒有額外自然邊界項
- 結尾直接做 `L₂_w` 和 `L₂_φ`

所以如果你要找「課題組最標準的 Mindlin plate mixed 測試模板」，應該先看這個。

它幾乎把標準骨架最乾淨地展現出來：

1. 製造解
2. domain 組裝
3. penalty 邊界
4. block solve
5. L2 驗證

---

## 5.2 `square.jl`

這個檔仍然是方板主線，但把製造解換成更高階多項式，難度明顯比 `patch_test.jl` 高。

它的特色是：

- 實際讀的是 `msh/patchtest_tri3_16.msh`
- 所以它雖然叫 `square.jl`，但不是用另一個獨立 square mesh，而是直接拿 patch test 系列中的三角形網格做測試
- `w`, `φ₁`, `φ₂`, `q` 都是高階多項式
- 仍然沿用四邊 penalty 的骨架
- 末尾有註解掉的 `WriteVTK` 後處理片段

這裡要特別提醒新手一點：

**`square.jl` 的名字容易讓人以為它對應另一份 square 專屬 mesh，但實際上它是讀 `patchtest_tri3_16.msh`。**

這種命名和實際檔案來源不完全一致的情況，在研究腳本中很常見，所以閱讀時一定要以 `gmsh.open(...)` 為準。

---

## 5.3 `circular.jl`

這是三個主腳本裡最複雜的一個。

它的特色是：

- 讀的是 `./msh/circular.msh`
- 除了 penalty 邊界之外，還用了自然邊界 operator：
  - `∫wVdΓ`
  - `∫φMdΓ`
- 邊界群組不是 `Γ¹~Γ⁴`，而是：
  - `Γᵉ`
  - `Γˡ`
  - `Γᵇ`
- 有 `USE_CIRCULAR_CLAMPED_REFERENCE` 開關
- 有 `ENABLE_MANUFACTURED_ERROR`、`PRINT_SANITY_CHECKS` 這類階段切換思路

### `circular.jl` 的特殊之處

它不是單純 patch test，而是同時保留了兩條思路：

1. 舊的製造解測試流程
2. 圓板 clamped reference 的思路

所以你會見到：

- 一邊定義 `w, φ, M, Q, q, m`
- 一邊又保留圓板精確解的註解區塊

這說明這個腳本已經從單純 manufactured solution 模板，擴展成更接近真實 benchmark 的研究版本。

### 圓板案例的邊界處理

`circular.jl` 最值得學的，是它如何把 penalty 與自然邊界混合起來：

- 在 `Γᵉ` 上做 penalty 約束
- 在 `Γˡ`、`Γᵇ` 上除了 penalty，還會額外組裝 `V` 和 `M`

因此它比 `patch_test.jl` 更接近「完整邊界物理條件」的寫法。

### `generateCircle.jl` 與 `circular.jl` 的對應

`generateCircle.jl` 會批量生成：

- `circular_tri3_3.msh`
- `circular_tri3_4.msh`
- ...

而 `circular.jl` 直接打開的是：

- `msh/circular.msh`

這裡不要誤會。repo 裏本身已經存了一份 `msh/circular.msh`，所以主腳本可以直接拿來跑；而批量生成腳本則負責產生一整系列不同密度的 mesh 供掃描或比較。

也就是說：

- `circular.msh` 是現成基準檔
- `circular_tri3_n.msh` 是批量生成的參數化版本

---

## 6. 這套代碼寫作方法為什麼是「標準寫法」

如果只看表面，你會覺得這幾個腳本只是普通 FEM 求解程式；但它們其實有一套很固定的研究寫法。

## 6.1 先定解析量，不是先定矩陣

它們全部先定義：

- `w`
- `φ₁`, `φ₂`
- 導數
- `M`
- `Q`
- `q`
- `m`

這樣做的意義是：

- 理論模型先清楚
- 弱式右端與邊界條件都能由解析量直接生成

這是非常標準的 manufactured solution / benchmark 驗證寫法。

## 6.2 mixed 系統一定拆 block

它們不把所有未知量混成一條不透明矩陣，而是固定拆成：

- `w-w`
- `φ-φ`
- `φ-w`

這樣的好處是：

- 看得出物理結構
- 方便逐塊除錯
- 方便研究不同 operator 對哪個子系統有影響

## 6.3 operator 一律寫成 `∫... => elements`

這是最典型的課題內風格。

寫法不是：

- 直接手寫 element stiffness
- 也不是包成一個總函數 `solveMindlinPlate(...)`

而是：

```julia
𝑎 = ∫... => elements
𝑎(K)
```

這種寫法的特點是：

- operator 名稱直接對應弱式
- 代碼與數學式關係很近
- 很適合研究型程式反覆替換 formulation

## 6.4 域內與邊界的 shape function 處理分開

這幾個腳本全部遵守：

- 域內：`set∇𝝭!`
- 邊界：`set𝝭!`
- 誤差：`set𝝭!`

這種分工非常清楚，也讓整個資料流更容易追。

## 6.5 不是直接保留解向量，而是回填節點欄位

```julia
push!(nodes, :d=>..., :d₁=>..., :d₂=>...)
```

這一步代表數值解之後的所有工作都重新回到節點資料結構。

好處是：

- `L₂` / `L₂φ` 可以直接從 nodes 取值
- VTK 或其他後處理也能用同一套資料
- 不需要每次再手動拆 DOF 向量

## 6.6 `TimerOutputs` 是標準研究腳本配件

這幾個主腳本都大量使用：

```julia
@timeit to "..."
```

這不是裝飾，而是研究型工作流的一部分。  
它代表這套代碼除了要算對，還要方便比較：

- 哪一段最耗時
- 域組裝和邊界組裝各花多少時間
- shape function 計算與 assemble 是否成為瓶頸

---

## 7. 你閱讀時最應該抓住的主線

如果把整個 repo 壓縮成一句話，它其實就是：

**先用 `BenchmarkExample` 產生或準備 `.msh`，再用 `ApproxOperator.GmshImport` 讀入，然後用 `ApproxOperator.MindlinPlate` 的 mixed operator 分塊組裝，最後回填節點並用 `L₂ / L₂φ` 驗證。**

更具體地說，主線是：

1. `generateCircle.jl` / `generatePatchtest.jl` 準備 mesh
2. `gmsh.open(...)` 開 mesh
3. `getPhysicalGroups()` 和 `get𝑿ᵢ()` 建資料結構
4. `getElements(...)` 抓域與邊界元素
5. `prescribe!` 餵材料、解析解、邊界資料
6. `set∇𝝭!` / `set𝝭!` 建立 shape function
7. `∫... => elements` 做 operator 組裝
8. block solve
9. `push!` 回寫節點
10. `L₂ / L₂φ` 驗證

這條主線就是你之後讀其他類似課題代碼時，最應該拿來比對的模板。

---

## 8. 新手推薦閱讀順序

如果你想最快讀懂這個 repo，我建議按這個順序：

1. `generatePatchtest.jl`  
   先看它怎樣批量生成方板 mesh，建立「mesh 不是 solver 做的」這個觀念。

2. `patch_test.jl`  
   這是最乾淨、最標準的 Mindlin plate mixed 測試模板。

3. `square.jl`  
   看同一套骨架如何換成更高階製造解與更實際的測試輸入。

4. `circular.jl`  
   最後再看圓板，因為它已經引入自然邊界、參考解切換、邊界分組差異，複雜度最高。

5. `BenchmarkExample` 的 `PatchTest` / `Circular`  
   如果你想回頭理解 `.msh` 是怎樣來的，再回去看 `D:\BenchmarkExample.jl` 中對應的生成函數。

---

## 9. 最後幫你總結成一句判斷

如果你問：「這個 repo 的標準寫法到底是什麼？」

最準確的回答是：

**它的標準寫法不是把所有流程包成黑箱，而是把理論場、mesh 匯入、domain operator、boundary operator、block matrix、節點回填、L2 驗證，全部拆成清楚且可逐段檢查的研究腳本流程。**

所以你之後寫類似課題代碼時，最應該模仿的不是某一條公式，而是這種固定順序：

**解析量定義 -> mesh 讀入 -> 分塊組裝 -> 邊界處理 -> block solve -> 回填節點 -> 誤差驗證。**
