# ApproxOperator.jl 新手導讀

這份筆記的目標不是把所有數學推導一次講完，而是先幫你建立一張「地圖」：

- 這個套件的檔案怎麼分工
- 資料在程式裡怎麼流動
- 幾個最常見的基礎函數各自代表什麼
- 寫新程式或看現有程式時，應該先看哪一層

如果你把這份文件先讀通，之後再去看 `ApproxOperator.jl` 的原始碼，會比較容易知道自己現在是在看「資料結構」、「shape function」、「網格匯入」，還是「物理算子」。

## 1. 一眼看懂整體結構

最重要的主線在 `src/`，可以先把它想成四層：

```text
ApproxOperator.jl
├─ node.jl / element.jl / operation.jl / shapefunction.jl
│  └─ 最核心的資料結構與組裝機制
├─ approximation/
│  └─ 各種元素或近似基底的 shape function 實作
├─ preprocession/
│  └─ 從 Gmsh 讀網格、產生元素、做型別轉換
└─ operation/
   └─ 各種物理問題的弱式積分算子與誤差函數
```

對新手來說，最值得先建立的流程是這條：

```text
Gmsh mesh
-> getPhysicalGroups()
-> get𝑿ᵢ()
-> getElements(...)
-> prescribe!(...)
-> set𝝭! / set∇𝝭! / set∇²𝝭!
-> 積分算子 => elements
-> 組裝 k, f
-> 求解
-> push!(nodes, :d => 解向量)
-> 誤差檢查，例如 L₂(...)
```

`test/runtests.jl` 就是在走這條主線，所以它是最值得先看的最小案例。

## 2. 這個資料夾裡的關係怎麼理解

### 主入口 `src/ApproxOperator.jl`

這個檔案本身不做複雜計算，它的主要角色是：

- 定義抽象型別：`AbstractElement`、`AbstractPiecewise`、`SpatialPartition`
- `include(...)` 其他檔案
- `export` 常用入口

可以把它理解成整個套件的總開關。

### 基礎層

| 檔案 | 作用 | 你可以怎麼理解 |
| --- | --- | --- |
| `src/node.jl` | 定義 `Node` | 所有資料最底層的容器 |
| `src/element.jl` | 定義 `Element` / `ReproducingKernel` / `prescribe!` | 把幾何節點和積分點打包成元素 |
| `src/operation.jl` | 定義 `Pair(form => elements)` 的呼叫方式 | 把局部積分算子組裝到全域矩陣/向量 |
| `src/shapefunction.jl` | 定義 `set𝝭!` 家族的入口 | 統一觸發 shape function 預計算 |

### approximation/

這一層的規律很明確：

- 一個檔案通常對應一種元素或一種近似方法
- 主要工作是定義 `set𝝭!(...)`、`set∇𝝭!(...)`、`set∇²𝝭!(...)`

例如：

- `approximation/tri3.jl`：三角形一次元素 `Tri3`
- `approximation/tri6.jl`：三角形二次元素 `Tri6`
- `approximation/seg2.jl`、`seg3.jl`：線元素
- `approximation/tet4.jl`、`tet10.jl`、`hex8.jl`：3D 元素
- `approximation/reproducingkernel.jl`：meshfree / reproducing kernel 類型
- `approximation/piecewise.jl`：piecewise polynomial / parametric 基底

你可以把這層理解成：「給定元素類型，告訴程式 shape function 應該怎麼算」。

### preprocession/

這一層負責把 Gmsh 或元素轉換接到主框架裡：

- `preprocession/importmsh.jl`
  - 定義子模組 `ApproxOperator.GmshImport`
  - 提供 `getPhysicalGroups()`、`get𝑿ᵢ()`、`getElements(...)`
- `preprocession/convert.jl`
  - 提供 `Tri3toTriHermite`、`Tri3toTriBell`、`Seg2toSegHermite` 等轉換函數

你可以把這層理解成：「把網格和元素資料整理成套件看得懂的形式」。

### operation/

這一層是物理問題本身，像是熱傳、彈性、板殼、Timoshenko 等。

目前主線中，很多檔案是**子模組**，例如：

- `ApproxOperator.Heat`
- `ApproxOperator.Elasticity`
- `ApproxOperator.Hyperelasticity`
- `ApproxOperator.Hamilton`
- `ApproxOperator.Timoshenko`
- `ApproxOperator.Stokes`
- `ApproxOperator.KirchhoffPlate`
- `ApproxOperator.MindlinPlate`
- `ApproxOperator.KirchhoffShell`

但要注意，不是每個 `operation/*.jl` 都長得一樣：

- `operation/curved_beam.jl` 是直接 include 到主模組，不是獨立子模組
- `operation/test.jl` 也是子模組

所以比較安全的理解方式是：

- `operation/` 主要放物理算子
- 其中「大多數」檔案是子模組
- 實際要怎麼 `import`，最好先看 `test/runtests.jl`

### 不是資料夾裡每個檔案都在主線啟用

這點很重要，因為新手很容易看到檔案就以為它一定會被載入。

依照 `src/ApproxOperator.jl` 目前的 `include(...)`：

- 已在主線載入的，是 `ApproxOperator.jl` 裡明確 `include` 的那些檔案
- `operation/phasefield.jl`、`operation/error_estimates.jl` 目前是被註解掉的
- `src/littletools.jl`、`approximation/grkgsi.jl`、`approximation/rkwave.jl`、`operation/plasticity.jl` 等檔案目前不在主入口主線裡

所以第一輪閱讀時，先以「主入口真的有 include 的檔案」為準，會最不容易混亂。

## 3. 核心資料結構

這個套件最核心的概念不是 class hierarchy，而是：

- `Node` 負責存資料
- `Element` 負責把一組幾何節點與一組積分點節點打包
- `operation` 負責把元素上的局部積分組裝成全域矩陣

### `Node`

`node.jl` 裡真正的底層型別是：

```julia
struct Node{T,N}
    index::NamedTuple{T,NTuple{N,Int}}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
```

這個設計的重點是：

- `Node` 不是直接把每個欄位存成 struct field
- 它是用 `index + data` 的方式，去對共享向量做「視圖式」存取

所以你看到：

```julia
xᵢ.x
ξ.𝝭
ξ[:∂𝝭∂x]
```

表面上像在讀一個欄位，實際上是在：

1. 先查這個 `Node` 的索引
2. 再去共享的 `data` 向量裡拿值

這也是為什麼這個專案裡常常看到 `push!`、`prescribe!` 在新增符號欄位，而不是去改 struct 定義。

### `𝑿ᵢ` 與 `𝑿ₛ`

| 名稱 | 定義位置 | 作用 |
| --- | --- | --- |
| `𝑿ᵢ` | `node.jl` | 幾何節點 / 網格節點 |
| `𝑿ₛ` | `node.jl` | 狀態節點 / 積分點節點 |

`𝑿ᵢ` 的 index 只有 `𝐼`，通常就是幾何節點編號。  
`𝑿ₛ` 的 index 是 `(:𝑔,:𝐺,:𝐶,:𝑠)`，主要用來管理積分點資料與 shape function 資料的共享向量。

如果你現在還看不懂 `𝑔 / 𝐺 / 𝐶 / 𝑠` 的全部差別，不用緊張。對新手來說，先知道它們是「內部索引用來共享資料」就夠了。

### `Element{T}`

`element.jl` 裡的基本元素型別是：

```julia
struct Element{T} <: AbstractElement
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end
```

可以直接把它理解成：

- `𝓒`：control/coordinate nodes，也就是這個元素的幾何節點
- `𝓖`：Gauss/integration nodes，也就是這個元素的積分點資料

這是整個專案最重要的資料配對。

很多後面的積分算子，都是在做這件事：

- 對每個積分點 `ξ in ap.𝓖`
- 讀出 `ξ[:𝝭]`、`ξ[:∂𝝭∂x]`、材料參數、邊界資料
- 再用 `ap.𝓒` 把局部量組裝到全域矩陣位置

### `ReproducingKernel{𝑝,𝑠,𝜙}`

這個型別的外型和 `Element` 很像，也同樣有：

- `𝓒`
- `𝓖`

差別在於：

- 普通 `Element{:Tri3}`、`Element{:Seg2}` 走的是標準元素 shape function
- `ReproducingKernel{...}` 走的是 meshfree / RK 型的近似方式
- 它會額外用到 `𝗠` 等 moment matrix 資料

所以可以把它理解成：「資料容器長得像元素，但 shape function 的計算方式不同」。

### `TRElement`

`TRElement` 來自 `approximation/CrouzeixRaviart.jl`，可以把它看成一個特殊元素族，用於 Crouzeix-Raviart 類型的三角形元素。

### `PiecewisePolynomial` / `PiecewiseParametric`

這兩個型別在 `approximation/piecewise.jl`。

你可以先粗略理解成：

- `PiecewisePolynomial`：局部分片多項式基底
- `PiecewiseParametric`：帶參數座標意義的局部分片基底

它們通常不是你第一眼要看的主線，但之後看混合格式、特殊基底或高階近似時會碰到。

### `RegularGrid`

`RegularGrid` 在 `approximation/meshfree.jl`，屬於 `SpatialPartition`。

它的作用不是做物理計算，而是做**鄰域搜尋**。對 meshfree / RK 類方法很重要，因為這類方法往往不是固定拿單個單元節點，而是要找一片鄰居點。

## 4. 幾個最常見基礎函數，分別在做什麼

### `getPhysicalGroups()`

來源：`ApproxOperator.GmshImport`

作用：

- 從 Gmsh 讀出 physical groups
- 回傳 `Dict{String, Pair{Int, Vector{Int}}}`

你可以把它理解成把 Gmsh 裡的 `"Ω"`、`"Γᵍ"` 這些物理區域名稱，轉成程式可以直接拿來抓元素的索引。

### `get𝑿ᵢ()`

來源：`ApproxOperator.GmshImport`

作用：

- 從 Gmsh 讀全部幾何節點
- 建立 `Vector{𝑿ᵢ}`

也就是把網格節點正式轉進 `ApproxOperator` 自己的資料格式。

### `getElements(...)`

來源：`ApproxOperator.GmshImport`

作用：

- 根據 physical group、元素型別、積分規則等資訊，產生元素向量
- 每個元素裡都會有 `𝓒` 和 `𝓖`

最常見的用法是：

```julia
elements = getElements(nodes, entities["Ω"])
```

這表示：

- 用 `nodes` 當幾何節點
- 從 `"Ω"` 那塊區域抓元素
- 用預設積分規則建立元素

它也有其他版本，可以額外指定：

- 積分階數
- 特定元素型別
- 手動積分點
- `SpatialPartition`，供 meshfree/RK 類型用

### `prescribe!(...)`

來源：主模組 `ApproxOperator`

作用：

- 把材料參數、邊界資料、精確解、罰參數等，寫入元素或節點資料欄位

最常見兩種寫法：

```julia
prescribe!(elements, :k => 1.0)
prescribe!(elements, :g => 𝑢)
```

它們分別表示：

- `:k => 1.0`：指定常數材料參數
- `:g => 𝑢`：指定函數型邊界條件或解析解

`prescribe!(ξ::Node, :g => f)` 會試著用不同參數形式呼叫 `f`，例如：

- `f(x, y, z)`
- `f(x, y, z, n₁)`
- `f(x, y, z, n₁, n₂)`

所以邊界函數如果需要法向資訊，也可以放進去。

大部分新手情況下，先用預設 `index=:𝐺` 即可，不必一開始就深究所有索引層級。

### `set𝝭!` / `set∇𝝭!` / `set∇²𝝭!` / `set∇̂³𝝭!`

來源：主模組 `ApproxOperator`，實作散在 `approximation/*.jl`

這組函數的功能，是先在元素的積分點上把 shape function 資料算好。

| 函數 | 代表的意思 |
| --- | --- |
| `set𝝭!` | 只算 shape function 值 |
| `set∇𝝭!` | 算 shape function 一階導數 |
| `set∇²𝝭!` | 算二階導數 |
| `set∇̂³𝝭!` | 算到三階相關資料 |

例如熱傳剛度矩陣常會用到 `∂𝝭∂x`、`∂𝝭∂y`，所以會先跑：

```julia
set∇𝝭!(elements)
```

而邊界上的 Nitsche 或罰項，如果只需要 `𝝭`，就可能只要：

```julia
set𝝭!(elements)
```

### `getDOFs(elements)`

來源：主模組 `ApproxOperator`

作用：

- 從元素集合裡把所有幾何節點編號 `𝐼` 收集出來

你可以把它理解成：「這批元素實際碰到了哪些自由度」。

### `Tri3toTriHermite` / `Tri3toTriBell` / `Seg2toSegHermite`

來源：`preprocession/convert.jl`

作用：

- 將一種元素資料轉成另一種元素資料
- 常用於高階元素、Hermite/Bell 類型的構造

可以把它理解成「幫你把幾何與積分資料重新包裝成另一種元素族」。

其中：

- `Tri3toTriHermite(...)` 會回傳元素、額外節點、edge 資訊
- `Tri3toTriBell(...)` 會回傳元素與額外節點
- `Seg2toSegHermite(...)` 會回傳新的 Hermite 線元素

## 5. 這份程式碼通常是怎麼寫的

這個專案的風格很有辨識度，看懂下面幾點之後，很多檔案會突然變得好讀很多。

### 1. 大量使用多重派發

例如：

```julia
set𝝭!(::Element{:Tri3}, x::Node)
set𝝭!(::Element{:Tri6}, x::Node)
set𝝭!(ap::ReproducingKernel, x::Node)
```

意思是：

- 同一個動作叫 `set𝝭!`
- 但不同元素型別，走不同實作

所以當你想知道某個元素的 shape function 怎麼算，不要先全域找一個大函數，而是直接去對應的 `approximation/*.jl` 找那個型別的方法。

### 2. 用 `Symbol` 當資料欄位名

這個專案不是把欄位全部寫死在 struct 裡，而是常常這樣做：

```julia
prescribe!(elements, :E => 210e9)
prescribe!(elements, :ν => 0.3)
prescribe!(elements, :u => exact_solution)
```

也就是：

- `:E`、`:ν`、`:u` 都只是 `Symbol`
- 真正的值是放在共享 `data` 裡

這樣的好處是靈活，但新手一開始會覺得「欄位從哪裡來」很不直觀。你可以把它理解成：這個專案是把欄位定義推遲到使用時才建立。

### 3. 很多看起來像 element 欄位的東西，其實是存在積分點資料裡

在 `element.jl` 中，`AbstractElement` 的 `getproperty` 會把未知欄位導向 `ap.𝓖[1]`。

所以：

```julia
ap.E
ap.𝐽
ap.n₁
```

很多時候其實不是 element 自己有這個欄位，而是去 element 第一個積分點的共享資料裡拿。

這是本專案閱讀上的一個關鍵觀念。

### 4. 用 `Pair(form => elements)` 組裝矩陣

`src/operation.jl` 把 `Pair` 做成可呼叫物件，所以你會看到這種寫法：

```julia
𝑎 = ∫∫∇v∇udxdy => elements
𝑎(k)

𝑓 = ∫vgdΓ => boundary_elements
𝑓(k, f)
```

它的意思非常直接：

- 左邊是局部積分算子
- 右邊是要套用的元素集合
- 呼叫時就逐元素組裝進 `k` 和 `f`

如果是混合格式，也可能是：

```julia
op = 某個雙場算子 => (elements_a, elements_b)
```

### 5. `operation/*.jl` 裡的函數名通常就是弱式積分的數學記號

例如：

- `∫∫∇v∇udxdy`
- `∫vgdΓ`
- `L₂`

這樣的命名很貼近公式，但對新手不一定友善。你可以先把它們讀成：

- 剛度項
- 邊界項
- 載重項
- 誤差範數

先抓角色，再回頭看公式細節。

## 6. 最小使用流程：以 `test/runtests.jl` 為例

目前最推薦的最小路徑，就是直接照著測試檔理解。

### 典型 import 寫法

```julia
using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy
import Gmsh: gmsh
```

這裡有一個很重要的實務提醒：

- 雖然主模組有 `export` 某些名稱
- 但依目前專案自己的測試檔，最穩妥的用法仍然是**直接從子模組匯入**

也就是說，新手先照 `test/runtests.jl` 的 import 方式走，最不容易出錯。

### Step 1. 讀入網格

```julia
gmsh.initialize()
gmsh.open("patchtest.msh")
```

### Step 2. 拿到 physical groups 與幾何節點

```julia
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
```

此時：

- `entities["Ω"]` 可能代表域內元素
- `entities["Γᵍ"]` 可能代表某段邊界
- `nodes` 是 `Vector{𝑿ᵢ}`

### Step 3. 建立全域矩陣與向量

```julia
nₚ = length(nodes)
k = zeros(nₚ, nₚ)
f = zeros(nₚ)
```

### Step 4. 建立域內元素，指定材料參數，預算 shape function

```julia
elements = getElements(nodes, entities["Ω"])
prescribe!(elements, :k => 1.0)
set∇𝝭!(elements)
```

這三行分別代表：

- 從 `"Ω"` 建立元素
- 指定熱傳導係數 `k`
- 預算一階導數，因為後面的熱傳剛度項會用到 `∂𝝭∂x`、`∂𝝭∂y`

### Step 5. 組裝域內剛度

```julia
𝑎 = ∫∫∇v∇udxdy => elements
𝑎(k)
```

### Step 6. 建立邊界元素，指定邊界條件，組裝邊界項

```julia
elements = getElements(nodes, entities["Γᵍ"])
prescribe!(elements, :g => 𝑢)
prescribe!(elements, :α => 1e7)
set𝝭!(elements)

𝑓Γ = ∫vgdΓ => elements
𝑓Γ(k, f)
```

這裡 `:g` 是邊界資料，`:α` 是罰參數。

### Step 7. 求解並把結果回填到節點

```julia
d = k \ f
push!(nodes, :d => d)
```

這步很值得記住：

- `d` 是解向量
- `push!(nodes, :d => d)` 之後，你就可以透過節點資料去存取離散解

### Step 8. 指定精確解並做誤差檢查

```julia
elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u => 𝑢)
set∇𝝭!(elements)
L₂error = L₂(elements)
```

這裡的思路是：

- 重新用較高積分階數建立元素
- 在積分點上放入精確解 `:u`
- 比較數值解與精確解，計算 `L₂` 誤差

## 7. 常見符號表

下面這張表是給新手看的「先能讀懂」版本，不是完整數學定義。

| 符號 | 直白意思 |
| --- | --- |
| `𝝭` | shape function 的值 |
| `∂𝝭∂x`, `∂𝝭∂y`, `∂𝝭∂z` | shape function 對座標的一階導數 |
| `∂²𝝭∂x²` 等 | 二階導數 |
| `∂³𝝭...` | 三階導數 |
| `𝑿ᵢ` | 幾何節點 |
| `𝑿ₛ` | 積分點 / 狀態節點 |
| `𝓒` | 元素的幾何節點集合 |
| `𝓖` | 元素的積分點集合 |
| `𝐼` | 幾何節點編號 |
| `𝐶` | 元素編號 |
| `𝐺` | 積分點資料的全域索引層 |
| `𝑔` | 積分點相關索引層 |
| `𝑠` | shape function 資料在共享向量中的偏移位置 |
| `ξ`, `η`, `ζ` | 參考元素座標 |
| `𝑤` | 積分權重乘上 Jacobian 後的權重資料 |
| `𝐽` | Jacobian 或幾何映射相關量 |
| `n₁`, `n₂`, `n₃` | 法向相關分量 |
| `𝗠` | reproducing kernel 用的 moment matrix |
| `L₂` | 常見的誤差範數 |

如果你剛開始看原始碼，只要先記住三件事就很夠用：

1. `𝝭` 相關通常是 shape function
2. `𝓒 / 𝓖` 通常是在分「幾何節點」和「積分點」
3. `∫...` 名稱通常是在寫弱式積分項

## 8. 建議閱讀順序

如果你接下來要自己去讀原始碼，我建議照這個順序：

1. `src/ApproxOperator.jl`
2. `src/node.jl`
3. `src/element.jl`
4. `src/operation.jl`
5. `src/shapefunction.jl`
6. `src/preprocession/importmsh.jl`
7. `src/approximation/tri3.jl`
8. `src/operation/heat.jl`
9. `test/runtests.jl`

這個順序的好處是：

- 先知道資料長什麼樣
- 再知道 shape function 怎麼掛進去
- 最後再看一個完整案例

## 9. 先記住的重點

如果你現在只想先抓主幹，先記下面這幾句就夠了：

- `Node` 是最底層資料容器，但它是靠 `index + data` 在共享向量裡取值
- `Element` = 幾何節點 `𝓒` + 積分點 `𝓖`
- `approximation/` 決定 shape function 怎麼算
- `preprocession/` 決定網格資料怎麼轉成元素
- `operation/` 決定弱式積分怎麼組裝
- `prescribe!` 是塞資料
- `set𝝭!` 家族是預算 shape function
- `算子 => elements` 是組裝的核心寫法
- `test/runtests.jl` 是目前最值得跟著走的最小路徑

等你接下來想往下鑽時，最自然的下一步通常會是下面三種之一：

- 我想看某個元素的 shape function 怎麼寫
- 我想看某個物理問題的剛度矩陣怎麼組裝
- 我想知道某個符號或欄位到底存在哪裡

到那時候，再針對你卡住的那一個點往下展開，就會比一開始直接硬啃整個資料夾輕鬆很多。
