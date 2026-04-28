# Timoshenko 1D/2D 標準寫法導讀

這組腳本是照 `mindlin_plate` 的研究腳本習慣搬過來的，目標不是做成 package，也不是 sweep framework，而是保留課題內常用的寫法節奏：

1. 先生成 `.msh`
2. `gmsh.open(...)`
3. `getPhysicalGroups()` / `get𝑿ᵢ()` / `getElements(...)`
4. `prescribe!`
5. `set∇𝝭!` 或 `set𝝭!`
6. `∫... => elements`
7. 組成 block matrix
8. `push!(nodes, ...)`
9. `L₂` / `L₂φ`

## 檔案對照

- `mindlin_plate/generateCircle.jl` 的角色，在這裡對應成 `src/generateBeam.jl`
- `mindlin_plate/generatePatchtest.jl` 的角色，在這裡對應成 `src/generatePatchtest.jl`
- `mindlin_plate/patch_test.jl` 的角色，在這裡對應成：
  - `src/beam.jl`：1D Timoshenko 梁 mixed block 主線
  - `src/patch_test.jl`：2D Mindlin 板 patch test
- `mindlin_plate/square.jl` 的角色，在這裡對應成 `src/square.jl`

## 為什麼 1D 要多補兩個地方

`mindlin_plate` 原本只有 2D 板，所以這次多補了兩個基礎件：

- `D:\BenchmarkExample.jl\src\TimoshenkoBeam.jl`
  - 新增 `BenchmarkExample.TimoshenkoBeam.generateMsh(...)`
  - 用來生成 1D 梁的 `Ω / Γ¹ / Γ²` mesh
- `D:\ApproxOperator.jl\src\operation\timoshenko.jl`
  - 保留原本的 combined-style operator
  - 另外新增 block-style operator：
    - `∫κκdΩ`
    - `∫wwdΩ`
    - `∫φφdΩ`
    - `∫φwdΩ`
    - `∫wqdΩ`
    - `∫φmdΩ`
    - `∫αwwdΓ`
    - `∫αφφdΓ`
    - `∫wVdΓ`
    - `∫φMdΓ`

這樣 `beam.jl` 就可以跟 `mindlin_plate/patch_test.jl` 一樣，用

```julia
kʷʷ, kᵠᵠ, kᵠʷ, fʷ, fᵠ
```

去顯式組裝 block system，而不是直接拼成單一 `2n × 2n` 矩陣。

## 1D 腳本主線

`src/beam.jl` 的流程是：

1. 定義材料與截面常數 `E, ν, κ, h, A, I, EI, kGA`
2. 定義解析量 `w, φ, q, m, M, V, g, g₁`
3. `gmsh.open("msh/beam.msh")`
4. `entities = getPhysicalGroups(); nodes = get𝑿ᵢ()`
5. 建立 `kʷʷ, kᵠᵠ, kᵠʷ, fʷ, fᵠ`
6. 對 `Ω` 做 domain 組裝：
   - `∫wwdΩ`
   - `∫φwdΩ`
   - `∫φφdΩ`
   - `∫κκdΩ`
   - `∫wqdΩ`
   - `∫φmdΩ`
7. 對 `Γ¹ ∪ Γ²` 做邊界組裝：
   - `∫αwwdΓ`
   - `∫wVdΓ`
   - `∫φMdΓ`
8. 解
   - `[kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]`
9. `push!(nodes, :d=>..., :φ=>...)`
10. 用 `L₂` / `L₂φ` 比對 `w_exact_ss` / `φ_exact_ss`

這裡的 1D unknown 仍然是位移 `w` 和轉角 `φ`，只是寫法上完全對齊 `mindlin_plate` 的 mixed block 習慣。

## 2D 腳本主線

`src/patch_test.jl` 和 `src/square.jl` 基本上就是照 `mindlin_plate` 原有風格整理：

- `patch_test.jl`
  - 讀 `msh/patchtest.msh`
  - 用製造解做最標準的 patch test
  - 四條邊 `Γ¹ ~ Γ⁴` 全部走 penalty
- `square.jl`
  - 預設讀 `msh/patchtest_tri3_16.msh`
  - 用較高階製造解測試方板
  - 保留註解式的 VTK 後處理風格

2D 仍然用 `ApproxOperator.MindlinPlate`，因為這就是 2D 板對應的 Timoshenko / Mindlin 類型模型。

## Mesh 生成入口

- `src/generateBeam.jl`
  - 生成 `msh/beam.msh`
  - 以及 `beam_seg2_n.msh`、`beam_seg3_n.msh`
- `src/generatePatchtest.jl`
  - 生成 `msh/patchtest.msh`
  - 以及 `patchtest_tri3_n.msh`、`patchtest_tri6_n.msh`、`patchtest_quad4_n.msh`、`patchtest_quad8_n.msh`

## 閱讀順序

建議新手按這個順序看：

1. `src/generateBeam.jl`
2. `src/beam.jl`
3. `src/generatePatchtest.jl`
4. `src/patch_test.jl`
5. `src/square.jl`

如果想對照課題組模板，再回去對照 `mindlin_plate` 的：

1. `generatePatchtest.jl`
2. `patch_test.jl`
3. `square.jl`

這樣會最容易看出兩邊其實只是物理量不同，但編寫節奏是一樣的。
