#!/usr/bin/env julia
# =====================================================================
#  fem_tri3_loop_shear_circle.jl  —  圓板（四分之一圓盤＋對稱）收斂迴圈（自由振動）
#  由 src/vibration/skew/fem_tri3_loop_shear_skew.jl 自動衍生（tools/gen_shape_loops.py）
#  ⚠ 不要手改本檔：要改就改 skew 版或產生器，然後重新產生
#
#  與 skew 版的差異：
#    1. 網格 msh/patchtest_circ_tri3_$(R)_n$(ndiv).msh
#    2. 邊界群組 3 個（Γ¹..Γ³）；Γ⁴ 已移除
#    3. 四分之一圓盤＋兩對稱邊（Γ² y=0、Γ³ x=0）；邊界條件 CC／SS
#    4. 輸出路徑 / 檔名含形狀名（circle）
#  所有路徑由 @__DIR__ 推導 → 任何 cwd 都能跑
#
#  執行：julia src/vibration/circle/fem_tri3_loop_shear_circle.jl
#  收斂範圍：改 PLATE_NDIVS；R 清單：改 Rs
# =====================================================================

const _ROOT = dirname(dirname(dirname(@__DIR__)))   # .../mindlin_code
const _SRC  = joinpath(_ROOT, "src")
const _MSH  = joinpath(_ROOT, "msh")
const _DATE = joinpath(_ROOT, "date")
const _VTK  = joinpath(_ROOT, "vtk")

const Rs = [1.0]    # 圓板半徑（四分之一圓盤；R = 1 與方板 a = b = 1 同口徑）

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ


using TimerOutputs, LinearAlgebra, WriteVTK, DelimitedFiles
import Gmsh: gmsh


E = 1.0
ν = 0.3
ρ = 1.0
# t/a：**不硬寫**（2026-09-23）；斜板振動參考解 = Liew 1993 JSV 168(1):39-69（掃描件待 OCR），
#        該文表值按 t/b 分檔 → OCR 確認後把預設改成論文那一檔（1e-3 / 1e-1 / 1e-2）
h = parse(Float64, get(ENV, "PLATE_H", "1e-3"))
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
σ₁₁ = 1.0
σ₂₂ = 0.0
σ₁₂ = 0.0
a = 1.0

# ============================================================
#  邊界條件設定
#  w_edges    : 參與 w 懲罰（∫αwwdΓ）的邊號
#  assemble_w : 是否執行 w 邊界組裝
#  φ_edges    : 參與 φ 懲罰（∫αφφdΓ）的邊號
#  assemble_φ : 是否執行 φ 邊界組裝
# ============================================================
# ============================================================
#  邊界條件設定（四分之一圓盤＋兩對稱邊）
#    Γ¹ = 圓弧：實體邊界
#         CC 夾支 → w 與 φ 皆鎖（n₁₁ = n₂₂ = 1）
#         SS 簡支 → 只鎖 w；φ 不鎖 ⇒ 彎矩 M_n 自然為 0
#    Γ² = y = 0 對稱邊 → 鎖 φ_y（dof 2I）      ⇒ n₂₂ = 1
#    Γ³ = x = 0 對稱邊 → 鎖 φ_x（dof 2I−1）    ⇒ n₁₁ = 1
#  對稱邊只需 φ_n = 0（必要條件）；w 的 ∂w/∂n = 0 由變分自然滿足，不施加 w 懲罰。
#  dof 對應：2I−1 = φ_x（∂w/∂x）、2I = φ_y（∂w/∂y）——由 ∫φwdΩ 的 B₁/B₂ 項確認。
# ============================================================
const NA = [
    (1.0, 0.0, 1.0),   # Γ¹ 圓弧（夾支用；SS 不組裝）
    (0.0, 0.0, 1.0),   # Γ² y = 0：φ_y = 0
    (1.0, 0.0, 0.0),   # Γ³ x = 0：φ_x = 0
]
bcs = [
    (name="CC", w_edges=[1], assemble_w=true, φ_edges=[1,2,3], assemble_φ=true),
    (name="SS", w_edges=[1], assemble_w=true, φ_edges=[2,3],   assemble_φ=true),
]

to = TimerOutput()

for R in Rs

    # 為每個邊界條件開啟獨立 CSV 檔
    ios = Dict{String,IOStream}()
    for bc in bcs
        mkpath(joinpath(_DATE, "fem/vibration/circle"))
        ios[bc.name] = open(joinpath(_DATE, "fem/vibration/circle/vib_fem_R$(R)_$(bc.name)_h$(h)_shear_km_circle.csv"), "w")
        write(ios[bc.name], "ndiv,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀\n")
    end

ndivs = [parse(Int, s) for s in split(get(ENV, "PLATE_NDIVS", "5,10,20"), ",")]   # 掃描範圍可用 PLATE_NDIVS=5,10 覆蓋（煙霧測試用）
for ndiv in ndivs
integrationOrder = 4
integrationOrder_shear = 1
gmsh.initialize()
@timeit to "open msh file" gmsh.open(joinpath(_MSH, "patchtest_circ_tri3_R$(R)_n$ndiv.msh"))
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()


nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᴳʷʷ = zeros(nʷ,nʷ)
kᴳᵠᵠ = zeros(2*nᵠ,2*nᵠ)
mʷʷ = zeros(nʷ,nʷ)
mᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)


@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"],integrationOrder)
    @timeit to "get elements" elements_s = getElements(nodes, entities["Ω"],integrationOrder_shear)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    prescribe!(elements_s, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set∇𝝭!(elements_s)
    𝑎ʷʷ = ∫∇w∇wdΩ=>elements_s
    𝑎ᵠʷ = ∫φwdΩ=>elements_s
    𝑎ᵠᵠ = [
        ∫φφdΩ=>elements_s,
        ∫κκdΩ=>elements,
    ]
    𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements
    𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements
    𝑎ᵐʷʷ = ∫ρwwdΩ=>elements_s
    𝑎ᵐᵠᵠ = ∫ρφφdΩ=>elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ᴳʷʷ(kᴳʷʷ)
    @timeit to "assemble" 𝑎ᴳᵠᵠ(kᴳᵠᵠ)
    @timeit to "assemble" 𝑎ᵐʷʷ(mʷʷ)
    @timeit to "assemble" 𝑎ᵐᵠᵠ(mᵠᵠ)
end


@timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    prescribe!(elements_1, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>NA[1][1], :n₁₂=>NA[1][2], :n₂₂=>NA[1][3])
    prescribe!(elements_2, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>NA[2][1], :n₁₂=>NA[2][2], :n₂₂=>NA[2][3])
    prescribe!(elements_3, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>NA[3][1], :n₁₂=>NA[3][2], :n₂₂=>NA[3][3])
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    # 邊界組裝移至下方 BC 迴圈中執行
end

gmsh.finalize()

# ----------------------------------------------------------
#  儲存域積分後的矩陣快照，供各 BC 還原
#  （kᵠʷ, kᴳ*, m* 與邊界條件無關，不需還原）
# ----------------------------------------------------------
kʷʷ₀ = copy(kʷʷ)
kᵠᵠ₀ = copy(kᵠᵠ)

# 各邊元素陣列，便於依 BC 動態組合
elements_edges = [elements_1, elements_2, elements_3]

# VTK 格點與 cells（與 BC 無關，只算一次）
xs = [node.x for node in nodes]'
ys = [node.y for node in nodes]'
zs = [node.z for node in nodes]'
points = [xs; ys; zs]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]

# ==============================================================
#  邊界條件迴圈
# ==============================================================
for bc in bcs
    @timeit to "BC $(bc.name)" begin

        # 還原至域積分後狀態
        kʷʷ_bc = copy(kʷʷ₀)
        kᵠᵠ_bc = copy(kᵠᵠ₀)

        # ---------- w 邊界懲罰 ----------
        w_elm = reduce(∪, [elements_edges[i] for i in bc.w_edges])
        𝑎ʷ = ∫αwwdΓ => w_elm
        if bc.assemble_w
            @timeit to "assemble" 𝑎ʷ(kʷʷ_bc)
        end

        # ---------- φ 邊界懲罰 ----------
        φ_elm = reduce(∪, [elements_edges[i] for i in bc.φ_edges])
        𝑎ᵠ = ∫αφφdΓ => φ_elm
        if bc.assemble_φ
            @timeit to "assemble" 𝑎ᵠ(kᵠᵠ_bc)
        end

        # ---------- 組裝整體矩陣 ----------
        k_mat  = [kᵠᵠ_bc kᵠʷ; kᵠʷ' kʷʷ_bc]
        kᴳ_mat = [kᴳᵠᵠ zeros(2*nᵠ,nʷ); zeros(nʷ,2*nᵠ) kᴳʷʷ]
        m_mat  = [mᵠᵠ  zeros(2*nᵠ,nʷ); zeros(nʷ,2*nᵠ) mʷʷ]

        # ---------- 特徵值求解 ----------
        # λ, v = eigen(k_mat, kᴳ_mat)
        λ, v = eigen(k_mat, m_mat)
        # LAPACK 對一般矩陣不保證特徵值排序 → 顯式排序，確保 index 對到最低階模態
        p = sortperm(real.(λ))
        λ = real.(λ[p]); v = real.(v[:, p])   # ← 取實部：eigen 對稠密一般矩陣會回複數噪聲
        index = findfirst(real.(λ).>1e-8)

        # ---------- 位移場（前 12 模態） ----------
        d₁  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index])
        d₂  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+1])
        d₃  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+2])
        d₄  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+3])
        d₅  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+4])
        d₆  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+5])
        d₇  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+6])
        d₈  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+7])
        d₉  = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+8])
        d₁₀ = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+9])
        d₁₁ = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+10])
        d₁₂ = real.(v[2*nᵠ+1:2*nᵠ+nʷ, index+11])
        push!(nodes,
            :d₁=>d₁,   :d₂=>d₂,   :d₃=>d₃,   :d₄=>d₄,
            :d₅=>d₅,   :d₆=>d₆,   :d₇=>d₇,   :d₈=>d₈,
            :d₉=>d₉,   :d₁₀=>d₁₀, :d₁₁=>d₁₁, :d₁₂=>d₁₂,
        )

        # ---------- VTK 輸出 ----------
        mkpath(joinpath(_VTK, "circle/fem/R$(R)/$(bc.name)"))
        vtk_grid(joinpath(_VTK, "circle/fem/R$(R)/$(bc.name)/vib_fem_R$(R)_$(bc.name)_h$(h)_shear_km_circle.vtu"), points, cells;
                 ascii=false, append=false, compress=false) do vtk
            vtk["v₁"]  = [node.d₁  for node in nodes]
            vtk["v₂"]  = [node.d₂  for node in nodes]
            vtk["v₃"]  = [node.d₃  for node in nodes]
            vtk["v₄"]  = [node.d₄  for node in nodes]
            vtk["v₅"]  = [node.d₅  for node in nodes]
            vtk["v₆"]  = [node.d₆  for node in nodes]
            vtk["v₇"]  = [node.d₇  for node in nodes]
            vtk["v₈"]  = [node.d₈  for node in nodes]
            vtk["v₉"]  = [node.d₉  for node in nodes]
            vtk["v₁₀"] = [node.d₁₀ for node in nodes]
            vtk["v₁₁"] = [node.d₁₁ for node in nodes]
            vtk["v₁₂"] = [node.d₁₂ for node in nodes]
        end

        # ---------- 結果輸出 ----------
        # k = (λ[index]*a^2/(π^2*Dᵇ)*h)
            # println("$(bc.name): ", k)
            # write(ios[bc.name], "$n,$k\n")
            # flush(ios[bc.name])

            n_modes = 10
            # 由 index 開始取前 10 個非零特徵值，轉成無量綱頻率參數 k
            k_vals = [sqrt(real(λ[index+j-1])) * R^2 * sqrt(ρ*h/Dᵇ) for j in 1:n_modes]
            println("$(bc.name): ", k_vals)
            write(ios[bc.name], "$ndiv,$(join(k_vals, ","))\n")
            flush(ios[bc.name])
    end
end
end

# 關閉所有 CSV 檔
for bc in bcs
    close(ios[bc.name])
end

end  # for R in Rs
