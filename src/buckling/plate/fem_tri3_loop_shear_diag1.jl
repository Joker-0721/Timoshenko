using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ


using TimerOutputs, LinearAlgebra, WriteVTK, DelimitedFiles
import Gmsh: gmsh


E = 1.0
ν = 0.3
ρ = 1.0
h = 1e-3
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
bcs = [
    (name="CCCC", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[1,2,3,4], assemble_φ=true),
    (name="SSSS", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[1,2,3,4], assemble_φ=false),
    (name="CSCS", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[1,3],     assemble_φ=true),
    (name="SCSC", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[2,4],     assemble_φ=true),
    (name="FSCS", w_edges=[2,3,4],   assemble_w=true,  φ_edges=[3],       assemble_φ=true),
    (name="FSSS", w_edges=[2,3,4],   assemble_w=true,  φ_edges=[3],       assemble_φ=false),
]

to = TimerOutput()

# 為每個邊界條件開啟獨立 CSV 檔
ios = Dict{String,IOStream}()
for bc in bcs
    mkpath("../diag/shear_diag1")
    ios[bc.name] = open("../diag/shear_diag1/vib_fem_$(bc.name)_shear_km_0.001.csv", "w")
    write(ios[bc.name], "ndiv,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀\n")
end

ndivs = 20:30
for ndiv in ndivs
integrationOrder = 2
integrationOrder_shear = 1
gmsh.initialize()
@timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv.msh")
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
    𝑎ᵐʷʷ = ∫ρwwdΩ=>elements  # B.3组1：质量保持2阶
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
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(elements_1, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
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
elements_edges = [elements_1, elements_2, elements_3, elements_4]

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
        mkpath("../vtk/fem/vibration/$(bc.name)")
        vtk_grid("../vtk/fem/vibration/$(bc.name)/vib_fem_$(bc.name)_shear_km_0.001.vtu", points, cells;
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
            k_vals = [sqrt(real(λ[index+j-1])) * a^2 * sqrt(ρ*h/Dᵇ) for j in 1:n_modes]
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
