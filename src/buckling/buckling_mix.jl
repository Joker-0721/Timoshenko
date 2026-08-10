const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫κκdΩ, ∫∇φσ∇φdΩ, ∫∇wσ∇wdΩ, ∫αwwdΓ, ∫αφφdΓ

using LinearAlgebra
using TimerOutputs, Printf
import Gmsh: gmsh

# 全局數據庫路徑規範
const DATA_DIR = "./data"
const UNIFIED_CSV = joinpath(DATA_DIR, "buckling_RKPM.csv")
mkpath(DATA_DIR)

# 全局幾何與材料物理常數
const E = 200e9
const ν = 0.3
const a = 1.0
const b = 1.0
const h = 1e-2
const Dᵇ = E * h^3 / (12 * (1 - ν^2))

const αʷ = 1.0e8 * Dᵇ
const αᵠ = 0.0
const σ₁₁ = 1.0
const σ₂₂ = 0.0
const σ₁₂ = 0.0

# 用公式計算解析解（不限制階數）
function compute_exact_modes(n_modes::Int)
    modes = Float64[]
    for m in 1:50
        for n in 1:50
            k = (m + n^2/m)^2
            push!(modes, k)
        end
    end
    sort!(modes)
    return modes[1:n_modes]
end

const integrationOrder = 2
const sʷ = 1.5
const sᵠ = 1.5
const to = TimerOutput()

const type_w = Element{:Quad}
const type_φ = Element{:Quad}
const n_div = 9
const s_size = 1.0 / (n_div - 1)

        mesh_file = "msh/st_q_$(n_div).msh"
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(mesh_file)

        # ======================================================================
        # (A) 載入網格並建立各場節點
        # ======================================================================
        # w 場節點
        nodes_w = get𝑿ᵢ()
        nʷ = length(nodes_w)
        push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ),
                       :s₂ => sʷ * s_size * ones(nʷ),
                       :s₃ => sʷ * s_size * ones(nʷ))
        sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

        # φ 場節點
        nodes_φ = get𝑿ᵢ()
        nᵠ = length(nodes_φ)
        push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ),
                       :s₂ => sᵠ * s_size * ones(nᵠ),
                       :s₃ => sᵠ * s_size * ones(nᵠ))
        sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

        # 物理群組（只需載入一次）
        entities = getPhysicalGroups()

        # ======================================================================
        # (B) 矩陣初始化
        # ======================================================================
        kʷʷ = zeros(nʷ, nʷ)
        kᵠᵠ = zeros(2nᵠ, 2nᵠ)
        kᵠʷ = zeros(2nᵠ, nʷ)
        kᵍʷ = zeros(nʷ, nʷ)
        kᵍᵠ = zeros(2nᵠ, 2nᵠ)

        # ======================================================================
        # (C) 定義元素、賦予物理常數、計算形函數
        # ======================================================================
        # 域內元素（使用 RKPM 節點）
        elements_w = getElements(nodes_w, entities["Ω"], type_w, integrationOrder, sp_w)
        elements_φ = getElements(nodes_φ, entities["Ω"], type_φ, integrationOrder, sp_φ)
        prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        set∇𝝭!(elements_w)
        set∇𝝭!(elements_φ)

        # 邊界懲罰用元素（使用 RKPM 節點）
        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_w_b = [getElements(nodes_w, entities[name], type_w, integrationOrder, sp_w, normal=true) 
                        for name in boundary_names]
        elements_φ_b = [getElements(nodes_φ, entities[name], type_φ, integrationOrder, sp_φ, normal=true) 
                        for name in boundary_names]
        for el in elements_w_b
            prescribe!(el, :α => αʷ, :g => (x, y, z) -> 0.0)
            set𝝭!(el)
        end
        for el in elements_φ_b
            prescribe!(el, :α => αᵠ, :g₁ => (x, y, z) -> 0.0, :g₂ => (x, y, z) -> 0.0,
                           :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
            set𝝭!(el)
        end
        
        @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
        @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
        @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
        @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
        @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
        prescribe!(elements_w_1, :α=>αʷ, :g=>w)
        prescribe!(elements_w_2, :α=>αʷ, :g=>w)
        prescribe!(elements_w_3, :α=>αʷ, :g=>w)
        prescribe!(elements_w_4, :α=>αʷ, :g=>w)
        @timeit to "calculate shape functions" set𝝭!(elements_q_1)
        @timeit to "calculate shape functions" set𝝭!(elements_q_2)
        @timeit to "calculate shape functions" set𝝭!(elements_q_3)
        @timeit to "calculate shape functions" set𝝭!(elements_q_4)
        @timeit to "calculate shape functions" set𝝭!(elements_w_1)
        @timeit to "calculate shape functions" set𝝭!(elements_w_2)
        @timeit to "calculate shape functions" set𝝭!(elements_w_3)
        @timeit to "calculate shape functions" set𝝭!(elements_w_4)

        # ======================================================================
        # (D) 組裝材料剛度
        # ======================================================================
        kʷʷ = ∫wwdΩ => elements_w
        kᵠʷ = ∫φwdΩ => (elements_φ, elements_w) # 注意：雙元素版本！
        kᵠᵠ = [∫φφdΩ => elements_φ, ∫κκdΩ => elements_φ]

        # ======================================================================
        # (E) 組裝幾何剛度
        # ======================================================================
        kᵍʷ = ∫∇wσ∇wdΩ => elements_w
        kᵍᵠ = ∫∇φσ∇φdΩ => elements_φ

        # ======================================================================
        # (F) 組裝邊界懲罰
        # ======================================================================
        for el in elements_w_b
            (∫αwwdΓ => el)(kʷʷ)
        end
        for el in elements_φ_b
            (∫αφφdΓ => el)(kᵠᵠ)
        end

        # ======================================================================
        # (G) 求解特徵值問題
        # ======================================================================
        K = [kʷʷ kᵠʷ'; kᵠʷ kᵠᵠ]
        Kᴳ = [kᵍʷ zeros(nʷ, 2nᵠ); zeros(2nᵠ, nʷ) kᵍᵠ]

        F = eigen(K, Kᴳ)
        λ = F.values

        # ======================================================================
        # (H) 輸出到 CSV
        # ======================================================================
        λ_sorted = sort(λ, by=real)
        exact_modes = compute_exact_modes(length(λ_sorted))

        file_exists = isfile(UNIFIED_CSV)
        open(UNIFIED_CSV, "a") do io
            if !file_exists
                println(io, "ndiv,lambda_numeric,lambda_exact")
            end
            for i in 1:length(λ_sorted)
                @printf(io, "%d,%.16e,%.6f\n", n_div, real(λ_sorted[i]), exact_modes[i])
            end
        end

        println("ndiv = $n_div, 共輸出 $(length(λ)) 個特徵值")

        gmsh.finalize()



# 图--------------------------------------------------------------------------------

# 坐标------------------------------------------------------------------------------
# nₚ = length(nodes)
# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d*4
#     # points[3,i] = us[i]*4
# end

nₚ = length(nodes)
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = node.d/15
    # points[3,i] = us[i]*4
end

# 二维------------------------------------------------------------------------------
# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
# vtk_grid("./vtk/hmd_2d/Tri3_d_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
#     # vtk["精确解"] = us
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# 三维------------------------------------------------------------------------------

# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/circular_tri3_"*string(ndiv), points, cells) do vtk
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]


# 三维误差-------------------------------------------------------------------------------------------------
# vtk_grid("./vtk/circular_tri3_" * string(ndiv) * ".vtu", points, cells;
#          ascii=true, append=false, compress=false) do vtk


#     vtk["L2_w", WriteVTK.VTKFieldData()] = [L₂_w]
#     vtk["L2_phi", WriteVTK.VTKFieldData()] = [L₂_φ]
# end
# -------------------------------------------------------------------------------------------------

# 三维变形------------------------------------------------------------------------------
# vtk_grid("./vtk/circular_Clamped_tri3_" * string(ndiv) * ".vtu", points, cells;
#          ascii=true, append=false, compress=false) do vtk
vtk_grid("./vtk/circular_Clamped_tri3_$n.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    # 挠度 w
    vtk["w"] = [node.d for node in nodes]
    # 转角 φ₁
    vtk["phi_1"] = [node.d₁ for node in nodes]
    # 转角 φ₂  
    vtk["phi_2"] = [node.d₂ for node in nodes]
end
# -------------------------------------------------------------------------------------------------
