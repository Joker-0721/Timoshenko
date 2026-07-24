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
const UNIFIED_CSV = joinpath(DATA_DIR, "buckling.csv")
mkpath(DATA_DIR)

# 全局幾何與材料物理常數
const E = 200e9
const ν = 0.3
const a = 1.0
const b = 1.0
const h = 1e-2
const Dᵇ = E * h^3 / (12 * (1 - ν^2))

const αʷ = 0.0
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

for n_div in 9:10
    @timeit to "RKPM Loop (ndiv=$n_div)" begin

        mesh_file = "msh/st_q_$(n_div).msh"

        type_w = :quad4
        type_φ = :quad4
        s_size = 1.0 / (n_div - 1)

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(mesh_file)

        # ======================================================================
        # (A) 載入網格並建立各場節點
        # ======================================================================
        # w 場節點
        # nodes_w = get𝑿ᵢ()
        # nʷ = length(nodes_w)
        # push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ),
        #                :s₂ => sʷ * s_size * ones(nʷ),
        #                :s₃ => sʷ * s_size * ones(nʷ))
        # sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

        # # φ 場節點
        # gmsh.clear()
        # gmsh.open(mesh_file)
        # nodes_φ = get𝑿ᵢ()
        # nᵠ = length(nodes_φ)
        # push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ),
        #                :s₂ => sᵠ * s_size * ones(nᵠ),
        #                :s₃ => sᵠ * s_size * ones(nᵠ))
        # sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

        gmsh.clear()
        gmsh.open(mesh_file)
        @timeit to "get nodes" nodes = get𝑿ᵢ()
        n = length(nodes)
        @timeit to "get entities" entities = getPhysicalGroups()
        elements_support = getElements(nodes, entities["Ω"], 1)
        # s_q, var_A = cal_area_support(elements_support)

        # 物理群組（只需載入一次）
        gmsh.clear()
        gmsh.open(mesh_file)
        entities = getPhysicalGroups()

        # ======================================================================
        # (B) 矩陣初始化
        # ======================================================================
        kʷʷ = zeros(n, n)
        kᵠᵠ = zeros(2n, 2n)
        kᵠʷ = zeros(2n, n)
        kᵍʷ = zeros(n, n)
        kᵍᵠ = zeros(2n, 2n)

        # ======================================================================
        # (C) 定義元素、賦予物理常數、計算形函數
        # ======================================================================
        # 域內元素（使用 RKPM 節點）
        elements_w = getElements(nodes, entities["Ω"], integrationOrder)
        elements_φ = getElements(nodes, entities["Ω"], integrationOrder)
        prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        set∇𝝭!(elements_w)
        set∇𝝭!(elements_φ)

        # 邊界懲罰用元素（使用 RKPM 節點）
        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_w_b = [getElements(nodes, entities[name], integrationOrder, normal=true) 
                        for name in boundary_names]
        elements_φ_b = [getElements(nodes, entities[name], integrationOrder, normal=true) 
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

        # ======================================================================
        # (D) 組裝材料剛度
        # ======================================================================
        𝑎ʷʷ = ∫wwdΩ => elements_w
        𝑎ᵠʷ = [∫φwdΩ => (elements_φ, elements_w)] # 注意：雙元素版本！
        𝑎ᵠᵠ = [∫φφdΩ => elements_φ, ∫κκdΩ => elements_φ]

        𝑎ʷʷ(kʷʷ)
        𝑎ᵠʷ(kᵠʷ)
        𝑎ᵠᵠ(kᵠᵠ)

        # ======================================================================
        # (E) 組裝幾何剛度
        # ======================================================================
        𝑎ᵍʷ = ∫∇wσ∇wdΩ => elements_w
        𝑎ᵍᵠ = ∫∇φσ∇φdΩ => elements_φ

        𝑎ᵍʷ(kᵍʷ)
        𝑎ᵍᵠ(kᵍᵠ)

        # ======================================================================
        # (G) 求解特徵值問題
        # ======================================================================
        K = [kʷʷ kᵠʷ'; kᵠʷ kᵠᵠ]
        Kᴳ = [kᵍʷ zeros(n, 2n); zeros(2n, n) kᵍᵠ]

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
    end
end
