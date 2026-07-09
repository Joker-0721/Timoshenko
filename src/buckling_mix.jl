const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ

using LinearAlgebra
using TimerOutputs, WriteVTK, Printf, Statistics
import Gmsh: gmsh

# 全局數據庫路徑規範
const DATA_DIR = "./data"
const UNIFIED_CSV = joinpath(DATA_DIR, "buckling_mix.csv")
mkpath(DATA_DIR)

# 全局幾何與材料物理常數
const E = 200e9
const ν = 0.3
const a = 1.0
const b = 1.0
const h = 1e-2
const Dᵇ = E * h^3 / (12 * (1 - ν^2))
const Dˢ = 5 / 6 * E * h / (2 * (1 + ν))

const αʷ = 1.0e8 * Dᵇ
const αᵠ = 1.0e8 * Dᵇ
const σ₁₁ = 1.0
const σ₂₂ = 0.0
const σ₁₂ = 0.0

const k_exact_modes = [
    (1 + 1^2/1)^2,  # Mode 1: (m=1, n=1)
    (2 + 1^2/2)^2,  # Mode 2: (m=2, n=1)
    (3 + 1^2/3)^2,  # Mode 3: (m=3, n=1)
    (2 + 2^2/2)^2,  # Mode 4: (m=2, n=2)
    (4 + 1^2/4)^2,  # Mode 5: (m=4, n=1)
    (3 + 2^2/3)^2   # Mode 6: (m=3, n=2)
]

const integrationOrder = 2
const sʷ = 1.5
const sᵠ = 1.5
const to = TimerOutput()


for n_div in 9:25
    @timeit to "Mix Loop (ndiv=$n_div)" begin

        mesh_file = "../msh/st_q_$(n_div).msh"

        type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
        type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
        type_M = :(PiecewisePolynomial{:Linear2D})
        s_size = 1.0 / (n_div - 1)

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)

        # ======================================================================
        # (A) 載入網格並建立各場節點
        # ======================================================================
        gmsh.open(mesh_file)
        nodes_w = get𝑿ᵢ()
        nʷ = length(nodes_w)
        push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ),
                       :s₂ => sʷ * s_size * ones(nʷ),
                       :s₃ => sʷ * s_size * ones(nʷ))
        sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

        gmsh.clear()
        gmsh.open(mesh_file)
        nodes_φ = get𝑿ᵢ()
        nᵠ = length(nodes_φ)
        push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ),
                       :s₂ => sᵠ * s_size * ones(nᵠ),
                       :s₃ => sᵠ * s_size * ones(nᵠ))
        sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

        gmsh.clear()
        gmsh.open(mesh_file)
        nodes = get𝑿ᵢ()
        entities = getPhysicalGroups()
        nˢ = length(nodes)

        # ======================================================================
        # (B) 矩陣初始化
        # ======================================================================
        nˢ = length(nodes)
        nₑ = length(elements_q)
        nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))
        nᵐ = nₑ * nᵖ

        # 位移子塊
        kʷʷ = zeros(nʷ, nʷ)
        kᵠᵠ = zeros(2nᵠ, 2nᵠ)
        kᵠʷ = zeros(2nᵠ, nʷ)   # 注意：混合形式中此項可能為零，取決於推導
        # 剪應力子塊
        kˢˢ = zeros(2nˢ, 2nˢ)
        kˢᵠ = zeros(2nˢ, 2nᵠ)
        kˢʷ = zeros(2nˢ, nʷ)
        # 彎矩子塊
        kᵐᵐ = zeros(3nᵐ, 3nᵐ)
        kᵐᵠ = zeros(3nᵐ, 2nᵠ)
        kᵐʷ = zeros(3nᵐ, nʷ)   # 可能為零
        kˢᵐ = zeros(2nˢ, 3nᵐ)  # 可能為零
        # 懲罰子塊（用於強制邊界條件）
        kᵅʷʷ = zeros(nʷ, nʷ)
        kᵅᵠᵠ = zeros(2nᵠ, 2nᵠ)

        # ======================================================================
        # (C) 定義元素、賦予物理常數、計算形函數
        # ======================================================================
        # ---- C1. 定義所有元素物件 ----
        # 剪切場 (Q)
        elements_q   = getElements(nodes, entities["Ω"], integrationOrder)
        elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)

        # w 場 (橫向位移)
        elements_w   = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
        elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)

        # 彎矩場 (M)
        elements_m   = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
        elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)

        # φ 場 (截面轉角)
        elements_φ   = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
        elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)

        # 邊界懲罰用元素 (w 與 φ 各自四個邊界)
        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_w_b = [getElements(nodes_w, entities[name], eval(type_w),integrationOrder, sp_w, normal=true) for name in boundary_names]
        elements_φ_b = [getElements(nodes_φ, entities[name], eval(type_φ),integrationOrder, sp_φ, normal=true) for name in boundary_names]

        # ---- C2. 統一賦予物理常數 ----
        # 區域元素：材料參數與應力 (w 場需要幾何剛度應力，提前給定)
        prescribe!(elements_q, :E => E, :ν => ν, :h => h)
        prescribe!(elements_w, :E => E, :ν => ν, :h => h, :σ₁₁ => σ₁₁, :σ₂₂ => σ₂₂, :σ₁₂ => σ₁₂)
        prescribe!(elements_m, :E => E, :ν => ν, :h => h)
        prescribe!(elements_φ, :E => E, :ν => ν, :h => h)

        # 邊界懲罰元素：懲罰參數與目標值
        for el in elements_w_b
            prescribe!(el, :α => αʷ, :g => (x, y, z) -> 0.0)
        end
        for el in elements_φ_b
            prescribe!(el, :α => αᵠ, :g₁ => (x, y, z) -> 0.0, :g₂ => (x, y, z) -> 0.0, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        end

        # ---- C3. 計算形函數 (值或梯度) ----
        set∇𝝭!(elements_q)  
        set𝝭!(elements_q_Γ)  
        set∇𝝭!(elements_w)  
        set𝝭!(elements_w_Γ) 
        set∇𝝭!(elements_m) 
        set𝝭!(elements_m_Γ) 
        set∇𝝭!(elements_φ)  
        set𝝭!(elements_φ_Γ) 


        # ======================================================================
        # (D) 組裝剛度矩陣 (不含邊界懲罰與幾何剛度)
        # ======================================================================
        # 材料剛度子塊（混合形式）
        kˢˢ = ∫QQdΩ => elements_q
        kˢʷ = -∫∇QwdΩ => (elements_q, elements_w)   # 域內部分
        kˢᵠ = ∫QφdΩ => (elements_q, elements_φ)
        kᵠᵠ = ∫κκdΩ => elements_φ                 # 彎曲貢獻
# 如果需要將剪切對 φ 的貢獻也放入 kᵠᵠ，可加上 ∫φφdΩ => elements_φ
# 但混合形式中，剪切貢獻已通過 Q 場的耦合體現，視推導而定

        # 邊界貢獻（自然條件）
        kˢʷ = [∫∇QwdΩ => (elements_q, elements_w),
                ∫QwdΓ  => (elements_q_Γ, elements_w_Γ)]  # 合併域內+邊界
        kᵐᵠ = [∫∇MφdΩ => (elements_m, elements_φ),
                ∫MφdΓ  => (elements_m_Γ, elements_φ_Γ)]

        # ======================================================================
        # (E) 幾何剛度矩陣 (初始應力效應)
        # ======================================================================
        kᵍʷ = ∫∇wσ∇wdΩ => elements_w
        kᵍᵠ = ∫∇φσ∇φdΩ => elements_φ

        # ======================================================================
        # (G) 靜態凝聚形成最終屈曲特徵值問題
        # ======================================================================
        k = -[kʷʷ kᵠʷ'; kᵠʷ kᵠᵠ]
        Kᴳ = [kᵍʷ kᵍʷᵠ ; kᵍʷᵠ' kᵍᵠ]

        F = eigen(k, Kᴳ)
        λ = F.values

        gmsh.finalize()

    end
end