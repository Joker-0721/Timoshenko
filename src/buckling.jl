const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ, ∫αφφdΓ

using LinearAlgebra
using TimerOutputs
using WriteVTK
using Printf
import Gmsh: gmsh

# 全局數據庫路徑規範
const DATA_DIR = "./data"
const UNIFIED_CSV = joinpath(DATA_DIR, "buckling_FEM.csv")
mkpath(DATA_DIR)

# 全局力學材料常數
const E = 200e9
const ν = 0.3
const h = 1e-2
const a = 1.0
const b = 1.0
const Dᵇ = E * h^3 / (12 * (1 - ν^2))

const α_base = 1.0e8 * Dᵇ
const σ₁₁ = 1.0
const σ₂₂ = 0.0
const σ₁₂ = 0.0

const k_exact_modes = [
    (1 + 1^2/1)^2,  # Mode 1: (m=1, n=1) -> 4.0000
    (2 + 1^2/2)^2,  # Mode 2: (m=2, n=1) -> 6.2500
    (3 + 1^2/3)^2,  # Mode 3: (m=3, n=1) -> 11.1111
    (2 + 2^2/2)^2,  # Mode 4: (m=2, n=2) -> 16.0000
    (4 + 1^2/4)^2,  # Mode 5: (m=4, n=1) -> 18.0625
    (3 + 2^2/3)^2   # Mode 6: (m=3, n=2) -> 18.7778
]

const integrationOrder = 2
const integrationOrder_shear = 1
const to = TimerOutput()

println("="^80)
println(" 執行 Pure FEM 模組：網格加密大循環 (ndiv = 9 ➔ 25) ")
println("="^80)

for n_div in 9:25
    @timeit to "Pure_FEM Loop (ndiv=$n_div)" begin

        mesh_path = "/home/a/Joker/msh/st_q_$(n_div).msh"
        if !isfile(mesh_path)
            println("  [WARN] 找不到指定的網格檔案，跳過此輪: ", mesh_path)
            continue
        end

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(mesh_path)

        entities = getPhysicalGroups()
        nodes = get𝑿ᵢ()

        nʷ = length(nodes)
        nᵠ = length(nodes)

        # ======================================================================
        # (A) 矩陣初始化
        # ======================================================================
        kʷʷ   = zeros(nʷ, nʷ)
        kᵠᵠ   = zeros(2 * nᵠ, 2 * nᵠ)
        kᵠʷ   = zeros(2 * nᵠ, nʷ)
        kᴳʷʷ = zeros(nʷ, nʷ)
        kᴳᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

        # ======================================================================
        # (B) 定義元素、賦予物理常數、計算形函數
        # ======================================================================
        # ---- B1. 定義所有元素物件 ----
        elements   = getElements(nodes, entities["Ω"], integrationOrder)        # 正常積分 (彎曲 + 幾何剛度)
        elements_s = getElements(nodes, entities["Ω"], integrationOrder_shear)  # 縮減積分 (剪切)

        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_w_b   = [getElements(nodes, entities[name], integrationOrder) for name in boundary_names]
        elements_φ_b   = [getElements(nodes, entities[name], integrationOrder) for name in boundary_names]

        # ---- B2. 統一賦予物理常數 ----
        # 區域元素：材料參數 + 幾何剛度所需應力
        prescribe!(elements,   :E => E, :ν => ν, :h => h,
                               :σ₁₁ => σ₁₁, :σ₂₂ => σ₂₂, :σ₁₂ => σ₁₂)
        prescribe!(elements_s, :E => E, :ν => ν, :h => h)

        # 邊界懲罰元素
        for el in elements_w_b
            prescribe!(el, :α => α_base, :g => 0.0)
        end
        for el in elements_φ_b
            prescribe!(el, :α => α_base, :g₁ => 0.0, :g₂ => 0.0,
                             :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        end

        # ---- B3. 計算形函數 (值或梯度) ----
        set∇𝝭!(elements)      # 彎曲、幾何剛度需要梯度
        set∇𝝭!(elements_s)    # 剪切需要梯度
        for el in elements_w_b
            set𝝭!(el)         # 邊界懲罰只需值
        end
        for el in elements_φ_b
            set𝝭!(el)
        end

        # ======================================================================
        # (C) 組裝剛度矩陣 (不含邊界懲罰)
        # ======================================================================
        (∫wwdΩ => elements_s)(kʷʷ)
        (∫φwdΩ => elements_s)(kᵠʷ)
        ([∫φφdΩ => elements_s, ∫κκdΩ => elements])(kᵠᵠ)

        # 幾何剛度 (初始應力效應)
        (∫∇wσ∇wdΩ => elements)(kᴳʷʷ)
        (∫∇φσ∇φdΩ => elements)(kᴳᵠᵠ)

        # ======================================================================
        # (D) 邊界懲罰項 (強制固支)
        # ======================================================================
        for el in elements_w_b
            (∫αwwdΓ => el)(kʷʷ)
        end
        for el in elements_φ_b
            (∫αφφdΓ => el)(kᵠᵠ)
        end

        # ======================================================================
        # (E) 形成總剛度與幾何剛度矩陣，求解特徵值
        # ======================================================================
        K   = [kᵠᵠ   kᵠʷ;   kᵠʷ'   kʷʷ]
        Kᴳ  = [kᴳᵠᵠ  zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ)  kᴳʷʷ]

        F = eigen(K, Kᴳ)
        λ = F.values

        mode_ids = sort!(
            collect(i for i in eachindex(λ)
                if isfinite(real(λ[i])) && abs(imag(λ[i])) < 1.0e-7 && real(λ[i]) > 0.0),
            by = i -> real(λ[i])
        )

        # ======================================================================
        # (F) 計算屈曲係數並記錄
        # ======================================================================
        h_size = 1.0 / n_div
        log10_h = log10(h_size)

        file_exists = isfile(UNIFIED_CSV)
        open(UNIFIED_CSV, "a") do io
            if !file_exists
                println(io, "method,ndiv,h,log10_h,mode_rank,lambda_cr,k_num,k_exact,error_y")
                file_exists = true
            end

            for rank in 1:min(6, length(mode_ids))
                m_id = mode_ids[rank]
                lam = real(λ[m_id])
                k_num = (lam * h * σ₁₁) * b^2 / (π^2 * Dᵇ)
                k_ex = k_exact_modes[rank]
                error_y = (k_num / k_ex) - 1.0

                @printf(io, "Pure_FEM,%d,%.6e,%.6f,%d,%.6e,%.6f,%.6f,%.6e\n",
                        n_div, h_size, log10_h, rank, lam, k_num, k_ex, error_y)
            end
        end

        gmsh.finalize()
        println("  [Pure_FEM] ndiv = $(n_div) 前 6 階數據導出成功。")
    end
end