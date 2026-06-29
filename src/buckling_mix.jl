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

println("="^80)
println(" 執行 Mix 模組：網格加密大循環 (ndiv = 9 ➔ 25) ")
println("="^80)

for n_div in 9:25
    @timeit to "Mix Loop (ndiv=$n_div)" begin

        mesh_file = "/home/a/Joker/msh/st_q_$(n_div).msh"
        if !isfile(mesh_file)
            println("  [WARN] 找不到指定的網格檔案，跳過此輪: ", mesh_file)
            continue
        end

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
        kˢˢ   = zeros(2 * nˢ, 2 * nˢ)
        kˢʷ   = zeros(2 * nˢ, nʷ)
        kˢᵠ   = zeros(2 * nˢ, 2 * nᵠ)
        kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
        kᵅʷʷ = zeros(nʷ, nʷ)
        kʷʷ   = zeros(nʷ, nʷ)      # 幾何剛度矩陣

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
        elements_w_b = [getElements(nodes_w, entities[name], eval(type_w),
                                    integrationOrder, sp_w, normal=true) for name in boundary_names]
        elements_φ_b = [getElements(nodes_φ, entities[name], eval(type_φ),
                                    integrationOrder, sp_φ, normal=true) for name in boundary_names]

        # ---- C2. 統一賦予物理常數 ----
        # 區域元素：材料參數與應力 (w 場需要幾何剛度應力，提前給定)
        prescribe!(elements_q, :E => E, :ν => ν, :h => h)
        prescribe!(elements_w, :E => E, :ν => ν, :h => h,
                              :σ₁₁ => σ₁₁, :σ₂₂ => σ₂₂, :σ₁₂ => σ₁₂)
        prescribe!(elements_m, :E => E, :ν => ν, :h => h)
        prescribe!(elements_φ, :E => E, :ν => ν, :h => h)

        # 邊界懲罰元素：懲罰參數與目標值
        for el in elements_w_b
            prescribe!(el, :α => αʷ, :g => (x, y, z) -> 0.0)
        end
        for el in elements_φ_b
            prescribe!(el, :α => αᵠ, :g₁ => (x, y, z) -> 0.0, :g₂ => (x, y, z) -> 0.0,
                             :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        end

        # ---- C3. 計算形函數 (值或梯度) ----
        set∇𝝭!(elements_q)     # Q 區域需梯度
        set𝝭!(elements_q_Γ)    # Q 邊界只需值

        set∇𝝭!(elements_w)     # w 區域需梯度 (用於剪切與幾何剛度)
        set𝝭!(elements_w_Γ)    # w 邊界只需值

        set∇𝝭!(elements_m)     # M 區域需梯度
        set𝝭!(elements_m_Γ)    # M 邊界只需值

        set∇𝝭!(elements_φ)     # φ 區域需梯度
        set𝝭!(elements_φ_Γ)    # φ 邊界只需值

        for el in elements_w_b
            set𝝭!(el)          # 邊界懲罰僅需值
        end
        for el in elements_φ_b
            set𝝭!(el)          # 邊界懲罰僅需值
        end

        # ======================================================================
        # (D) 組裝剛度矩陣 (不含邊界懲罰與幾何剛度)
        # ======================================================================
        (∫QQdΩ => elements_q)(kˢˢ)
        ([∫∇QwdΩ => (elements_q, elements_w),
          ∫QwdΓ  => (elements_q_Γ, elements_w_Γ)])(kˢʷ)
        (∫QφdΩ => (elements_q, elements_φ))(kˢᵠ)

        nₑ = length(elements_q)
        nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))
        nᵐ = nₑ * nᵖ
        kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ)
        kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)

        (∫MMdΩ => elements_m)(kᵐᵐ)
        ([∫∇MφdΩ => (elements_m, elements_φ),
          ∫MφdΓ  => (elements_m_Γ, elements_φ_Γ)])(kᵐᵠ)

        # ======================================================================
        # (E) 幾何剛度矩陣 (初始應力效應)
        # ======================================================================
        (∫∇wσ∇wdΩ => elements_w)(kʷʷ)

        # ======================================================================
        # (F) 邊界懲罰項 (強制固支)
        # ======================================================================
        dummy_fʷ = zeros(nʷ)
        dummy_fᵠ = zeros(2 * nᵠ)
        (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kᵅʷʷ, dummy_fʷ)
        (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵅᵠᵠ, dummy_fᵠ)

        # ======================================================================
        # (G) 靜態凝聚形成最終屈曲特徵值問題
        # ======================================================================
        k_cond = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ   kˢᵠ'*(kˢˢ\kˢʷ)
                   kˢʷ'*(kˢˢ\kˢᵠ)                            kˢʷ'*(kˢˢ\kˢʷ)-kᵅʷʷ]
        ks = Symmetric(0.5 * (k_cond + k_cond'))
        Kᴳ_cond = [zeros(2*nᵠ, 2*nᵠ)  zeros(2*nᵠ, nʷ)
                   zeros(nʷ, 2*nᵠ)    kʷʷ]

        F = eigen(ks, Kᴳ_cond)
        λ = F.values

        # ======================================================================
        # (H) 屈曲模態篩選與數據記錄
        # ======================================================================
        mode_ids = sort!([i for i in eachindex(λ)
                          if isfinite(real(λ[i])) && abs(imag(λ[i])) < 1.0e-7 && real(λ[i]) > 0.0],
                         by = i -> real(λ[i]))

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

                @printf(io, "Mix,%d,%.6e,%.6f,%d,%.6e,%.6f,%.6f,%.6e\n",
                        n_div, h_size, log10_h, rank, lam, k_num, k_ex, error_y)
            end
        end

        gmsh.finalize()
        println("  [Mix] ndiv = $(n_div) 前 6 階數據導出成功。")
    end
end