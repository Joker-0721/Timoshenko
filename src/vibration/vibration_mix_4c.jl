const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ  

using LinearAlgebra
using TimerOutputs, WriteVTK, Printf
import Gmsh: gmsh

# ==================== 1. 全域物理參數與基準解設定 ====================
E = 1.0e6
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.001
h = h_over_b * b
ρ = 1.0
I_moment = h^3 / 12
Dᵇ = E * h^3 / (12 * (1 - ν^2))
Dˢ = 5 / 6 * E * h / (2 * (1 + ν))
αʷ = 1.0e8 * Dᵇ
αᵠ = 1.0e8 * Dᵇ          # 四邊固支 (CCCC)

eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12
integrationOrder = 2
sʷ = 1.5
sᵠ = 1.5

case_prefix = "vibration_mix_4c"
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
if !ispath(data_dir)
    mkpath(data_dir)
end

# 🌟 四邊固支 (CCCC) 前六階無因次頻率參數 λ² = ω a² √(ρh/D)
# 取自 Leissa (1969) 表 4.1 (Kirchhoff 板)，因 h/b=0.001 極薄，Mindlin 修正可忽略。
# 視為 Ritz 法收斂基準解，用於誤差評估。
const CCCC_ref_w = [35.992,   # Mode 1  (1,1)
                    73.413,   # Mode 2  (2,1)
                    73.413,   # Mode 3  (1,2)
                    108.27,   # Mode 4  (2,2)
                    131.64,   # Mode 5  (3,1)
                    132.24]   # Mode 6  (1,3)

# 主數據容器
ndiv_series = 9:25
master_results = []

# ==================== 2. 網格收斂迴圈 ====================
for n_div in ndiv_series
    println("🔷 [網格級別] 正在計算 n_div = $n_div ...")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()

    mesh_file = normpath(joinpath(@__DIR__, "..", "msh", "st_q_$(n_div).msh"))
    if !isfile(mesh_file)
        gmsh.finalize()
        continue
    end

    # --- 定義場類型與輔助尺度 ---
    type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_Q = :tri3
    type_M = :(PiecewisePolynomial{:Linear2D})
    s_size = 1.0 / (n_div - 1)

    # ==========================================================================
    # (A) 載入網格並建立各場的節點與離散點
    # ==========================================================================
    # --- w 場 (橫向位移) ---
    gmsh.clear()
    gmsh.open(mesh_file)
    nodes_w = get𝑿ᵢ()
    nʷ = length(nodes_w)
    push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ),
                   :s₂ => sʷ * s_size * ones(nʷ),
                   :s₃ => sʷ * s_size * ones(nʷ))
    sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

    # --- φ 場 (截面轉角) ---
    gmsh.clear()
    gmsh.open(mesh_file)
    nodes_φ = get𝑿ᵢ()
    nᵠ = length(nodes_φ)
    push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ),
                   :s₂ => sᵠ * s_size * ones(nᵠ),
                   :s₃ => sᵠ * s_size * ones(nᵠ))
    sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

    # --- Q 場與 M 場共用的底層網格 ---
    gmsh.clear()
    gmsh.open(mesh_file)
    nodes = get𝑿ᵢ()
    entities = getPhysicalGroups()
    nˢ = length(nodes)

    # ==========================================================================
    # (B) 矩陣初始化
    # ==========================================================================
    kˢˢ   = zeros(2 * nˢ, 2 * nˢ)
    kˢʷ   = zeros(2 * nˢ, nʷ)
    kˢᵠ   = zeros(2 * nˢ, 2 * nᵠ)
    kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
    kᵅʷʷ = zeros(nʷ, nʷ)
    mʷʷ   = zeros(nʷ, nʷ)
    mᵠᵠ   = zeros(2 * nᵠ, 2 * nᵠ)

    # ==========================================================================
    # (C) 定義元素、賦予物理常數、計算形函數
    # ==========================================================================
    # ---- C1. 定義所有元素物件 ----
    elements_q    = getElements(nodes, entities["Ω"], integrationOrder)
    elements_w    = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    elements_w_Γ  = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    elements_q_Γ  = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    elements_m    = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
    elements_φ    = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    elements_φ_Γ  = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    elements_m_Γ  = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)

    # ---- C2. 統一賦予物理常數 ----
    prescribe!(elements_q, :E => E, :ν => ν, :h => h)
    prescribe!(elements_w, :E => E, :ν => ν, :h => h)
    prescribe!(elements_m, :E => E, :ν => ν, :h => h)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h)

    # ---- C3. 計算形函數（值與梯度） ----
    set∇𝝭!(elements_q)     # Q 場需要梯度
    set𝝭!(elements_w)      # w 場區域內只需值
    set𝝭!(elements_w_Γ)    # w 邊界只需值
    set𝝭!(elements_q_Γ)    # Q 邊界只需值
    set∇𝝭!(elements_m)     # M 場需要梯度
    set∇𝝭!(elements_φ)     # φ 場需要梯度
    set𝝭!(elements_φ_Γ)    # φ 邊界只需值
    set𝝭!(elements_m_Γ)    # M 邊界只需值

    # ==========================================================================
    # (D) 組裝剛度矩陣 (不含邊界懲罰項)
    # ==========================================================================
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

    # ==========================================================================
    # (E) 邊界懲罰項 (強制固支：w=0, ψ=0)
    # ==========================================================================
    boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]

    elements_w_b = [getElements(nodes_w, entities[name], eval(type_w),
                                integrationOrder, sp_w, normal=true) for name in boundary_names]
    for el in elements_w_b
        prescribe!(el, :α => αʷ, :g => (x, y, z) -> 0.0)
        set𝝭!(el)
    end
    dummy_fʷ = zeros(nʷ)
    (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kᵅʷʷ, dummy_fʷ)

    elements_φ_b = [getElements(nodes_φ, entities[name], eval(type_φ),
                                integrationOrder, sp_φ, normal=true) for name in boundary_names]
    for el in elements_φ_b
        prescribe!(el, :α => αᵠ, :g₁ => (x, y, z) -> 0.0, :g₂ => (x, y, z) -> 0.0,
                         :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        set𝝭!(el)
    end
    dummy_fᵠ = zeros(2 * nᵠ)
    (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵅᵠᵠ, dummy_fᵠ)

    # ==========================================================================
    # (F) 質量矩陣
    # ==========================================================================
    prescribe!(elements_w, :ρ => ρ, :h => h)
    (∫ρhwwdΩ => elements_w)(mʷʷ)

    prescribe!(elements_φ, :ρ => ρ, :h => h)
    (∫ρIφφdΩ => elements_φ)(mᵠᵠ)

    M_total = [mᵠᵠ                  zeros(2 * nᵠ, nʷ)
               zeros(nʷ, 2 * nᵠ)   mʷʷ]

    # ==========================================================================
    # (G) 靜態凝聚形成最終特徵值問題
    # ==========================================================================
    k_cond = -[kˢᵠ' * (kˢˢ \ kˢᵠ) + kᵐᵠ' * (kᵐᵐ \ kᵐᵠ) - kᵅᵠᵠ   kˢᵠ' * (kˢˢ \ kˢʷ)
               kˢʷ' * (kˢˢ \ kˢᵠ)                                    kˢʷ' * (kˢˢ \ kˢʷ) - kᵅʷʷ]
    ks = Symmetric(0.5 * (k_cond + k_cond'))
    ms = Symmetric(0.5 * (M_total + M_total'))

    # ==========================================================================
    # (H) 求解廣義特徵值
    # ==========================================================================
    F = eigen(ks, ms)
    ω² = F.values

    mode_ids = sort!([i for i in eachindex(ω²)
                      if real(ω²[i]) > omega_sq_tol && abs(imag(ω²[i])) < eigen_imag_tol],
                     by = i -> real(ω²[i]))
    n_modes_output = min(length(mode_ids), 20)

    # ==========================================================================
    # (I) 頻率比對與數據記錄（使用 Ritz 基準解）
    # ==========================================================================
    n_ref = length(CCCC_ref_w)
    n_compare = min(n_modes_output, n_ref)
    for r in 1:n_compare
        mode_id = mode_ids[r]
        ω_real_m = sqrt(real(ω²[mode_id]))
        w_h_FEM_m = ω_real_m * a^2 * sqrt(ρ * h / Dᵇ)
        w_h_exact_m = CCCC_ref_w[r]          # 基準解
        log10_Error_w_h = log10(max(abs(w_h_FEM_m - w_h_exact_m) / w_h_exact_m, eps(Float64)))

        # 輸出時 matched_m 設為模態排序編號，matched_n 設為 0
        push!(master_results, (n_div, r, r, 0,
                               w_h_FEM_m, w_h_exact_m, log10_Error_w_h))
    end

    gmsh.finalize()
    println("   -> n_div = $n_div 計算完成。")
end

# ==================== 3. 匯出 CSV 總表 ====================
csv_master_path = joinpath(data_dir, "$(case_prefix).csv")
open(csv_master_path, "w") do io
    println(io, join(["n_div", "mode_rank", "matched_m", "matched_n",
                      "w_h_FEM", "w_h_exact", "log10_Error_w_h"], ","))
    for data in master_results
        n_div, rank, m, n, w_FEM, w_exact, err_wh = data
        println(io, join([n_div, rank, m, n,
                          @sprintf("%.4f", w_FEM),
                          @sprintf("%.4f", w_exact),
                          @sprintf("%.6f", err_wh)], ","))
    end
end
println("📊 成功導出 Mix (CCCC) 總表至: $csv_master_path")