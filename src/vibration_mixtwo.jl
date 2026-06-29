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

# ==================== 1. 全域物理參數與理論解設定 ====================
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
αᵠ = 0.0

eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12
integrationOrder = 2
sʷ = 1.5
sᵠ = 1.5

case_prefix = "vibration_mixtwo"
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
if !ispath(data_dir)
    mkpath(data_dir)
end

function exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
    ν_ex = 0.3
    E_ex = Dᵇ * 12 * (1 - ν_ex^2) / h^3
    Dˢ_ex = (5 / 6) * E_ex * h / (2 * (1 + ν_ex))

    α_m = m * π / a
    β_n = n * π / b
    λ² = α_m^2 + β_n^2

    I₂ = ρ * h^3 / 12
    A_coef = ρ * h * I₂
    B_coef = ρ * h * Dᵇ * λ² + Dˢ_ex * (ρ * h + I₂ * λ²)
    C_coef = Dᵇ * Dˢ_ex * (λ²^2)

    ω_1 = sqrt((B_coef - sqrt(B_coef^2 - 4 * A_coef * C_coef)) / (2 * A_coef))
    w_h_exact = ω_1 * a^2 * sqrt(ρ * h / Dᵇ)

    return w_h_exact, w_h_exact
end

ndiv_series = 9:25
master_results = []

println("🚀 開始執行全網格自動化特徵值計算 (無 L2 積分優化版)...")

# ==================== 2. 網格收斂迴圈 ====================
for n_div in ndiv_series
    println("🔷 [網格級別] 正在計算 n_div = $n_div ...")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()

    ndiv_w = n_div - 2
    ndiv_φ = n_div
    ndiv = n_div
    mesh_file_w = normpath(joinpath(@__DIR__, "..", "msh", "st_q_$(ndiv_w).msh"))
    mesh_file_φ = normpath(joinpath(@__DIR__, "..", "msh", "st_q_$(ndiv_φ).msh"))
    mesh_file_s = normpath(joinpath(@__DIR__, "..", "msh", "st_q_$(ndiv).msh"))

    if !isfile(mesh_file_w) || !isfile(mesh_file_φ) || !isfile(mesh_file_s)
        gmsh.finalize()
        continue
    end

    type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_Q = :tri3
    type_M = :(PiecewisePolynomial{:Linear2D})

    # ==========================================================================
    # (A) 載入多場網格數據
    # ==========================================================================
    gmsh.open(mesh_file_w)
    nodes_w = get𝑿ᵢ()
    nʷ = length(nodes_w)
    s_w = 1.0 / (ndiv_w - 1)
    push!(nodes_w, :s₁ => sʷ * s_w * ones(nʷ),
                   :s₂ => sʷ * s_w * ones(nʷ),
                   :s₃ => sʷ * s_w * ones(nʷ))
    sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

    gmsh.clear()
    gmsh.open(mesh_file_φ)
    nodes_φ = get𝑿ᵢ()
    nᵠ = length(nodes_φ)
    s_𝜙 = 1.0 / (ndiv_φ - 1)
    push!(nodes_φ, :s₁ => sᵠ * s_𝜙 * ones(nᵠ),
                   :s₂ => sᵠ * s_𝜙 * ones(nᵠ),
                   :s₃ => sᵠ * s_𝜙 * ones(nᵠ))
    sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

    gmsh.clear()
    gmsh.open(mesh_file_s)
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
    # 區域元素：材料參數與質量密度
    prescribe!(elements_q, :E => E, :ν => ν, :h => h)
    prescribe!(elements_w, :E => E, :ν => ν, :h => h, :ρ => ρ)
    prescribe!(elements_m, :E => E, :ν => ν, :h => h)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :ρ => ρ)

    # 邊界懲罰元素：懲罰參數與目標值
    for el in elements_w_b
        prescribe!(el, :α => αʷ, :g => (x, y, z) -> 0.0)
    end
    for el in elements_φ_b
        prescribe!(el, :α => αᵠ, :g₁ => (x, y, z) -> 0.0, :g₂ => (x, y, z) -> 0.0,
                         :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    end

    # ---- C3. 計算形函數 (值或梯度) ----
    set∇𝝭!(elements_q)     # Q 場需要梯度
    set𝝭!(elements_q_Γ)    # Q 邊界只需值

    set𝝭!(elements_w)      # w 區域內積分僅需形函數值
    set𝝭!(elements_w_Γ)    # w 邊界積分僅需形函數值

    set∇𝝭!(elements_m)     # M 場需要梯度
    set𝝭!(elements_m_Γ)    # M 邊界只需值

    set∇𝝭!(elements_φ)     # φ 場需要梯度
    set𝝭!(elements_φ_Γ)    # φ 邊界只需值

    # 邊界懲罰元素僅需形函數值
    for el in elements_w_b
        set𝝭!(el)
    end
    for el in elements_φ_b
        set𝝭!(el)
    end

    # ==========================================================================
    # (D) 組裝剛度矩陣 (不含邊界懲罰)
    # ==========================================================================
    # 剪切相關
    (∫QQdΩ => elements_q)(kˢˢ)
    ([∫∇QwdΩ => (elements_q, elements_w),
      ∫QwdΓ  => (elements_q_Γ, elements_w_Γ)])(kˢʷ)
    (∫QφdΩ => (elements_q, elements_φ))(kˢᵠ)

    # 彎矩相關 (需先計算 M 場自由度數)
    nₑ = length(elements_q)
    nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))
    nᵐ = nₑ * nᵖ
    kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ)
    kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)

    (∫MMdΩ => elements_m)(kᵐᵐ)
    ([∫∇MφdΩ => (elements_m, elements_φ),
      ∫MφdΓ  => (elements_m_Γ, elements_φ_Γ)])(kᵐᵠ)

    # ==========================================================================
    # (E) 邊界懲罰項 (強制固支)
    # ==========================================================================
    dummy_fʷ = zeros(nʷ)
    dummy_fᵠ = zeros(2 * nᵠ)
    (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kᵅʷʷ, dummy_fʷ)
    (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵅᵠᵠ, dummy_fᵠ)

    # ==========================================================================
    # (F) 質量矩陣
    # ==========================================================================
    (∫ρhwwdΩ => elements_w)(mʷʷ)
    (∫ρIφφdΩ => elements_φ)(mᵠᵠ)
    M_total = [mᵠᵠ zeros(2 * nᵠ, nʷ); zeros(nʷ, 2 * nᵠ) mʷʷ]

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
    # (I) 頻率比對與數據記錄
    # ==========================================================================
    for r in 1:n_modes_output
        mode_id = mode_ids[r]
        ω_real_m = sqrt(real(ω²[mode_id]))
        w_h_FEM_m = ω_real_m * a^2 * sqrt(ρ * h / Dᵇ)

        best_err = Inf
        best_mn = (1, 1)
        for m_harm in 1:15, n_harm in 1:15
            _, w_h_ex = exact_omega_sq(m_harm, n_harm, a, b, Dᵇ, ρ, h)
            if abs(w_h_FEM_m - w_h_ex) < best_err
                best_err = abs(w_h_FEM_m - w_h_ex)
                best_mn = (m_harm, n_harm)
            end
        end
        m_match, n_match = best_mn
        _, w_h_exact_m = exact_omega_sq(m_match, n_match, a, b, Dᵇ, ρ, h)
        log10_Error_w_h = log10(max(abs(w_h_FEM_m - w_h_exact_m) / w_h_exact_m, eps(Float64)))

        push!(master_results, (n_div, r, m_match, n_match,
                               w_h_FEM_m, w_h_exact_m, log10_Error_w_h))
    end

    gmsh.finalize()
end

# ==================== 3. 匯出單一匯總 CSV ====================
csv_master_path = joinpath(data_dir, "$(case_prefix)_master_all_mesh.csv")
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
println("成功！")