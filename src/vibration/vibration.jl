const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
# 🌟 核心修正：完全移除了 L₂w, L₂φ 函數引進
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫αwwdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ

using LinearAlgebra
using Printf
import Gmsh: gmsh

# ==================== 1. 全域物理參數與理論解設定 ====================
E = 1.0e6; ν = 0.3; a = 1.0; b = 1.0; h_over_b = 0.001; h = h_over_b * b; ρ = 1.0
Dᵇ = E * h^3 / 12 / (1 - ν^2)
α = 1.0e8 * Dᵇ
eigen_imag_tol = 1.0e-7; omega_sq_tol = 1.0e-12

case_prefix = "vibration_fem"
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
if !ispath(data_dir); mkpath(data_dir); end

function exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
    ν_ex = 0.3
    E_ex = Dᵇ * 12 * (1 - ν_ex^2) / h^3
    Dˢ_ex = (5/6) * E_ex * h / (2 * (1 + ν_ex))
    
    α_m = m * π / a
    β_n = n * π / b
    λ² = α_m^2 + β_n^2
    
    # 1. 定義顯式代數係數
    I₂ = ρ * h^3 / 12
    A_coef = ρ * h * I₂
    B_coef = ρ * h * Dᵇ * λ² + Dˢ_ex * (ρ * h + I₂ * λ²)
    C_coef = Dᵇ * Dˢ_ex * (λ²^2)  # 即 λ⁴
    
    # 2. 直接帶入公式求出最低頻彎曲角頻率 ω_1
    ω_1 = sqrt((B_coef - sqrt(B_coef^2 - 4 * A_coef * C_coef)) / (2 * A_coef))
    
    # 3. 轉換為無因次頻率
    w_h_exact = ω_1 * a^2 * sqrt(ρ * h / Dᵇ)
    
    # 🌟 保持原本代碼的元組回傳格式，確保解包相容性
    return w_h_exact, w_h_exact
end

ndiv_series = 9:25
master_results = []

# ==================== 2. 自動化網格大迴圈 ====================

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
    
    gmsh.open(mesh_file)
    entities = getPhysicalGroups(); nodes = get𝑿ᵢ()
    nʷ = length(nodes); nᵠ = length(nodes)

    # ==========================================================================
    # ---- (B) 矩陣初始化與標準有限元組裝 ----
    # ==========================================================================
    kʷʷ = zeros(nʷ, nʷ); kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ); kᵠʷ = zeros(2 * nᵠ, nʷ)
    mʷʷ = zeros(nʷ, nʷ); mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    elements = getElements(nodes, entities["Ω"], 2)          # 正常積分(彎曲)
    elements_s = getElements(nodes, entities["Ω"], 1)        # 縮減積分(剪切)
    
    prescribe!(elements, :E => E, :ν => ν, :h => h); set∇𝝭!(elements)
    prescribe!(elements_s, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_s)

    (∫wwdΩ => elements_s)(kʷʷ)
    (∫φwdΩ => elements_s)(kᵠʷ)
    ([∫φφdΩ => elements_s, ∫κκdΩ => elements])(kᵠᵠ)
    
    prescribe!(elements, :ρ => ρ); (∫ρhwwdΩ => elements)(mʷʷ)
    (∫ρIφφdΩ => elements)(mᵠᵠ)

    boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
    elements_b = [getElements(nodes, entities[name], 2) for name in boundary_names]
    for el in elements_b; prescribe!(el, :α => α, :g => (x,y,z)->0.0); set𝝭!(el); end
    dummy_f = zeros(nʷ)
    (∫αwwdΓ => elements_b[1] ∪ elements_b[2] ∪ elements_b[3] ∪ elements_b[4])(kʷʷ, dummy_f)

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    M = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
    Ks = Symmetric(0.5*(K + K'))
    Ms = Symmetric(0.5*(M + M'))

    # ==========================================================================
    # ---- (C) 求解特徵值問題 ----
    # ==========================================================================
    F = eigen(Ks, Ms); ω² = F.values

    mode_ids = sort!([i for i in eachindex(ω²) if real(ω²[i]) > omega_sq_tol && abs(imag(ω²[i])) < eigen_imag_tol], by = i -> real(ω²[i]))
    n_modes_output = min(length(mode_ids), 20)

    # ==========================================================================
    # ---- (D) 數據比對與提取 (不包含任何空間 L2 殘差運算) ----
    # ==========================================================================
    for r in 1:n_modes_output
        mode_id = mode_ids[r]
        ω_real_m = sqrt(real(ω²[mode_id]))
        w_h_FEM_m = ω_real_m * a^2 * sqrt(ρ * h / Dᵇ)

        best_err = Inf; best_mn = (1, 1)
        for m_harm in 1:15, n_harm in 1:15
            _, w_h_ex = exact_omega_sq(m_harm, n_harm, a, b, Dᵇ, ρ, h)
            if abs(w_h_FEM_m - w_h_ex) < best_err
                best_err = abs(w_h_FEM_m - w_h_ex)
                best_mn = (m_harm, n_harm)
            end
        end
        m_match, n_match = best_mn
        _,w_h_exact_m = exact_omega_sq(m_match, n_match, a, b, Dᵇ, ρ, h)
        log10_Error_w_h = log10(max(abs(w_h_FEM_m - w_h_exact_m)/w_h_exact_m, eps(Float64)))

        push!(master_results, (n_div, r, m_match, n_match, w_h_FEM_m, w_h_exact_m, log10_Error_w_h))
    end

    gmsh.finalize()
    println("   -> n_div = $n_div 計算完成。")
end

# ==================== 3. 匯出單一匯總 CSV ====================
csv_master_path = joinpath(data_dir, "$(case_prefix)_master_all_mesh.csv")
open(csv_master_path, "w") do io
    println(io, join(["n_div", "mode_rank", "matched_m", "matched_n", "w_h_FEM", "w_h_exact", "log10_Error_w_h"], ","))
    for data in master_results
        n_div, rank, m, n, w_FEM, w_exact, err_wh = data
        println(io, join([n_div, rank, m, n, @sprintf("%.4f", w_FEM), @sprintf("%.4f", w_exact), @sprintf("%.6f", err_wh)], ","))
        # println("$w_FEM, $w_exact")
    end
end
println("📊 成功導出 FEM 總表至: $csv_master_path")