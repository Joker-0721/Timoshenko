const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ

using LinearAlgebra
using Printf
import Gmsh: gmsh

# ==================== 1. 全域物理參數與基準解設定 ====================
E = 1.0e6
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.001
h = h_over_b * b
ρ = 1.0
Dᵇ = E * h^3 / 12 / (1 - ν^2)
αʷ = 1.0e8 * Dᵇ          # w 邊界懲罰剛度
αᵠ = 1.0e8 * Dᵇ          # φ 邊界懲罰剛度
eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12

case_prefix = "vibration_fem_4c"
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
if !ispath(data_dir)
    mkpath(data_dir)
end

# 🌟 四邊固支 (CCCC) 前六階無因次頻率參數 λ² = ω a² √(ρh/D)
# 數據取自 Leissa (1969) 表 4.1 (Kirchhoff 板)，因 h/b=0.001 極薄，Mindlin 修正可忽略。
# 亦可視為 Ritz 法收斂數值，用於誤差評估。
const CCCC_ref_w = [35.992,   # Mode 1  (1,1)
                    73.413,   # Mode 2  (2,1)
                    73.413,   # Mode 3  (1,2)
                    108.27,   # Mode 4  (2,2)
                    131.64,   # Mode 5  (3,1)
                    132.24]   # Mode 6  (1,3)

# 原簡支解析解函數保留但不再使用（若需要可移除）
# function exact_omega_sq(...) ... end

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
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    nʷ = length(nodes)
    nᵠ = length(nodes)

    # ==========================================================================
    # (B) 矩陣初始化
    # ==========================================================================
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
    kᵠʷ = zeros(2 * nᵠ, nʷ)
    mʷʷ = zeros(nʷ, nʷ)
    mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    # ==========================================================================
    # (C) 定義元素、賦予物理常數、計算形函數
    # ==========================================================================
    elements     = getElements(nodes, entities["Ω"], 2)    # 彎曲 (完全積分)
    elements_s   = getElements(nodes, entities["Ω"], 1)    # 剪切 (縮減積分)

    boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
    elements_w_b = [getElements(nodes, entities[name], 2) for name in boundary_names]
    elements_φ_b = [getElements(nodes, entities[name], 2) for name in boundary_names]

    # 區域物理常數
    prescribe!(elements,   :E => E, :ν => ν, :h => h, :ρ => ρ)
    prescribe!(elements_s, :E => E, :ν => ν, :h => h)

    # 邊界懲罰條件
    for el in elements_w_b
        prescribe!(el, :α => αʷ, :g => (x, y, z) -> 0.0)
    end
    for el in elements_φ_b
        prescribe!(el, :α => αᵠ, :g₁ => (x, y, z) -> 0.0, :g₂ => (x, y, z) -> 0.0,
                         :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    end

    # 形函數計算
    set∇𝝭!(elements)
    set∇𝝭!(elements_s)
    for el in elements_w_b
        set𝝭!(el)
    end
    for el in elements_φ_b
        set𝝭!(el)
    end

    # ==========================================================================
    # (D) 組裝剛度矩陣 (不含邊界懲罰)
    # ==========================================================================
    (∫wwdΩ => elements_s)(kʷʷ)
    (∫φwdΩ => elements_s)(kᵠʷ)
    ([∫φφdΩ => elements_s, ∫κκdΩ => elements])(kᵠᵠ)

    # ==========================================================================
    # (E) 組裝質量矩陣
    # ==========================================================================
    (∫ρhwwdΩ => elements)(mʷʷ)
    (∫ρIφφdΩ => elements)(mᵠᵠ)

    # ==========================================================================
    # (F) 邊界懲罰項 (四邊固支：w=0, ψ=0)
    # ==========================================================================
    dummy_fʷ = zeros(nʷ)
    dummy_fᵠ = zeros(2 * nᵠ)
    (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kʷʷ, dummy_fʷ)
    (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵠᵠ, dummy_fᵠ)

    # ==========================================================================
    # (G) 組裝全域矩陣並對稱化
    # ==========================================================================
    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    M = [mᵠᵠ zeros(2 * nᵠ, nʷ); zeros(nʷ, 2 * nᵠ) mʷʷ]
    Ks = Symmetric(0.5 * (K + K'))
    Ms = Symmetric(0.5 * (M + M'))

    # ==========================================================================
    # (H) 求解廣義特徵值
    # ==========================================================================
    F = eigen(Ks, Ms)
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
        w_h_exact_m = CCCC_ref_w[r]    # 取自 Ritz 法基準值
        log10_Error_w_h = log10(max(abs(w_h_FEM_m - w_h_exact_m) / w_h_exact_m, eps(Float64)))

        # 輸出時，matched_m 儲存模態排序編號，matched_n 設為 0（無波數含義）
        push!(master_results, (n_div, r, r, 0,
                               w_h_FEM_m, w_h_exact_m, log10_Error_w_h))
    end

    gmsh.finalize()
    println("   -> n_div = $n_div 計算完成。")
end

# ==================== 3. 匯出單一匯總 CSV ====================
csv_master_path = joinpath(data_dir, "$(case_prefix)_FEM.csv")
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
println("📊 成功導出 FEM 總表至: $csv_master_path")