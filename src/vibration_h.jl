# ==============================================================================
#  2D Mindlin Plate Free Vibration Analysis - Q4 Element Baseline (Thickness Test)
#  流派：Standard Displacement FEM Q4 (標準有限元縮減積分基準線)
#  🌟 最新修改：網格死鎖固定為 17，厚度 h 從 0.1 跑到 0.00001 幾何大循環
#  特性：完美移植智慧幅值與正負號對齊，產出具備厚度鎖死對比意義的 L2 直線總表
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫αwwdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ,
                                     L₂w, L₂φ  # 🌟 修正點 1：改為引入底層真正完工的標準幾何算子

using LinearAlgebra, Printf, TimerOutputs, WriteVTK, DelimitedFiles
import Gmsh: gmsh

# ==================== 1. 材料與基礎幾何設定 ====================
E = 1.0e6; ν = 0.3; a = 1.0; b = 1.0; ρ = 1.0
eigen_imag_tol = 1.0e-7; omega_sq_tol = 1.0e-12; residual_warn_tol = 1.0e-6

# 🌟 滿足要求 1：網格死鎖固定使用 17 號網格
mesh_file = normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh"))
mesh_name = splitext(basename(mesh_file))[1]

# 🌟 滿足要求 2：H 自動設為幾何數列循環從 0.1 跑到 0.00001
h_series = [0.1, 0.01, 0.001, 0.0001, 0.00001]

# 🌟 滿足要求 4：宣告分流資料夾路徑
vtk_dir  = normpath(joinpath(@__DIR__, "..", "VTK"))
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
mkpath(vtk_dir); mkpath(data_dir)

# 建立用來繪製厚度鎖死圖的容器
h_plot_data = Float64[]
err_omega_plot_data = Float64[]
err_w_plot_data = Float64[]
err_phi_plot_data = Float64[]

const to = TimerOutput()

# ---------------------- 後處理輔助函數 ----------------------
function selected_eigenpair(ω²ᵢ, vᵢ)
    return isfinite(real(ω²ᵢ)) && isfinite(imag(ω²ᵢ)) && abs(imag(ω²ᵢ)) < eigen_imag_tol && real(ω²ᵢ) > omega_sq_tol
end
function center_node_id(nodes)
    return nodes[argmin([(node.x - a / 2)^2 + (node.y - b / 2)^2 for node in nodes])].𝐼
end
function normalize_mode(vᵢ, nodes, nᵠ)
    dm = real.(vᵢ); w_vec = dm[2*nᵠ+1:end]; dm ./= maximum(abs.(w_vec))
    if dm[2*nᵠ + center_node_id(nodes)] < -sqrt(eps(Float64)); dm .*= -1.0; end
    return dm
end
function relative_residual(K, M, ω²ᵢ, vᵢ)
    Kv = K * vᵢ; Mv = M * vᵢ
    return norm(Kv - ω²ᵢ * Mv) / max(norm(Kv), abs(ω²ᵢ) * norm(Mv), eps(Float64))
end
function vtk_cell(elm)
    return length(elm.𝓒) == 4 ? MeshCell(VTKCellTypes.VTK_QUAD, [x.𝐼 for x in elm.𝓒]) : MeshCell(VTKCellTypes.VTK_TRIANGLE, [x.𝐼 for x in elm.𝓒])
end
function exact_omega_sq(m, n, a, b, D, ρ, h)
    Ω_exact = π^2 * (m^2 + n^2); return Ω_exact * sqrt(D / (ρ * h)) / a^2, Ω_exact
end

# ==================== 2. 主程式厚度大循環 ====================
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)

for h in h_series
    println("\n" * "="^70)
    println("固定網格 17 ➔ 執行有限元厚度變更分析：h = $(h)")
    println("="^70)
    
    case_prefix = "vibration_fem_st_q_17_h_$(h)"
    integrationOrder = 2; integrationOrder_shear = 1

    # 🌟 隨著厚度改變，重新計算抗彎剛度與罰邊界參數
    Dᵇ = E * h^3 / 12 / (1 - ν^2)
    α = 1.0e8 * Dᵇ

    gmsh.clear(); gmsh.open(mesh_file); entities = getPhysicalGroups(); nodes = get𝑿ᵢ()
    nʷ = length(nodes); nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ); kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ); kᵠʷ = zeros(2 * nᵠ, nʷ)
    mʷʷ = zeros(nʷ, nʷ); mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    @timeit to "assemble stiffness and mass" begin
        elements = getElements(nodes, entities["Ω"], integrationOrder)
        elements_s = getElements(nodes, entities["Ω"], integrationOrder_shear)
        prescribe!(elements, :E => E, :ν => ν, :h => h); prescribe!(elements_s, :E => E, :ν => ν, :h => h)
        set∇𝝭!(elements); set∇𝝭!(elements_s)

        A_ww = ∫wwdΩ => elements_s; A_fw = ∫φwdΩ => elements_s; A_ff = [∫φφdΩ => elements_s, ∫κκdΩ => elements]
        A_ww(kʷʷ); A_fw(kᵠʷ); A_ff(kᵠᵠ)

        prescribe!(elements, :ρ => ρ); (∫ρhwwdΩ => elements)(mʷʷ); (∫ρIφφdΩ => elements)(mᵠᵠ)
        global elements_Ω = elements
    end

    @timeit to "apply SSSS boundary" begin
        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_b = [getElements(nodes, entities[name], integrationOrder) for name in boundary_names]
        for el in elements_b; prescribe!(el, :α => α, :g => (x,y,z)->0.0); set𝝭!(el); end
        A_b = ∫αwwdΓ => elements_b[1] ∪ elements_b[2] ∪ elements_b[3] ∪ elements_b[4]
        dummy_f = zeros(nʷ); A_b(kʷʷ, dummy_f)
    end

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]; M = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
    Ks = Symmetric(0.5*(K + K')); Ms = Symmetric(0.5*(M + M'))
    F = eigen(Ks, Ms); ω² = F.values; V = F.vectors

    global mode_ids = sort!([i for i in eachindex(ω²) if selected_eigenpair(ω²[i], V[:, i])], by = i -> real(ω²[i]))
    dm_modes = [normalize_mode(V[:, mode_id], nodes, nᵠ) for mode_id in mode_ids]
    residuals = [relative_residual(K, M, real(ω²[mode_id]), dm_modes[mode_rank]) for (mode_rank, mode_id) in enumerate(mode_ids)]

    # 🌟 滿足要求 2：提取第一階有效模態進行智慧幅值歸一化與方向對齊
    mode_id_1 = mode_ids[1]; dm_1 = dm_modes[1]; vᵠ_1 = dm_1[1:2*nᵠ]; vʷ_1 = dm_1[2*nᵠ+1:end]
    w_exact_func  = (x,y,z) -> sin(1 * π * x / a) * sin(1 * π * y / b)
    φx_exact_func = (x,y,z) -> (1 * π / a) * cos(1 * π * x / a) * sin(1 * π * y / b)
    φy_exact_func = (x,y,z) -> (1 * π / b) * sin(1 * π * x / a) * cos(1 * π * y / b)

    max_idx = argmax(abs.(vʷ_1))
    scale_factor = abs(vʷ_1[max_idx]) > 1e-12 ? w_exact_func(nodes[max_idx].x, nodes[max_idx].y, 0.0) / vʷ_1[max_idx] : 1.0
    vʷ_scaled = vʷ_1 .* scale_factor; vᵠ_scaled = vᵠ_1 .* scale_factor

    # 預配置數據抽屜，原地覆寫注入
    push!(nodes, :d => zeros(nʷ)); push!(nodes, :d₁ => zeros(nᵠ), :d₂ => zeros(nᵠ))
    for (i, node) in enumerate(nodes)
        node.d = vʷ_scaled[i]; nodes[i].d₁ = vᵠ_scaled[2*i-1]; nodes[i].d₂ = vᵠ_scaled[2*i]
    end

    prescribe!(elements_Ω, :w => w_exact_func); prescribe!(elements_Ω, :φ₁ => φx_exact_func, :φ₂ => φy_exact_func)
    L2_w_val   = L₂w(elements_Ω)
    L2_phi_val = L₂φ(elements_Ω)

    ω_real = sqrt(real(ω²[mode_id_1])); Ω_FEM = ω_real * a^2 * sqrt(ρ * h / Dᵇ)
    _, Ω_exact = exact_omega_sq(1, 1, a, b, Dᵇ, ρ, h)
    log10_Error_Omega = log10(max(abs(Ω_FEM - Ω_exact)/Ω_exact, eps(Float64)))

    push!(h_plot_data, h)
    push!(err_omega_plot_data, log10_Error_Omega)
    push!(err_w_plot_data, log10(max(L2_w_val, eps(Float64))))
    push!(err_phi_plot_data, log10(max(L2_phi_val, eps(Float64))))

    # ---- 🌟 滿足要求 4：VTK 雲圖導出至 vtk_dir ----
    vtu_path = joinpath(vtk_dir, "$(case_prefix).vtu")
    cells = [vtk_cell(elm) for elm in elements_Ω]; points = zeros(3, length(nodes))
    for node in nodes; points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0; end
    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        for (mode_rank, dm) in enumerate(dm_modes[1:min(length(dm_modes), 20)])
            vtk["Mode_$(mode_rank)_w"] = dm[2*nᵠ+1:end]
        end
    end

    # ---- 🌟 滿足要求 5：全頻主數據整合至指定的單一 CSV 命名形式中 ----
    csv_summary_path = joinpath(data_dir, "vibration_st_q_17_h_$(h)_modes.csv")
    open(csv_summary_path, "w") do io
        println(io, join(["mode_rank", "matched_m", "matched_n", "Omega_FEM", "Omega_exact", "log10_Error_Omega", "relative_residual"], ","))
        for r in 1:length(mode_ids)
            ω_r = sqrt(real(ω²[mode_ids[r]])); Ω_F = ω_r * a^2 * sqrt(ρ * h / Dᵇ)
            println(io, join([r, 1, 1, @sprintf("%.4f", Ω_F), @sprintf("%.4f", Ω_exact), @sprintf("%.6f", log10_Error_Omega), @sprintf("%.2e", residuals[r])], ","))
        end
    end
end

# 🌟 滿足要求 3：整個厚度大循環結束後，單獨導出「有限元基準線專用厚度收斂表」
h_convergence_csv = joinpath(data_dir, "vibration_thickness_fem_line.csv")
open(h_convergence_csv, "w") do io
    println(io, "h,log10_h,log10_Error_Omega,log10_Error_L2_w,log10_Error_L2_phi")
    for i in 1:length(h_plot_data)
        println(io, join([h_plot_data[i], log10(h_plot_data[i]), err_omega_plot_data[i], err_w_plot_data[i], err_phi_plot_data[i]], ","))
    end
end
println("\n[完全成功] 有限元基準線厚度逼近與鎖死大測試完工！總表導出至: ", h_convergence_csv)
gmsh.finalize()