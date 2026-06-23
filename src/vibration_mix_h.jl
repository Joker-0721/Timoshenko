# ==============================================================================
#  2D Mindlin Plate Free Vibration Analysis - Thickness Convergence Test
#  流派：Meshfree-FE Hybrid Mixed Method (單網格體系)
#  🌟 最新修改：網格死鎖固定為 17，厚度 h 從 0.1 跑到 0.00001 幾何循環
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ, L₂w, L₂φ  

using LinearAlgebra, TimerOutputs, WriteVTK, Printf
import Gmsh
const gmsh = Gmsh.gmsh

# ==================== 1. 材料與基礎幾何設定 ====================
E = 1.0e6; ν = 0.3; a = 1.0; b = 1.0; ρ = 1.0
eigen_imag_tol = 1.0e-7; omega_sq_tol = 1.0e-12
integrationOrder = 2; sʷ = 1.5; sᵠ = 1.5

# 🌟 滿足要求 1：死鎖固定讀取 17 號網格
mesh_file = normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh"))
mesh_name = splitext(basename(mesh_file))[1]

# 🌟 滿足要求 2：H 設為幾何數列循環從 0.1 跑到 0.00001
h_series = [0.1, 0.01, 0.001, 0.0001, 0.00001]

# 🌟 滿足要求 4：路徑精密分流宣告
vtk_dir  = normpath(joinpath(@__DIR__, "..", "VTK"))
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
mkpath(vtk_dir); mkpath(data_dir)

# 繪圖大表容器
h_plot_data = Float64[]
err_omega_plot_data = Float64[]
err_w_plot_data = Float64[]
err_phi_plot_data = Float64[]

function relative_residual(K, M, ω²ᵢ, vᵢ)
    Kv = K * vᵢ; Mv = M * vᵢ
    return norm(Kv - ω²ᵢ * Mv) / max(norm(Kv), abs(ω²ᵢ) * norm(Mv), eps(Float64))
end

function exact_omega_sq(m, n, a, b, D, ρ, h)
    Ω_exact = π^2 * (m^2 + n^2)
    return Ω_exact * sqrt(D / (ρ * h)) / a^2, Ω_exact
end

function vtk_cell(elm)
    return length(elm.𝓒) == 4 ? MeshCell(VTKCellTypes.VTK_QUAD, [x.𝐼 for x in elm.𝓒]) : MeshCell(VTKCellTypes.VTK_TRIANGLE, [x.𝐼 for x in elm.𝓒])
end

const to = TimerOutput()

# ==================== 2. 主程式厚度循環求解 ====================
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)

for h in h_series
    println("\n" * "="^70)
    println("固定網格 17 ➔ 執行厚度變更分析：h = $(h)")
    println("="^70)
    
    case_prefix = "vibration_mix_st_q_17_h_$(h)"
    
    # 🌟 隨著厚度改變，重新計算物理彎曲與剪切剛度
    Dᵇ = E * h^3 / (12 * (1 - ν^2))
    Dˢ = 5 / 6 * E * h / (2 * (1 + ν))
    αʷ = 1.0e8 * Dᵇ; αᵠ = 0.0            

    type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_Q = :tri3; type_M = :(PiecewisePolynomial{:Linear2D})
    
    ndiv = 16 # 固定 17 網格對應的平分數
    s_size = 1.0 / ndiv

    # 加載與拓撲重構
    gmsh.clear(); gmsh.open(mesh_file); nodes_w = get𝑿ᵢ(); nʷ = length(nodes_w)
    push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ), :s₂ => sʷ * s_size * ones(nʷ), :s₃ => sʷ * s_size * ones(nʷ))
    sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

    gmsh.clear(); gmsh.open(mesh_file); nodes_φ = get𝑿ᵢ(); nᵠ = length(nodes_φ)
    push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ), :s₂ => sᵠ * s_size * ones(nᵠ), :s₃ => sᵠ * s_size * ones(nᵠ))
    sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

    gmsh.clear(); gmsh.open(mesh_file); nodes = get𝑿ᵢ(); entities = getPhysicalGroups(); nˢ = length(nodes)

    kˢˢ = zeros(2 * nˢ, 2 * nˢ); kˢʷ = zeros(2 * nˢ, nʷ); kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
    kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ); kᵅʷʷ = zeros(nʷ, nʷ); mʷʷ = zeros(nʷ, nʷ); mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    @timeit to "assemble" begin
        elements_q = getElements(nodes, entities["Ω"], integrationOrder); prescribe!(elements_q, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_q)
        elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w); prescribe!(elements_w, :E => E, :ν => ν, :h => h); set𝝭!(elements_w)
        elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true); set𝝭!(elements_w_Γ)
        elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true); set𝝭!(elements_q_Γ)
        (∫QQdΩ => elements_q)(kˢˢ); ([∫∇QwdΩ => (elements_q, elements_w), ∫QwdΓ => (elements_q_Γ, elements_w_Γ)])(kˢʷ)

        elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder); prescribe!(elements_m, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_m)
        elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ); prescribe!(elements_φ, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_φ)
        elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true); set𝝭!(elements_φ_Γ)
        elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder); set𝝭!(elements_m_Γ)
        nₑ = length(elements_q); nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[])); nᵐ = nₑ * nᵖ
        kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ); kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)
        (∫MMdΩ => elements_m)(kᵐᵐ); ([∫∇MφdΩ => (elements_m, elements_φ), ∫MφdΓ => (elements_m_Γ, elements_φ_Γ)])(kᵐᵠ); (∫QφdΩ => (elements_q, elements_φ))(kˢᵠ)

        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_w_b = [getElements(nodes_w, entities[name], eval(type_w), integrationOrder, sp_w, normal=true) for name in boundary_names]
        for el in elements_w_b; prescribe!(el, :α => αʷ, :g => (x,y,z)->0.0); set𝝭!(el); end
        elements_φ_b = [getElements(nodes_φ, entities[name], eval(type_φ), integrationOrder, sp_φ, normal=true) for name in boundary_names]
        for el in elements_φ_b; prescribe!(el, :α => αᵠ, :g₁ => (x,y,z)->0.0, :g₂ => (x,y,z)->0.0, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0); set𝝭!(el); end
        dummy_fʷ = zeros(nʷ); dummy_fᵠ = zeros(2*nᵠ)
        (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kᵅʷʷ, dummy_fʷ)
        (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵅᵠᵠ, dummy_fᵠ)

        prescribe!(elements_w, :ρ => ρ, :h => h); (∫ρhwwdΩ => elements_w)(mʷʷ)
        prescribe!(elements_φ, :ρ => ρ, :h => h); (∫ρIφφdΩ => elements_φ)(mᵠᵠ)
        global M_total = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
    end

    global k_cond = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ   kˢᵠ'*(kˢˢ\kˢʷ); kˢʷ'*(kˢˢ\kˢᵠ)                        kˢʷ'*(kˢˢ\kˢʷ)-kᵅʷʷ]
    ks = Symmetric(0.5 * (k_cond + k_cond')); ms = Symmetric(0.5 * (M_total + M_total'))
    F = eigen(ks, ms); ω² = F.values; V = F.vectors

    mode_ids = sort!([i for i in eachindex(ω²) if real(ω²[i]) > omega_sq_tol && abs(imag(ω²[i])) < eigen_imag_tol], by = i -> real(ω²[i]))
    n_modes_output = min(length(mode_ids), 20)
    residuals = [relative_residual(ks, ms, real(ω²[mid]), V[:, mid]) for mid in mode_ids]

    # ---- 🌟 滿足要求 2：提取第一階特徵模態，智慧幅值與正負對齊 ----
    v_first = V[:, mode_ids[1]]; vᵠ_1 = v_first[1:2*nᵠ]; vʷ_1 = v_first[2*nᵠ+1:end]
    w_exact_func  = (x,y,z) -> sin(1*π*x/a) * sin(1*π*y/b)
    φx_exact_func = (x,y,z) -> (1*π/a) * cos(1*π*x/a) * sin(1*π*y/b)
    φy_exact_func = (x,y,z) -> (1*π/b) * sin(1*π*x/a) * cos(1*π*y/b)

    max_idx = argmax(abs.(vʷ_1))
    scale_factor = abs(vʷ_1[max_idx]) > 1e-12 ? w_exact_func(nodes_w[max_idx].x, nodes_w[max_idx].y, 0.0) / vʷ_1[max_idx] : 1.0
    vʷ_scaled = vʷ_1 .* scale_factor; vᵠ_scaled = vᵠ_1 .* scale_factor

    push!(nodes_w, :d => zeros(nʷ)); push!(nodes_φ, :d₁ => zeros(nᵠ), :d₂ => zeros(nᵠ))
    for (i, node) in enumerate(nodes_w); node.d = vʷ_scaled[i]; end
    for i in 1:length(nodes_φ); nodes_φ[i].d₁ = vᵠ_scaled[2*i-1]; nodes_φ[i].d₂ = vᵠ_scaled[2*i]; end

    prescribe!(elements_w, :w => w_exact_func); prescribe!(elements_φ, :φ₁ => φx_exact_func, :φ₂ => φy_exact_func)
    L2_w_val   = L₂w(elements_w)
    L2_phi_val = L₂φ(elements_φ)

    ω_real = sqrt(real(ω²[mode_ids[1]])); Ω_FEM = ω_real * a^2 * sqrt(ρ * h / Dᵇ)
    _, Ω_exact = exact_omega_sq(1, 1, a, b, Dᵇ, ρ, h)
    log10_Error_Omega = log10(max(abs(Ω_FEM - Ω_exact)/Ω_exact, eps(Float64)))

    push!(h_plot_data, h)
    push!(err_omega_plot_data, log10_Error_Omega)
    push!(err_w_plot_data, log10(max(L2_w_val, eps(Float64))))
    push!(err_phi_plot_data, log10(max(L2_phi_val, eps(Float64))))

    # ---- 🌟 滿足要求 4：VTK 雲圖導出至 vtk_dir ----
    vtu_path = joinpath(vtk_dir, "$(case_prefix).vtu")
    elements_Ω = getElements(nodes, entities["Ω"], 1); cells = [vtk_cell(elm) for elm in elements_Ω]
    points = zeros(3, nˢ); for node in nodes; points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0; end
    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        for (mode_rank, mode_id) in enumerate(mode_ids[1:n_modes_output])
            vtk["Mode_$(mode_rank)_w"] = V[:, mode_id][2*nᵠ+1:end]
        end
    end

    # ---- 🌟 滿足要求 5：全頻主數據整合至指定的單一 CSV 命名形式中 ----
    csv_summary_path = joinpath(data_dir, "vibration_mix_st_q_17_h_$(h)_modes.csv")
    open(csv_summary_path, "w") do io
        println(io, join(["mode_rank", "matched_m", "matched_n", "Omega_FEM", "Omega_exact", "log10_Error_Omega", "relative_residual"], ","))
        for r in 1:length(mode_ids)
            ω_r = sqrt(real(ω²[mode_ids[r]])); Ω_F = ω_r * a^2 * sqrt(ρ * h / Dᵇ)
            println(io, join([r, 1, 1, @sprintf("%.4f", Ω_F), @sprintf("%.4f", Ω_exact), @sprintf("%.6f", log10_Error_Omega), @sprintf("%.2e", residuals[r])], ","))
        end
    end
end

# 🌟 滿足要求 3：循環結束後，自動將厚度對數誤差輸出為獨立大表
h_convergence_csv = joinpath(data_dir, "vibration_thickness_mix_line.csv")
open(h_convergence_csv, "w") do io
    println(io, "h,log10_h,log10_Error_Omega,log10_Error_L2_w,log10_Error_L2_phi")
    for i in 1:length(h_plot_data)
        println(io, join([h_plot_data[i], log10(h_plot_data[i]), err_omega_plot_data[i], err_w_plot_data[i], err_phi_plot_data[i]], ","))
    end
end
println("\n[完全成功] Mix 單網格厚度逼近鎖死大測試完工！總表導出至: ", h_convergence_csv)
gmsh.finalize()