# ==============================================================================
#  2D Mindlin Plate Free Vibration Analysis - Pure Clean Log-Error Convergence
#  流派：Standard Displacement FEM Q4 (標準有限元縮減積分基準線)
#  🌟 最新終極修正：
#    1. 引入 3x3 Mindlin 連續域特徵矩陣，確保理論解頻譜起點完美對齊 1.0
#    2. 恢復 scale_factor 幾何對齊，確保 Log-Log 圖準確度
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫αwwdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ, L₂w, L₂φ

using LinearAlgebra
using Printf
using TimerOutputs
using WriteVTK
import Gmsh: gmsh

E = 1.0e6; ν = 0.3; a = 1.0; b = 1.0; h_over_b = 0.001; h = h_over_b * b; ρ = 1.0
Dᵇ = E * h^3 / 12 / (1 - ν^2)
α = 1.0e8 * Dᵇ

eigen_imag_tol = 1.0e-7; omega_sq_tol = 1.0e-12

mesh_files = [normpath(joinpath(@__DIR__, "..", "msh", "st_q_$i.msh")) for i in 9:25]

vtk_dir  = normpath(joinpath(@__DIR__, "..", "VTK")); data_dir = normpath(joinpath(@__DIR__, "..", "date"))
mkpath(vtk_dir); mkpath(data_dir)

h_plot_data = Float64[]; err_w_h_plot_data = Float64[]; err_w_plot_data = Float64[]; err_phi_plot_data = Float64[]

function selected_eigenpair(ω²ᵢ, vᵢ)
    return isfinite(real(ω²ᵢ)) && isfinite(imag(ω²ᵢ)) && abs(imag(ω²ᵢ)) < eigen_imag_tol && real(ω²ᵢ) > omega_sq_tol
end

# 🌟 真正的 Mindlin 3x3 矩陣解析解
function exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
    ν_ex = 0.3
    E_ex = Dᵇ * 12 * (1 - ν_ex^2) / h^3
    Dˢ_ex = (5/6) * E_ex * h / (2 * (1 + ν_ex))
    
    α_m = m * π / a; β_n = n * π / b; λ² = α_m^2 + β_n^2
    
    K_ex = zeros(3,3)
    K_ex[1,1] = Dˢ_ex * λ²; K_ex[1,2] = Dˢ_ex * α_m; K_ex[1,3] = Dˢ_ex * β_n
    K_ex[2,1] = Dˢ_ex * α_m; K_ex[2,2] = Dᵇ * α_m^2 + Dᵇ * ((1-ν_ex)/2) * β_n^2 + Dˢ_ex; K_ex[2,3] = Dᵇ * ((1+ν_ex)/2) * α_m * β_n
    K_ex[3,1] = Dˢ_ex * β_n; K_ex[3,2] = K_ex[2,3]; K_ex[3,3] = Dᵇ * β_n^2 + Dᵇ * ((1-ν_ex)/2) * α_m^2 + Dˢ_ex
    
    M_ex = zeros(3,3)
    M_ex[1,1] = ρ * h; M_ex[2,2] = ρ * h^3 / 12; M_ex[3,3] = ρ * h^3 / 12
    
    vals = eigvals(K_ex, M_ex)
    ω_exact = sqrt(minimum(real.(vals)))
    w_h_exact = ω_exact * a^2 * sqrt(ρ * h / Dᵇ)
    return w_h_exact, w_h_exact
end

function vtk_cell(elm)
    node_ids = [xᵢ.𝐼 for xᵢ in elm.𝓒]
    return length(node_ids) == 4 ? MeshCell(VTKCellTypes.VTK_QUAD, node_ids) : MeshCell(VTKCellTypes.VTK_TRIANGLE, node_ids)
end

const to = TimerOutput()
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)

for mesh_file in mesh_files
    mesh_name = splitext(basename(mesh_file))[1]
    println("\n" ^ 2 * "=" ^ 70 * "\n傳統有限元全譜對數分析：$(mesh_name)\n" * "=" ^ 70)
    
    case_prefix = "vibration_fem_$(mesh_name)"
    integrationOrder = 2; integrationOrder_shear = 1

    ndiv = match(r"(\d+)", mesh_name) === nothing ? 16 : parse(Int, match(r"(\d+)", mesh_name).captures[1]) - 1
    s_size = 1.0 / ndiv

    gmsh.clear(); gmsh.open(mesh_file); entities = getPhysicalGroups(); nodes = get𝑿ᵢ()
    nʷ = length(nodes); nᵠ = length(nodes)

    kʷʷ = zeros(nʷ, nʷ); kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ); kᵠʷ = zeros(2 * nᵠ, nʷ); mʷʷ = zeros(nʷ, nʷ); mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    @timeit to "assemble stiffness and mass" begin
        elements = getElements(nodes, entities["Ω"], integrationOrder)
        elements_s = getElements(nodes, entities["Ω"], integrationOrder_shear)
        prescribe!(elements, :E => E, :ν => ν, :h => h); prescribe!(elements_s, :E => E, :ν => ν, :h => h)
        set∇𝝭!(elements); set∇𝝭!(elements_s)

        (∫wwdΩ => elements_s)(kʷʷ); (∫φwdΩ => elements_s)(kᵠʷ); ([∫φφdΩ => elements_s, ∫κκdΩ => elements])(kᵠᵠ)
        prescribe!(elements, :ρ => ρ); (∫ρhwwdΩ => elements)(mʷʷ); (∫ρIφφdΩ => elements)(mᵠᵠ)
    end

    @timeit to "apply SSSS boundary" begin
        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        elements_b = [getElements(nodes, entities[name], integrationOrder) for name in boundary_names]
        for el in elements_b; prescribe!(el, :α => α, :g => (x,y,z)->0.0); set𝝭!(el); end
        dummy_f = zeros(nʷ); (∫αwwdΓ => elements_b[1] ∪ elements_b[2] ∪ elements_b[3] ∪ elements_b[4])(kʷʷ, dummy_f)
    end

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]; M = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
    Ks = Symmetric(0.5*(K + K')); Ms = Symmetric(0.5*(M + M'))
    F = eigen(Ks, Ms); ω² = F.values; V = F.vectors

    mode_ids = sort!([i for i in eachindex(ω²) if selected_eigenpair(ω²[i], V[:, i])], by = i -> real(ω²[i]))
    n_modes_output = min(length(mode_ids), 20)

    modal_post_records = []
    push!(nodes, :d => zeros(nʷ), :d₁ => zeros(nᵠ), :d₂ => zeros(nᵠ))

    for r in 1:length(mode_ids)
        mode_id = mode_ids[r]; v_full = V[:, mode_id]
        vᵠ = v_full[1:2*nᵠ]; vʷ = v_full[2*nᵠ+1:end]

        ω_real_m = sqrt(real(ω²[mode_id])); w_h_FEM_m = ω_real_m * a^2 * sqrt(ρ * h / Dᵇ)

        best_err = Inf; best_mn = (1, 1)
        for m_harm in 1:15, n_harm in 1:15
            _, w_h_ex = exact_omega_sq(m_harm, n_harm, a, b, Dᵇ, ρ, h)
            if abs(w_h_FEM_m - w_h_ex) < best_err
                best_err = abs(w_h_FEM_m - w_h_ex); best_mn = (m_harm, n_harm)
            end
        end
        m_match, n_match = best_mn
        _, w_h_exact_m = exact_omega_sq(m_match, n_match, a, b, Dᵇ, ρ, h)
        log10_Error_w_h = log10(max(abs(w_h_FEM_m - w_h_exact_m)/w_h_exact_m, eps(Float64)))

        # 🌟 恢復幾何振幅對齊 (scale_factor)
        w_ex_func  = (x,y,z) -> sin(m_match * π * x / a) * sin(n_match * π * y / b)
        φx_ex_func = (x,y,z) -> (m_match * π / a) * cos(m_match * π * x / a) * sin(n_match * π * y / b)
        φy_ex_func = (x,y,z) -> (n_match * π / b) * sin(m_match * π * x / a) * cos(n_match * π * y / b)

        max_idx = argmax(abs.(vʷ))
        scale_factor = abs(vʷ[max_idx]) > 1e-12 ? w_ex_func(nodes[max_idx].x, nodes[max_idx].y, 0.0) / vʷ[max_idx] : 1.0
        vʷ_scaled = vʷ .* scale_factor; vᵠ_scaled = vᵠ .* scale_factor

        for (i, node) in enumerate(nodes)
            node.d = vʷ_scaled[i]; node.d₁ = vᵠ_scaled[2*i-1]; node.d₂ = vᵠ_scaled[2*i]
        end

        prescribe!(elements, :w => w_ex_func); prescribe!(elements, :φ₁ => φx_ex_func); prescribe!(elements, :φ₂ => φy_ex_func)

        log10_Error_L2_w   = log10(max(L₂w(elements), eps(Float64)))
        log10_Error_L2_phi = log10(max(L₂φ(elements), eps(Float64)))

        push!(modal_post_records, (m_match, n_match, w_h_exact_m, log10_Error_w_h, log10_Error_L2_w, log10_Error_L2_phi))
    end

    vtu_path = joinpath(vtk_dir, "$(case_prefix).vtu")
    cells = [vtk_cell(elm) for elm in elements]; points = zeros(3, length(nodes))
    for node in nodes; points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0; end
    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        for (mode_rank, mode_id) in enumerate(mode_ids[1:n_modes_output])
            vtk["Mode_$(mode_rank)_w"] = V[:, mode_id][2*nᵠ+1:end]
        end
    end

    csv_summary_path = joinpath(data_dir, "$(case_prefix)_log_convergence.csv")
    open(csv_summary_path, "w") do io
        println(io, join(["mode_rank", "matched_m", "matched_n", "w_h_FEM", "w_h_exact", "log10_Error_w_h", "log10_Error_L2_w", "log10_Error_L2_phi"], ","))
        for r in 1:length(mode_ids)
            m, n, w_h_ex, log10_err_wh, log10_l2_w, log10_l2_φ = modal_post_records[r]
            ω_real = sqrt(real(ω²[mode_ids[r]])); w_h_F = ω_real * a^2 * sqrt(ρ * h / Dᵇ)
            println(io, join([r, m, n, @sprintf("%.4f", w_h_F), @sprintf("%.4f", w_h_ex), @sprintf("%.6f", log10_err_wh), @sprintf("%.6f", log10_l2_w), @sprintf("%.6f", log10_l2_φ)], ","))
        end
    end

    m_1, n_1, w_h_ex_1, log10_err_wh_1, log10_l2_w_1, log10_l2_φ_1 = modal_post_records[1]
    push!(h_plot_data, s_size); push!(err_w_h_plot_data, log10_err_wh_1); push!(err_w_plot_data, log10_l2_w_1); push!(err_phi_plot_data, log10_l2_φ_1)
end

h_convergence_csv = joinpath(data_dir, "vibration_FEM_h_line.csv")
open(h_convergence_csv, "w") do io
    println(io, "h,log10_h,log10_Error_w_h,log10_Error_L2_w,log10_Error_L2_phi")
    for i in 1:length(h_plot_data)
        println(io, join([h_plot_data[i], log10(h_plot_data[i]), err_w_h_plot_data[i], err_w_plot_data[i], err_phi_plot_data[i]], ","))
    end
end
println("\n[完全成功] 傳統有限元(FEM)收斂率大閉環！純淨版收斂折線已成功導出至: ", h_convergence_csv)
gmsh.finalize()