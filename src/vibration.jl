# 2D Mindlin plate vibration analysis with SSSS boundary conditions using Q4 elements
# 升級：全譜 CSV 填滿、VTK 保持前 20 階，並加入 w, phi1, phi2 能量定量分流拆解
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫αwwdΓ, ∫ρhwwdΩ, ∫ρIφφdΩ

using LinearAlgebra
using Printf
using TimerOutputs
using WriteVTK
using DelimitedFiles
import Gmsh: gmsh

# ==================== 參數設定 ====================
E = 1.0e6
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.001
h = h_over_b * b
ρ = 1.0

# 薄板彎曲剛度
Dᵇ = E * h^3 / 12 / (1 - ν^2)
α = 1.0e8 * Dᵇ

eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12
residual_warn_tol = 1.0e-6

# 要分析的網格檔案
mesh_files = [
    normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh")),
]

output_dir = normpath(joinpath(@__DIR__, "..", "vtk"))
mkpath(output_dir)
output_dir_csv = normpath(joinpath(@__DIR__, "..", "2d4s_vi_q4_FEM"))
mkpath(output_dir_csv)

const to = TimerOutput()

# ---------------------- 篩選與後處理函數 ----------------------
function selected_eigenpair(ω²ᵢ, vᵢ)
    return isfinite(real(ω²ᵢ)) &&
           isfinite(imag(ω²ᵢ)) &&
           abs(imag(ω²ᵢ)) < eigen_imag_tol &&
           real(ω²ᵢ) > omega_sq_tol &&
           all(isfinite, real.(vᵢ)) &&
           all(isfinite, imag.(vᵢ)) &&
           norm(imag.(vᵢ)) <= eigen_imag_tol * max(norm(real.(vᵢ)), eps(Float64))
end

function center_node_id(nodes)
    distances = [(node.x - a / 2)^2 + (node.y - b / 2)^2 for node in nodes]
    return nodes[argmin(distances)].𝐼
end

function normalize_mode(vᵢ, nodes, nᵠ)
    dm = real.(vᵢ)
    w = dm[2*nᵠ+1:end]
    max_abs_w = maximum(abs.(w))
    max_abs_w > 0.0 || error("mode has zero transverse displacement.")
    dm ./= max_abs_w

    center_id = center_node_id(nodes)
    if dm[2*nᵠ + center_id] < -sqrt(eps(Float64))
        dm .*= -1.0
    end
    return dm
end

function relative_residual(K, M, ω²ᵢ, vᵢ)
    Kv = K * vᵢ
    Mv = M * vᵢ
    residual = Kv - ω²ᵢ * Mv
    denominator = max(norm(Kv), abs(ω²ᵢ) * norm(Mv), eps(Float64))
    return norm(residual) / denominator
end

function vtk_cell(elm)
    node_ids = [xᵢ.𝐼 for xᵢ in elm.𝓒]
    if length(node_ids) == 4
        return MeshCell(VTKCellTypes.VTK_QUAD, node_ids)
    elseif length(node_ids) == 3
        return MeshCell(VTKCellTypes.VTK_TRIANGLE, node_ids)
    else
        error("unsupported cell with $(length(node_ids)) nodes")
    end
end

# ---------------------- 解析解讀取與 MAC 計算 ----------------------
struct ExactMode
    m::Int
    n::Int
    w::Vector{Float64}
    φx::Vector{Float64}
    φy::Vector{Float64}
end

function read_exact_csv(csv_path)
    data = readdlm(csv_path, ',', header=true)
    headers = data[2][:]
    values = data[1]
    n_rows = size(values, 1)
    if !any(startswith.(headers, "m"))
        @warn "CSV 缺少 (m,n) 資訊，將使用頻率匹配代替 MAC 配對。"
        return nothing, nothing
    end
    mode_indices = Set{Int}()
    for h_field in headers
        if startswith(h_field, "w")
            idx = parse(Int, h_field[2:end])
            push!(mode_indices, idx)
        end
    end
    n_modes = length(mode_indices)
    modes = ExactMode[]
    nodes = []
    for i in 1:n_rows
        x = values[i, 2]
        y = values[i, 3]
        push!(nodes, (id=i, x=x, y=y))
    end
    for k in 1:n_modes
        m_col = "m$k"
        n_col = "n$k"
        if !(m_col in headers) || !(n_col in headers)
            error("CSV 缺少 $m_col 或 $n_col 欄位")
        end
        m_val = Int(values[1, findfirst(==(m_col), headers)])
        n_val = Int(values[1, findfirst(==(n_col), headers)])
        w = Float64[]
        φx = Float64[]
        φy = Float64[]
        w_col = "w$k"
        φx_col = "φx$k"
        φy_col = "φy$k"
        for i in 1:n_rows
            push!(w, values[i, findfirst(==(w_col), headers)])
            push!(φx, values[i, findfirst(==(φx_col), headers)])
            push!(φy, values[i, findfirst(==(φy_col), headers)])
        end
        push!(modes, ExactMode(m_val, n_val, w, φx, φy))
    end
    return nodes, modes
end

function compute_mac(u1::Vector{Float64}, u2::Vector{Float64}, nᵠ::Int)
    dot_ij = dot(u1, u2)
    norm_i = dot(u1, u1)
    norm_j = dot(u2, u2)
    if norm_i == 0.0 || norm_j == 0.0
        return 0.0
    end
    return (dot_ij^2) / (norm_i * norm_j)
end

function align_exact_to_numerical_nodes(exact_nodes, exact_modes, numerical_nodes)
    coord_to_idx = Dict{Tuple{Float64,Float64}, Int}()
    for (i, node) in enumerate(numerical_nodes)
        key = (round(node.x, digits=12), round(node.y, digits=12))
        coord_to_idx[key] = i
    end
    new_modes = []
    for mode in exact_modes
        w_new = zeros(length(numerical_nodes))
        φx_new = zeros(length(numerical_nodes))
        φy_new = zeros(length(numerical_nodes))
        for (j, exact_node) in enumerate(exact_nodes)
            key = (round(exact_node.x, digits=12), round(exact_node.y, digits=12))
            if haskey(coord_to_idx, key)
                idx = coord_to_idx[key]
                w_new[idx] = mode.w[j]
                φx_new[idx] = mode.φx[j]
                φy_new[idx] = mode.φy[j]
            else
                error("節點座標 ($(exact_node.x), $(exact_node.y)) 在數值網格中找不到")
            end
        end
        push!(new_modes, ExactMode(mode.m, mode.n, w_new, φx_new, φy_new))
    end
    return new_modes
end

function exact_omega_sq(m, n, a, b, D, ρ, h)
    Ω_exact = π^2 * (m^2 + n^2)
    ω_exact = Ω_exact * sqrt(D / (ρ * h)) / a^2
    return ω_exact, Ω_exact
end

# ==================== 主程式 ====================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
all_results = Dict{String, Any}()

for mesh_file in mesh_files
    mesh_name = splitext(basename(mesh_file))[1]
    println("\n" * "="^70)
    println("分析：$(mesh_name) 網格")
    println("="^70)

    case_prefix = "vibration_0001$(mesh_name)"
    integrationOrder = 2
    integrationOrder_shear = 1

    gmsh.clear()
    @timeit to "open msh file" gmsh.open(mesh_file)
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
    kᵠʷ = zeros(2 * nᵠ, nʷ)
    mʷʷ = zeros(nʷ, nʷ)
    mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    @timeit to "assemble stiffness and mass" begin
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"], integrationOrder)
        @timeit to "get shear elements" elements_s = getElements(nodes, entities["Ω"], integrationOrder_shear)
        prescribe!(elements, :E => E, :ν => ν, :h => h)
        prescribe!(elements_s, :E => E, :ν => ν, :h => h)
        @timeit to "calculate shape functions" set∇𝝭!(elements)
        @timeit to "calculate shear shape functions" set∇𝝭!(elements_s)

        𝑎ʷʷ = ∫wwdΩ => elements_s
        𝑎ᵠʷ = ∫φwdΩ => elements_s
        𝑎ᵠᵠ = [∫φφdΩ => elements_s, ∫κκdΩ => elements]
        @timeit to "assemble Kww" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble Kφw" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble Kφφ" 𝑎ᵠᵠ(kᵠᵠ)

        @timeit to "assemble consistent mass" begin
            prescribe!(elements, :ρ => ρ)
            (∫ρhwwdΩ => elements)(mʷʷ)
            (∫ρIφφdΩ => elements)(mᵠᵠ)
        end
        global elements_Ω = elements
    end

    @timeit to "apply SSSS boundary (penalty)" begin
        elements_left   = getElements(nodes, entities["Γ¹"], integrationOrder)
        elements_right  = getElements(nodes, entities["Γ²"], integrationOrder)
        elements_top    = getElements(nodes, entities["Γ³"], integrationOrder)
        elements_bottom = getElements(nodes, entities["Γ⁴"], integrationOrder)

        w_boundary(x, y, z) = 0.0
        prescribe!(elements_left,   :α => α, :g => w_boundary)
        prescribe!(elements_right,  :α => α, :g => w_boundary)
        prescribe!(elements_top,    :α => α, :g => w_boundary)
        prescribe!(elements_bottom, :α => α, :g => w_boundary)
        set𝝭!(elements_left); set𝝭!(elements_right); set𝝭!(elements_top); set𝝭!(elements_bottom)

        𝑎ʷ = ∫αwwdΓ => elements_left ∪ elements_right ∪ elements_top ∪ elements_bottom
        fᵅ = zeros(nʷ)
        𝑎ʷ(kʷʷ, fᵅ)
    end

    @timeit to "solve vibration eigenvalue" begin
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        M = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
        Ks = Symmetric(0.5*(K + K'))
        Ms = Symmetric(0.5*(M + M'))
        F = eigen(Ks, Ms)
        ω² = F.values
        V = F.vectors

        # 篩選有效模態
        global mode_ids = sort!(
            collect(i for i in eachindex(ω²) if selected_eigenpair(ω²[i], V[:, i])),
            by = i -> real(ω²[i]),
        )
        dm_modes = [normalize_mode(V[:, mode_id], nodes, nᵠ) for mode_id in mode_ids]
        residuals = [relative_residual(K, M, real(ω²[mode_id]), dm_modes[mode_rank]) for (mode_rank, mode_id) in enumerate(mode_ids)]

        ω²₁ = real(ω²[first(mode_ids)])
        ω₁ = sqrt(ω²₁)

        all_results[mesh_name] = (nodes=nodes, elements_Ω=elements_Ω, ω²=ω², V=V, mode_ids=mode_ids, dm_modes=dm_modes, residuals=residuals, nᵠ=nᵠ)
    end

    # ==================== 讀取與對齊解析解 ====================
    m_regex = match(r"(\d+)", mesh_name)
    n_side = m_regex === nothing ? "unknown" : m_regex.captures[1]
    exact_csv_path = joinpath(output_dir_csv, "2d4s_vi_q4_ex_$(n_side)x$(n_side).csv")
    
    mac_matches = []
    if isfile(exact_csv_path)
        exact_nodes, exact_modes = read_exact_csv(exact_csv_path)
        if exact_modes !== nothing
            aligned_modes = align_exact_to_numerical_nodes(exact_nodes, exact_modes, nodes)
            for (r, dm) in enumerate(dm_modes)
                best_mac = -1.0
                best_mn = (0, 0)
                for em in aligned_modes
                    v_ex = [em.φx; em.φy; em.w]  # 對齊自由度排序
                    mac_val = compute_mac(dm, v_ex, nᵠ)
                    if mac_val > best_mac
                        best_mac = mac_val
                        best_mn = (em.m, em.n)
                    end
                end
                push!(mac_matches, (r, 0, best_mac, best_mn[1], best_mn[2]))
            end
        end
    else
        # 若無解析解 CSV，使用頻率逼近匹配
        for r in 1:length(dm_modes)
            ωᵢ = sqrt(real(ω²[mode_ids[r]]))
            Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)
            best_err = Inf; best_mn = (0,0)
            for m_harm in 1:15, n_harm in 1:15
                _, Ω_ex = exact_omega_sq(m_harm, n_harm, a, b, Dᵇ, ρ, h)
                err = abs(Ω_FEM - Ω_ex) / Ω_ex * 100
                if err < best_err; best_err = err; best_mn = (m_harm, n_harm); end
            end
            push!(mac_matches, (r, 0, 1.0, best_mn[1], best_mn[2]))
        end
    end

    # 🌟 圖片保持：僅將前 20 階導出為 VTK 雲圖
    n_modes_output = min(length(dm_modes), 20)
    vtu_path = joinpath(output_dir, "$(case_prefix).vtu")
    cells = [vtk_cell(elm) for elm in elements_Ω]
    points = zeros(3, length(nodes))
    for node in nodes; points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0; end

    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        ω²_selected = [real(ω²[mid]) for mid in mode_ids[1:n_modes_output]]
        vtk["Omega_dimensionless", WriteVTK.VTKFieldData()] = [sqrt(w2) * a^2 * sqrt(ρ * h / Dᵇ) for w2 in ω²_selected]
        for (mode_rank, dm) in enumerate(dm_modes[1:n_modes_output])
            vtk["Mode_$(mode_rank)_w"]    = dm[2*nᵠ+1:end]
            vtk["Mode_$(mode_rank)_phi1"] = dm[1:2:2*nᵠ]
            vtk["Mode_$(mode_rank)_phi2"] = dm[2:2:2*nᵠ]
        end
    end
    println("  [VTK] 圖形場輸出成功 (前 $n_modes_output 階): ", vtu_path)

    # 🌟 數據滿載：全頻譜所有有效模態完整寫入 CSV 摘要，並附帶量化能量拆解
    csv_summary_path = joinpath(output_dir, "$(case_prefix)_modes.csv")
    open(csv_summary_path, "w") do io
        println(io, join(["mode_rank", "eigen_index", "omega_sq", "omega",
                          "Omega_FEM", "Omega_exact", "matched_m", "matched_n",
                          "MAC", "W_energy_pct", "Phi1_energy_pct", "Phi2_energy_pct",
                          "relative_residual", "residual_ok"], ","))
        for r in 1:length(dm_modes)
            mode_id = mode_ids[r]
            ω²ᵢ = real(ω²[mode_id])
            ωᵢ = sqrt(ω²ᵢ)
            Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)
            _, _, mac_val, m, n = mac_matches[r]
            Ω_exact = (m != 0 && n != 0) ? π^2 * (m^2 + n^2) : NaN
            
            # 🌟 【有限元動能物理拆解核心】
            dm = dm_modes[r]
            vᵠ = dm[1:2*nᵠ]
            vʷ = dm[2*nᵠ+1:end]
            E_kinetic_w = vʷ' * mʷʷ * vʷ
            vᵠ₁ = vᵠ[1:2:end]; vᵠ₂ = vᵠ[2:2:end]
            mᵠ₁ᵠ₁ = mᵠᵠ[1:2:end, 1:2:end]; mᵠ₂ᵠ₂ = mᵠᵠ[2:2:end, 2:2:end]
            E_kinetic_𝜙1 = vᵠ₁' * mᵠ₁ᵠ₁ * vᵠ₁
            E_kinetic_𝜙2 = vᵠ₂' * mᵠ₂ᵠ₂ * vᵠ₂
            
            E_total = E_kinetic_w + E_kinetic_𝜙1 + E_kinetic_𝜙2
            W_pct = (E_kinetic_w / E_total) * 100
            Phi1_pct = (E_kinetic_𝜙1 / E_total) * 100
            Phi2_pct = (E_kinetic_𝜙2 / E_total) * 100

            println(io, join([r, mode_id, ω²ᵢ, ωᵢ, Ω_FEM, Ω_exact, m, n, mac_val,
                              @sprintf("%.4f", W_pct), @sprintf("%.4f", Phi1_pct), @sprintf("%.4f", Phi2_pct),
                              residuals[r], residuals[r] < residual_warn_tol], ","))
        end
    end
    println("  [CSV] 有限元全資料數據摘要（含變量能量拆解）成功導出: ", csv_summary_path)
end
gmsh.finalize()
println("\n所有分析完成！")