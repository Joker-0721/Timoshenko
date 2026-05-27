# 2D Mindlin plate vibration analysis with SSSS boundary conditions using Q4 elements
# 新增：讀取解析解 CSV，透過 MAC 自動配對模態

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
E = 1.0e8
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.1
h = h_over_b * b
ρ = 1.0
α = 1.0e8 * E

# 薄板彎曲剛度
Dᵇ = E * h^3 / 12 / (1 - ν^2)

eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12
residual_warn_tol = 1.0e-6

# 要分析的網格檔案
mesh_files = [
    normpath(joinpath(@__DIR__, "..", "msh", "struct_quad_17.msh")),
    # normpath(joinpath(@__DIR__, "..", "msh", "struct_tri_17.msh"))
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
"""
     read_exact_csv(csv_path)
讀取解析解 CSV 檔案，回傳:
  - exact_nodes::Vector{Node} (自訂結構，含 x, y, id)
  - modes::Vector{ExactMode}，其中 ExactMode 包含 (m,n, w, φx, φy)
注意：假設 CSV 欄位為 node_id, x, y, w1, φx1, φy1, w2, φx2, φy2, ...
且模態數量由欄位數推得。
"""
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
    # 解析標題，取得模態資訊
    mode_cols = Dict{Int,Tuple{Int,Int}}()  # mode_index -> (m,n)
    # 標題格式範例: "w1", "φx1", "φy1", "w2", ...
    # 我們需要從 "w1" 中提取 1，但沒有儲存 m,n。實際 (m,n) 並未寫在 CSV 中。
    # 因此我們必須從外部知道 (m,n) 順序。替代方案：假設 CSV 中的模態順序與 generate_mn_pairs 一致，
    # 但由於沒有儲存 (m,n)，我們無法在此取得。因此建議修改 exact_vibration.jl 將 (m,n) 寫入 CSV 的額外欄位。
    # 為了讓此函數通用，我們改為期望 CSV 包含 "m1","n1","m2","n2",... 欄位。
    # 如果沒有，我們只能回傳模式索引，由外部提供 (m,n) 列表。
    # 這裡我們假設 CSV 中已經存在 "m1","n1","m2","n2",... 欄位（需修改 exact_vibration.jl）。
    # 若無，則回傳空的 modes，呼叫者改用頻率匹配。
    # 為了簡化，我們先檢查是否有 "m1" 欄位。
    if !any(startswith.(headers, "m"))
        @warn "CSV 缺少 (m,n) 資訊，將使用頻率匹配代替 MAC 配對。"
        return nothing, nothing
    end
    # 解析模態數量
    mode_indices = Set{Int}()
    for h in headers
        if startswith(h, "w")
            idx = parse(Int, h[2:end])
            push!(mode_indices, idx)
        end
    end
    n_modes = length(mode_indices)
    modes = ExactMode[]
    # 提取節點座標
    nodes = []
    for i in 1:n_rows
        x = values[i, 2]
        y = values[i, 3]
        push!(nodes, (id=i, x=x, y=y))
    end
    # 提取每個模態的 (m,n) 及向量
    for k in 1:n_modes
        # 找 m_k, n_k
        m_col = "m$k"
        n_col = "n$k"
        if !(m_col in headers) || !(n_col in headers)
            error("CSV 缺少 m$k 或 n$k 欄位")
        end
        m = Int(values[1, findfirst(==(m_col), headers)])
        n = Int(values[1, findfirst(==(n_col), headers)])
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
        push!(modes, ExactMode(m, n, w, φx, φy))
    end
    return nodes, modes
end

"""
    compute_mac(u1, u2, nᵠ)
計算兩個振型的 MAC 值。u1, u2 為完整自由度向量 [φx; φy; w]。
MAC = |u1' * u2|^2 / ((u1'*u1)*(u2'*u2))
"""
function compute_mac(u1::Vector{Float64}, u2::Vector{Float64}, nᵠ::Int)
    # u1, u2 長度應為 2*nᵠ + nʷ，但 nʷ = nᵠ (因為節點數相同)
    dot_ij = dot(u1, u2)
    norm_i = dot(u1, u1)
    norm_j = dot(u2, u2)
    if norm_i == 0.0 || norm_j == 0.0
        return 0.0
    end
    return (dot_ij^2) / (norm_i * norm_j)
end

"""
    align_exact_to_numerical_nodes(exact_nodes, exact_modes, numerical_nodes)
將解析解的模態向量重新排序，使其對應到數值節點的順序（基於座標匹配）。
回傳新的 ExactMode 列表，其中 w, φx, φy 已按 numerical_nodes 順序排列。
"""
function align_exact_to_numerical_nodes(exact_nodes, exact_modes, numerical_nodes)
    # 建立座標到索引的映射 (精確到小數點後 12 位)
    coord_to_idx = Dict{Tuple{Float64,Float64}, Int}()
    for (i, node) in enumerate(numerical_nodes)
        key = (round(node.x, digits=12), round(node.y, digits=12))
        coord_to_idx[key] = i
    end
    # 檢查每個 exact node 是否都能找到
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

# ---------------------- 主程式 ----------------------
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

all_results = Dict{String, Any}()

for mesh_file in mesh_files
    mesh_name = splitext(basename(mesh_file))[1]
    println("\n" * "="^70)
    println("分析：$(mesh_name) 網格")
    println("="^70)

    case_prefix = "mindlin_2d_ssss_$(mesh_name)_vibration"
    vtu_path = joinpath(output_dir, "$(case_prefix)_modes.vtu")

    integrationOrder = 3
    integrationOrder_shear = 2

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
            # 注意：∫ρhwwdΩ 和 ∫ρIφφdΩ 應返回矩陣，這裡假設它們可直接累加到 mʷʷ, mᵠᵠ
            # 如果它們是算子，需用 (∫ρhwwdΩ(ρ, h) => elements)(mʷʷ) 形式
            # 根據用戶之前提供的程式碼，應該是算子形式，但用戶寫成了賦值，需修正。
            # 這裡更正為正確的算子呼叫：
            (∫ρhwwdΩ(ρ, h) => elements)(mʷʷ)
            (∫ρIφφdΩ(ρ, h) => elements)(mᵠᵠ)
        end

        global elements_Ω = elements
    end

    @timeit to "apply SSSS boundary (penalty)" begin
        @timeit to "get boundary elements" begin
            elements_left   = getElements(nodes, entities["Γ¹"], integrationOrder)
            elements_right  = getElements(nodes, entities["Γ²"], integrationOrder)
            elements_top    = getElements(nodes, entities["Γ³"], integrationOrder)
            elements_bottom = getElements(nodes, entities["Γ⁴"], integrationOrder)
        end

        w_boundary(x, y) = 0.0
        prescribe!(elements_left,   :α => α, :g => w_boundary)
        prescribe!(elements_right,  :α => α, :g => w_boundary)
        prescribe!(elements_top,    :α => α, :g => w_boundary)
        prescribe!(elements_bottom, :α => α, :g => w_boundary)
        @timeit to "calculate boundary shape functions" begin
            set𝝭!(elements_left)
            set𝝭!(elements_right)
            set𝝭!(elements_top)
            set𝝭!(elements_bottom)
        end

        𝑎ʷ = ∫αwwdΓ => elements_left ∪ elements_right ∪ elements_top ∪ elements_bottom
        fᵅ = zeros(nʷ)
        @timeit to "assemble boundary penalty" 𝑎ʷ(kʷʷ, fᵅ)
    end

    @timeit to "solve vibration eigenvalue" begin
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        M = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
        Ks = Symmetric(0.5*(K + K'))
        Ms = Symmetric(0.5*(M + M'))
        F = eigen(Ks, Ms)
        ω² = F.values
        V = F.vectors

        mode_ids = sort!(
            collect(i for i in eachindex(ω²) if selected_eigenpair(ω²[i], V[:, i])),
            by=i -> real(ω²[i]),
        )
        isempty(mode_ids) && error("no positive finite vibration eigenvalue found")

        dm_modes = [normalize_mode(V[:, mode_id], nodes, nᵠ) for mode_id in mode_ids]
        residuals = [
            relative_residual(K, M, real(ω²[mode_id]), dm_modes[mode_rank])
            for (mode_rank, mode_id) in enumerate(mode_ids)
        ]

        ω²₁ = real(ω²[first(mode_ids)])
        ω₁ = sqrt(ω²₁)

        all_results[mesh_name] = (nodes=nodes, elements_Ω=elements_Ω, ω²=ω², V=V,
                                  mode_ids=mode_ids, dm_modes=dm_modes,
                                  residuals=residuals, nᵠ=nᵠ)
    end

    # ==================== 讀取解析解 CSV ====================
    # 根據 mesh_name 推測解析解檔案名稱（假設解析解已生成）
    # 例如 "struct_quad_17" 對應 "2d4s_vi_q4_ex_17x17.csv" (因為 exact_vibration.jl 使用 nodes_per_side)
    # 但 mesh_name 是 "struct_quad_17"，我們需要提取 17 -> "17x17"
    # 更通用：從 mesh_file 中提取數字部分。假設檔名為 struct_quad_17.msh，則數字為 17。
    m = match(r"(\d+)", mesh_name)
    if m === nothing
        @warn "無法從網格名稱提取節點數，跳過解析解比對"
        exact_nodes = nothing
        exact_modes = nothing
    else
        n_side = parse(Int, m.captures[1])
        exact_csv_path = joinpath(output_dir_csv, "2d4s_vi_q4_ex_$(n_side)x$(n_side).csv")
        if !isfile(exact_csv_path)
            @warn "解析解 CSV 檔案不存在: $exact_csv_path，跳過比對"
            exact_nodes, exact_modes = nothing, nothing
        else
            exact_nodes, exact_modes = read_exact_csv(exact_csv_path)
            if exact_modes !== nothing
                # 將解析解振型對齊到數值節點順序
                exact_modes = align_exact_to_numerical_nodes(exact_nodes, exact_modes, nodes)
                println("成功讀取解析解，共 $(length(exact_modes)) 個模態")
            end
        end
    end

    # ==================== 模態配對（MAC） ====================
    mac_matches = []  # 每個元素為 (mode_rank, matched_mode_index, mac_value, m, n)
    if exact_modes !== nothing && length(exact_modes) > 0
        n_modes_num = length(dm_modes)
        n_modes_exact = length(exact_modes)
        mac_matrix = zeros(n_modes_num, n_modes_exact)
        # 組裝數值模態完整向量 [φx; φy; w]
        num_vecs = []
        for dm in dm_modes
            φx = dm[1:2:2*nᵠ]
            φy = dm[2:2:2*nᵠ]
            w = dm[2*nᵠ+1:end]
            vec = vcat(φx, φy, w)
            push!(num_vecs, vec)
        end
        # 組裝解析模態完整向量
        exact_vecs = []
        for mode in exact_modes
            vec = vcat(mode.φx, mode.φy, mode.w)
            push!(exact_vecs, vec)
        end
        # 計算 MAC 矩陣
        for i in 1:n_modes_num
            for j in 1:n_modes_exact
                mac_matrix[i,j] = compute_mac(num_vecs[i], exact_vecs[j], nᵠ)
            end
        end
        # 為每個數值模態找到最佳匹配的解析模態
        for i in 1:n_modes_num
            best_j = argmax(mac_matrix[i,:])
            best_mac = mac_matrix[i, best_j]
            matched_mode = exact_modes[best_j]
            push!(mac_matches, (i, best_j, best_mac, matched_mode.m, matched_mode.n))
        end
    else
        # 若無解析解，則使用頻率匹配（原有方式，但需要 exact_omega_sq 函數）
        # 注意：原程式已有 exact_omega_sq，這裡保留作為備用
        for r in 1:length(dm_modes)
            ω²ᵢ = real(ω²[mode_ids[r]])
            ωᵢ = sqrt(ω²ᵢ)
            Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)
            best_err = Inf
            best_mn = (0,0)
            for m in 1:10
                for n in 1:10
                    _, Ω_ex = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
                    err = abs(Ω_FEM - Ω_ex) / Ω_ex * 100
                    if err < best_err
                        best_err = err
                        best_mn = (m,n)
                    end
                end
            end
            push!(mac_matches, (r, 0, 0.0, best_mn[1], best_mn[2]))  # MAC 設為 0
        end
    end

    # ==================== 輸出結果 ====================
    println("\n結果摘要 ($(mesh_name)):")
    println("  節點數: ", length(nodes))
    println("  選取的有效模態數: ", length(mode_ids))
    println("  ω²₁ (基頻平方): ", ω²₁)
    println("  ω₁  (基頻): ", ω₁)
    println("  最大相對殘差: ", maximum(residuals))
    if maximum(residuals) >= residual_warn_tol
        @printf("  警告: 部分殘差超過 %.1e\n", residual_warn_tol)
    end

    # ---- 模態比較表（使用 MAC 配對結果） ----
    println("\n  模態比較 (MAC 配對):")
    println("  " * @sprintf("%-8s %-8s %-12s %-12s %-12s", "mode", "(m,n)", "MAC", "Ω_FEM", "residual"))
    n_modes_show = min(length(dm_modes), 15)
    for r in 1:n_modes_show
        _, _, mac_val, m, n = mac_matches[r]
        ω²ᵢ = real(ω²[mode_ids[r]])
        ωᵢ = sqrt(ω²ᵢ)
        Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)
        println("  " * @sprintf("%-8d (%-d,%-d)     %-10.4e  %-10.4f  %-10.2e",
                                r, m, n, mac_val, Ω_FEM, residuals[r]))
    end

    # ---- 寫入 VTK（數值模態）----
    cells = [vtk_cell(elm) for elm in elements_Ω]
    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for node in nodes
        points[1, node.𝐼] = node.x
        points[2, node.𝐼] = node.y
        points[3, node.𝐼] = 0.0
    end
    n_modes_output = min(length(dm_modes), nₚ)

    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        ω²_selected = [real(ω²[mode_id]) for mode_id in mode_ids[1:n_modes_output]]
        Ω_selected = [sqrt(ω²ᵢ) * a^2 * sqrt(ρ * h / Dᵇ) for ω²ᵢ in ω²_selected]
        vtk["omega_sq", WriteVTK.VTKFieldData()] = ω²_selected
        vtk["omega", WriteVTK.VTKFieldData()] = sqrt.(ω²_selected)
        vtk["Omega_dimensionless", WriteVTK.VTKFieldData()] = Ω_selected
        vtk["relative_residual", WriteVTK.VTKFieldData()] = residuals[1:n_modes_output]
        vtk["rho", WriteVTK.VTKFieldData()] = [ρ]
        vtk["alpha_boundary_penalty", WriteVTK.VTKFieldData()] = [α]
        vtk["color_range_min", WriteVTK.VTKFieldData()] = [-1.0]
        vtk["color_range_max", WriteVTK.VTKFieldData()] = [1.0]

        for (mode_rank, dm) in enumerate(dm_modes[1:n_modes_output])
            w = dm[2*nᵠ+1:end]
            φx = dm[1:2:2*nᵠ]
            φy = dm[2:2:2*nᵠ]
            vtk["w$(mode_rank)"] = w
            vtk["φx$(mode_rank)"] = φx
            vtk["φy$(mode_rank)"] = φy
        end
        # 可選：將解析振型也寫入 VTK (需對齊節點順序)
        # 這裡省略以保持檔案簡潔
    end
    println("  VTU 輸出模態數: $n_modes_output")
    println("  VTU 輸出: ", vtu_path)

    # ---- 寫入 CSV 摘要（含配對資訊）----
    csv_summary_path = joinpath(output_dir, "$(case_prefix)_modes.csv")
    open(csv_summary_path, "w") do io
        println(io, join(["mode_rank", "eigen_index", "omega_sq", "omega", "Omega_dimensionless",
                          "matched_m", "matched_n", "MAC", "relative_residual", "residual_ok"], ","))
        for r in 1:length(dm_modes)
            mode_id = mode_ids[r]
            ω²ᵢ = real(ω²[mode_id])
            ωᵢ = sqrt(ω²ᵢ)
            Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)
            _, _, mac_val, m, n = mac_matches[r]
            println(io, join([r, mode_id, ω²ᵢ, ωᵢ, Ω_FEM, m, n, mac_val, residuals[r],
                              residuals[r] < residual_warn_tol], ","))
        end
    end
    println("  CSV 摘要輸出: ", csv_summary_path)

    # ---- 輸出節點數據 CSV（數值模態）----
    csv_nodal_path = joinpath(output_dir_csv, "2d4s_vi_FEM_$(mesh_name).csv")
    open(csv_nodal_path, "w") do io
        header = ["node_id", "x", "y"]
        for k in 1:n_modes_output
            push!(header, "w$k", "φx$k", "φy$k")
        end
        println(io, join(header, ","))
        for node in nodes
            i = node.𝐼
            row = [string(i), string(node.x), string(node.y)]
            for mode_rank in 1:n_modes_output
                dm = dm_modes[mode_rank]
                w = dm[2*nᵠ+1:end][i]
                φx = dm[1:2:2*nᵠ][i]
                φy = dm[2:2:2*nᵠ][i]
                push!(row, string(w), string(φx), string(φy))
            end
            println(io, join(row, ","))
        end
    end
    println("  節點數據 CSV 輸出: ", csv_nodal_path)
end

gmsh.finalize()

println("\n" * "="^70)
println("所有分析完成！")
println("="^70)
println(to)