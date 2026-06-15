# 2D Mindlin plate exact vibration analysis with SSSS boundary conditions
# 升級：全譜 CSV 滿載填入、加入 Navier 連續域純理論 w, phi1, phi2 能量分析
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements

using WriteVTK, Printf
import Gmsh: gmsh

# ==================== 參數設定 ====================
a = 1.0   # 板長 (x 方向)
b = 1.0   # 板寬 (y 方向)
h_over_b = 0.001
h = h_over_b * b
ρ = 1.0
I_moment = h^3 / 12  # 板截面二次矩 (轉動慣量基底)

# 要分析的网格文件
mesh_files = [
    normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh")),
]

output_dir = normpath(joinpath(@__DIR__, "..", "vtk"))
mkpath(output_dir)
output_dir_csv = normpath(joinpath(@__DIR__, "..", "2d4s_vi_q4_FEM"))
mkpath(output_dir_csv)

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

# SSSS 薄板的 exact 模态空间波形
function exact_mode_shapes(x, y, m, n)
    wx = sin(m * π * x / a) * sin(n * π * y / b)
    φx = (m * π / a) * cos(m * π * x / a) * sin(n * π * y / b)
    φy = (n * π / b) * sin(m * π * x / a) * cos(n * π * y / b)
    return wx, φx, φy
end

# 生成 (m,n) 配对，按理论频率从小到大严格排序
function generate_mn_pairs(max_modes)
    pairs = Tuple{Int,Int}[]
    max_mn = ceil(Int, sqrt(max_modes)) + 2
    for total in 2:(2*max_mn)
        for m_harm in 1:max_mn
            n_harm = total - m_harm
            if n_harm >= 1 && n_harm <= max_mn
                push!(pairs, (m_harm, n_harm))
            end
        end
    end
    sort!(pairs, by = p -> p[1]^2 + p[2]^2)
    return pairs[1:min(max_modes, length(pairs))]
end

# ==================== 主程式 ====================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

for mesh_file in mesh_files
    mesh_name = splitext(basename(mesh_file))[1]
    println("\n" * "="^60)
    println("生成 exact 解與理論全譜能量解析：$(mesh_name)")
    println("="^60)

    gmsh.clear()
    gmsh.open(mesh_file)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = getElements(nodes, entities["Ω"], 1)

    cells = [vtk_cell(elm) for elm in elements]
    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for node in nodes; points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0; end

    # 輸出模態總數等於全域節點自由度基底
    mn_pairs = generate_mn_pairs(nₚ)
    println("  節點數: $nₚ, 理論計算全譜總模態數: $(length(mn_pairs))")

    # ---- 1. 輸出 VTU 圖形場（保持僅前 20 階，防止體積過大） ----
    n_vtu_output = min(length(mn_pairs), 20)
    vtu_path = joinpath(output_dir, "vibration_$(mesh_name)_ex.vtu")
    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        for mode_idx in 1:n_vtu_output
            m_val, n_val = mn_pairs[mode_idx]
            w_exact = zeros(nₚ); φx_exact = zeros(nₚ); φy_exact = zeros(nₚ)
            for node in nodes
                w_exact[node.𝐼], φx_exact[node.𝐼], φy_exact[node.𝐼] = exact_mode_shapes(node.x, node.y, m_val, n_val)
            end
            vtk["Mode_$(mode_idx)_w"]    = w_exact
            vtk["Mode_$(mode_idx)_phi1"] = φx_exact
            vtk["Mode_$(mode_idx)_phi2"] = φy_exact
        end
    end
    println("  [VTU] 理論圖形場輸出成功 (前 $n_vtu_output 階): ", vtu_path)

    # ---- 2. 🌟 新增：輸出理論全譜能量摘要 CSV (滿載全資料) ----
    exact_summary_csv = joinpath(output_dir, "vibration_$(mesh_name)_ex_summary.csv")
    open(exact_summary_csv, "w") do io
        println(io, join(["mode_rank", "m", "n", "Omega_exact", "W_energy_pct", "Phi1_energy_pct", "Phi2_energy_pct"], ","))
        for (mode_idx, (m_val, n_val)) in enumerate(mn_pairs)
            Omega_ex = π^2 * (m_val^2 + n_val^2)
            
            # 🌟 【Navier 全域連續體動能閉合積分公式】
            T_w = ρ * h
            T_𝜙x = ρ * I_moment * (m_val * π / a)^2
            T_𝜙y = ρ * I_moment * (n_val * π / b)^2
            T_total = T_w + T_𝜙x + T_𝜙y
            
            W_pct = (T_w / T_total) * 100
            Phi1_pct = (T_𝜙x / T_total) * 100
            Phi2_pct = (T_𝜙y / T_total) * 100
            
            println(io, join([mode_idx, m_val, n_val, @sprintf("%.4f", Omega_ex),
                              @sprintf("%.6f", W_pct), @sprintf("%.6f", Phi1_pct), @sprintf("%.6f", Phi2_pct)], ","))
        end
    end
    println("  [CSV] 理論全譜能量摘要 CSV 成功導出: ", exact_summary_csv)

    # ---- 3. 原本的空間節點場數據大 CSV (滿載) ----
    m_regex = match(r"(\d+)", mesh_name)
    n_side = m_regex === nothing ? "unknown" : m_regex.captures[1]
    csv_path = joinpath(output_dir_csv, "2d4s_vi_q4_ex_$(n_side)x$(n_side).csv")
    open(csv_path, "w") do io
        header = ["node_id", "x", "y"]
        for k in 1:length(mn_pairs); push!(header, "m$k", "n$k", "w$k", "φx$k", "φy$k"); end
        println(io, join(header, ","))
        for node in nodes
            row = [string(node.𝐼), string(node.x), string(node.y)]
            for (m_val, n_val) in mn_pairs
                w, φx, φy = exact_mode_shapes(node.x, node.y, m_val, n_val)
                push!(row, string(m_val), string(n_val), string(w), string(φx), string(φy))
            end
            println(io, join(row, ","))
        end
    end
    println("  [CSV] 空間節點場 CSV 輸出: ", csv_path)
end
gmsh.finalize()
println("\n全部理論解析解計算完成！")