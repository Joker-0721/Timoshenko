# 2D Mindlin plate exact vibration analysis with SSSS boundary conditions
# 生成解析解模态，输出 VTU 和 CSV (包含 (m,n) 信息)

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements

using WriteVTK
import Gmsh: gmsh

# ==================== 參數設定 ====================
a = 1.0   # 板長 (x 方向)
b = 1.0   # 板寬 (y 方向)

# 要分析的网格文件（与数值计算使用相同的 mesh 文件）
mesh_files = [
    normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh")),
    # normpath(joinpath(@__DIR__, "..", "msh", "st_t_17.msh"))
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

# SSSS 薄板的 exact 模态
function exact_mode_shapes(x, y, m, n)
    wx = sin(m * π * x / a) * sin(n * π * y / b)
    φx = (m * π / a) * cos(m * π * x / a) * sin(n * π * y / b)
    φy = (n * π / b) * sin(m * π * x / a) * cos(n * π * y / b)
    return wx, φx, φy
end

# 生成 (m,n) 配对，按 Ω = π²(m²+n²) 排序
function generate_mn_pairs(max_modes)
    pairs = Tuple{Int,Int}[]
    max_mn = ceil(Int, sqrt(max_modes)) + 2
    for total in 2:(2*max_mn)
        for m in 1:max_mn
            n = total - m
            if n >= 1 && n <= max_mn
                push!(pairs, (m, n))
                length(pairs) >= max_modes && break
            end
        end
        length(pairs) >= max_modes && break
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
    println("生成 exact 解：$(mesh_name)")
    println("="^60)

    # 读取网格（确保与数值计算使用相同的节点坐标）
    gmsh.clear()
    gmsh.open(mesh_file)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    # 获取域单元（仅用于输出 VTK 网格，形状无关紧要）
    elements = getElements(nodes, entities["Ω"], 1)

    # 准备 VTK 点与单元
    cells = [vtk_cell(elm) for elm in elements]
    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for node in nodes
        points[1, node.𝐼] = node.x
        points[2, node.𝐼] = node.y
        points[3, node.𝐼] = 0.0
    end

    # 模态数量 = 节点数（每个节点一个模态）
    n_modes = nₚ
    mn_pairs = generate_mn_pairs(n_modes)

    println("  节点数: $nₚ, 输出模态数: $(length(mn_pairs))")

    # ---- 输出 VTU ----
    vtu_path = joinpath(output_dir,
                        "vibration_$(mesh_name)_ex.vtu")
    vtk_grid(vtu_path, points, cells; ascii = true, append = false, compress = false) do vtk
        for (mode_idx, (m, n)) in enumerate(mn_pairs)
            w_exact = zeros(nₚ)
            φx_exact = zeros(nₚ)
            φy_exact = zeros(nₚ)
            for node in nodes
                i = node.𝐼
                w_exact[i], φx_exact[i], φy_exact[i] = exact_mode_shapes(node.x, node.y, m, n)
            end
            vtk["w$(mode_idx)"] = w_exact
            vtk["φx$(mode_idx)"] = φx_exact
            vtk["φy$(mode_idx)"] = φy_exact

            if mode_idx <= 10 || mode_idx == length(mn_pairs)
                println("  模态 #$mode_idx: (m,n)=($m,$n)")
            end
        end
        println("  ... 共 $(length(mn_pairs)) 个模态")
    end
    println("  VTU 输出: ", vtu_path)

    # ---- 输出 CSV（节点坐标 + 所有模态数据 + (m,n) 信息）----
    # 从 mesh_name 中提取网格尺寸（例如 struct_quad_17 -> 17）
    m = match(r"(\d+)", mesh_name)
    n_side = m === nothing ? "unknown" : m.captures[1]
    csv_path = joinpath(output_dir_csv,
                        "2d4s_vi_q4_ex_$(n_side)x$(n_side).csv")
    open(csv_path, "w") do io
        # 表头：node_id, x, y, m1, n1, w1, φx1, φy1, m2, n2, w2, φx2, φy2, ...
        header = ["node_id", "x", "y"]
        for k in 1:length(mn_pairs)
            push!(header, "m$k", "n$k", "w$k", "φx$k", "φy$k")
        end
        println(io, join(header, ","))

        # 每个节点的数据
        for node in nodes
            i = node.𝐼
            row = [string(i), string(node.x), string(node.y)]
            for (mode_idx, (m, n)) in enumerate(mn_pairs)
                w, φx, φy = exact_mode_shapes(node.x, node.y, m, n)
                push!(row, string(m), string(n), string(w), string(φx), string(φy))
            end
            println(io, join(row, ","))
        end
    end
    println("  CSV 输出: ", csv_path)
end

gmsh.finalize()
println("\n完成！")