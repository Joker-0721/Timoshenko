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

# 要分析的網格尺寸
mesh_sizes = [5, 17]

output_dir_vtk = normpath(joinpath(@__DIR__, "..", "vtk"))
mkpath(output_dir_vtk)
output_dir_csv = normpath(joinpath(@__DIR__, "..", "2d4s_vi_q4_ex"))
mkpath(output_dir_csv)

# ==================== 網格生成 ====================
function generate_square_mesh!(nodes_per_side)
    gmsh.clear()
    gmsh.model.add("square_mesh")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(a, 0.0, 0.0)
    p3 = gmsh.model.geo.addPoint(a, b, 0.0)
    p4 = gmsh.model.geo.addPoint(0.0, b, 0.0)

    bottom = gmsh.model.geo.addLine(p1, p2)
    right = gmsh.model.geo.addLine(p2, p3)
    top = gmsh.model.geo.addLine(p3, p4)
    left = gmsh.model.geo.addLine(p4, p1)
    loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(left, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(right, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteSurface(surface)
    gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()

    domain_group = gmsh.model.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, domain_group, "domain")

    gmsh.model.mesh.generate(2)
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

# SSSS 薄板的 exact 模態
# w_mn(x,y)   = sin(mπx/a) * sin(nπy/b)
# φx_mn(x,y)  = (mπ/a) * cos(mπx/a) * sin(nπy/b)
# φy_mn(x,y)  = (nπ/b) * sin(mπx/a) * cos(nπy/b)
function exact_mode_shapes(x, y, m, n)
    wx = sin(m * π * x / a) * sin(n * π * y / b)
    φx = (m * π / a) * cos(m * π * x / a) * sin(n * π * y / b)
    φy = (n * π / b) * sin(m * π * x / a) * cos(n * π * y / b)
    return wx, φx, φy
end

# 生成 (m,n) 配對，按 Ω = π²(m²+n²) 排序
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

for nodes_per_side in mesh_sizes
    println("\n" * "="^60)
    println("生成 exact VTU：$(nodes_per_side)×$(nodes_per_side) 網格")
    println("="^60)

    # 生成網格
    generate_square_mesh!(nodes_per_side)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = getElements(nodes, entities["domain"], 1)

    # 準備 VTK 點與單元
    cells = [vtk_cell(elm) for elm in elements]
    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for node in nodes
        points[1, node.𝐼] = node.x
        points[2, node.𝐼] = node.y
        points[3, node.𝐼] = 0.0
    end

    # 模態數量 = 節點數（每個節點一個模態）
    n_modes = nₚ
    mn_pairs = generate_mn_pairs(n_modes)

    println("  節點數: $nₚ, 輸出模態數: $(length(mn_pairs))")

    # ---- 輸出 VTU ----
    vtu_path = joinpath(output_dir_vtk,
                        "mindlin_2d_ssss_q4_$(nodes_per_side)x$(nodes_per_side)_exact_modes.vtu")
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
                println("  模態 #$mode_idx: (m,n)=($m,$n)")
            end
        end
        println("  ... 共 $(length(mn_pairs)) 個模態")
    end
    println("  VTU 輸出: ", vtu_path)

    # ---- 輸出 CSV（節點座標 + 所有模態數據）----
    csv_path = joinpath(output_dir_csv,
                        "2d4s_vi_q4_ex_$(nodes_per_side)x$(nodes_per_side).csv")
    open(csv_path, "w") do io
        # 表頭
        header = ["node_id", "x", "y"]
        for k in 1:length(mn_pairs)
            push!(header, "w$k", "φx$k", "φy$k")
        end
        println(io, join(header, ","))

        # 每個節點的數據
        for node in nodes
            i = node.𝐼
            row = [string(i), string(node.x), string(node.y)]
            for (mode_idx, (m, n)) in enumerate(mn_pairs)
                w, φx, φy = exact_mode_shapes(node.x, node.y, m, n)
                push!(row, string(w), string(φx), string(φy))
            end
            println(io, join(row, ","))
        end
    end
    println("  CSV 輸出: ", csv_path)
end

gmsh.finalize()
println("\n完成！")
