# ==============================================================================
#  2D Mindlin Plate Exact Vibration Analysis with SSSS Boundary Conditions
#  流派：Navier Analytical Exact Solution (理論解析解連續域閉合特徵場)
#  🌟 最新修改：網格死鎖固定為 17，厚度 h 從 0.1 跑到 0.00001 幾何大循環
#  特性：分流儲存 CSV 至 date 資料夾，VTK 雲圖至 VTK 資料夾
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements

using WriteVTK, Printf
import Gmsh: gmsh

# ==================== 1. 參數與循環序列設定 ====================
a = 1.0   # 板長 (x 方向)
b = 1.0   # 板寬 (y 方向)
ρ = 1.0

# 🌟 滿足要求 1：網格死鎖固定使用 17 號網格
mesh_file = normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh"))
mesh_name = splitext(basename(mesh_file))[1]

# 🌟 滿足要求 2：H 自動設為幾何數列循環從 0.1 跑到 0.00001
h_series = [0.1, 0.01, 0.001, 0.0001, 0.00001]

# 🌟 滿足要求 4：宣告分流資料夾路徑
vtk_dir  = normpath(joinpath(@__DIR__, "..", "VTK"))
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
mkpath(vtk_dir); mkpath(data_dir)

function vtk_cell(elm)
    node_ids = [xᵢ.𝐼 for xᵢ in elm.𝓒]
    return length(node_ids) == 4 ? MeshCell(VTKCellTypes.VTK_QUAD, node_ids) : MeshCell(VTKCellTypes.VTK_TRIANGLE, node_ids)
end

# SSSS 薄板的 exact 模態空間波形
function exact_mode_shapes(x, y, m, n)
    wx = sin(m * π * x / a) * sin(n * π * y / b)
    φx = (m * π / a) * cos(m * π * x / a) * sin(n * π * y / b)
    φy = (n * π / b) * sin(m * π * x / a) * cos(n * π * y / b)
    return wx, φx, φy
end

# 生成 (m,n) 配對，按理論頻率從小到大嚴格排序
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

# ==================== 2. 主程式厚度大循環 ====================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

for h in h_series
    println("\n" * "="^60)
    println("固定網格 17 ➔ 生成理論解厚度分析：h = $(h)")
    println("="^60)

    # 🌟 隨著厚度改變，重新計算截面二次矩
    I_moment = h^3 / 12  

    gmsh.clear()
    gmsh.open(mesh_file)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = getElements(nodes, entities["Ω"], 1)

    cells = [vtk_cell(elm) for elm in elements]
    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for node in nodes
        points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0
    end

    mn_pairs = generate_mn_pairs(nₚ)

    # ---- A. 輸出不同厚度下的 VTU 圖形場（分流至 VTK 資料夾） ----
    n_vtu_output = min(length(mn_pairs), 20)
    vtu_path = joinpath(vtk_dir, "vibration_exact_st_q_17_h_$(h)_ex.vtu")
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

    # ---- B. 🌟 數據大整合：統一輸出至指定的對稱單一 CSV 文檔（分流至 date 資料夾） ----
    exact_summary_csv = joinpath(data_dir, "vibration_ex_st_q_17_h_$(h)_modes.csv")
    open(exact_summary_csv, "w") do io
        println(io, join(["mode_rank", "matched_m", "matched_n", "Omega_exact", "W_energy_pct", "Phi1_energy_pct", "Phi2_energy_pct"], ","))
        for (mode_idx, (m_val, n_val)) in enumerate(mn_pairs)
            Omega_ex = π^2 * (m_val^2 + n_val^2)
    
            # 【Navier 全域連續體動能閉合積分公式】
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
    println("  [CSV] 數據整合成功！理論全譜厚度文檔導出至: ", exact_summary_csv)
end
gmsh.finalize()
println("\n全部理論解析解厚度極限分析計算完成！")