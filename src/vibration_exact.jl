# ==============================================================================
# 2D Mindlin plate exact vibration analysis with SSSS boundary conditions
# 🌟 最新修改：
#   1. CSV 數據強制分流儲存至 ../date 資料夾
#   2. VTU 雲圖強制分流儲存至 ../VTK 資料夾
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements

using LinearAlgebra
using WriteVTK, Printf
import Gmsh: gmsh

# ==================== 參數設定 ====================
a = 1.0
E = 1.0e6
ν = 0.3
b = 1.0   
h_over_b = 0.001
h = h_over_b * b
ρ = 1.0
I_moment = h^3 / 12  
Dᵇ = E * h^3 / (12 * (1 - ν^2))


mesh_files = [normpath(joinpath(@__DIR__, "..", "msh", "st_q_$i.msh")) for i in 15:25]

# 🌟 實施路徑精密分流
vtk_dir  = normpath(joinpath(@__DIR__, "..", "VTK"))
data_dir = normpath(joinpath(@__DIR__, "..", "date"))
mkpath(vtk_dir)
mkpath(data_dir)

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

function exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
    # 反推連續域的材料參數
    ν = 0.3
    E = Dᵇ * 12 * (1 - ν^2) / h^3
    Dˢ = (5/6) * E * h / (2 * (1 + ν))
    
    α_m = m * π / a
    β_n = n * π / b
    λ² = α_m^2 + β_n^2
    
    # 建立 Mindlin 理論精確的 3x3 剛度矩陣
    K_ex = zeros(3,3)
    K_ex[1,1] = Dˢ * λ²
    K_ex[1,2] = Dˢ * α_m
    K_ex[1,3] = Dˢ * β_n
    K_ex[2,1] = Dˢ * α_m
    K_ex[2,2] = Dᵇ * α_m^2 + Dᵇ * ((1-ν)/2) * β_n^2 + Dˢ
    K_ex[2,3] = Dᵇ * ((1+ν)/2) * α_m * β_n
    K_ex[3,1] = Dˢ * β_n
    K_ex[3,2] = K_ex[2,3]
    K_ex[3,3] = Dᵇ * β_n^2 + Dᵇ * ((1-ν)/2) * α_m^2 + Dˢ
    
    # 建立 Mindlin 理論精確的 3x3 質量矩陣 (含轉動慣量)
    M_ex = zeros(3,3)
    M_ex[1,1] = ρ * h
    M_ex[2,2] = ρ * h^3 / 12
    M_ex[3,3] = ρ * h^3 / 12
    
    # 直接求解這組 3x3 矩陣，得到真實的 Mindlin 理論頻率
    vals = eigvals(K_ex, M_ex)
    ω_exact = sqrt(minimum(real.(vals)))
    
    # 轉化為無因次頻率
    w_h_exact = ω_exact * a^2 * sqrt(ρ * h / Dᵇ)
    return w_h_exact, w_h_exact
end

function exact_mode_shapes(x, y, m, n)
    wx = sin(m * π * x / a) * sin(n * π * y / b)
    φx = (m * π / a) * cos(m * π * x / a) * sin(n * π * y / b)
    φy = (n * π / b) * sin(m * π * x / a) * cos(n * π * y / b)
    return wx, φx, φy
end

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
    println("生成 exact 處理解析與理論全譜數據：$(mesh_name)")
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

    mn_pairs = generate_mn_pairs(nₚ)
    println("  節點數: $nₚ, 理論計算全譜總模態數: $(length(mn_pairs))")

    # ---- 1. 輸出 VTU 圖形場 至 ../VTK ----
    n_vtu_output = min(length(mn_pairs), 20)
    vtu_path = joinpath(vtk_dir, "vibration_$(mesh_name)_ex.vtu")
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

    # ---- 2. 🌟 輸出純淨的理論全譜 CSV 至 ../date ----
    exact_csv = joinpath(data_dir, "vibration_$(mesh_name)_ex.csv")
    open(exact_csv, "w") do io
        println(io, join(["mode_rank", "m", "n", "w_h_exact"], ","))
        for (mode_idx, (m_val, n_val)) in enumerate(mn_pairs)
            _, w_h_ex = exact_omega_sq(m_val, n_val, a, b, Dᵇ, ρ, h)
            println(io, join([mode_idx, m_val, n_val, @sprintf("%.4f", w_h_ex)], ","))
        end
    end
    println("  [CSV] 純淨版理論全譜 CSV 成功導出: ", exact_csv)
end

gmsh.finalize()
println("\n全部純淨版理論解析解計算完成！")