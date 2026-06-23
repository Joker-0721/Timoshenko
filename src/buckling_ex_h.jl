# ==============================================================================
#  高級計算固體力學研究平台：Mindlin 屈曲分析理論解析解批量生成器
#  功能：固定 17 網格 + 厚度 h 自動化循環 (0.1 ➔ 0.00001) + 真實幾何座標對齊
#  對齊路徑規範：將 VTK 檔案完全投遞至指定的 ./VTK 資料夾下
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using WriteVTK
import Gmsh: gmsh

# 將檔案精確投遞至指定名稱的 VTK 資料夾
const VTK_DIR = "./VTK"
mkpath(VTK_DIR)

# 全局幾何邊長常數
const a = 1.0  
const b = 1.0  

# 定義解析解的位移與轉角連續函數場
w_exact_func(x, y, z) = sin(π * x / a) * sin(π * y / b)
phi1_exact_func(x, y, z) = (π / a) * cos(π * x / a) * sin(π * y / b)
phi2_exact_func(x, y, z) = (π / b) * sin(π * x / a) * cos(π * y / b)

println("="^80)
println(" 啟動 exact_buckling.jl 厚度系列解析解生成器 ")
println("  固定使用網格: /home/a/Joker/msh/st_q_17.msh")
println("="^80)

# 🚀 同步定義與主程式完全相同的厚度掃描序列
h_series = [0.1, 0.01, 0.001, 0.0001, 0.00001]

# 固定讀取 17 網格
mesh_path = "/home/a/Joker/msh/st_q_17.msh"

if !isfile(mesh_path)
    error("【錯誤】找不到指定的 17 階基礎網格檔案: ", mesh_path)
end

for h_val in h_series
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(mesh_path)
    
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = getElements(nodes, entities["Ω"], 2)
    
    # 準備 VTK 曲面資料結構
    n_points = length(nodes)
    points = zeros(3, n_points)
    w_analytical = zeros(n_points)
    phi1_analytical = zeros(n_points)
    phi2_analytical = zeros(n_points)
    
    # 填寫真實網格節點座標，並同時計算該點的精確理論解析解
    for (i, node) in enumerate(nodes)
        points[1, i] = node.x
        points[2, i] = node.y
        points[3, i] = 0.0
        
        # 帶入連續解析公式
        w_analytical[i]    = w_exact_func(node.x, node.y, 0.0)
        phi1_analytical[i] = phi1_exact_func(node.x, node.y, 0.0)
        phi2_analytical[i] = phi2_exact_func(node.x, node.y, 0.0)
    end
    
    # 建立結構化四邊形單元拓撲連接關係
    cells = [MeshCell(VTKCellTypes.VTK_QUAD, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]
    
    # 將理論場精確儲存在 VTK 資料夾下，並帶上厚度後綴
    vtu_path = joinpath(VTK_DIR, "exact_buckling_mesh_17_h_$(h_val).vtu")
    
    vtk_grid(vtu_path, points, cells; ascii=true) do vtk
        vtk["w_exact", VTKPointData()] = w_analytical
        vtk["phi1_exact", VTKPointData()] = phi1_analytical
        vtk["phi2_exact", VTKPointData()] = phi2_analytical
    end
    
    gmsh.finalize()
    println("  [成功] 厚度 h = $(h_val) 的解析解 VTU 已導出: ", vtu_path)
end

println("="^80)
println(" 【大成功】17網格下所有解析解對照場檔案（h 變量系列）已生成完畢 ")
println("="^80)