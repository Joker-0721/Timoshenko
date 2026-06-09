# 引入寫入 VTK 檔案所需的套件
using WriteVTK

# 1. 定義幾何與網格參數
a = 1.0  # 板長 (x 方向)
b = 1.0  # 板寬 (y 方向)
nx = 20  # x 方向的元素數量 (網格密度)
ny = 20  # y 方向的元素數量

# 產生網格節點座標 (x, y)
x_coords = range(0, a, length=nx+1)
y_coords = range(0, b, length=ny+1)

# 2. 定義解析解的位移函數
function analytical_mode_shape(x, y, a, b, m=1, n=1)
    A = 1.0 
    return A * sin(m * π * x / a) * sin(n * π * y / b)
end

# 3. 準備 VTK 網格資料結構
n_points = (nx + 1) * (ny + 1)
points = zeros(3, n_points)
w_analytical = zeros(n_points)

# 填寫節點座標，並同時計算該節點的解析解撓度 w
for j in 1:(ny + 1)
    for i in 1:(nx + 1)
        # 【導師修正】：直接利用 i, j 索引計算 point_id，避開全域變數作用域的問題
        point_id = (j - 1) * (nx + 1) + i
        
        x = x_coords[i]
        y = y_coords[j]
        
        points[1, point_id] = x
        points[2, point_id] = y
        points[3, point_id] = 0.0
        
        # 計算解析解的撓度
        w_analytical[point_id] = analytical_mode_shape(x, y, a, b, 1, 1)
    end
end

# 定義網格元素 (Cells) 的連接關係
cells = MeshCell[]
for j in 1:ny
    for i in 1:nx
        n1 = (j - 1) * (nx + 1) + i
        n2 = n1 + 1
        n3 = n2 + (nx + 1)
        n4 = n1 + (nx + 1)
        
        push!(cells, MeshCell(VTKCellTypes.VTK_QUAD, [n1, n2, n3, n4]))
    end
end

# 4. 輸出成 VTK 檔案
output_filename = "buckling_exact"
vtk_grid(output_filename, points, cells) do vtk
    vtk["w_exact", VTKPointData()] = w_analytical
end

println("解析解 VTU 檔案已成功生成: $(output_filename).vtu")