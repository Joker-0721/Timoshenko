using Gmsh

function generate_cylinder_mesh(R=1.0, L=2.0, t=0.01; mesh_size=0.05)
    gmsh.initialize()
    gmsh.model.add("cylinder_shell")
    
    # 引入 OpenCASCADE 幾何引擎
    occ = gmsh.model.occ
    
    # 1. 【核心修正】：先建立 3D 空間中的幾何點，取得它們的 Int32 標籤 (Tags)
    p1 = occ.addPoint(R, 0, 0)
    p_center = occ.addPoint(0, 0, 0)
    p2 = occ.addPoint(-R, 0, 0)
    
    # 2. 利用點的標籤建立圓弧線段（起點、圓心、終點）
    c1 = occ.addCircleArc(p1, p_center, p2)
    c2 = occ.addCircleArc(p2, p_center, p1)
    
    # 沿著 Z 軸（軸向）擠壓長度 L，生成兩個半圓柱面
    extru1 = occ.extrude([(1, c1)], 0, 0, L)
    extru2 = occ.extrude([(1, c2)], 0, 0, L)
    
    occ.synchronize()
    
    # 3. 【優化安全機制】：從擠壓返回的元件中，過濾出 dim == 2 的拓撲結構，精確抓取曲面 ID
    surf1_id = [pair[2] for pair in extru1 if pair[1] == 2][1]
    surf2_id = [pair[2] for pair in extru2 if pair[1] == 2][1]
    
    # 4. 強制網格重新組裝成四邊形 (Q4 Elements)
    gmsh.model.mesh.setRecombine(2, surf1_id)
    gmsh.model.mesh.setRecombine(2, surf2_id)
    
    # 設定網格特徵尺寸
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
    
    # 5. 定義邊界條件與物理群組 (Physical Groups) —— 網格定義域
    gmsh.model.addPhysicalGroup(2, [surf1_id, surf2_id], 1, "Ω")
    
    # 提取邊界：找出位於 z=0 和 z=L 的兩端全固支閉合曲線
    boundary_entities = gmsh.model.getBoundary([(2, surf1_id), (2, surf2_id)], true, false, false)
    boundary_curves = [entity[2] for entity in boundary_entities]
    
    gmsh.model.addPhysicalGroup(1, boundary_curves, 2, "Γ_clamped")
    
    # 6. 生成 2D 網格並儲存
    gmsh.model.mesh.generate(2)
    mkpath("./msh")
    gmsh.write("./msh/cylinder_shell_q4.msh")
    
    gmsh.finalize()
    println("【成功】3D 圓柱殼 Q4 網格模型已生成：./msh/cylinder_shell_q4.msh")
end

# 執行生成（半徑 1.0, 長度 2.0, 網格密度 0.05）
generate_cylinder_mesh(1.0, 2.0, 0.01, mesh_size=0.05)