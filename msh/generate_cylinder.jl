using Gmsh

function generate_cylinder_structured(R=1.0, L=2.0; n_circumference=20, n_axial=40)
    gmsh.initialize()
    gmsh.model.add("cylinder_structured")
    
    occ = gmsh.model.occ
    
    # 1. 建立幾何控制點
    p1 = occ.addPoint(R, 0, 0)
    p_center = occ.addPoint(0, 0, 0)
    p2 = occ.addPoint(-R, 0, 0)
    
    # 2. 建立兩個半圓弧
    c1 = occ.addCircleArc(p1, p_center, p2)
    c2 = occ.addCircleArc(p2, p_center, p1)
    
    # 3. 沿著 Z 軸擠壓生成兩個半圓柱面
    extru1 = occ.extrude([(1, c1)], 0, 0, L)
    extru2 = occ.extrude([(1, c2)], 0, 0, L)
    
    occ.synchronize()
    
    # 安全抓取兩個曲面的 ID
    surf1_id = [pair[2] for pair in extru1 if pair[1] == 2][1]
    surf2_id = [pair[2] for pair in extru2 if pair[1] == 2][1]
    
    # ==========================================================================
    # 【導師核心升級】：動態超限結構化網格設定
    # ==========================================================================
    # 遍歷模型中所有的「線段 (Dim=1)」
    for (dim, tag) in gmsh.model.getEntities(1)
        type = gmsh.model.getType(dim, tag)
        
        if type == "Line"
            # 如果是直線（軸向長度方向），強制均勻切成 n_axial 等分
            gmsh.model.mesh.setTransfiniteCurve(tag, n_axial + 1)
        elseif type == "Circle"
            # 如果是圓弧（周向圓周方向），強制均勻切成 n_circumference 等分
            gmsh.model.mesh.setTransfiniteCurve(tag, n_circumference + 1)
        end
    end
    
    # 宣告這兩個曲面啟用超限映射演算法（強迫網格橫平豎直）
    gmsh.model.mesh.setTransfiniteSurface(surf1_id)
    gmsh.model.mesh.setTransfiniteSurface(surf2_id)
    
    # 結合 Recombine，產出 100% 純淨的結構化矩形單元
    gmsh.model.mesh.setRecombine(2, surf1_id)
    gmsh.model.mesh.setRecombine(2, surf2_id)
    # ==========================================================================
    
    # 4. 生成網格並儲存
    gmsh.model.mesh.generate(2)
    mkpath("./msh")
    gmsh.write("./msh/cylinder_shell_q4.msh")
    
    gmsh.finalize()
    println("【大成功】絕對結構化的 3D 圓柱殼 Q4 網格已順利生成！")
    println("配置：圓周方向單側 $(n_circumference) 等分，軸向 $(n_axial) 等分。")
end

# 執行生成（你可以自由調整這兩個控制參數來改變網格密度，這對收斂性分析超方便！）
generate_cylinder_structured(1.0, 2.0, n_circumference=20, n_axial=40)