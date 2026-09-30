# =====================================================================
#  gen_circle_msh.jl — 圓板（四分之一圓盤）網格生成器
#  2026-09-23 新增（S7；B5 定案：四分之一圓盤 + 兩對稱邊）
#
#  幾何：半徑 R 的四分之一圓盤（x ≥ 0、y ≥ 0），3 塊 transfinite 拓樸
#        （拓樸沿用 BenchmarkExample.jl/src/Circular.jl，但座標由 R=5 硬寫改為參數 R）
#
#  物理組（solver 只認這套命名；見 src/*/circle/*.jl）：
#        Γ¹ = 圓弧（兩段）        → 實體邊界：CC 夾支 / SS 簡支
#        Γ² = y = 0 直邊          → 對稱邊，鎖 φ_y
#        Γ³ = x = 0 直邊          → 對稱邊，鎖 φ_x
#        A  = 圓心點
#        Γ  = 邊元素群（1-D）      → MF 迴圈的 normal=true 與 elements_m_* 需要
#        Ω  = 3 個面
#
#  為什麼需要 Γ²/Γ³ 而不是方板式 Γ¹..Γ⁴：圓板邊界只有圓弧一條實體邊，
#  其餘兩條是對稱線；對稱線上只需 φ_n = 0（必要條件），w 的 ∂w/∂n = 0
#  由變分自然滿足 → 不需要 w 懲罰，故不能套方板 BC 命名。
#
#  檔名：msh/patchtest_circ_tri3_R$(R)_n$(n).msh
#        ⚠ $(R) 用 Julia 預設 Float64 格式（1.0 → "1.0"），必須與
#          src/*/circle/*.jl 內的 gmsh.open 組法完全一致
#
#  執行：cd mindlin_code && julia --startup-file=no src/make_gmsh/gen_circle_msh.jl
#        已存在的檔案會跳過（不覆蓋既有網格）；面積自檢在檔尾
# =====================================================================
import BenchmarkExample
import Gmsh: gmsh

const RSR = [1.0]                # 圓板半徑（R = 1 與方板 a = b = 1 同口徑）
# 細分數：5/10/15/20 為收斂序列；4/9/14/19 是 MF「最優比」的 Q 場網格
# （mf_opt_ndivs 讓 n_q ≤ n_s ⇒ ndiv_q = ndiv−1，見 src/mf_ratio.jl）
const NS  = [4, 5, 9, 10, 14, 15, 19, 20]      # transfinite：每條曲線細分數 n → h ≈ R/(2n)

# ---------------------------------------------------------------- 生成
function gen_quarter_circle(file::String, R::Float64, n::Int; order::Int = 1)
    gmsh.initialize()
    gmsh.model.add("QuarterCircle")

    lc = R / n
    gmsh.model.geo.addPoint(0.0,      0.0,     0.0, lc, 1)   # 1 圓心
    gmsh.model.geo.addPoint(R/2,      0.0,     0.0, lc, 2)   # 2 (R/2, 0)
    gmsh.model.geo.addPoint(R,        0.0,     0.0, lc, 3)   # 3 (R, 0)
    gmsh.model.geo.addPoint(R/√2,     R/√2,    0.0, lc, 4)   # 4 45°
    gmsh.model.geo.addPoint(0.0,      R,       0.0, lc, 5)   # 5 (0, R)
    gmsh.model.geo.addPoint(0.0,      R/2,     0.0, lc, 6)   # 6 (0, R/2)
    gmsh.model.geo.addPoint(2.0*R/5,  2.0*R/5, 0.0, lc, 7)   # 7 Circular.jl 的 (2,2) 按 R/5 縮放

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addCircleArc(3, 1, 4, 3)
    gmsh.model.geo.addCircleArc(4, 1, 5, 4)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(6, 1, 6)
    gmsh.model.geo.addLine(7, 2, 7)
    gmsh.model.geo.addLine(7, 6, 8)
    gmsh.model.geo.addLine(7, 4, 9)

    gmsh.model.geo.addCurveLoop([1, -7, 8, 6], 1)
    gmsh.model.geo.addCurveLoop([2, 3, -9, 7], 2)
    gmsh.model.geo.addCurveLoop([9, 4, 5, -8], 3)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [3, 4], -1, "Γ¹")   # 圓弧
    gmsh.model.addPhysicalGroup(1, [1, 2], -1, "Γ²")   # y = 0 對稱邊
    gmsh.model.addPhysicalGroup(1, [5, 6], -1, "Γ³")   # x = 0 對稱邊
    gmsh.model.addPhysicalGroup(0, [1],    -1, "A")    # 圓心
    gmsh.model.addPhysicalGroup(2, [1, 2, 3], -1, "Ω")

    for i in 1:9
        gmsh.model.mesh.setTransfiniteCurve(i, n + 1)
    end
    for i in 1:3
        i == 2 ? gmsh.model.mesh.setTransfiniteSurface(i, "Right") :
                 gmsh.model.mesh.setTransfiniteSurface(i)
    end

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    tag = BenchmarkExample.addEdgeElements((2, 1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()

    gmsh.write(file)
    gmsh.finalize()
    return nothing
end

# ---------------------------------------------------------------- 面積自檢
"""
    check_area(file, Aref) -> (A, ntri, nnode)

讀回 .msh，用直線三角形面積和核算（面積守恆檢查，見執行計劃 §檢查清單）。
"""
function check_area(file::String, Aref::Float64)
    gmsh.initialize()
    gmsh.open(file)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    idx = Dict{Int,Int}()
    for (i, t) in enumerate(tags)
        idx[Int(t)] = i
    end
    px(i) = coords[3*(i-1) + 1]
    py(i) = coords[3*(i-1) + 2]

    etypes, _, enodes = gmsh.model.mesh.getElements(2)
    A = 0.0
    ntri = 0
    for (t, conn) in zip(etypes, enodes)
        t == 2 || continue                      # 3 節點三角形
        for k in 1:3:length(conn)
            a1 = idx[Int(conn[k])]; a2 = idx[Int(conn[k+1])]; a3 = idx[Int(conn[k+2])]
            x1 = px(a1); y1 = py(a1)
            x2 = px(a2); y2 = py(a2)
            x3 = px(a3); y3 = py(a3)
            A += abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) / 2
            ntri += 1
        end
    end
    gmsh.finalize()
    return A, ntri, length(tags)
end

# ---------------------------------------------------------------- 主流程
println("=== 圓板（四分之一圓盤） R = $(RSR)，逐條曲線細分 n = $(NS) ===")
for R in RSR
    for n in NS
        f = "msh/patchtest_circ_tri3_R$(R)_n$n.msh"
        if isfile(f)
            println("skip      $f  ($(round(filesize(f)/1024, digits=1)) KB)")
            continue
        end
        gen_quarter_circle(f, R, n)
        println("generated $f  ($(round(filesize(f)/1024, digits=1)) KB)")
    end
end

println("\n=== 面積自檢（πR²/4 對比三角形面積和） ===")
for R in RSR
    Aref = π * R^2 / 4
    for n in NS
        f = "msh/patchtest_circ_tri3_R$(R)_n$n.msh"
        isfile(f) || continue
        A, ntri, nnode = check_area(f, Aref)
        println("  R=$R n=$n  nodes=$(lpad(nnode,5))  tri=$(lpad(ntri,5))  " *
                "A=$(round(A, digits=8))  Aref=$(round(Aref, digits=8))  " *
                "ratio=$(round(A/Aref, digits=5))")
    end
end

println("\n=== 邊界節點數（Γ¹ 應為 2n+1、Γ²/Γ³ 各 n+1）===")
for R in RSR
    n = NS[end]
    f = "msh/patchtest_circ_tri3_R$(R)_n$n.msh"
    isfile(f) || continue
    gmsh.initialize()
    gmsh.open(f)
    for name in ["Γ¹", "Γ²", "Γ³", "Γ", "Ω"]
        try
            ents = gmsh.model.getEntitiesForPhysicalName(name)
            println("  $name -> $(length(ents)) 實體")
        catch e
            println("  $name -> (讀取失敗: $(sprint(showerror, e)[1:min(end,60)]))")
        end
    end
    gmsh.finalize()
end

println("\n=== done ===")
