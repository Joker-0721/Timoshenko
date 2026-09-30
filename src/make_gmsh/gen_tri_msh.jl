# =====================================================================
#  gen_tri_msh.jl — 三角形板（等腰）網格生成器
#  2026-09-23 新增（S9 / C6；Kitipornchai 1993 等腰三角形 Mindlin 板）
#
#  幾何：底邊 b（y = 0，A(0,0) → B(b,0)），頂點 C(b/2, H)
#        頂角 θ（度）→ 高 H = (b/2)/tan(θ/2)
#        等邊三角形 θ = 60°；θ < 60° 為銳角、θ > 60° 為鈍角等腰
#
#  物理組（solver 只認這套命名）：
#        Γ¹ = 底邊 A–B
#        Γ² = 右腰 B–C
#        Γ³ = 左腰 C–A
#        Γ  = 邊元素群（1-D）  → MF 迴圈的 normal=true 與 elements_m_* 需要
#        Ω  = 面
#
#  檔名：msh/patchtest_tri_tri3_θ$(θ)_n$(n).msh
#        ⚠ $(θ) 用 Julia 預設 Float64 格式（60.0 → "60.0"），必須與
#          src/*/tri/*.jl 內的 gmsh.open 組法完全一致
#
#  執行：cd mindlin_code && julia --startup-file=no src/make_gmsh/gen_tri_msh.jl
#        已存在的檔案會跳過（不覆蓋既有網格）；面積自檢在檔尾
# =====================================================================
import BenchmarkExample
import Gmsh: gmsh

const THETAS = [60.0]                 # 頂角（度）；文獻有多個頂角時在此加值
const B      = 1.0                    # 底邊長（與方板 a = b = 1 同口徑）
# 細分數：5/10/20/40 為收斂序列；4/9/19/39 是 MF「最優比」的 Q 場網格
# （mf_opt_ndivs 讓 n_q ≤ n_s ⇒ ndiv_q = ndiv−1，見 src/mf_ratio.jl）
const NS     = [4, 5, 9, 10, 19, 20, 39, 40]   # 每條邊細分數 n → h = b/n

height(θdeg::Float64) = (B/2) / tand(θdeg/2)

# ---------------------------------------------------------------- 生成
function gen_iso_tri(file::String, θdeg::Float64, n::Int; order::Int = 1)
    H  = height(θdeg)
    lc = max(B, H) / n
    gmsh.initialize()
    gmsh.model.add("IsoTriangle")

    gmsh.model.geo.addPoint(0.0,   0.0, 0.0, lc, 1)   # A 左下
    gmsh.model.geo.addPoint(B,     0.0, 0.0, lc, 2)   # B 右下
    gmsh.model.geo.addPoint(B/2,   H,   0.0, lc, 3)   # C 頂點

    gmsh.model.geo.addLine(1, 2, 1)      # Γ¹ 底邊
    gmsh.model.geo.addLine(2, 3, 2)      # Γ² 右腰
    gmsh.model.geo.addLine(3, 1, 3)      # Γ³ 左腰
    gmsh.model.geo.addCurveLoop([1, 2, 3], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [1], -1, "Γ¹")
    gmsh.model.addPhysicalGroup(1, [2], -1, "Γ²")
    gmsh.model.addPhysicalGroup(1, [3], -1, "Γ³")
    gmsh.model.addPhysicalGroup(2, [1], -1, "Ω")

    for i in 1:3
        gmsh.model.mesh.setTransfiniteCurve(i, n + 1)
    end
    gmsh.model.mesh.setTransfiniteSurface(1)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    tag = BenchmarkExample.addEdgeElements((2, 1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()

    gmsh.write(file)
    gmsh.finalize()
    return H
end

# ---------------------------------------------------------------- 面積自檢
"""
    check_area(file) -> (A, ntri, nnode)

讀回 .msh，用直線三角形面積和核算（面積守恆檢查）。
"""
function check_area(file::String)
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
        t == 2 || continue
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
println("=== 三角形板（等腰） b = $(B)，頂角 θ = $(THETAS)，每邊細分 n = $(NS) ===")
for θ in THETAS
    for n in NS
        f = "msh/patchtest_tri_tri3_θ$(θ)_n$n.msh"
        if isfile(f)
            println("skip      $f  ($(round(filesize(f)/1024, digits=1)) KB)")
            continue
        end
        H = gen_iso_tri(f, θ, n)
        println("generated $f  ($(round(filesize(f)/1024, digits=1)) KB)  H = $(round(H, digits=6))")
    end
end

println("\n=== 面積自檢（½bH 對比三角形面積和） ===")
for θ in THETAS
    Aref = B * height(θ) / 2
    for n in NS
        f = "msh/patchtest_tri_tri3_θ$(θ)_n$n.msh"
        isfile(f) || continue
        A, ntri, nnode = check_area(f)
        println("  θ=$θ n=$n  nodes=$(lpad(nnode,5))  tri=$(lpad(ntri,5))  " *
                "A=$(round(A, digits=8))  Aref=$(round(Aref, digits=8))  " *
                "ratio=$(round(A/Aref, digits=5))")
    end
end

println("\n=== 邊界節點數（Γ¹/Γ²/Γ³ 應各為 n+1）===")
for θ in THETAS
    n = NS[end]
    f = "msh/patchtest_tri_tri3_θ$(θ)_n$n.msh"
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
