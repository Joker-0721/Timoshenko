import BenchmarkExample: SkewPlate

# =====================================================================
#  generate_skew_b75.jl — 補齊 β=75° 結構化三角網格
#  背景：msh/ 原本只有 β=30/45/60（generate_skew.jl 的 angles 雖含 75 但當時未生成）
#        論文表格新增 β=75° 一組（SSSS/SCSC/CCCC/FSFS × 5 方法）→ 需要這批網格
#  需求 ndiv：loop 預設 5,10,20；MF「最優比」另用 9（=10−1）、19（=20−1）
#  這裡比照 β=30 既有檔名集合一次補齊（n=4,5,9,10..35,39,40），已存在者跳過
#  執行（在 /home/a/Joker 下）：julia --startup-file=no src/make_gmsh/generate_skew_b75.jl
# =====================================================================

const NS = vcat([4, 5, 9], collect(10:35), [39, 40])

n_done = 0
for n in NS
    f = "msh/patchtest_skew_tri3_β75_n$n.msh"
    if isfile(f)
        println("skip (exists) ", f)
        continue
    end
    SkewPlate.generateMsh(f, β_deg=75, order=1, quad=false, transfinite=n + 1)
    println("written ", f)
    global n_done += 1
end

println("β=75 meshes: $(n_done) new / $(length(NS)) requested, dir = msh/")
