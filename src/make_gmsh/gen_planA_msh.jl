# =====================================================================
#  gen_planA_msh.jl — 方案 A 專用網格補生成（2026-09-23 定案）
#
#  方案 A：不改求解器，共同序列 n = 5, 10, 20, 40（h = a/n 倍增）
#    φ/q 場：n = ndiv     = 5, 10, 20, 40
#    w   場：n = ndiv - 1 = 4, 9, 19, 39   （MF 的 w 場刻意粗一級）
#  方板 MF 的 w 場同樣要 n = 4（patchtest_tri3_4.msh）
#
#  執行：cd /home/a/Joker && julia --startup-file=no src/make_gmsh/gen_planA_msh.jl
#  已存在的檔案會跳過（不覆蓋既有網格）
# =====================================================================
import BenchmarkExample: SkewPlate, PatchTest

const NS = [4, 5, 9, 10, 19, 20, 39, 40]

println("=== 斜板結構化三角 β = 30/45/60 ===")
for β_deg in [30, 45, 60]
    for n in NS
        f = "msh/patchtest_skew_tri3_β$(β_deg)_n$n.msh"
        if isfile(f)
            println("skip      $f")
            continue
        end
        SkewPlate.generateMsh(f, β_deg=β_deg, order=1, quad=false, transfinite=n+1)
        println("generated $f  ($(round(filesize(f)/1024, digits=1)) KB)")
    end
end

println("=== 方板（w 場 n = 4）===")
for n in [4]
    f = "msh/patchtest_tri3_$n.msh"
    if isfile(f)
        println("skip      $f")
        continue
    end
    PatchTest.generateMsh(f, transfinite=n+1, order=1, quad=false)
    println("generated $f  ($(round(filesize(f)/1024, digits=1)) KB)")
end

println("=== done ===")
