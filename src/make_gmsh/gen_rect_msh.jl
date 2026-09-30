# =====================================================================
#  gen_rect_msh.jl — 矩形板網格（方板 loop 開出 α = a/b 之後用）
#  2026-09-23 新增（使用者：方板的 α 都沒開出來用）
#
#  幾何：x ∈ [0, a]，y ∈ [0, b]，a = 1.0 固定，b = a/α
#        ⇒ α = a/b，正方板 α = 1（沿用既有 patchtest_tri3_$n.msh，不重複生成）
#
#  α 取值 = Liew 1993 Part I 矩形板表的 a/b：
#        {0.4, 0.6, 0.8, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5}
#        （來源 deepseek_harness/data/liew1993/liew1993_lambda_table.csv，
#          142 行 × t/b = 0.001/0.1/0.2，5 種邊界）
#
#  檔名：msh/patchtest_rect_tri3_a$(α)_n$n.msh
#        與 loop 內 _msh(n) 的組法**必須完全一致**（α 用 Julia 預設輸出格式）
#
#  執行：cd ~/Joker && julia --startup-file=no src/make_gmsh/gen_rect_msh.jl
#  已存在的檔案會跳過（不覆蓋既有網格）
# =====================================================================
import BenchmarkExample: PatchTest

const ALPHAS = [0.4, 0.6, 0.8, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]   # a/b（不含 1.0）
const NS = [4, 5, 9, 10, 19, 20, 39, 40]                          # 方案 A 序列

println("=== 矩形板 α = a/b（a = 1.0, b = a/α）===")
for α in ALPHAS
    for n in NS
        f = "msh/patchtest_rect_tri3_a$(α)_n$n.msh"
        if isfile(f)
            println("skip      $f")
            continue
        end
        PatchTest.generateMsh(f, transfinite = n + 1, order = 1, quad = false, a = 1.0, b = 1.0/α)
        println("generated $f  ($(round(filesize(f)/1024, digits=1)) KB)")
    end
end

println("=== 產出清單（節點數自檢用 src/make_gmsh/check_rect_msh.jl）===")
for α in ALPHAS
    for n in [10, 40]
        f = "msh/patchtest_rect_tri3_a$(α)_n$n.msh"
        isfile(f) || continue
        println("  α=$α n=$n  $(basename(f))  $(filesize(f)) B")
    end
end

println("=== done ===")
