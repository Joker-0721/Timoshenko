import BenchmarkExample: SkewPlate

# ================================================================
#  斜板網格生成脚本
#  三種模式：
#    1. 結構化三角：transfinite=n+1
#    2. 非結構化三角（n_boundary 方法）：n_boundary=n, lc=1/n
#       （邊界節點數精確固定，內部 Delaunay 非結構化，節點數平滑增長）
#    3. 結構化四邊形：transfinite=n+1, quad=true
# ================================================================

# ---- 1. 結構化三角（多個角度）----
for β_deg in [30, 45, 60]
    for n in 10:35
        SkewPlate.generateMsh("msh/patchtest_skew_tri3_β$(β_deg)_n$n.msh",
                               β_deg=β_deg, order=1, quad=false, transfinite=n+1)
    end
end

# ---- 2. 非結構化三角（n_boundary 方法，推薦用嚟做收斂性分析）----
# n_boundary：每條邊精確 n 個節點（唔會跳）
# lc：內部網格密度（lc=1/n 令內部同邊界密度一致）
# 總節點數會平滑增長，每次 n+1 大約增加 10-20 個節點
for β_deg in [30, 45, 60]
    for n in 10:35
        SkewPlate.generateMsh("msh/patchtest_skew_un_tri3_nb_β$(β_deg)_n$n.msh",
                               β_deg=β_deg,
                               n_boundary=n,    # 邊界固定 n 節點
                               lc=1/n,           # 內部密度
                               order=1,
                               quad=false)
    end
end

# ---- （對比）純非結構化三角（純 lc，會跳變，唔推薦）----
# for n in 10:35
#     SkewPlate.generateMsh("msh/patchtest_skew_un_tri3_lc_β30_n$n.msh",
#                            lc=1/n, β_deg=30, order=1, quad=false)
# end

# ---- 3. 結構化四邊形 ----
for n in 10:35
    SkewPlate.generateMsh("msh/patchtest_skew_quad4_β30_n$n.msh",
                           β_deg=30, order=1, quad=true, transfinite=n+1)
end
