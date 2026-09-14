import BenchmarkExample: SkewPlate

# ---- 1. 結構化三角（多個角度）----
for β_deg in [0, 15, 30, 45, 60, 75]
    for n in 10:35
        SkewPlate.generateMsh("msh/patchtest_skew_tri3_β$(β_deg)_n$n.msh",
                               β_deg=β_deg, order=1, quad=false, transfinite=n+1)
    end
end

# ---- 2. 非結構化三角（n_boundary 方法：邊界固定 + 內部可控）----
for β_deg in [0, 15, 30, 45, 60, 75]
    for n in 10:35
        SkewPlate.generateMsh("msh/patchtest_skew_un_tri3_nb_β$(β_deg)_n$n.msh",
                               β_deg=β_deg,
                               n_boundary=n,    # 每條邊精確 n 個節點
                               lc=1/n,           # 內部密度
                               order=1,
                               quad=false)
    end
end

# ---- 3. 結構化四邊形 ----
for n in 10:35
    SkewPlate.generateMsh("msh/patchtest_skew_quad4_β30_n$n.msh",
                           β_deg=30, order=1, quad=true, transfinite=n+1)
end
