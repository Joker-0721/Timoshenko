import BenchmarkExample: PatchTest

# ---- 1. 結構化三角 ----
# for n in 20:30
#     PatchTest.generateMsh("msh/patchtest_tri3_$n.msh", transfinite=n+1, order=1, quad=false)
# end

# ---- 2. 非結構化三角（n_boundary 方法：邊界固定 + 內部可控）----
for n in 10:35
    PatchTest.generateMsh("msh/patchtest_un_tri3_nb_$n.msh",
                           n_boundary=n,    # 每條邊精確 n 個節點
                           lc=1/n,           # 內部密度
                           order=1,
                           quad=false)
end

# ---- 3. 非結構化四邊形 ----
# for n in 10:35
#     PatchTest.generateMsh("msh/patchtest_un_quad4_$n.msh", lc=1/n, order=1, quad=true)
# end
