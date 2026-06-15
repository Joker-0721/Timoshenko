import BenchmarkExample: PatchTest

# for n in 2:25
#     PatchTest.generateMsh("msh/patchtest_tri3_$n.msh", transfinite=n+1, order=1, quad=false)
#     PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+1, order=2, quad=false)
# end

# for n in 2:3
#     PatchTest.generateMsh("msh/patchtest_quad4_$n.msh", transfinite=n+1, order=1, quad=true)
#     PatchTest.generateMsh("msh/patchtest_quad8_$n.msh", transfinite=n+1, order=2, quad=true)
# end

# for n in 2:32
#     PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = 1/(n-1), order=1, quad=false)
# end

for n in 3:32
    PatchTest.generateMsh("msh/st_q_$n.msh", lc = 1/(n-1), order=1, quad=true)
    # PatchTest.generateMsh("msh/patchtest_quad8_$n.msh", transfinite=n+1, order=2, quad=true)
end

# ---------------------------------------------------------------------------------------------------------

# for n in (0.04166)
#     PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = n, order=1, quad=false)
#     # PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+1, order=2, quad=false)
# end

# for nw in (0.1,0.2)
#     for nq in (0.3,0.4)
#     PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = nw*nq, order=1, quad=false)
#    # PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+0.001, order=2, quad=false)
#     end
# end