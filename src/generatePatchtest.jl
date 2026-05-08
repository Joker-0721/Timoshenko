import BenchmarkExample: PatchTest

# PatchTest.generateMsh("msh/patchtest.msh", transfinite = 3, order = 1, quad = false)

# for n in 2:25
#     PatchTest.generateMsh("msh/patchtest_tri3_$n.msh", transfinite = n + 1, order = 1, quad = false)
#     PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite = n + 1, order = 2, quad = false)
# end

# for n in 2:25
#     PatchTest.generateMsh("msh/patchtest_quad4_$n.msh", transfinite = n + 1, order = 1, quad = true)
#     PatchTest.generateMsh("msh/patchtest_quad8_$n.msh", transfinite = n + 1, order = 2, quad = true)
# end


n = 12

PatchTest.generateMsh("D:\\Joker\\Timoshenko\\msh\\mindlin_ssss.msh", transfinite = n + 1, order = 1, quad = false)
# PatchTest.generateMsh("msh/mindlin_biaxial_ssss.msh", transfinite = n + 1, order = 1, quad = false)
# PatchTest.generateMsh("msh/mindlin_shear_ssss.msh", transfinite = n + 1, order = 1, quad = false)
# PatchTest.generateMsh("msh/mindlin_combined_load_ssss.msh", transfinite = n + 1, order = 1, quad = false)