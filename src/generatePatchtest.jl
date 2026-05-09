import BenchmarkExample: PatchTest
using Printf

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
b = 1.0
msh_dir = normpath(joinpath(@__DIR__, "..", "msh"))

bui_2011_node_patterns = [
    7,
    9,
    11,
    13,
    15,
    17,
]

table_9_1_aspect_ratios = [
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    0.90,
    1.00,
    1.10,
    1.20,
    1.30,
    1.40,
    1.41,
]

table_9_12_aspect_ratios = [
    1.00,
    1.50,
    2.00,
    2.50,
]

function aspect_tag(r)
    return replace(@sprintf("%.2f", r), "."=>"p")
end

function generate_rect_msh(prefix, r)
    a = r*b
    nx = max(2, round(Int, n*r)) + 1
    ny = n + 1
    PatchTest.generateMsh(
        joinpath(msh_dir, "$(prefix)_ab_$(aspect_tag(r)).msh"),
        a = a,
        b = b,
        transfinite_x = nx,
        transfinite_y = ny,
        order = 1,
        quad = false,
    )
end

function generate_square_msh(prefix, nodes_per_side)
    PatchTest.generateMsh(
        joinpath(msh_dir, "$(prefix)_$(nodes_per_side)x$(nodes_per_side).msh"),
        a = b,
        b = b,
        transfinite_x = nodes_per_side,
        transfinite_y = nodes_per_side,
        order = 1,
        quad = false,
    )
end

for nodes_per_side in bui_2011_node_patterns
    generate_square_msh("bui_2011_square", nodes_per_side)
end

for r in table_9_1_aspect_ratios
    generate_rect_msh("mindlin_uniaxial_ssss", r)
end

for r in table_9_12_aspect_ratios
    generate_rect_msh("mindlin_shear_ssss", r)
end

generate_rect_msh("mindlin_biaxial_ssss", 1.0)
generate_rect_msh("mindlin_combined_load_ssss", 1.0)
