#!/usr/bin/env julia
#=
    run_all.jl
    按顺序跑晒 skew/ 入面所有 case
    用法: julia run_all.jl
    注意: 本档案放在 skew/ 入面，但会自动 cd 去上一层 src/，
          令 case 入面嘅 ../msh/、../date/ 都解析到 project root
=#
using Dates

# ---- 按顺序运行的档案列表（相对 src/ 目录）----
const FILES = [
    # "./buckling/skew/CCCC_fem_mix_tri3.jl",
    # "./buckling/skew/CCCC_fem_tri3_shear.jl",
    # "./buckling/skew/CCCC_fem_tri3.jl",
    # "./buckling/skew/CCCC_mf_w_mix_tri3.jl",
    # "./buckling/skew/CCCC_mf1_w_mix_tri3.jl",

    # "./buckling/skew/CSCS_fem_mix_tri3.jl",
    # "./buckling/skew/CSCS_fem_tri3_shear.jl",
    # "./buckling/skew/CSCS_fem_tri3.jl",
    # "./buckling/skew/CSCS_mf_w_mix_tri3.jl",
    # "./buckling/skew/CSCS_mf1_w_mix_tri3.jl",

    # "./buckling/skew/FSCS_fem_mix_tri3.jl",
    # "./buckling/skew/FSCS_fem_tri3_shear.jl",
    # "./buckling/skew/FSCS_fem_tri3.jl",
    # "./buckling/skew/FSCS_mf_w_mix_tri3.jl",
    # "./buckling/skew/FSCS_mf1_w_mix_tri3.jl",

    # "./buckling/skew/FSSS_fem_mix_tri3.jl",
    # "./buckling/skew/FSSS_fem_tri3_shear.jl",
    # "./buckling/skew/FSSS_fem_tri3.jl",
    # "./buckling/skew/FSSS_mf_w_mix_tri3.jl",
    # "./buckling/skew/FSSS_mf1_w_mix_tri3.jl",

    # "./buckling/skew/FCFC_fem_mix_tri3.jl",
    # "./buckling/skew/FCFC_fem_tri3_shear.jl",
    "./buckling/skew/FCFC_fem_tri3.jl",
    # "./buckling/skew/FCFC_mf_w_mix_tri3.jl",
    # "./buckling/skew/FCFC_mf1_w_mix_tri3.jl",

    # "./buckling/skew/SSSS_fem_mix_tri3.jl",
    # "./buckling/skew/SSSS_fem_tri3_shear.jl",
    # "./buckling/skew/SSSS_fem_tri3.jl",
    # "./buckling/skew/SSSS_mf_w_mix_tri3.jl",
    # "./buckling/skew/SSSS_mf1_w_mix_tri3.jl",
]

# cd 去 src/（skew/ 上兩層），令 case 入面嘅 ../msh/、../date/ 解析正確
cd(joinpath(@__DIR__, "..", ".."))

println("="^70)
println("批量运行开始: $(now())")
println("cwd: $(pwd())")
println("共 $(length(FILES)) 个档案")
println("="^70)

for (i, f) in enumerate(FILES)
    println("\n", "─"^70)
    println("[$i/$(length(FILES))] $(now())  开始运行: $f")
    println("─"^70)
    if !isfile(f)
        error("档案不存在，终止运行: $f")
    end
    run(`julia $f`)
    println("[$i/$(length(FILES))] $(now())  完成: $f")
end

println("\n", "="^70)
println("全部完成: $(now())")
println("="^70)
