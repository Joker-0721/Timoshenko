#!/usr/bin/env julia
#=
    run_all.jl
    按顺序跑晒 plate/ 入面所有 case
    用法: julia run_all.jl
    注意: 本档案需放在 src/ 目录下（同 buckling/ 同级）
=#
using Dates

# ---- 按顺序运行的档案列表 ----
const FILES = [
    "./buckling/plate/CCCC_fem_mix_tri3.jl",
    "./buckling/plate/CCCC_fem_tri3_shear.jl",
    "./buckling/plate/CCCC_fem_tri3.jl",
    "./buckling/plate/CCCC_mf_w_mix_tri3.jl",
    "./buckling/plate/CCCC_mf1_w_mix_tri3.jl",

    # "./buckling/plate/CSCS_fem_mix_tri3.jl",
    # "./buckling/plate/CSCS_fem_tri3_shear.jl",
    # "./buckling/plate/CSCS_fem_tri3.jl",
    # "./buckling/plate/CSCS_mf_w_mix_tri3.jl",
    # "./buckling/plate/CSCS_mf1_w_mix_tri3.jl",

    # "./buckling/plate/FSCS_fem_mix_tri3.jl",
    # "./buckling/plate/FSCS_fem_tri3_shear.jl",
    # "./buckling/plate/FSCS_fem_tri3.jl",
    # "./buckling/plate/FSCS_mf_w_mix_tri3.jl",
    # "./buckling/plate/FSCS_mf1_w_mix_tri3.jl",

    # "./buckling/plate/FSSS_fem_mix_tri3.jl",
    # "./buckling/plate/FSSS_fem_tri3_shear.jl",
    # "./buckling/plate/FSSS_fem_tri3.jl",
    # "./buckling/plate/FSSS_mf_w_mix_tri3.jl",
    # "./buckling/plate/FSSS_mf1_w_mix_tri3.jl",

    # "./buckling/plate/SCSC_fem_mix_tri3.jl",
    # "./buckling/plate/SCSC_fem_tri3_shear.jl",
    # "./buckling/plate/SCSC_fem_tri3.jl",
    # "./buckling/plate/SCSC_mf_w_mix_tri3.jl",
    # "./buckling/plate/SCSC_mf1_w_mix_tri3.jl",

    "./buckling/plate/SSSS_fem_mix_tri3.jl",
    "./buckling/plate/SSSS_fem_tri3_shear.jl",
    "./buckling/plate/SSSS_fem_tri3.jl",
    "./buckling/plate/SSSS_mf_w_mix_tri3.jl",
    "./buckling/plate/SSSS_mf1_w_mix_tri3.jl",
]

# 切换到本脚本所在目录，确保相对路径正确（../msh/、../date/ 都解析到 project root）
cd(@__DIR__)

println("="^70)
println("批量运行开始: $(now())")
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
