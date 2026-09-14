#!/usr/bin/env julia
#=
  run_all.jl
  按顺序运行三个有限元分析脚本
  用法: julia run_all.jl
  注意: 本档案需与三个脚本放在同一目录下
=#

using Dates

# ---- 按顺序运行的档案列表 ----
const FILES = [
    # "fem_tri3_loop.jl",
    "fem_tri3_loop_shear.jl",
    # "fem_mix_tri3_loop.jl",
    # "mf_w_loop.jl",
    # "mf_w_loop2.jl",
#     "mf_w_loop_un.jl",
#     "mf_w_loop_un_2.jl",
    # "./vibration/skew/fem_tri3_loop_skew.jl",
    # "./vibration/skew/fem_tri3_loop_shear_skew.jl",
    # "./vibration/skew/fem_mix_tri3_loop_skew.jl",
    # "./vibration/skew/mf_w_loop_skew.jl",
    # "./vibration/skew/mf_w_loop2_skew.jl",
]

# 切换到本脚本所在目录，确保相对路径正确
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

    # include: 同进程运行，包只需加载一次，速度快
    # 如需完全隔离（避免变量冲突），可改为: run(`julia $f`)
    run(`julia $f`)

    println("[$i/$(length(FILES))] $(now())  ✓ 完成: $f")
end

println("\n", "="^70)
println("全部完成: $(now())")
println("="^70)
