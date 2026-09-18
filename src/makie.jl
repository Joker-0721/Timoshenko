using GLMakie, ReadVTK
using GeometryBasics: TriangleFace, Point2f
using Base64: base64decode

cd(@__DIR__)   # 確保相對路徑以 src/ 為基準

# 讀點座標 + 全部 PointData 欄位
function read_vtu_points(filename)
    txt  = read(filename, String)
    npts = parse(Int, match(r"NumberOfPoints=\"(\d+)\"", txt).captures[1])
    points = nothing
    fields = Dict{String,Vector{Float64}}()
    pat = r"(?s)<DataArray\b(?<attrs>[^>]*)>(?<body>.*?)</DataArray>"
    for mm in eachmatch(pat, txt)
        attrs = mm[:attrs]
        body  = mm[:body]
        fmt   = match(r"format=\"(\w+)\"", attrs)
        fmt === nothing && continue
        fmt = fmt.captures[1]
        name = match(r"Name=\"([^\"]*)\"", attrs)
        name = name === nothing ? "" : name.captures[1]
        typ  = match(r"type=\"([^\"]*)\"", attrs)
        typ  = typ === nothing ? "" : typ.captures[1]
        nc   = match(r"NumberOfComponents=\"(\d+)\"", attrs)
        nc   = nc === nothing ? 1 : parse(Int, nc.captures[1])
        typ == "Float64" || continue
        if fmt == "ascii"
            vals = [parse(Float64, s) for s in split(strip(body)) if !isempty(s)]
            if name == "Points"
                points = reshape(vals, nc, :)
            else
                fields[name] = vals
            end
        else
            b = base64decode(strip(replace(body, "\n" => "")))
            d = collect(reinterpret(Float64, b[9:end]))
            expected = npts * nc
            if length(d) == expected - 1
                d = [d; d[end]]
            elseif length(d) > expected
                d = d[1:expected]
            elseif length(d) < expected - 1
                d = [d; d[end] * ones(Float64, expected - length(d))]
            end
            if name == "Points"
                points = reshape(d, nc, :)
            else
                fields[name] = d
            end
        end
    end
    points === nothing && error("找不到 Points DataArray：$filename")
    return points, fields
end

# 讀三角形 faces
function read_vtu_faces(filename)
    vtk    = ReadVTK.VTKFile(filename)
    mcells = ReadVTK.to_meshcells(ReadVTK.get_cells(vtk))
    return [TriangleFace{Int}(c.connectivity...) for c in mcells]
end

# 三角形頂點值線性插值到規則網格
function tri_interp_grid(xs, ys, vals, faces; nx=200, ny=200, pad=1e-10)
    xmin, xmax = extrema(xs); ymin, ymax = extrema(ys)
    xg = range(xmin, xmax, length=nx)
    yg = range(ymin, ymax, length=ny)
    zg = fill(NaN, (nx, ny))
    for f in faces
        i1, i2, i3 = f
        x1,x2,x3 = xs[i1], xs[i2], xs[i3]
        y1,y2,y3 = ys[i1], ys[i2], ys[i3]
        z1,z2,z3 = vals[i1], vals[i2], vals[i3]
        xlo, xhi = min(x1,x2,x3), max(x1,x2,x3)
        ylo, yhi = min(y1,y2,y3), max(y1,y2,y3)
        ix1 = searchsortedfirst(xg, xlo); ix2 = searchsortedlast(xg, xhi)
        iy1 = searchsortedfirst(yg, ylo); iy2 = searchsortedlast(yg, yhi)
        for ix in ix1:ix2, iy in iy1:iy2
            x, y = xg[ix], yg[iy]
            det = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
            abs(det) < eps() && continue
            λ1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / det
            λ2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / det
            λ3 = 1 - λ1 - λ2
            if λ1 >= -pad && λ2 >= -pad && λ3 >= -pad
                zg[ix, iy] = λ1*z1 + λ2*z2 + λ3*z3
            end
        end
    end
    return xg, yg, zg
end

names12 = ["v₁","v₂","v₃","v₄","v₅","v₆","v₇","v₈","v₉","v₁₀","v₁₁","v₁₂"]
BCs = ["CCCC","SSSS","CSCS","SCSC","FSCS","FSSS","FSFS","FCFC","CFFF"]

# 通用繪圖函數（2D 雲圖，無 colorbar、無標題、無坐標軸、無邊框）
function plot_mode(pts, faces, col, out_path)
    pt2 = [Point2f(pts[1, i], pts[2, i]) for i in axes(pts, 2)]
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1], aspect=DataAspect())
    mesh!(ax, pt2, faces, color=col, colormap=:RdBu)
    xg, yg, zg = tri_interp_grid(pts[1, :], pts[2, :], col, faces)
    contour!(ax, xg, yg, zg, levels=12, color=:black, linewidth=0.7)
    hidedecorations!(ax)
    hidespines!(ax)
    mkpath(dirname(out_path))
    save(out_path, fig)
    println("saved: ", out_path)
end

# 處理單個 VTK 檔案（v₁..v₁₂）
function process_vtk(filename, tag, bc, out_root)
    if !isfile(filename)
        println("skip: ", filename)
        return
    end
    pts, fields = read_vtu_points(filename)
    faces = read_vtu_faces(filename)
    for m in 1:12
        if !haskey(fields, names12[m])
            println("skip missing field: ", filename, " :: ", names12[m])
            continue
        end
        plot_mode(pts, faces, fields[names12[m]], "$out_root/$bc/$(tag)_$(bc)_$(names12[m]).png")
    end
end

# # ================================================================
# #  第一段：方板振動 - FEM
# #  VTK 路徑: vtk/fem/vibration/plate/$bc/vib_fem_$bc.vtu
# # ================================================================
# println("="^50)
# println("第一段：方板振動 - FEM")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/fem/vibration/plate/$bc/vib_fem_$bc.vtu",
#                 "fem", bc, "../Fig/fem/vibration/plate/fem_vibration")
# end

# # ================================================================
# #  第二段：方板振動 - FEM shear
# #  VTK 路徑: vtk/fem/vibration/plate/$bc/vib_fem_$bc_shear.vtu
# # ================================================================
# println("="^50)
# println("第二段：方板振動 - FEM shear")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/fem/vibration/plate/$bc/vib_fem_$(bc)_shear.vtu",
#                 "fem_shear", bc, "../Fig/fem/vibration/plate/fem_shear")
# end

# # ================================================================
# #  第三段：方板振動 - FEM mix
# #  VTK 路徑: vtk/fem/vibration/plate/$bc/vib_fem_mix_$bc.vtu
# # ================================================================
# println("="^50)
# println("第三段：方板振動 - FEM mix")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/fem/vibration/plate/$bc/vib_fem_mix_$bc.vtu",
#                 "fem_mix", bc, "../Fig/fem/vibration/plate/fem_mix")
# end

# # ================================================================
# #  第四段：方板振動 - MF
# #  VTK 路徑: vtk/mf/vibration/plate/$bc/vib_mf_w_φ_$bc.vtu
# # ================================================================
# println("="^50)
# println("第四段：方板振動 - MF")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/mf/vibration/plate/$bc/vib_mf_w_φ_$bc.vtu",
#                 "mf", bc, "../Fig/mf/vibration/plate/mf_vibration")
# end

# # ================================================================
# #  第五段：方板振動 - MF2
# #  VTK 路徑: vtk/mf/vibration/plate/$bc/vib_mf_w_φ_${bc}2.vtu
# # ================================================================
# println("="^50)
# println("第五段：方板振動 - MF2")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/mf/vibration/plate/$bc/vib_mf_w_φ_$(bc)2.vtu",
#                 "mf2", bc, "../Fig/mf/vibration/plate/mf_vibration2")
# end

# # ================================================================
# #  第六段：方板屈曲 - FEM
# #  VTK 路徑: vtk/fem/buckling/$bc/fem_$bc.vtu
# # ================================================================
# println("="^50)
# println("第六段：方板屈曲 - FEM")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/fem/buckling/$bc/fem_$bc.vtu",
#                 "fem_buckling", bc, "../Fig/fem/buckling/plate/fem_buckling")
# end

# # ================================================================
# #  第七段：方板屈曲 - FEM mix
# #  VTK 路徑: vtk/fem/buckling/$bc/fem_mix_$bc.vtu
# # ================================================================
# println("="^50)
# println("第七段：方板屈曲 - FEM mix")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/fem/buckling/$bc/fem_mix_$bc.vtu",
#                 "fem_mix_buckling", bc, "../Fig/fem/buckling/plate/fem_mix_buckling")
# end

# # ================================================================
# #  第八段：方板屈曲 - MF
# #  VTK 路徑: vtk/mf/buckling/$bc/mf_w_φ_$bc.vtu
# # ================================================================
# println("="^50)
# println("第八段：方板屈曲 - MF")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/mf/buckling/$bc/mf_w_φ_$bc.vtu",
#                 "mf_buckling", bc, "../Fig/mf/buckling/plate/mf_buckling")
# end

# # ================================================================
# #  第九段：方板屈曲 - MF2
# #  VTK 路徑: vtk/mf/buckling/$bc/mf_w_φ_${bc}2.vtu
# # ================================================================
# println("="^50)
# println("第九段：方板屈曲 - MF2")
# println("="^50)
# for bc in BCs
#     process_vtk("../vtk/mf/buckling/$bc/mf_w_φ_$(bc)2.vtu",
#                 "mf2_buckling", bc, "../Fig/mf/buckling/plate/mf_buckling2")
# end

# # ================================================================
# #  第十段：斜板振動 - FEM（β=30/45/60）
# #  VTK 路徑: vtk/fem/vibration/skew/$bc/vib_fem_β${β}_${bc}_skew.vtu
# # ================================================================
# println("="^50)
# println("第十段：斜板振動 - FEM")
# println("="^50)
# for β in [30, 45, 60]
#     println("--- 斜角 β = $(β)° ---")
#     for bc in BCs
#         process_vtk("../vtk/fem/vibration/skew/$bc/vib_fem_β$(β)_$(bc)_skew.vtu",
#                     "fem", bc, "../Fig/fem/vibration/skew/fem_vibration_β$(β)")
#     end
# end

# # ================================================================
# #  第十一段：斜板振動 - MF（β=30/45/60）
# #  VTK 路徑: vtk/mf/vibration/skew/$bc/vib_mf_w_φ_β${β}_${bc}_0.001.vtu
# # ================================================================
# println("="^50)
# println("第十一段：斜板振動 - MF")
# println("="^50)
# for β in [30, 45, 60]
#     println("--- 斜角 β = $(β)° ---")
#     for bc in BCs
#         process_vtk("../vtk/mf/vibration/skew/$bc/vib_mf_w_φ_β$(β)_$(bc)_0.001.vtu",
#                     "mf", bc, "../Fig/mf/vibration/skew/mf_vibration_β$(β)")
#     end
# end

# # ================================================================
# #  第十二段：斜板屈曲 - FEM（β=30/45/60）
# #  VTK 路徑: vtk/fem/buckling/skew/$bc/fem_β${β}_${bc}_skew.vtu
# # ================================================================
# println("="^50)
# println("第十二段：斜板屈曲 - FEM")
# println("="^50)
# for β in [30, 45, 60]
#     println("--- 斜角 β = $(β)° ---")
#     for bc in BCs
#         process_vtk("../vtk/fem/buckling/skew/$bc/fem_β$(β)_$(bc)_skew.vtu",
#                     "fem_buckling", bc, "../Fig/fem/buckling/skew/fem_buckling_β$(β)")
#     end
# end

# # ================================================================
# #  第十三段：斜板屈曲 - MF（β=30/45/60）
# #  VTK 路徑: vtk/mf/buckling/skew/$bc/mf_w_φ_β${β}_${bc}.vtu
# # ================================================================
# println("="^50)
# println("第十三段：斜板屈曲 - MF")
# println("="^50)
# for β in [30, 45, 60]
#     println("--- 斜角 β = $(β)° ---")
#     for bc in BCs
#         process_vtk("../vtk/mf/buckling/skew/$bc/mf_w_φ_β$(β)_$(bc).vtu",
#                     "mf_buckling", bc, "../Fig/mf/buckling/skew/mf_buckling_β$(β)")
#     end
# end

# ================================================================
#  第十四段：new_vib（u₁..u₉ + v₁..v₁₂，新格式）
#  VTK 路徑: vtk/new_vib/*.vtu
#  輸出:     Fig/new_vib/<檔名>/<檔名>_<場名>.png
# ================================================================
const U_FIELDS = ["u₁","u₂","u₃","u₄","u₅","u₆","u₇","u₈","u₉"]
function process_new_vib(filename, out_root)
    isfile(filename) || (println("skip: ", filename); return)
    tag = splitext(basename(filename))[1]
    pts, fields = read_vtu_points(filename)
    faces = read_vtu_faces(filename)
    for name in [U_FIELDS; names12]
        if !haskey(fields, name)
            println("  skip missing: $tag :: $name")
            continue
        end
        plot_mode(pts, faces, fields[name], "$out_root/$tag/$(tag)_$(name).png")
    end
end

println("="^50)
println("第十四段：new_vib")
println("="^50)
vib_dir  = joinpath(@__DIR__, "..", "vtk", "new_vib")
out_root = joinpath(@__DIR__, "..", "Fig", "new_vib")
for f in sort([x for x in readdir(vib_dir; join=true) if endswith(x, ".vtu")])
    println("處理: ", basename(f))
    process_new_vib(f, out_root)
end

println("\n全部完成。")
