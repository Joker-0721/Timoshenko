using GLMakie, ReadVTK
using GeometryBasics: TriangleFace
using Base64: base64decode

cd(@__DIR__)   # 確保相對路徑以 src/ 為基準

# 讀點座標 + 全部 PointData 欄位（自製 base64 解碼，繞過 ReadVTK 讀 256 值時的 off-by-one）
function read_vtu_points(filename)
    txt  = read(filename, String)
    npts = parse(Int, match(r"NumberOfPoints=\"(\d+)\"", txt).captures[1])

    points = nothing
    fields = Dict{String,Vector{Float64}}()

    # 抓取所有 DataArray 區塊
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

    typ == "Float64" || continue      # ★ 只處理 Float64；cells 的 UInt8/Int64 全部跳過

    if fmt == "ascii"
        vals = [parse(Float64, s) for s in split(strip(body)) if !isempty(s)]
        if name == "Points"
            points = reshape(vals, nc, :)
        else
            fields[name] = vals
        end
    else   # fmt == "binary"
        b = base64decode(strip(replace(body, "\n" => "")))
        d = collect(reinterpret(Float64, b[9:end]))
        expected = npts * nc
        if length(d) == expected - 1
            # 檔案 base64 尾端被截斷 8 bytes，補齊最後一個 Float64
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


# 讀三角形 faces（cells 不受 off-by-one 影響，用 ReadVTK）
function read_vtu_faces(filename)
    vtk    = ReadVTK.VTKFile(filename)
    mcells = ReadVTK.to_meshcells(ReadVTK.get_cells(vtk))
    return [TriangleFace{Int}(c.connectivity...) for c in mcells]
end

# 將三角形頂點值插值到規則網格（回傳 xg, yg, zg）
function tri_interp_grid(xs, ys, vals, faces; nx=200, ny=200, pad=1e-10)
    xmin, xmax = extrema(xs); ymin, ymax = extrema(ys)
    xg = range(xmin, xmax, length=nx)
    yg = range(ymin, ymax, length=ny)
    zg = fill(NaN, (nx, ny))

    # 對每個三角形，只在其包圍盒內的網格點做重心插值（快速）
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
            # 重心座標
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
BCs = ["CCCC","CSCS","FSCS","FSSS","SCSC","SSSS"]

# ================================================================
#  開關：想跑邊段就將對應改成 true，唔想跑就 false
# ================================================================
DO_FEM  = true    # 純有限元位移法（fem + fem_shear）
DO_MIX  = false    # 混合有限元（fem_mix）
DO_MF   = false    # 無網格法（mf + mf2）
PLOT_3D = false   # true = 3D 曲面圖, false = 2D 雲圖+等高線

# ================================================================
#  通用繪圖函數
# ================================================================
function plot_mode(pts, faces, col, tag, bc, mode_name, out_dir)
    if PLOT_3D
        # ---- 3D 曲面圖 ----
        fig = Figure(size=(800, 600))
        ax = Axis3(fig[1, 1])
        plt = mesh!(ax, pts, faces, color=col, colormap=:RdBu, shading=false)
        ax.title = "$tag  $bc  $mode_name"
        Colorbar(fig[1, 2], plt, label=mode_name)
    else
        # ---- 2D 雲圖 + 等高線 ----
        pt2 = [Point2f(pts[1, i], pts[2, i]) for i in axes(pts, 2)]
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], aspect=DataAspect())
        plt = mesh!(ax, pt2, faces, color=col, colormap=:RdBu)
        # 等高線
        xg, yg, zg = tri_interp_grid(pts[1, :], pts[2, :], col, faces)
        contour!(ax, xg, yg, zg, levels=12, color=:black, linewidth=0.7)
        Colorbar(fig[1, 2], plt, label=mode_name)
        ax.title = "$tag  $bc  $mode_name"
    end
    mkpath(out_dir)
    out_path = "$out_dir/$(tag)_$(bc)_$(mode_name).png"
    save(out_path, fig)
    println("saved: ", out_path)
end

# 處理單個 VTK 檔案的 12 個模態
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
        col = fields[names12[m]]
        plot_mode(pts, faces, col, tag, bc, names12[m], "$out_root/$bc")
    end
end

# ================================================================
#  第一段：純有限元（fem 位移法 + shear 降階積分）
# ================================================================
if DO_FEM
    println("="^50)
    println("DO_FEM：純有限元")
    println("="^50)
    for bc in BCs
        process_vtk("../vtk/fem/vibration/$bc/vib_fem_$(bc)_0.001.vtu",
                    "fem", bc, "../Fig/fem/vibration/fem_vibration_0.001")
        process_vtk("../vtk/fem/vibration/$bc/vib_fem_$(bc)_shear.vtu",
                    "fem_shear", bc, "../Fig/fem/vibration/fem_shear")
        process_vtk("../vtk/fem/vibration/$bc/vib_fem_$(bc)_shear_km_0.001.vtu",
                    "fem_shear_km", bc, "../Fig/fem/vibration/fem_shear_km")    
        process_vtk("../vtk/fem/buckling/$bc/fem_$bc.vtu",
                    "fem", bc, "../Fig/fem/buckling/fem_buckling")
        process_vtk("../vtk/fem/buckling/$bc/fem2_$(bc)_shear.vtu",
                    "fem_shear", bc, "../Fig/fem/buckling/fem_shear")
    end
end

# ================================================================
#  第二段：混合有限元（fem_mix）
# ================================================================
if DO_MIX
    println("="^50)
    println("DO_MIX：混合有限元")
    println("="^50)
    for bc in BCs
        process_vtk("../vtk/fem/vibration/$bc/vib_fem_mix_$(bc)_0.001.vtu",
                    "fem_mix", bc, "../Fig/fem/vibration/fem_mix_0.001")
        process_vtk("../vtk/fem/buckling/$bc/fem_mix_$bc.vtu",
                    "fem_mix", bc, "../Fig/fem/buckling/fem_mix")
    end
end

# ================================================================
#  第三段：無網格法（mf + mf2）
# ================================================================
if DO_MF
    println("="^50)
    println("DO_MF：無網格法")
    println("="^50)
    for bc in BCs
        process_vtk("../vtk/mf/vibration/$bc/vib_mf_w_φ_$(bc)_0.001.vtu",
                    "mf", bc, "../Fig/mf/vibration/mf_vibration_0.001")
        process_vtk("../vtk/mf/vibration/$bc/vib_mf_w_φ_$(bc)2_0.001.vtu",
                    "mf2", bc, "../Fig/mf/vibration/mf_vibration2_0.001")
        process_vtk("../vtk/mf/buckling/$bc/mf_w_φ_$bc.vtu",
                    "mf", bc, "../Fig/mf/buckling/mf_buckling")
        process_vtk("../vtk/mf/buckling/$bc/mf1_w_φ_$(bc).vtu",
                    "mf2", bc, "../Fig/mf/buckling/mf_buckling2")
    end
end

println("\n全部完成。")
