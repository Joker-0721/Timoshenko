using GLMakie, ReadVTK
using GeometryBasics: TriangleFace
using Base64: base64decode

# 讀點座標 + 全部 PointData 欄位（自製 base64 解碼，繞過 ReadVTK 讀 256 值時的 off-by-one）
function read_vtu_points(filename)
    txt  = read(filename, String)
    npts = parse(Int, match(r"NumberOfPoints=\"(\d+)\"", txt).captures[1])

    # Points（3 分量 Float64，跳前 8 bytes header）
    m = match(r"<DataArray[^>]*Name=\"Points\"[^>]*format=\"binary\"[^>]*>\s*([A-Za-z0-9+/=\n]+?)\s*</DataArray>", txt)
    pb = base64decode(strip(replace(m.captures[1], "\n" => "")))
    pts = reshape(reinterpret(Float64, pb[9:end]), 3, :)

    # PointData 欄位（Float64、1 分量、binary）
    fields = Dict{String,Vector{Float64}}()
    pat = r"<DataArray[^>]*type=\"Float64\"[^>]*Name=\"([^\"]+)\"[^>]*NumberOfComponents=\"1\"[^>]*format=\"binary\"[^>]*>\s*([A-Za-z0-9+/=\n]+?)\s*</DataArray>"
    for mm in eachmatch(pat, txt)
        name = mm.captures[1]
        name == "Points" && continue
        b = base64decode(strip(replace(mm.captures[2], "\n" => "")))
        d = reinterpret(Float64, b[9:end])
        if length(d) == npts - 1
            # 檔案 base64 尾端被截斷 8 bytes（最後一個 Float64），補齊
            d = [d; d[end]]
        elseif length(d) == npts + 1
            d = d[1:npts]
        end
        fields[name] = collect(d)
    end
    return pts, fields
end

# 讀三角形 faces（cells 不受 off-by-one 影響，用 ReadVTK）
function read_vtu_faces(filename)
    vtk    = ReadVTK.VTKFile(filename)
    mcells = ReadVTK.to_meshcells(ReadVTK.get_cells(vtk))
    return [TriangleFace{Int}(c.connectivity...) for c in mcells]
end

names12 = ["v₁","v₂","v₃","v₄","v₅","v₆","v₇","v₈","v₉","v₁₀","v₁₁","v₁₂"]

DO_FEM = false   # fem： binary、每 mode 一檔
DO_MIX = true   # mix： binary、每 mode 一檔（檔名 n 後有空格）

if DO_FEM
    BCs  = ["CCCC","CSCS","FSCS","FSSS","SCSC","SSSS"]
    nses = [7,9,11,13,15,17]
    for bc in BCs, n in nses, m in 1:6
        filename = "./vtk/fem_$(n)_$(bc)_mode_$(m).vtu"
        if !isfile(filename)
            println("skip: ", filename)
            continue
        end
        pts, fields = read_vtu_points(filename)
        faces = read_vtu_faces(filename)
        col = fields[names12[m]]
        fig = Figure(size=(800,600)); ax = Axis3(fig[1,1])
        plt = mesh!(ax, pts, faces, color=col, colormap=:jet, shading=false)
        ax.title = "fem_$(bc)_$(n) : $(names12[m])"
        Colorbar(fig[1,2], plt, label=names12[m])
        mkpath("Fig/fem/$(bc)")
        save("Fig/fem/$(bc)/fig_fem_$(bc)_$(n)_mode_$(m).png", fig)
    end
end

if DO_MIX
    BCs  = ["CCCC","CSCS","FSCS","FSSS","SCSC","SSSS"]
    nses = 7:17
    for bc in BCs, n in nses, m in 1:6
        # filename = "./vtk/mix/$(bc)/mix_$(n) _$(bc)_mode_$(m).vtu"   # 注意 n 後有空格
        filename = "./vtk/mix/$(bc)/mix_$(n)_$(bc)_roit1_mode_$(m).vtu" 
        if !isfile(filename)
            println("skip: ", filename)
            continue
        end
        pts, fields = read_vtu_points(filename)
        faces = read_vtu_faces(filename)
        col = fields[names12[m]]
        fig = Figure(size=(800,600)); ax = Axis3(fig[1,1])
        plt = mesh!(ax, pts, faces, color=col, colormap=:jet, shading=false)
        ax.title = "mix_$(bc)_roit_$(n) : $(names12[m])"
        Colorbar(fig[1,2], plt, label=names12[m])
        mkpath("Fig/mix/$(bc)")
        # save("Fig/mix/$(bc)/fig_mix_$(bc)_$(n)_mode_$(m).png", fig)
        save("Fig/mix1/$(bc)/fig_mix_$(bc)_roit1_$(n)_mode_$(m).png", fig)
    end
end