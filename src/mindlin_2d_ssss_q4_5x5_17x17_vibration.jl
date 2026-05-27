const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫αwwdΓ

using LinearAlgebra
using Printf
using TimerOutputs
using WriteVTK
import Gmsh: gmsh

# ==================== 參數設定 ====================
E = 1.0e8
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.1
h = h_over_b * b
ρ = 1.0
α = 1.0e8 * E

# 薄板彎曲剛度
Dᵇ = E * h^3 / 12 / (1 - ν^2)

eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12
residual_warn_tol = 1.0e-6

# 要分析的網格尺寸
mesh_sizes = [5, 17]

output_dir = normpath(joinpath(@__DIR__, "..", "vtk"))
mkpath(output_dir)
output_dir_csv = normpath(joinpath(@__DIR__, "..", "2d4s_vi_q4_FEM"))
mkpath(output_dir_csv)

const to = TimerOutput()

# ==================== 輔助函數 ====================
function generate_square_mesh!(nodes_per_side)
    gmsh.clear()
    gmsh.model.add("mindlin_2d_ssss_q4_vibration")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(a, 0.0, 0.0)
    p3 = gmsh.model.geo.addPoint(a, b, 0.0)
    p4 = gmsh.model.geo.addPoint(0.0, b, 0.0)

    bottom = gmsh.model.geo.addLine(p1, p2)
    right = gmsh.model.geo.addLine(p2, p3)
    top = gmsh.model.geo.addLine(p3, p4)
    left = gmsh.model.geo.addLine(p4, p1)
    loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(left, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(right, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteSurface(surface)
    gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()

    left_group = gmsh.model.addPhysicalGroup(1, [left])
    gmsh.model.setPhysicalName(1, left_group, "left")
    right_group = gmsh.model.addPhysicalGroup(1, [right])
    gmsh.model.setPhysicalName(1, right_group, "right")
    top_group = gmsh.model.addPhysicalGroup(1, [top])
    gmsh.model.setPhysicalName(1, top_group, "top")
    bottom_group = gmsh.model.addPhysicalGroup(1, [bottom])
    gmsh.model.setPhysicalName(1, bottom_group, "bottom")
    domain_group = gmsh.model.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, domain_group, "domain")

    gmsh.model.mesh.generate(2)
end

function assemble_consistent_mass!(mᵠᵠ, mʷʷ, elements, ρ, h)
    ρh = ρ * h
    ρI = ρ * h^3 / 12.0
    for elm in elements
        𝓒 = elm.𝓒
        𝓖 = elm.𝓖
        for ξ in 𝓖
            N = ξ[:𝝭]
            𝑤 = ξ.𝑤
            for (i, xᵢ) in enumerate(𝓒)
                I = xᵢ.𝐼
                for (j, xⱼ) in enumerate(𝓒)
                    J = xⱼ.𝐼
                    mʷʷ[I, J] += ρh * N[i] * N[j] * 𝑤
                    mᵠᵠ[2*I-1, 2*J-1] += ρI * N[i] * N[j] * 𝑤
                    mᵠᵠ[2*I, 2*J] += ρI * N[i] * N[j] * 𝑤
                end
            end
        end
    end
end

function selected_eigenpair(ω²ᵢ, vᵢ)
    return isfinite(real(ω²ᵢ)) &&
           isfinite(imag(ω²ᵢ)) &&
           abs(imag(ω²ᵢ)) < eigen_imag_tol &&
           real(ω²ᵢ) > omega_sq_tol &&
           all(isfinite, real.(vᵢ)) &&
           all(isfinite, imag.(vᵢ)) &&
           norm(imag.(vᵢ)) <= eigen_imag_tol * max(norm(real.(vᵢ)), eps(Float64))
end

function center_node_id(nodes)
    distances = [(node.x - a / 2)^2 + (node.y - b / 2)^2 for node in nodes]
    return nodes[argmin(distances)].𝐼
end

function normalize_mode(vᵢ, nodes, nᵠ)
    dm = real.(vᵢ)
    w = dm[2*nᵠ+1:end]
    max_abs_w = maximum(abs.(w))
    max_abs_w > 0.0 || error("mode has zero transverse displacement.")
    dm ./= max_abs_w

    center_id = center_node_id(nodes)
    if dm[2*nᵠ + center_id] < -sqrt(eps(Float64))
        dm .*= -1.0
    end
    return dm
end

function relative_residual(K, M, ω²ᵢ, vᵢ)
    Kv = K * vᵢ
    Mv = M * vᵢ
    residual = Kv - ω²ᵢ * Mv
    denominator = max(norm(Kv), abs(ω²ᵢ) * norm(Mv), eps(Float64))
    return norm(residual) / denominator
end

function vtk_cell(elm)
    node_ids = [xᵢ.𝐼 for xᵢ in elm.𝓒]
    if length(node_ids) == 4
        return MeshCell(VTKCellTypes.VTK_QUAD, node_ids)
    elseif length(node_ids) == 3
        return MeshCell(VTKCellTypes.VTK_TRIANGLE, node_ids)
    else
        error("unsupported cell with $(length(node_ids)) nodes")
    end
end

# 薄板 SSSS 的無因次頻率解析解
function exact_omega_sq(m, n, a, b, D, ρ, h)
    # Ω_exact = π²(m² + n²)  for a square plate (a=b)
    # ω = Ω * sqrt(D / (ρ*h)) / a²
    Ω_exact = π^2 * (m^2 + n^2)
    ω_exact = Ω_exact * sqrt(D / (ρ * h)) / a^2
    return ω_exact, Ω_exact
end

# ==================== 主程式：對每種網格進行分析 ====================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

all_results = Dict{Int, Any}()

for nodes_per_side in mesh_sizes
    println("\n" * "="^70)
    println("分析：$(nodes_per_side)×$(nodes_per_side) Q4 網格")
    println("="^70)

    case_prefix = "mindlin_2d_ssss_q4_$(nodes_per_side)x$(nodes_per_side)_vibration"
    vtu_path = joinpath(output_dir, "$(case_prefix)_modes.vtu")

    integrationOrder = 2
    integrationOrder_shear = 1

    @timeit to "generate Q4 square mesh" generate_square_mesh!(nodes_per_side)
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
    kᵠʷ = zeros(2 * nᵠ, nʷ)
    mʷʷ = zeros(nʷ, nʷ)
    mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    @timeit to "assemble stiffness and mass" begin
        @timeit to "get elements" elements = getElements(nodes, entities["domain"], integrationOrder)
        @timeit to "get shear elements" elements_s = getElements(nodes, entities["domain"], integrationOrder_shear)
        prescribe!(elements, :E => E, :ν => ν, :h => h)
        prescribe!(elements_s, :E => E, :ν => ν, :h => h)
        @timeit to "calculate shape functions" set∇𝝭!(elements)
        @timeit to "calculate shear shape functions" set∇𝝭!(elements_s)

        𝑎ʷʷ = ∫wwdΩ => elements_s
        𝑎ᵠʷ = ∫φwdΩ => elements_s
        𝑎ᵠᵠ = [
            ∫φφdΩ => elements_s,
            ∫κκdΩ => elements,
        ]
        @timeit to "assemble Kww" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble Kφw" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble Kφφ" 𝑎ᵠᵠ(kᵠᵠ)
        @timeit to "assemble consistent mass" assemble_consistent_mass!(mᵠᵠ, mʷʷ, elements, ρ, h)

        global elements_domain = elements
    end

    @timeit to "apply SSSS boundary (penalty)" begin
        @timeit to "get boundary elements" begin
            elements_left = getElements(nodes, entities["left"], integrationOrder)
            elements_right = getElements(nodes, entities["right"], integrationOrder)
            elements_top = getElements(nodes, entities["top"], integrationOrder)
            elements_bottom = getElements(nodes, entities["bottom"], integrationOrder)
        end

        w_boundary(x, y, z) = 0.0
        prescribe!(elements_left, :α => α, :g => w_boundary)
        prescribe!(elements_right, :α => α, :g => w_boundary)
        prescribe!(elements_top, :α => α, :g => w_boundary)
        prescribe!(elements_bottom, :α => α, :g => w_boundary)
        @timeit to "calculate boundary shape functions" begin
            set𝝭!(elements_left)
            set𝝭!(elements_right)
            set𝝭!(elements_top)
            set𝝭!(elements_bottom)
        end

        𝑎ʷ = ∫αwwdΓ => elements_left ∪ elements_right ∪ elements_top ∪ elements_bottom
        fᵅ = zeros(nʷ)
        @timeit to "assemble boundary penalty" 𝑎ʷ(kʷʷ, fᵅ)
    end

    @timeit to "solve vibration eigenvalue" begin
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        M = [
            mᵠᵠ zeros(2 * nᵠ, nʷ)
            zeros(nʷ, 2 * nᵠ) mʷʷ
        ]
        Ks = Symmetric(0.5 * (K + K'))
        Ms = Symmetric(0.5 * (M + M'))
        F = eigen(Ks, Ms)
        ω² = F.values
        V = F.vectors

        mode_ids = sort!(
            collect(i for i in eachindex(ω²) if selected_eigenpair(ω²[i], V[:, i])),
            by = i -> real(ω²[i]),
        )
        isempty(mode_ids) && error("no positive finite vibration eigenvalue found")

        dm_modes = [normalize_mode(V[:, mode_id], nodes, nᵠ) for mode_id in mode_ids]
        residuals = [
            relative_residual(K, M, real(ω²[mode_id]), dm_modes[mode_rank])
            for (mode_rank, mode_id) in enumerate(mode_ids)
        ]

        ω²₁ = real(ω²[first(mode_ids)])
        ω₁ = sqrt(ω²₁)

        # 儲存結果
        all_results[nodes_per_side] = (
            nodes = nodes,
            elements_domain = elements_domain,
            ω² = ω²,
            V = V,
            mode_ids = mode_ids,
            dm_modes = dm_modes,
            residuals = residuals,
            nᵠ = nᵠ,
        )
    end

    # ==================== 輸出結果 ====================
    println("\n結果摘要 ($(nodes_per_side)×$(nodes_per_side)):")
    println("  節點數: ", length(nodes))
    println("  選取的有效模態數: ", length(mode_ids))
    println("  ω²₁ (基頻平方): ", ω²₁)
    println("  ω₁  (基頻): ", ω₁)
    println("  最大相對殘差: ", maximum(residuals))
    if maximum(residuals) >= residual_warn_tol
        @printf("  警告: 部分殘差超過 %.1e\n", residual_warn_tol)
    end

    # ---- 模態比較表 ----
    println("\n  模態比較 (無因次頻率 Ω = ω·a²·√(ρ·h/Dᵇ)):")
    println("  " * @sprintf("%-8s %-8s %-14s %-14s %-14s %-12s",
                            "mode", "(m,n)", "Ω_exact", "Ω_FEM", "error(%)", "residual"))
    n_modes_show = min(length(mode_ids), 15)
    for r in 1:n_modes_show
        mode_id = mode_ids[r]
        ω²ᵢ = real(ω²[mode_id])
        ωᵢ = sqrt(ω²ᵢ)
        Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)

        # 找出最接近的 (m,n) 配對
        best_err = Inf
        best_mn = (0, 0)
        for m in 1:10
            for n in 1:10
                _, Ω_ex = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
                err = abs(Ω_FEM - Ω_ex) / Ω_ex * 100
                if err < best_err
                    best_err = err
                    best_mn = (m, n)
                end
            end
        end
        m, n = best_mn
        _, Ω_ex = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
        println("  " * @sprintf("%-8d (%-d,%-d)     %-10.4f  %-10.4f  %-10.4f  %-10.2e",
                                r, m, n, Ω_ex, Ω_FEM, best_err, residuals[r]))
    end

    # ---- 寫入 VTK ----
    cells = [vtk_cell(elm) for elm in elements_domain]
    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for node in nodes
        points[1, node.𝐼] = node.x
        points[2, node.𝐼] = node.y
        points[3, node.𝐼] = 0.0
    end

    # 限制輸出模態數量 = 節點數（與 exact 一致）
    n_modes_output = min(length(dm_modes), nₚ)

    vtk_grid(vtu_path, points, cells; ascii = true, append = false, compress = false) do vtk
        ω²_selected = [real(ω²[mode_id]) for mode_id in mode_ids[1:n_modes_output]]
        Ω_selected = [sqrt(ω²ᵢ) * a^2 * sqrt(ρ * h / Dᵇ) for ω²ᵢ in ω²_selected]
        vtk["omega_sq", WriteVTK.VTKFieldData()] = ω²_selected
        vtk["omega", WriteVTK.VTKFieldData()] = sqrt.(ω²_selected)
        vtk["Omega_dimensionless", WriteVTK.VTKFieldData()] = Ω_selected
        vtk["relative_residual", WriteVTK.VTKFieldData()] = residuals[1:n_modes_output]
        vtk["rho", WriteVTK.VTKFieldData()] = [ρ]
        vtk["alpha_boundary_penalty", WriteVTK.VTKFieldData()] = [α]
        vtk["color_range_min", WriteVTK.VTKFieldData()] = [-1.0]
        vtk["color_range_max", WriteVTK.VTKFieldData()] = [1.0]

        for (mode_rank, dm) in enumerate(dm_modes[1:n_modes_output])
            w = dm[2*nᵠ+1:end]
            φx = dm[1:2:2*nᵠ]
            φy = dm[2:2:2*nᵠ]
            vtk["w$(mode_rank)"] = w
            vtk["φx$(mode_rank)"] = φx
            vtk["φy$(mode_rank)"] = φy
        end
    end
    println("  VTU 輸出模態數: $n_modes_output")
    println("  VTU 輸出: ", vtu_path)

    # ---- 寫入 CSV 摘要 ----
    csv_path = joinpath(output_dir,
                        "mindlin_2d_ssss_q4_$(nodes_per_side)x$(nodes_per_side)_vibration_modes.csv")
    open(csv_path, "w") do io
        println(io, join([
            "mode_rank", "eigen_index",
            "omega_sq", "omega",
            "Omega_dimensionless",
            "m_match", "n_match", "Omega_exact",
            "rel_error_percent",
            "relative_residual",
            "residual_ok",
        ], ","))
        for r in 1:n_modes_show
            mode_id = mode_ids[r]
            ω²ᵢ = real(ω²[mode_id])
            ωᵢ = sqrt(ω²ᵢ)
            Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)

            best_err = Inf
            best_mn = (0, 0)
            best_Ω_ex = 0.0
            for m in 1:10
                for n in 1:10
                    _, Ω_ex = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
                    err = abs(Ω_FEM - Ω_ex) / Ω_ex * 100
                    if err < best_err
                        best_err = err
                        best_mn = (m, n)
                        best_Ω_ex = Ω_ex
                    end
                end
            end
            m, n = best_mn

            println(io, join([
                r, mode_id,
                ω²ᵢ, ωᵢ,
                Ω_FEM,
                m, n, best_Ω_ex,
                best_err,
                residuals[r],
                residuals[r] < residual_warn_tol,
            ], ","))
        end
    end
    println("  CSV 輸出: ", csv_path)

    # ---- 輸出節點數據 CSV（節點座標 + 所有模態的 w, φx, φy）----
    csv_nodal_path = joinpath(output_dir_csv,
                              "2d4s_vi_q4_FEM_$(nodes_per_side)x$(nodes_per_side).csv")
    open(csv_nodal_path, "w") do io
        # 表頭
        header = ["node_id", "x", "y"]
        for k in 1:n_modes_output
            push!(header, "w$k", "φx$k", "φy$k")
        end
        println(io, join(header, ","))

        # 每個節點的數據
        for node in nodes
            i = node.𝐼
            row = [string(i), string(node.x), string(node.y)]
            for mode_rank in 1:n_modes_output
                dm = dm_modes[mode_rank]
                w = dm[2*nᵠ+1:end][i]
                φx = dm[1:2:2*nᵠ][i]
                φy = dm[2:2:2*nᵠ][i]
                push!(row, string(w), string(φx), string(φy))
            end
            println(io, join(row, ","))
        end
    end
    println("  節點數據 CSV 輸出: ", csv_nodal_path)
end

gmsh.finalize()

println("\n" * "="^70)
println("所有分析完成！")
println("="^70)
println(to)
