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

case_prefix = "mindlin_2d_ssss_q4_17x17_vibration"
output_dir = normpath(joinpath(@__DIR__, "..", "vtk"))

E = 1.0e8
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.1
h = h_over_b*b
ρ = 1.0
α = 1.0e8*E

nodes_per_side = 17
eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12
residual_warn_tol = 1.0e-6

eigen_check_path = joinpath(output_dir, "$(case_prefix)_eigen_check.csv")
mode_summary_path = joinpath(output_dir, "$(case_prefix)_mode_summary.csv")
mode_node_path = joinpath(output_dir, "$(case_prefix)_mode_node_values.csv")
vtu_path = joinpath(output_dir, "$(case_prefix)_modes.vtu")

const to = TimerOutput()

function generate_square_mesh!(nodes_per_side)
    gmsh.clear()
    gmsh.model.add("mindlin_2d_ssss_q4_17x17_vibration")

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

function assemble_consistent_mass!(mᵠᵠ, mʷʷ, elements)
    ρh = ρ*h
    ρI = ρ*h^3/12.0
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
                    mʷʷ[I, J] += ρh*N[i]*N[j]*𝑤
                    mᵠᵠ[2*I-1, 2*J-1] += ρI*N[i]*N[j]*𝑤
                    mᵠᵠ[2*I, 2*J] += ρI*N[i]*N[j]*𝑤
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
           norm(imag.(vᵢ)) <= eigen_imag_tol*max(norm(real.(vᵢ)), eps(Float64))
end

function center_node_id(nodes)
    distances = [(node.x - a/2)^2 + (node.y - b/2)^2 for node in nodes]
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
    Kv = K*vᵢ
    Mv = M*vᵢ
    residual = Kv - ω²ᵢ*Mv
    denominator = max(norm(Kv), abs(ω²ᵢ)*norm(Mv), eps(Float64))
    return norm(residual)/denominator
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

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

integrationOrder = 2
integrationOrder_shear = 1

@timeit to "generate Q4 square mesh" generate_square_mesh!(nodes_per_side)
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ, nʷ)
kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
kᵠʷ = zeros(2*nᵠ, nʷ)
mʷʷ = zeros(nʷ, nʷ)
mᵠᵠ = zeros(2*nᵠ, 2*nᵠ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, M" begin
    @timeit to "get elements" elements = getElements(nodes, entities["domain"], integrationOrder)
    @timeit to "get shear elements" elements_s = getElements(nodes, entities["domain"], integrationOrder_shear)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_s, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shear shape functions" set∇𝝭!(elements_s)

    𝑎ʷʷ = ∫wwdΩ=>elements_s
    𝑎ᵠʷ = ∫φwdΩ=>elements_s
    𝑎ᵠᵠ = [
        ∫φφdΩ=>elements_s,
        ∫κκdΩ=>elements,
    ]
    @timeit to "assemble Kww" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble Kφw" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble Kφφ" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble consistent mass" assemble_consistent_mass!(mᵠᵠ, mʷʷ, elements)

    global elements_domain = elements
end

@timeit to "calculate ∫αwwdΓ" begin
    @timeit to "get left elements" elements_left = getElements(nodes, entities["left"], integrationOrder)
    @timeit to "get right elements" elements_right = getElements(nodes, entities["right"], integrationOrder)
    @timeit to "get top elements" elements_top = getElements(nodes, entities["top"], integrationOrder)
    @timeit to "get bottom elements" elements_bottom = getElements(nodes, entities["bottom"], integrationOrder)

    w_boundary(x, y, z) = 0.0
    prescribe!(elements_left, :α=>α, :g=>w_boundary)
    prescribe!(elements_right, :α=>α, :g=>w_boundary)
    prescribe!(elements_top, :α=>α, :g=>w_boundary)
    prescribe!(elements_bottom, :α=>α, :g=>w_boundary)
    @timeit to "calculate left shape functions" set𝝭!(elements_left)
    @timeit to "calculate right shape functions" set𝝭!(elements_right)
    @timeit to "calculate top shape functions" set𝝭!(elements_top)
    @timeit to "calculate bottom shape functions" set𝝭!(elements_bottom)

    𝑎ʷ = ∫αwwdΓ=>elements_left∪elements_right∪elements_top∪elements_bottom
    fᵅ = zeros(nʷ)
    @timeit to "assemble boundary penalty" 𝑎ʷ(kʷʷ, fᵅ)
end

@timeit to "solve vibration eigenvalue" begin
    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    M = [
        mᵠᵠ zeros(2*nᵠ, nʷ)
        zeros(nʷ, 2*nᵠ) mʷʷ
    ]
    Ks = Symmetric(0.5*(K + K'))
    Ms = Symmetric(0.5*(M + M'))
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
end

gmsh.finalize()

println(to)

println("2D Mindlin SSSS Q4 5x5 vibration case")
println("nodes: ", length(nodes))
println("selected positive finite modes: ", length(mode_ids))
println("ρ: ", ρ)
println("ω² first: ", ω²₁)
println("ω first: ", ω₁)
println("max relative residual: ", maximum(residuals))
if maximum(residuals) >= residual_warn_tol
    @printf("warning: some residuals exceed %.1e; check %s\n", residual_warn_tol, eigen_check_path)
end

function write_eigen_check(filepath, K, M, ω², V, mode_ids, dm_modes, residuals)
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        println(io, join([
            "mode_rank",
            "eigen_index",
            "omega_sq",
            "omega",
            "omega_sq_isfinite",
            "omega_sq_is_positive",
            "vector_isfinite",
            "norm_v",
            "norm_Kv",
            "norm_Mv",
            "residual_norm",
            "relative_residual",
            "residual_within_tolerance",
            "Mv_over_v",
            "selected_positive_mode",
        ], ","))

        selected = Set(mode_ids)
        for (mode_rank, mode_id) in enumerate(mode_ids)
            ω²ᵢ = real(ω²[mode_id])
            dm = dm_modes[mode_rank]
            Kv = K*dm
            Mv = M*dm
            residual = Kv - ω²ᵢ*Mv
            residual_norm = norm(residual)
            norm_v = norm(dm)
            norm_Kv = norm(Kv)
            norm_Mv = norm(Mv)
            relative_residual = residuals[mode_rank]
            Mv_over_v = norm_Mv/max(norm_v, eps(Float64))

            println(io, join([
                mode_rank,
                mode_id,
                ω²ᵢ,
                sqrt(ω²ᵢ),
                isfinite(ω²ᵢ),
                ω²ᵢ > 0.0,
                all(isfinite, dm),
                norm_v,
                norm_Kv,
                norm_Mv,
                residual_norm,
                relative_residual,
                relative_residual < residual_warn_tol,
                Mv_over_v,
                mode_id in selected,
            ], ","))
        end
    end
end

write_eigen_check(eigen_check_path, K, M, ω², V, mode_ids, dm_modes, residuals)
println("eigen check csv: ", eigen_check_path)

function write_mode_data(summary_path, node_path, nodes, ω², nᵠ, mode_ids, dm_modes)
    mkpath(dirname(summary_path))

    open(summary_path, "w") do summary_io
        println(summary_io, join([
            "mode_rank",
            "eigen_index",
            "vtk_w_field",
            "omega_sq",
            "omega",
            "selected_positive_mode",
            "vector_isfinite",
            "vector_norm",
            "w_norm",
            "w_min",
            "w_min_node",
            "w_max",
            "w_max_node",
            "max_abs_w",
            "normalized_max_abs_w",
            "rebuild_error",
            "relative_rebuild_error",
        ], ","))

        open(node_path, "w") do node_io
            println(node_io, join([
                "mode_rank",
                "eigen_index",
                "vtk_w_field",
                "omega_sq",
                "omega",
                "node_row",
                "node_id",
                "x",
                "y",
                "w_normalized",
            ], ","))

            for (mode_rank, mode_id) in enumerate(mode_ids)
                ω²ᵢ = real(ω²[mode_id])
                dm = dm_modes[mode_rank]
                w = dm[2*nᵠ+1:end]

                rebuilt = similar(dm)
                rebuilt[1:2*nᵠ] .= dm[1:2*nᵠ]
                rebuilt[2*nᵠ+1:end] .= w
                rebuild_error = norm(rebuilt - dm)
                relative_rebuild_error = rebuild_error/max(norm(dm), eps(Float64))

                w_min, w_min_node = findmin(w)
                w_max, w_max_node = findmax(w)
                max_abs_w = maximum(abs.(w))

                println(summary_io, join([
                    mode_rank,
                    mode_id,
                    "w$(mode_rank)",
                    ω²ᵢ,
                    sqrt(ω²ᵢ),
                    true,
                    all(isfinite, dm),
                    norm(dm),
                    norm(w),
                    w_min,
                    w_min_node,
                    w_max,
                    w_max_node,
                    max_abs_w,
                    isapprox(max_abs_w, 1.0; atol=1.0e-12, rtol=1.0e-12),
                    rebuild_error,
                    relative_rebuild_error,
                ], ","))

                for (node_row, node) in enumerate(nodes)
                    node_id = node.𝐼
                    println(node_io, join([
                        mode_rank,
                        mode_id,
                        "w$(mode_rank)",
                        ω²ᵢ,
                        sqrt(ω²ᵢ),
                        node_row,
                        node_id,
                        node.x,
                        node.y,
                        w[node_id],
                    ], ","))
                end
            end
        end
    end
end

write_mode_data(mode_summary_path, mode_node_path, nodes, ω², nᵠ, mode_ids, dm_modes)
println("mode summary csv: ", mode_summary_path)
println("mode node values csv: ", mode_node_path)

cells = [vtk_cell(elm) for elm in elements_domain]
nₚ = length(nodes)
points = zeros(3, nₚ)
for node in nodes
    points[1, node.𝐼] = node.x
    points[2, node.𝐼] = node.y
    points[3, node.𝐼] = 0.0
end

vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
    ω²_selected = [real(ω²[mode_id]) for mode_id in mode_ids]
    vtk["omega_sq", WriteVTK.VTKFieldData()] = ω²_selected
    vtk["omega", WriteVTK.VTKFieldData()] = sqrt.(ω²_selected)
    vtk["relative_residual", WriteVTK.VTKFieldData()] = residuals
    vtk["rho", WriteVTK.VTKFieldData()] = [ρ]
    vtk["alpha_boundary_penalty", WriteVTK.VTKFieldData()] = [α]
    vtk["color_range_min", WriteVTK.VTKFieldData()] = [-1.0]
    vtk["color_range_max", WriteVTK.VTKFieldData()] = [1.0]

    for (mode_rank, dm) in enumerate(dm_modes)
        w = dm[2*nᵠ+1:end]
        vtk["w$(mode_rank)"] = w
    end
end
println("VTU output: ", vtu_path)
