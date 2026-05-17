const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
using LinearAlgebra
using Printf
using WriteVTK

import Gmsh: gmsh
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ

const CASE_PREFIX = "mindlin_2d_ssss_q4_5x5_buckling"
const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "vtk"))

const E = 1.0e8
const ν = 0.3
const a = 1.0
const b = 1.0
const h_over_b = 0.1
const h = h_over_b*b
const Dᵇ = E*h^3/(12.0*(1.0 - ν^2))
const α = 1.0e8*E

const NODES_PER_SIDE = 5
const K_REF_FIRST_MODE = 4.0
const EIGEN_IMAG_TOL = 1.0e-7
const RESIDUAL_TOL = 1.0e-6

struct BucklingSolution
    nodes
    elements
    K::Matrix{Float64}
    KG::Matrix{Float64}
    lambdas::Vector{Float64}
    full_modes::Vector{Vector{Float64}}
    full_residuals::Vector{Float64}
end

function bending_k(lambda::Float64)
    return lambda*b^2/(π^2*Dᵇ)
end

function generate_square_mesh!()
    gmsh.clear()
    gmsh.model.add(CASE_PREFIX)

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

    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, NODES_PER_SIDE)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, NODES_PER_SIDE)
    gmsh.model.geo.mesh.setTransfiniteCurve(left, NODES_PER_SIDE)
    gmsh.model.geo.mesh.setTransfiniteCurve(right, NODES_PER_SIDE)
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

function assemble_boundary_penalty!(kʷʷ, nodes, entities)
    fᵅ = zeros(size(kʷʷ, 1))
    w_boundary(x, y, z) = 0.0
    for side in ("left", "right", "top", "bottom")
        elements = getElements(nodes, entities[side])
        prescribe!(elements, :α => α, :g => w_boundary)
        set𝝭!(elements)
        (∫αwwdΓ => elements)(kʷʷ, fᵅ)
    end
    return fᵅ
end

function assemble_case()
    generate_square_mesh!()

    entities = getPhysicalGroups()
    for name in ("domain", "left", "right", "top", "bottom")
        haskey(entities, name) || error("missing physical group: $name")
    end

    nodes = get𝑿ᵢ()
    nʷ = length(nodes)
    nᵠ = length(nodes)

    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
    kᵠʷ = zeros(2*nᵠ, nʷ)
    kᴳʷʷ = zeros(nʷ, nʷ)
    kᴳᵠᵠ = zeros(2*nᵠ, 2*nᵠ)

    elements = getElements(nodes, entities["domain"])
    prescribe!(elements, :E => E, :ν => ν, :h => h)
    set∇𝝭!(elements)
    (∫wwdΩ => elements)(kʷʷ)
    (∫φwdΩ => elements)(kᵠʷ)
    ([∫φφdΩ => elements, ∫κκdΩ => elements])(kᵠᵠ)

    σ₁₁(x, y, z) = 1.0
    σ₂₂(x, y, z) = 0.0
    σ₁₂(x, y, z) = 0.0
    prescribe!(elements, :σ₁₁ => σ₁₁, :σ₂₂ => σ₂₂, :σ₁₂ => σ₁₂)
    (∫∇wσ∇wdΩ => elements)(kᴳʷʷ)
    (∫∇φσ∇φdΩ => elements)(kᴳᵠᵠ)

    assemble_boundary_penalty!(kʷʷ, nodes, entities)

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    KG = [
        kᴳᵠᵠ zeros(2*nᵠ, nʷ)
        zeros(nʷ, 2*nᵠ) kᴳʷʷ
    ]

    return nodes, elements, K, KG
end

function is_real_finite_value(lambda)
    return isfinite(real(lambda)) &&
           isfinite(imag(lambda)) &&
           abs(imag(lambda)) < EIGEN_IMAG_TOL &&
           real(lambda) > 0.0
end

function is_real_finite_vector(v)
    all(isfinite, real.(v)) || return false
    all(isfinite, imag.(v)) || return false
    return norm(imag.(v)) <= EIGEN_IMAG_TOL*max(norm(real.(v)), eps(Float64))
end

function center_node_id(nodes)
    distances = [(node.x - a/2)^2 + (node.y - b/2)^2 for node in nodes]
    return nodes[argmin(distances)].𝐼
end

function normalize_mode(v, nodes)
    nᵠ = length(nodes)
    mode = real.(v)
    w = mode[2*nᵠ + 1:end]
    max_abs_w = maximum(abs.(w))
    max_abs_w > 0.0 || error("mode has zero transverse displacement.")
    mode ./= max_abs_w

    center_id = center_node_id(nodes)
    if mode[2*nᵠ + center_id] < -sqrt(eps(Float64))
        mode .*= -1.0
    end
    return mode
end

function relative_residual(K, KG, lambda::Float64, v::Vector{Float64})
    Kv = K*v
    KGv = KG*v
    residual = Kv - lambda*KGv
    denominator = max(norm(Kv), abs(lambda)*norm(KGv), eps(Float64))
    return norm(residual)/denominator
end

function solve_full_eigen(nodes, K, KG)
    F = eigen(K, KG)
    ids = [
        i for i in eachindex(F.values)
        if is_real_finite_value(F.values[i]) && is_real_finite_vector(F.vectors[:, i])
    ]
    sort!(ids, by = i -> real(F.values[i]))
    isempty(ids) && error("no positive finite buckling eigenvalue found.")

    lambdas = Float64[real(F.values[i]) for i in ids]
    full_modes = [normalize_mode(F.vectors[:, i], nodes) for i in ids]
    full_residuals = [
        relative_residual(K, KG, lambdas[i], full_modes[i])
        for i in eachindex(lambdas)
    ]
    return lambdas, full_modes, full_residuals
end

function solve_case()
    nodes, elements, K, KG = assemble_case()
    lambdas, full_modes, full_residuals = solve_full_eigen(nodes, K, KG)
    return BucklingSolution(nodes, elements, K, KG, lambdas, full_modes, full_residuals)
end

csv_line(values) = join(values, ",")

function write_eigen_check(path, solution::BucklingSolution)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, csv_line((
            "mode",
            "lambda",
            "k_num",
            "k_ref",
            "relative_error_vs_ref",
            "full_relative_residual",
            "lambda_isfinite",
            "lambda_is_positive",
            "vector_isfinite",
            "selected_positive_mode",
        )))

        for mode in eachindex(solution.lambdas)
            lambda = solution.lambdas[mode]
            k_num = bending_k(lambda)
            rel_error = abs(k_num - K_REF_FIRST_MODE)/abs(K_REF_FIRST_MODE)
            vector_isfinite = all(isfinite, solution.full_modes[mode])
            println(io, csv_line((
                mode,
                lambda,
                k_num,
                K_REF_FIRST_MODE,
                rel_error,
                solution.full_residuals[mode],
                isfinite(lambda),
                lambda > 0.0,
                vector_isfinite,
                true,
            )))
        end
    end
end

function write_mode_data(summary_path, node_path, solution::BucklingSolution)
    nodes = solution.nodes
    nᵠ = length(nodes)
    mkpath(dirname(summary_path))

    open(summary_path, "w") do summary_io
        println(summary_io, csv_line((
            "mode",
            "vtk_w_field",
            "lambda",
            "k_num",
            "k_ref",
            "w_min",
            "w_min_node",
            "w_max",
            "w_max_node",
            "max_abs_w",
            "w_norm",
            "phi_1_norm",
            "phi_2_norm",
            "normalized_max_abs_w",
            "rebuild_error",
            "relative_rebuild_error",
        )))

        open(node_path, "w") do node_io
            println(node_io, csv_line((
                "mode",
                "vtk_w_field",
                "lambda",
                "k_num",
                "node_row",
                "node_id",
                "x",
                "y",
                "w_normalized",
                "phi_1",
                "phi_2",
            )))

            for mode in eachindex(solution.full_modes)
                dm = solution.full_modes[mode]
                phi_1 = dm[1:2:2*nᵠ]
                phi_2 = dm[2:2:2*nᵠ]
                w = dm[2*nᵠ + 1:end]
                rebuilt = zeros(3*nᵠ)
                rebuilt[1:2:2*nᵠ] .= phi_1
                rebuilt[2:2:2*nᵠ] .= phi_2
                rebuilt[2*nᵠ + 1:end] .= w
                rebuild_error = norm(rebuilt - dm)
                relative_rebuild_error = rebuild_error/max(norm(dm), eps(Float64))

                w_min, w_min_node = findmin(w)
                w_max, w_max_node = findmax(w)
                max_abs_w = maximum(abs.(w))
                lambda = solution.lambdas[mode]
                k_num = bending_k(lambda)

                println(summary_io, csv_line((
                    mode,
                    "w$(mode)",
                    lambda,
                    k_num,
                    K_REF_FIRST_MODE,
                    w_min,
                    w_min_node,
                    w_max,
                    w_max_node,
                    max_abs_w,
                    norm(w),
                    norm(phi_1),
                    norm(phi_2),
                    isapprox(max_abs_w, 1.0; atol=1.0e-12, rtol=1.0e-12),
                    rebuild_error,
                    relative_rebuild_error,
                )))

                for (node_row, node) in enumerate(nodes)
                    node_id = node.𝐼
                    println(node_io, csv_line((
                        mode,
                        "w$(mode)",
                        lambda,
                        k_num,
                        node_row,
                        node_id,
                        node.x,
                        node.y,
                        w[node_id],
                        phi_1[node_id],
                        phi_2[node_id],
                    )))
                end
            end
        end
    end
end

function vtk_cell(element)
    node_ids = [node.𝐼 for node in element.𝓒]
    if length(node_ids) == 4
        return MeshCell(VTKCellTypes.VTK_QUAD, node_ids)
    elseif length(node_ids) == 3
        return MeshCell(VTKCellTypes.VTK_TRIANGLE, node_ids)
    else
        error("unsupported cell with $(length(node_ids)) nodes")
    end
end

function write_vtu(path, solution::BucklingSolution)
    nodes = solution.nodes
    nᵠ = length(nodes)
    points = zeros(3, nᵠ)
    for node in nodes
        points[1, node.𝐼] = node.x
        points[2, node.𝐼] = node.y
        points[3, node.𝐼] = 0.0
    end
    cells = [vtk_cell(element) for element in solution.elements]

    mkpath(dirname(path))
    vtk_grid(path, points, cells; ascii=true, append=false, compress=false) do vtk
        vtk["lambda", WriteVTK.VTKFieldData()] = solution.lambdas
        vtk["k_num", WriteVTK.VTKFieldData()] = bending_k.(solution.lambdas)
        vtk["k_ref_first_mode", WriteVTK.VTKFieldData()] = fill(K_REF_FIRST_MODE, length(solution.lambdas))
        vtk["full_relative_residual", WriteVTK.VTKFieldData()] = solution.full_residuals
        vtk["alpha_boundary_penalty", WriteVTK.VTKFieldData()] = [α]
        vtk["color_range_min", WriteVTK.VTKFieldData()] = [-1.0]
        vtk["color_range_max", WriteVTK.VTKFieldData()] = [1.0]

        for mode in eachindex(solution.full_modes)
            dm = solution.full_modes[mode]
            w = dm[2*nᵠ + 1:end]
            phi_1 = dm[1:2:2*nᵠ]
            phi_2 = dm[2:2:2*nᵠ]
            vtk["w$(mode)"] = w
            vtk["phi_1_$(mode)"] = phi_1
            vtk["phi_2_$(mode)"] = phi_2
        end
    end
end

function write_outputs(solution::BucklingSolution)
    eigen_path = joinpath(OUTPUT_DIR, "$(CASE_PREFIX)_eigen_check.csv")
    summary_path = joinpath(OUTPUT_DIR, "$(CASE_PREFIX)_mode_summary.csv")
    node_path = joinpath(OUTPUT_DIR, "$(CASE_PREFIX)_mode_node_values.csv")
    vtu_path = joinpath(OUTPUT_DIR, "$(CASE_PREFIX)_modes.vtu")

    write_eigen_check(eigen_path, solution)
    write_mode_data(summary_path, node_path, solution)
    write_vtu(vtu_path, solution)
    return eigen_path, summary_path, node_path, vtu_path
end

function assert_solution(solution::BucklingSolution)
    !isempty(solution.lambdas) || error("no selected finite positive modes")
    all(isfinite, solution.lambdas) || error("non-finite eigenvalue detected")
    all(>(0.0), solution.lambdas) || error("non-positive eigenvalue detected")
    maximum(solution.full_residuals) < RESIDUAL_TOL ||
        error("full residual exceeds tolerance: $(maximum(solution.full_residuals))")

    nᵠ = length(solution.nodes)
    for (mode, dm) in enumerate(solution.full_modes)
        w = dm[2*nᵠ + 1:end]
        isapprox(maximum(abs.(w)), 1.0; atol=1.0e-12, rtol=1.0e-12) ||
            error("mode $mode is not normalized by max(abs(w))")
    end
end

function main()
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        solution = solve_case()
        assert_solution(solution)
        eigen_path, summary_path, node_path, vtu_path = write_outputs(solution)

        @printf("2D Mindlin SSSS Q4 5x5 buckling case\n")
        @printf("nodes = %d, selected positive finite modes = %d\n", length(solution.nodes), length(solution.lambdas))
        @printf("h/b = %.6f, E = %.6e, nu = %.6f\n", h_over_b, E, ν)
        @printf("first lambda = %.12e\n", solution.lambdas[1])
        @printf("first k_num  = %.12f (thin-plate reference %.6f)\n", bending_k(solution.lambdas[1]), K_REF_FIRST_MODE)
        @printf("max full residual = %.6e\n", maximum(solution.full_residuals))
        println("eigen check csv: $eigen_path")
        println("mode summary csv: $summary_path")
        println("mode node values csv: $node_path")
        println("VTU output: $vtu_path")
    finally
        gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__ || !isinteractive()
    main()
end
