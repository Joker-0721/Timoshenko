using ApproxOperator
using LinearAlgebra

import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Timoshenko: ∫κEIκds, ∫γkGAγds, ∫vqds, L₂, L₂φ, w_exact_ss
import Gmsh: gmsh

const MODEL_NAME = "beam1d"
const RAW_CSV_PATH = normpath(joinpath(@__DIR__, "..", "data", "static_bending_1d.csv"))
const OUTPUT_COLUMNS = (
    :model,
    :element_type,
    :integration_strategy,
    :gauss_bending,
    :gauss_shear,
    :gauss_L2,
    :h_over_L,
    :mesh_level,
    :bc,
    :reference_kind,
    :width_over_L,
    :E,
    :ν,
    :κ,
    :L,
    :q,
    :w_ref,
    :w_fem,
    :L2_w_pct,
    :L2_φ_pct,
    :dofs,
    :runtime_s,
    :status,
    :note,
)

const GAUSS_SET = (1, 2, 3, 5, 10, 15)
const H_OVER_L_SET = (0.1, 0.2, 0.4)
const MESH_LEVEL_SET = (10, 20, 40)
const BC_SET = (:SS,)

const E0 = 1.0e8
const ν0 = 0.3
const κ0 = 5.0 / 6.0
const L0 = 1.0
const q0 = 1.0
const WIDTH_OVER_L0 = 1.0

section_width() = WIDTH_OVER_L0 * L0
beam_thickness(h_over_L::Float64) = h_over_L * L0
observation_x(bc::Symbol) = bc == :SS ? L0 / 2.0 : L0

function section_data(h_over_L::Float64)
    width = section_width()
    h = beam_thickness(h_over_L)
    G = E0 / (2.0 * (1.0 + ν0))
    A = width * h
    I = width * h^3 / 12.0
    EI = E0 * I
    kGA = κ0 * G * A
    return (; width, h, G, A, I, EI, kGA)
end

function w_exact_cf(x::Float64, E::Float64, I::Float64, κ::Float64, G::Float64, A::Float64, L::Float64, q::Float64)
    w_b = q * x^2 * (6.0 * L^2 - 4.0 * L * x + x^2) / (24.0 * E * I)
    w_s = q * x * (2.0 * L - x) / (2.0 * κ * G * A)
    return w_b + w_s
end

function φ_exact_ss(x::Float64, E::Float64, I::Float64, L::Float64, q::Float64)
    return q * (L^3 - 6.0 * L * x^2 + 4.0 * x^3) / (24.0 * E * I)
end

function φ_exact_cf(x::Float64, E::Float64, I::Float64, L::Float64, q::Float64)
    return q * x * (3.0 * L^2 - 3.0 * L * x + x^2) / (6.0 * E * I)
end

function w_exact(bc::Symbol, x::Float64, h_over_L::Float64)
    s = section_data(h_over_L)
    if bc == :SS
        return w_exact_ss(x, E0, s.I, κ0, s.G, s.A, L0, q0)
    elseif bc == :CF
        return w_exact_cf(x, E0, s.I, κ0, s.G, s.A, L0, q0)
    end
    error("unsupported boundary condition $bc")
end

function φ_exact(bc::Symbol, x::Float64, h_over_L::Float64)
    s = section_data(h_over_L)
    if bc == :SS
        return φ_exact_ss(x, E0, s.I, L0, q0)
    elseif bc == :CF
        return φ_exact_cf(x, E0, s.I, L0, q0)
    end
    error("unsupported boundary condition $bc")
end

function line_gauss_rule(n::Int)
    n >= 1 || error("Gauss point count must be positive.")
    if n == 1
        return [0.0], [2.0]
    end
    β = [k / sqrt(4.0 * k^2 - 1.0) for k in 1:(n - 1)]
    J = SymTridiagonal(zeros(n), β)
    λ, V = eigen(J)
    return collect(λ), [2.0 * V[1, i]^2 for i in 1:n]
end

function line_integration_tuple(n::Int)
    ξs, ws = line_gauss_rule(n)
    local_coordinates = Float64[]
    for ξ in ξs
        append!(local_coordinates, (ξ, 0.0, 0.0))
    end
    return local_coordinates, ws
end

function closest_node(nodes, x::Float64)
    best = nodes[1]
    best_distance = abs(nodes[1].x - x)
    for node in nodes
        distance = abs(node.x - x)
        if distance < best_distance
            best = node
            best_distance = distance
        end
    end
    return best
end

function solve_zero_dirichlet(K::Matrix{Float64}, F::Vector{Float64}, fixed::Vector{Int})
    fixed = sort(unique(fixed))
    free = setdiff(collect(eachindex(F)), fixed)
    U = zeros(length(F))
    U[free] = K[free, free] \ F[free]
    return U
end

function csv_cell(value)
    if value isa AbstractFloat
        return isfinite(value) ? string(value) : ""
    elseif value isa Integer
        return string(value)
    end
    text = string(value)
    if occursin('"', text)
        text = replace(text, "\"" => "\"\"")
    end
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * text * "\""
    end
    return text
end

function write_rows_csv(path::AbstractString, rows::Vector{<:NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(OUTPUT_COLUMNS), ","))
        for row in rows
            values = [csv_cell(getproperty(row, column)) for column in OUTPUT_COLUMNS]
            println(io, join(values, ","))
        end
    end
    return path
end

function build_beam_model!(mesh_level::Int)
    gmsh.model.add("beam1d-ss-$(mesh_level)")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(L0, 0.0, 0.0)
    Ω = gmsh.model.geo.addLine(p1, p2)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setTransfiniteCurve(Ω, mesh_level + 1)
    gmsh.model.addPhysicalGroup(1, [Ω], 1)
    gmsh.model.setPhysicalName(1, 1, "Ω")
    gmsh.model.mesh.generate(1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    return entities, nodes
end

function integration_strategy_name(gauss_bending::Int, gauss_shear::Int)
    return gauss_bending == gauss_shear ? "full" : "reduced"
end

function success_row(;
    gauss_bending::Int,
    gauss_shear::Int,
    gauss_L2::Int,
    h_over_L::Float64,
    mesh_level::Int,
    bc::Symbol,
    w_ref,
    w_fem,
    L2_w_pct,
    L2_φ_pct,
    dofs::Int,
    runtime_s::Float64,
    note::String,
)
    return (
        model = MODEL_NAME,
        element_type = "Seg2",
        integration_strategy = integration_strategy_name(gauss_bending, gauss_shear),
        gauss_bending = gauss_bending,
        gauss_shear = gauss_shear,
        gauss_L2 = gauss_L2,
        h_over_L = h_over_L,
        mesh_level = mesh_level,
        bc = String(bc),
        reference_kind = "beam_analytic",
        width_over_L = WIDTH_OVER_L0,
        E = E0,
        ν = ν0,
        κ = κ0,
        L = L0,
        q = q0,
        w_ref = w_ref,
        w_fem = w_fem,
        L2_w_pct = L2_w_pct,
        L2_φ_pct = L2_φ_pct,
        dofs = dofs,
        runtime_s = runtime_s,
        status = "ok",
        note = note,
    )
end

function failure_row(;
    gauss_bending::Int,
    gauss_shear::Int,
    gauss_L2::Int,
    h_over_L::Float64,
    mesh_level::Int,
    bc::Symbol,
    runtime_s::Float64,
    note::String,
)
    return (
        model = MODEL_NAME,
        element_type = "Seg2",
        integration_strategy = integration_strategy_name(gauss_bending, gauss_shear),
        gauss_bending = gauss_bending,
        gauss_shear = gauss_shear,
        gauss_L2 = gauss_L2,
        h_over_L = h_over_L,
        mesh_level = mesh_level,
        bc = String(bc),
        reference_kind = "beam_analytic",
        width_over_L = WIDTH_OVER_L0,
        E = E0,
        ν = ν0,
        κ = κ0,
        L = L0,
        q = q0,
        w_ref = "",
        w_fem = "",
        L2_w_pct = "",
        L2_φ_pct = "",
        dofs = 0,
        runtime_s = runtime_s,
        status = "failed",
        note = note,
    )
end

function run_case(gauss_bending::Int, gauss_shear::Int, gauss_L2::Int, h_over_L::Float64, mesh_level::Int, bc::Symbol)
    println("[beam1d] bc=$(bc) h/L=$(h_over_L) mesh=$(mesh_level) gb=$(gauss_bending) gs=$(gauss_shear) gL2=$(gauss_L2)")
    started = time_ns()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    try
        entities, nodes = build_beam_model!(mesh_level)
        s = section_data(h_over_L)
        load_order = max(gauss_bending, gauss_shear)

        elements_b = getElements(nodes, entities["Ω"], gauss_bending)
        prescribe!(elements_b, :EI => s.EI)
        set∇𝝭!(elements_b)

        elements_s = getElements(nodes, entities["Ω"], gauss_shear)
        prescribe!(elements_s, :kGA => s.kGA)
        set∇𝝭!(elements_s)

        elements_q = getElements(nodes, entities["Ω"], load_order)
        prescribe!(elements_q, :q => q0)
        set∇𝝭!(elements_q)

        ndofs = 2 * length(nodes)
        k = zeros(ndofs, ndofs)
        f = zeros(ndofs)

        𝑎 = [
            ∫κEIκds => elements_b,
            ∫γkGAγds => elements_s,
        ]
        𝑎(k)

        𝑓 = ∫vqds => elements_q
        𝑓(f)

        atol = max(1.0e-12, L0 * 1.0e-10)
        fixed = Int[]
        if bc == :SS
            for node in nodes
                if isapprox(node.x, 0.0; atol = atol) || isapprox(node.x, L0; atol = atol)
                    push!(fixed, 2 * node.𝐼 - 1)
                end
            end
        elseif bc == :CF
            for node in nodes
                if isapprox(node.x, 0.0; atol = atol)
                    push!(fixed, 2 * node.𝐼 - 1)
                    push!(fixed, 2 * node.𝐼)
                end
            end
        else
            error("unsupported boundary condition $bc")
        end

        U = solve_zero_dirichlet(k, f, fixed)
        push!(nodes, :d => U[1:2:end], :φ => U[2:2:end])

        elements_l2 = getElements(nodes, entities["Ω"], line_integration_tuple(gauss_L2))
        prescribe!(
            elements_l2,
            :u => (x, y, z) -> w_exact(bc, x, h_over_L),
            :φ => (x, y, z) -> φ_exact(bc, x, h_over_L),
        )
        set𝝭!(elements_l2)

        xobs = observation_x(bc)
        obs = closest_node(nodes, xobs)
        runtime_s = (time_ns() - started) / 1.0e9

        return success_row(
            gauss_bending = gauss_bending,
            gauss_shear = gauss_shear,
            gauss_L2 = gauss_L2,
            h_over_L = h_over_L,
            mesh_level = mesh_level,
            bc = bc,
            w_ref = w_exact(bc, xobs, h_over_L),
            w_fem = obs.d,
            L2_w_pct = 100.0 * L₂(elements_l2),
            L2_φ_pct = 100.0 * L₂φ(elements_l2),
            dofs = ndofs,
            runtime_s = runtime_s,
            note = "Test-style 1D Timoshenko sweep",
        )
    catch err
        runtime_s = (time_ns() - started) / 1.0e9
        return failure_row(
            gauss_bending = gauss_bending,
            gauss_shear = gauss_shear,
            gauss_L2 = gauss_L2,
            h_over_L = h_over_L,
            mesh_level = mesh_level,
            bc = bc,
            runtime_s = runtime_s,
            note = string("case execution failed: ", sprint(showerror, err)),
        )
    finally
        gmsh.finalize()
    end
end

benchmark_case() = (gauss_bending = 2, gauss_shear = 1, gauss_L2 = 5, h_over_L = 0.001, mesh_level = 100, bc = :SS)

function smoke_cases()
    return [
        (gauss_bending = 2, gauss_shear = 1, gauss_L2 = 5, h_over_L = 0.1, mesh_level = 40, bc = :SS),
    ]
end

function baseline_cases()
    cases = NamedTuple[]
    for bc in BC_SET, h_over_L in H_OVER_L_SET, mesh_level in MESH_LEVEL_SET,
        gauss_bending in GAUSS_SET, gauss_shear in GAUSS_SET, gauss_L2 in GAUSS_SET
        push!(
            cases,
            (
                gauss_bending = gauss_bending,
                gauss_shear = gauss_shear,
                gauss_L2 = gauss_L2,
                h_over_L = h_over_L,
                mesh_level = mesh_level,
                bc = bc,
            ),
        )
    end
    return cases
end

function run_cases(cases)
    rows = NamedTuple[]
    for case in cases
        push!(rows, run_case(case.gauss_bending, case.gauss_shear, case.gauss_L2, case.h_over_L, case.mesh_level, case.bc))
    end
    return rows
end

function main(args = ARGS)
    cases = if isempty(args) || args[1] == "smoke"
        smoke_cases()
    elseif args[1] == "baseline"
        baseline_cases()
    else
        error("usage: julia static_bending_1d_sweep.jl [smoke|baseline]")
    end

    rows = run_cases(cases)
    output_path = write_rows_csv(RAW_CSV_PATH, rows)
    println("wrote $(length(rows)) rows to $output_path ($MODEL_NAME)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
