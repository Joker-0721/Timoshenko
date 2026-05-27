const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    push!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ

using LinearAlgebra
using Printf
using TimerOutputs
using XLSX
import Gmsh: gmsh

const BUI_E = 200.0e9
const BUI_ν = 0.3
const BUI_B = 1.0
const BUI_H_OVER_B = 0.01
const BUI_H = BUI_H_OVER_B*BUI_B
const BUI_D = BUI_E*BUI_H^3/(12.0*(1.0 - BUI_ν^2))
const BUI_OUTPUT_XLSX = normpath(joinpath(@__DIR__, "..", "results", "bui_2011_paper_mindlin.xlsx"))
const BUI_MSH_DIR = normpath(joinpath(@__DIR__, "..", "msh"))
const BUI_NODE_PATTERNS = (7, 9, 11, 13, 15, 17)
const BUI_BOUNDARIES = ("CCCC", "CSCS", "SCSC", "SSSS", "FSCS", "FSSS")
const EIGEN_IMAG_TOL = 1.0e-7

const to = TimerOutput()

const BUI_EXACT = Dict(
    "CCCC" => 10.070,
    "CSCS" => 7.690,
    "SCSC" => 6.750,
    "SSSS" => 4.000,
    "FSCS" => 1.700,
    "FSSS" => 1.440,
)

const BUI_TABLE2_REFERENCES = Dict(
    "SSSS" => (4.000, 4.011, 4.041, 3.999, 4.000, 4.00, 4.017, 3.9977, 4.011),
    "CCCC" => (10.070, 10.392, 10.387, 10.142, 10.109, 10.08, 10.308, 10.052, 10.310),
    "CSCS" => (7.690, 7.796, 7.757, 7.683, 7.712, 7.70, missing, missing, 7.681),
    "SCSC" => (6.750, 6.882, 6.972, 6.781, missing, missing, missing, missing, missing),
    "FSCS" => (1.700, 1.718, 1.724, 1.712, missing, missing, missing, missing, missing),
    "FSSS" => (1.440, 1.422, 1.417, 1.428, missing, missing, missing, missing, missing),
)

struct BuiCase
    boundary::String
    nodes_per_side::Int
end

mesh_file(nodes_per_side::Int) = joinpath(BUI_MSH_DIR, "bui_2011_square_$(nodes_per_side)x$(nodes_per_side).msh")

function positive_finite_eigenvalues(k, kg)
    λ = eigvals(k, kg)
    values = Float64[]
    for λᵢ in λ
        if isfinite(real(λᵢ)) && isfinite(imag(λᵢ)) && abs(imag(λᵢ)) < EIGEN_IMAG_TOL && real(λᵢ) > 0.0
            push!(values, real(λᵢ))
        end
    end
    return sort!(values)
end

function constrained_nodes(nodes, entities, boundary::String)
    length(boundary) == 4 || error("Bui boundary must have four letters.")
    side_names = ("Γ¹", "Γ²", "Γ³", "Γ⁴")
    fixed_w = Set{Int}()
    fixed_φ = Set{Int}()

    # Boundary order follows the generated Bui square meshes:
    # Γ¹ = bottom, Γ² = right, Γ³ = top, Γ⁴ = left.
    for (bc, side) in zip(collect(boundary), side_names)
        bc == 'F' && continue
        elements = getElements(nodes, entities[side])
        side_nodes = (node.𝐼 for element in elements for node in element.𝓒)
        if bc == 'S'
            union!(fixed_w, side_nodes)
        elseif bc == 'C'
            union!(fixed_w, side_nodes)
            union!(fixed_φ, side_nodes)
        else
            error("unsupported boundary code '$bc' in $boundary")
        end
    end

    return sort!(collect(fixed_w)), sort!(collect(fixed_φ))
end

function assemble_mindlin_matrices(nodes, entities)
    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
    kᵠʷ = zeros(2*nᵠ, nʷ)
    kgʷʷ = zeros(nʷ, nʷ)

    @timeit to "calculate domain operators" begin
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :E=>BUI_E, :ν=>BUI_ν, :h=>BUI_H)
        @timeit to "calculate shape functions" set∇𝝭!(elements)

        𝑎ʷʷ = ∫wwdΩ=>elements
        𝑎ᵠʷ = ∫φwdΩ=>elements
        𝑎ᵠᵠ = [
            ∫φφdΩ=>elements,
            ∫κκdΩ=>elements,
        ]

        @timeit to "assemble Kww" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble Kφw" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble Kφφ" 𝑎ᵠᵠ(kᵠᵠ)

        N₁₁(x, y, z) = 1.0
        N₂₂(x, y, z) = 0.0
        N₁₂(x, y, z) = 0.0
        prescribe!(elements, :σ₁₁=>N₁₁, :σ₂₂=>N₂₂, :σ₁₂=>N₁₂)

        𝑎ᴳʷʷ = ∫wwGdΩ2D=>elements
        @timeit to "assemble KGww" 𝑎ᴳʷʷ(kgʷʷ)
    end

    return kᵠᵠ, kᵠʷ, kʷʷ, kgʷʷ
end

function solve_bui_case(case::BuiCase)
    filepath = mesh_file(case.nodes_per_side)
    isfile(filepath) || error("mesh file not found: $filepath. Run Timoshenko/src/generatePatchtest.jl first.")

    gmsh.clear()
    @timeit to "open msh file" gmsh.open(filepath)
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    kᵠᵠ, kᵠʷ, kʷʷ, kgʷʷ = assemble_mindlin_matrices(nodes, entities)
    nᵠ = length(nodes)

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    KG = [
        zeros(2*nᵠ, 2*nᵠ) zeros(2*nᵠ, length(nodes))
        zeros(length(nodes), 2*nᵠ) kgʷʷ
    ]

    fixed_w_nodes, fixed_φ_nodes = constrained_nodes(nodes, entities, case.boundary)
    fixed_φ_dofs = reduce(vcat, ([2*i - 1, 2*i] for i in fixed_φ_nodes), init=Int[])
    fixed_w_dofs = 2*nᵠ .+ fixed_w_nodes
    fixed_dofs = sort!(unique!(vcat(fixed_φ_dofs, fixed_w_dofs)))
    free_dofs = setdiff(1:size(K, 1), fixed_dofs)
    free_φ = [dof for dof in free_dofs if dof <= 2*nᵠ]
    free_w = [dof for dof in free_dofs if dof > 2*nᵠ]

    Kpp = K[free_φ, free_φ]
    Kpw = K[free_φ, free_w]
    Kwp = K[free_w, free_φ]
    Kww = K[free_w, free_w]
    KGww = KG[free_w, free_w]
    Keff = isempty(free_φ) ? Kww : Kww - Kwp*(Kpp\Kpw)

    λ = positive_finite_eigenvalues(Keff, KGww)
    isempty(λ) && error("no positive finite buckling eigenvalue found for $(case.boundary), $(case.nodes_per_side)x$(case.nodes_per_side)")

    λcr = first(λ)
    k_num = λcr*BUI_B^2/(π^2*BUI_D)
    σcr = λcr/BUI_H

    return (;
        λcr,
        k_num,
        σcr,
        nnodes = length(nodes),
        nfixed = length(fixed_dofs),
        mesh = basename(filepath),
    )
end

function excel_column_name(j::Int)
    name = ""
    while j > 0
        j, r = divrem(j - 1, 26)
        name = string(Char('A' + r), name)
    end
    return name
end

excel_cell(row::Int, col::Int) = string(excel_column_name(col), row)

function write_sheet!(sheet, headers, rows)
    for (j, header) in enumerate(headers)
        sheet[excel_cell(1, j)] = header
    end
    for (i, row) in enumerate(rows)
        for (j, value) in enumerate(row)
            sheet[excel_cell(i + 1, j)] = ismissing(value) ? nothing : value
        end
    end
end

function table1_rows()
    rows = Vector{Any}[]
    values_by_boundary = Dict{String, Dict{Int, Float64}}()

    for boundary in BUI_BOUNDARIES
        values_by_boundary[boundary] = Dict{Int, Float64}()
        row_values = Any[boundary]

        for nodes_per_side in BUI_NODE_PATTERNS
            case = BuiCase(boundary, nodes_per_side)
            result = solve_bui_case(case)
            values_by_boundary[boundary][nodes_per_side] = result.k_num
            push!(row_values, result.k_num)
            @printf("%-4s %2dx%-2d mesh=%s k=%12.6f λ=%12.6e nodes=%d fixed=%d\n",
                boundary,
                nodes_per_side,
                nodes_per_side,
                result.mesh,
                result.k_num,
                result.λcr,
                result.nnodes,
                result.nfixed)
        end

        push!(row_values, BUI_EXACT[boundary])
        push!(rows, row_values)
    end

    return rows, values_by_boundary
end

function table2_rows(values_by_boundary)
    rows = Vector{Any}[]
    for boundary in BUI_BOUNDARIES
        refs = BUI_TABLE2_REFERENCES[boundary]
        present = values_by_boundary[boundary][13]
        push!(rows, Any[boundary, refs..., present])
    end
    return rows
end

function write_workbook(table1, table2)
    mkpath(dirname(BUI_OUTPUT_XLSX))

    table1_headers = Any["Boundary"]
    append!(table1_headers, ["$(n)x$(n)" for n in BUI_NODE_PATTERNS])
    push!(table1_headers, "Exact [1]")

    table2_headers = [
        "Boundary", "Exact [1]", "FEM [45]", "BEM [45]", "DRM [45]",
        "SFSM [34]", "SFSM [33]", "RPIM [70]", "DQEM [42]", "DSC [58]", "Present",
    ]

    XLSX.openxlsx(BUI_OUTPUT_XLSX, mode="w") do xf
        sheet1 = xf[1]
        XLSX.rename!(sheet1, "Table 1")
        write_sheet!(sheet1, table1_headers, table1)

        sheet2 = XLSX.addsheet!(xf, "Table 2")
        write_sheet!(sheet2, table2_headers, table2)
    end
end

function main()
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        table1, values_by_boundary = table1_rows()
        table2 = table2_rows(values_by_boundary)
        write_workbook(table1, table2)
        println("Excel output: ", BUI_OUTPUT_XLSX)
        println(to)
    finally
        gmsh.finalize()
    end
end

main()
