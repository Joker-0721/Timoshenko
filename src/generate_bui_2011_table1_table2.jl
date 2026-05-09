const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    push!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ

using LinearAlgebra
using Printf
using XLSX
import Gmsh: gmsh

const BUI_H_OVER_B = 0.01
const BUI_E = 200.0e9
const BUI_NU = 0.3
const BUI_B = 1.0
const BUI_OUTPUT_XLSX = normpath(joinpath(@__DIR__, "..", "results", "bui_2011_table1_table2.xlsx"))
const BUI_NODE_PATTERNS = (7, 9, 11, 13, 15, 17)

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

const BUI_BOUNDARIES = ("CCCC", "CSCS", "SCSC", "SSSS", "FSCS", "FSSS")

struct BuiCase
    boundary::String
    internal_boundary::String
    mesh_n::Int
end

function bui_boundary_to_internal(boundary::String)
    length(boundary) == 4 || error("Bui boundary must have four letters.")
    b = collect(boundary)
    # Bui: Γ1=bottom, Γ2=right, Γ3=top, Γ4=left.
    # Internal order: left, right, top, bottom.
    return String([b[4], b[2], b[3], b[1]])
end

function bending_rigidity(E, ν, h)
    return E*h^3/(12.0*(1.0 - ν^2))
end

function generate_square_mesh!(nodes_per_side::Int)
    gmsh.clear()
    gmsh.model.add("bui_2011_square")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(BUI_B, 0.0, 0.0)
    p3 = gmsh.model.geo.addPoint(BUI_B, BUI_B, 0.0)
    p4 = gmsh.model.geo.addPoint(0.0, BUI_B, 0.0)

    bottom = gmsh.model.geo.addLine(p1, p2)
    right = gmsh.model.geo.addLine(p2, p3)
    top = gmsh.model.geo.addLine(p3, p4)
    left = gmsh.model.geo.addLine(p4, p1)
    loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(right, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteCurve(left, nodes_per_side)
    gmsh.model.geo.mesh.setTransfiniteSurface(surface)
    gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()

    for (name, entity) in (("left", left), ("right", right), ("top", top), ("bottom", bottom))
        group = gmsh.model.addPhysicalGroup(1, [entity])
        gmsh.model.setPhysicalName(1, group, name)
    end
    domain = gmsh.model.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, domain, "domain")

    gmsh.model.mesh.generate(2)
end

function constrained_nodes(nodes, entities, internal_boundary::String)
    side_names = ("left", "right", "top", "bottom")
    fixed_w = Set{Int}()
    fixed_phi = Set{Int}()

    for (bc, side) in zip(collect(internal_boundary), side_names)
        bc == 'F' && continue
        elements = getElements(nodes, entities[side])
        side_nodes = collect(node.𝐼 for element in elements for node in element.𝓒)
        if bc == 'S'
            union!(fixed_w, side_nodes)
        elseif bc == 'C'
            union!(fixed_w, side_nodes)
            union!(fixed_phi, side_nodes)
        else
            error("unsupported boundary code '$bc' in $internal_boundary")
        end
    end

    return sort!(collect(fixed_w)), sort!(collect(fixed_phi))
end

function positive_finite_eigenvalues(k, kg)
    λ = eigvals(k, kg)
    values = Float64[]
    for λᵢ in λ
        if isfinite(real(λᵢ)) && isfinite(imag(λᵢ)) && abs(imag(λᵢ)) < 1.0e-7 && real(λᵢ) > 0.0
            push!(values, real(λᵢ))
        end
    end
    return sort!(values)
end

function solve_bui_case(case::BuiCase)
    generate_square_mesh!(case.mesh_n + 1)

    h = BUI_H_OVER_B*BUI_B
    D = bending_rigidity(BUI_E, BUI_NU, h)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
    kᵠʷ = zeros(2*nᵠ, nʷ)
    kgʷʷ = zeros(nʷ, nʷ)

    elements = getElements(nodes, entities["domain"])
    prescribe!(elements, :E=>BUI_E, :ν=>BUI_NU, :h=>h)
    set∇𝝭!(elements)
    (∫wwdΩ => elements)(kʷʷ)
    (∫φwdΩ => elements)(kᵠʷ)
    ([∫φφdΩ => elements, ∫κκdΩ => elements])(kᵠᵠ)

    N₁₁(x,y,z) = 1.0
    N₂₂(x,y,z) = 0.0
    N₁₂(x,y,z) = 0.0
    prescribe!(elements, :σ₁₁=>N₁₁, :σ₂₂=>N₂₂, :σ₁₂=>N₁₂)
    (∫wwGdΩ2D => elements)(kgʷʷ)

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    KG = [
        zeros(2*nᵠ, 2*nᵠ) zeros(2*nᵠ, nʷ)
        zeros(nʷ, 2*nᵠ) kgʷʷ
    ]

    fixed_w_nodes, fixed_phi_nodes = constrained_nodes(nodes, entities, case.internal_boundary)
    fixed_phi_dofs = reduce(vcat, ([2*i - 1, 2*i] for i in fixed_phi_nodes), init=Int[])
    fixed_w_dofs = 2*nᵠ .+ fixed_w_nodes
    fixed_dofs = sort!(unique!(vcat(fixed_phi_dofs, fixed_w_dofs)))
    free_dofs = setdiff(1:size(K, 1), fixed_dofs)
    free_phi = [dof for dof in free_dofs if dof <= 2*nᵠ]
    free_w = [dof for dof in free_dofs if dof > 2*nᵠ]

    Kpp = K[free_phi, free_phi]
    Kpw = K[free_phi, free_w]
    Kwp = K[free_w, free_phi]
    Kww = K[free_w, free_w]
    KGww = KG[free_w, free_w]
    Keff = isempty(free_phi) ? Kww : Kww - Kwp*(Kpp\Kpw)

    λ = positive_finite_eigenvalues(Keff, KGww)
    isempty(λ) && error("no positive finite buckling eigenvalue found for $(case.boundary)")

    λcr = first(λ)
    k_num = λcr*BUI_B^2/(π^2*D)
    return (; λcr, k_num, nnodes=length(nodes), nfixed=length(fixed_dofs))
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

function table1_rows()
    rows = Vector{Any}[]
    values_by_boundary = Dict{String, Dict{Int, Float64}}()

    for boundary in BUI_BOUNDARIES
        internal_boundary = bui_boundary_to_internal(boundary)
        row_values = Any[boundary]
        values_by_boundary[boundary] = Dict{Int, Float64}()
        for nodes_per_side in BUI_NODE_PATTERNS
            case = BuiCase(boundary, internal_boundary, nodes_per_side - 1)
            result = solve_bui_case(case)
            values_by_boundary[boundary][nodes_per_side] = result.k_num
            push!(row_values, result.k_num)
            @printf("%-4s %2dx%-2d internal=%-4s k=%12.6f nodes=%d\n",
                boundary, nodes_per_side, nodes_per_side, internal_boundary, result.k_num, result.nnodes)
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
    finally
        gmsh.finalize()
    end
end

main()
