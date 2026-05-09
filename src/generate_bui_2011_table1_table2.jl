include("mindlin_square_buckling.jl")

using Printf
using XLSX

const BUI_H_OVER_B = 0.01
const BUI_E = 200.0e9
const BUI_NU = 0.3
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

function bui_boundary_to_internal(boundary::String)
    length(boundary) == 4 || error("Bui boundary must have four letters.")
    b = collect(boundary)
    # Bui: gamma1=bottom, gamma2=right, gamma3=top, gamma4=left.
    # Solver: left, right, top, bottom.
    return String([b[4], b[2], b[3], b[1]])
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
            case = BucklingCase(
                "Bui 2011 Table 1 $(boundary) $(nodes_per_side)x$(nodes_per_side)",
                :uniaxial,
                internal_boundary,
                1.0,
                0.0,
                nodes_per_side - 1,
            )
            result = solve_buckling(case, BUI_H_OVER_B; E=BUI_E, ν=BUI_NU)
            values_by_boundary[boundary][nodes_per_side] = result.k_num
            push!(row_values, result.k_num)
            @printf("%-4s %2dx%-2d  internal=%-4s  k=%12.6f  nodes=%d\n",
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
