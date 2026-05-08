using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫αwwdΓ

using LinearAlgebra
using Printf
using TimerOutputs
using XLSX
import Gmsh: gmsh

E = 30.0e6
ν = 0.3
b = 1.0
h_over_b = 0.01
h = h_over_b*b
Dᵇ = E*h^3/12/(1-ν^2)

# T2 表 9-12：四边简支矩形板纯剪切屈曲因子 k。
table_9_12 = [
    (1.0, 14.71),
    (1.5, 11.50),
    (2.0, 10.34),
    (2.5, 10.85),
]

w(x,y,z) = 0.0
σ₁₁(x,y,z) = 0.0
σ₂₂(x,y,z) = 0.0
σ₁₂(x,y,z) = 1.0

const to = TimerOutput()
const output_xlsx = normpath(joinpath(@__DIR__, "..", "results", "mindlin_shear_ssss_table_9_12.xlsx"))

function aspect_tag(r)
    return replace(@sprintf("%.2f", r), "."=>"p")
end

function mesh_file(r)
    return normpath(joinpath(@__DIR__, "..", "msh", "mindlin_shear_ssss_ab_$(aspect_tag(r)).msh"))
end

function excel_column_name(j)
    name = ""
    while j > 0
        j, r = divrem(j - 1, 26)
        name = string(Char('A' + r), name)
    end
    return name
end

excel_cell(row, col) = string(excel_column_name(col), row)

function write_xlsx(filepath, rows)
    mkpath(dirname(filepath))
    headers = ["a/b", "mesh", "lambda_cr", "k_num", "k_book", "k_error", "tau_cr_num", "tau_cr_book", "nodes"]
    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "T2 9-12 shear")
        for (j, header) in enumerate(headers)
            sheet[excel_cell(1, j)] = header
        end
        for (i, row) in enumerate(rows)
            values = [row.r, row.mesh, row.λcr, row.k_num, row.k_book, row.k_err, row.τcr_num, row.τcr_book, row.nnodes]
            for (j, value) in enumerate(values)
                sheet[excel_cell(i + 1, j)] = value
            end
        end
    end
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

function solve_case(r)
    gmsh.clear()
    @timeit to "open msh file" gmsh.open(mesh_file(r))
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ,nʷ)
    kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kᵠʷ = zeros(2*nᵠ,nʷ)
    kgʷʷ = zeros(nʷ,nʷ)

    @timeit to "calculate stiffness and geometric stiffness" begin
        elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
        set∇𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ=>elements
        𝑎ᵠʷ = ∫φwdΩ=>elements
        𝑎ᵠᵠ = [∫φφdΩ=>elements, ∫κκdΩ=>elements]
        𝑎ʷʷ(kʷʷ)
        𝑎ᵠʷ(kᵠʷ)
        𝑎ᵠᵠ(kᵠᵠ)

        prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        𝑎ᴳʷʷ = ∫wwGdΩ2D=>elements
        𝑎ᴳʷʷ(kgʷʷ)
    end

    @timeit to "calculate ∫αwwdΓ" begin
        elements_1 = getElements(nodes, entities["Γ¹"])
        elements_2 = getElements(nodes, entities["Γ²"])
        elements_3 = getElements(nodes, entities["Γ³"])
        elements_4 = getElements(nodes, entities["Γ⁴"])
        prescribe!(elements_1, :α=>1e8*E, :g=>w)
        prescribe!(elements_2, :α=>1e8*E, :g=>w)
        prescribe!(elements_3, :α=>1e8*E, :g=>w)
        prescribe!(elements_4, :α=>1e8*E, :g=>w)
        set𝝭!(elements_1)
        set𝝭!(elements_2)
        set𝝭!(elements_3)
        set𝝭!(elements_4)
        𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
        fᵅ = zeros(nʷ)
        𝑎ʷ(kʷʷ,fᵅ)
    end

    @timeit to "solve buckling eigenvalue" begin
        Keff = kʷʷ - kᵠʷ'*(kᵠᵠ\kᵠʷ)
        λ_positive = positive_finite_eigenvalues(Keff, kgʷʷ)
        isempty(λ_positive) && error("no positive finite buckling eigenvalue found for a/b = $r")
        λcr = first(λ_positive)
        k_num = λcr*b^2/(π^2*Dᵇ)
        τcr_num = λcr/h
    end

    return λcr, k_num, τcr_num, nʷ
end

gmsh.initialize()
try
    println("T2 table 9-12: SSSS rectangular plate under pure in-plane shear")
    results = []
    @printf("%7s %18s %12s %14s %10s %12s %14s %14s %8s\n",
        "a/b", "msh", "λcr", "k_num", "k_book", "k_err", "τcr_num", "τcr_book", "nodes")
    println("-"^115)
    for (r, k_book) in table_9_12
        λcr, k_num, τcr_num, nnodes = solve_case(r)
        k_err = abs(k_num-k_book)/abs(k_book)
        τcr_book = k_book*π^2*Dᵇ/(b^2*h)
        @printf("%7.2f %18s %12.6e %14.10f %10.2f %12.4e %14.6f %14.6f %8d\n",
            r, basename(mesh_file(r)), λcr, k_num, k_book, k_err, τcr_num, τcr_book, nnodes)
        push!(results, (
            r = r, mesh = basename(mesh_file(r)), λcr = λcr, k_num = k_num,
            k_book = k_book, k_err = k_err, τcr_num = τcr_num,
            τcr_book = τcr_book, nnodes = nnodes,
        ))
    end
    write_xlsx(output_xlsx, results)
    println("Excel output: ", output_xlsx)
finally
    gmsh.finalize()
end

println(to)
