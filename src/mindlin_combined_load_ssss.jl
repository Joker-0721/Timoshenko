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

# T2 表 9-13：正方形四边简支板，压缩 σ 与剪切 τ 同时作用。
table_9_13 = [
    (0.0, 14.71),
    (0.5,  7.09),
    (1.0,  4.50),
    (1.5,  3.24),
    (2.0,  2.51),
]

w(x,y,z) = 0.0

const to = TimerOutput()
const output_xlsx = normpath(joinpath(@__DIR__, "..", "results", "mindlin_combined_load_ssss_table_9_13.xlsx"))

function mesh_file()
    return normpath(joinpath(@__DIR__, "..", "msh", "mindlin_combined_load_ssss_ab_1p00.msh"))
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
    headers = [
        "sigma_over_tau", "mesh", "lambda_tau_cr", "k_num", "k_book", "k_error",
        "tau_cr_num", "tau_cr_book", "sigma_cr_num", "sigma_cr_book", "nodes",
    ]
    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "T2 9-13 combined")
        for (j, header) in enumerate(headers)
            sheet[excel_cell(1, j)] = header
        end
        for (i, row) in enumerate(rows)
            values = [
                row.ratio, row.mesh, row.λcr, row.k_num, row.k_book, row.k_err,
                row.τcr_num, row.τcr_book, row.σcr_num, row.σcr_book, row.nnodes,
            ]
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

function solve_case(ratio)
    σ₁₁(x,y,z) = ratio
    σ₂₂(x,y,z) = 0.0
    σ₁₂(x,y,z) = 1.0

    gmsh.clear()
    @timeit to "open msh file" gmsh.open(mesh_file())
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
        isempty(λ_positive) && error("no positive finite buckling eigenvalue found for σ/τ = $ratio")
        λcr = first(λ_positive)
        k_num = λcr*b^2/(π^2*Dᵇ)
        τcr_num = λcr/h
        σcr_num = ratio*τcr_num
    end

    return λcr, k_num, τcr_num, σcr_num, nʷ
end

gmsh.initialize()
try
    println("T2 table 9-13: SSSS square plate under combined compression and shear")
    results = []
    @printf("%8s %18s %12s %14s %10s %12s %14s %14s %14s %14s %8s\n",
        "σ/τ", "msh", "λτcr", "k_num", "k_book", "k_err",
        "τcr_num", "τcr_book", "σcr_num", "σcr_book", "nodes")
    println("-"^138)
    for (ratio, k_book) in table_9_13
        λcr, k_num, τcr_num, σcr_num, nnodes = solve_case(ratio)
        k_err = abs(k_num-k_book)/abs(k_book)
        τcr_book = k_book*π^2*Dᵇ/(b^2*h)
        σcr_book = ratio*τcr_book
        @printf("%8.3f %18s %12.6e %14.10f %10.2f %12.4e %14.6f %14.6f %14.6f %14.6f %8d\n",
            ratio, basename(mesh_file()), λcr, k_num, k_book, k_err,
            τcr_num, τcr_book, σcr_num, σcr_book, nnodes)
        push!(results, (
            ratio = ratio, mesh = basename(mesh_file()), λcr = λcr, k_num = k_num,
            k_book = k_book, k_err = k_err, τcr_num = τcr_num,
            τcr_book = τcr_book, σcr_num = σcr_num, σcr_book = σcr_book,
            nnodes = nnodes,
        ))
    end
    write_xlsx(output_xlsx, results)
    println("Excel output: ", output_xlsx)
finally
    gmsh.finalize()
end

println(to)
