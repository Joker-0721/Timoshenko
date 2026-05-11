using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫ψxψxGdΩ2D, ∫ψyψyGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫αwwdΓ

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

# T1 表 9-1：均勻受壓簡支矩形板，E = 30e6 psi, ν = 0.3, h/b = 0.01。
# 後兩欄是書上列印值，計算時另用解析公式得到更多位數。
table_9_1 = [
    (0.20, 27.00, 73200.0),
    (0.30, 13.20, 35800.0),
    (0.40,  8.41, 22800.0),
    (0.50,  6.25, 16900.0),
    (0.60,  5.14, 13900.0),
    (0.70,  4.53, 12300.0),
    (0.80,  4.20, 11400.0),
    (0.90,  4.04, 11000.0),
    (1.00,  4.00, 10800.0),
    (1.10,  4.04, 11000.0),
    (1.20,  4.13, 11200.0),
    (1.30,  4.28, 11600.0),
    (1.40,  4.47, 12100.0),
    (1.41,  4.49, 12200.0),
]

function exact_uniaxial_k(r; mmax=80)
    k_exact = Inf
    m_exact = 0
    for m in 1:mmax
        k = (m/r + r/m)^2
        if k < k_exact
            k_exact = k
            m_exact = m
        end
    end
    return k_exact, m_exact
end

exact_λcr(k_exact) = k_exact*π^2*Dᵇ/b^2
exact_σcr(k_exact) = exact_λcr(k_exact)/h

w(x,y,z) = 0.0
σ₁₁(x,y,z) = 1.0
σ₂₂(x,y,z) = 0.0
σ₁₂(x,y,z) = 0.0

const to = TimerOutput()
const output_xlsx = normpath(joinpath(@__DIR__, "..", "results", "mindlin_uniaxial_ssss_table_9_1.xlsx"))

function aspect_tag(r)
    return replace(@sprintf("%.2f", r), "."=>"p")
end

function mesh_file(r)
    return normpath(joinpath(@__DIR__, "..", "msh", "mindlin_uniaxial_ssss_ab_$(aspect_tag(r)).msh"))
end

function excel_column_name(j)
    name = ""
    while j > 0
        j, r = divrem(j - 1, 26)
        name = string(Char('A' + r), name)
    end
    return name
end

function excel_cell(row, col)
    return string(excel_column_name(col), row)
end

function write_table_9_1_xlsx(filepath, rows)
    mkpath(dirname(filepath))
    headers = [
        "a/b",
        "m",
        "mesh",
        "lambda_cr",
        "k_num",
        "k_exact",
        "k_book",
        "k_error",
        "sigma_cr_num",
        "sigma_cr_exact",
        "sigma_cr_book",
        "nodes",
    ]

    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Table 9-1")
        for (j, header) in enumerate(headers)
            sheet[excel_cell(1, j)] = header
        end
        for (i, row) in enumerate(rows)
            values = [
                row.r,
                row.m_exact,
                row.mesh,
                row.λcr,
                row.k_num,
                row.k_exact,
                row.k_book,
                row.k_err,
                row.σcr_num,
                row.σcr_exact,
                row.σcr_book,
                row.nnodes,
            ]
            for (j, value) in enumerate(values)
                sheet[excel_cell(i + 1, j)] = value
            end
        end
    end
end

function solve_table_9_1_case(r)
    gmsh.clear()
    @timeit to "open msh file" gmsh.open(mesh_file(r))
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ,nʷ)
    kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kᵠʷ = zeros(2*nᵠ,nʷ)
    kgᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kgʷʷ = zeros(nʷ,nʷ)

    @timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ, ∫wwGdΩ2D, ∫ψψGdΩ2D" begin
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
        @timeit to "calculate shape functions" set∇𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ=>elements
        𝑎ᵠʷ = ∫φwdΩ=>elements
        𝑎ᵠᵠ = [
            ∫φφdΩ=>elements,
            ∫κκdΩ=>elements,
        ]
        @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)

        prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        𝑎ᴳʷʷ = ∫wwGdΩ2D=>elements
        𝑎ᴳᵠˣᵠˣ = ∫ψxψxGdΩ2D=>elements
        𝑎ᴳᵠʸᵠʸ = ∫ψyψyGdΩ2D=>elements
        @timeit to "assemble" 𝑎ᴳʷʷ(kgʷʷ)
        @timeit to "assemble" 𝑎ᴳᵠˣᵠˣ(view(kgᵠᵠ, 1:2:2*nᵠ, 1:2:2*nᵠ))
        @timeit to "assemble" 𝑎ᴳᵠʸᵠʸ(view(kgᵠᵠ, 2:2:2*nᵠ, 2:2:2*nᵠ))
    end

    @timeit to "calculate ∫αwwdΓ" begin
        @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
        @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
        @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"])
        @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"])
        prescribe!(elements_1, :α=>1e8*E, :g=>w)
        prescribe!(elements_2, :α=>1e8*E, :g=>w)
        prescribe!(elements_3, :α=>1e8*E, :g=>w)
        prescribe!(elements_4, :α=>1e8*E, :g=>w)
        @timeit to "calculate shape functions" set𝝭!(elements_1)
        @timeit to "calculate shape functions" set𝝭!(elements_2)
        @timeit to "calculate shape functions" set𝝭!(elements_3)
        @timeit to "calculate shape functions" set𝝭!(elements_4)
        𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
        fᵅ = zeros(nʷ)
        @timeit to "assemble" 𝑎ʷ(kʷʷ,fᵅ)
    end

    @timeit to "solve buckling eigenvalue" begin
        K = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
        KG = [
            kgᵠᵠ zeros(2*nᵠ,nʷ)
            zeros(nʷ,2*nᵠ) kgʷʷ
        ]
        λ = eigvals(K, KG)
        λ_positive = sort!(collect(real(λᵢ) for λᵢ in λ if isfinite(real(λᵢ)) && isfinite(imag(λᵢ)) && abs(imag(λᵢ)) < 1.0e-7 && real(λᵢ) > 0.0))
        isempty(λ_positive) && error("no positive finite buckling eigenvalue found for a/b = $r")
        λcr = first(λ_positive)
        k_num = λcr*b^2/(π^2*Dᵇ)
        σcr_num = λcr/h
    end

    return λcr, k_num, σcr_num, nʷ
end

gmsh.initialize()
try
    println("T1 表 9-1: SSSS rectangular plate under uniaxial in-plane compression")
    println("E = 30e6 psi, ν = 0.3, h/b = 0.01")
    println("Required mesh names: msh/mindlin_uniaxial_ssss_ab_<a_over_b>.msh, e.g. ab_1p00.")
    results = []
    @printf("%7s %3s %18s %12s %14s %14s %10s %12s %14s %14s %10s %8s\n",
        "a/b", "m", "msh", "λcr", "k_num", "k_exact",
        "k_book", "k_err", "σcr_num", "σcr_exact", "σ_book", "nodes")
    println("-"^153)

    for (r, k_book, σcr_book) in table_9_1
        k_exact, m_exact = exact_uniaxial_k(r)
        σcr_exact = exact_σcr(k_exact)
        λcr, k_num, σcr_num, nnodes = solve_table_9_1_case(r)
        k_err = abs(k_num-k_exact)/abs(k_exact)
        @printf("%7.2f %3d %18s %12.6e %14.10f %14.10f %10.2f %12.4e %14.6f %14.6f %10.1f %8d\n",
            r, m_exact, basename(mesh_file(r)), λcr, k_num, k_exact,
            k_book, k_err, σcr_num, σcr_exact, σcr_book, nnodes)
        push!(results, (
            r = r,
            m_exact = m_exact,
            mesh = basename(mesh_file(r)),
            λcr = λcr,
            k_num = k_num,
            k_exact = k_exact,
            k_book = k_book,
            k_err = k_err,
            σcr_num = σcr_num,
            σcr_exact = σcr_exact,
            σcr_book = σcr_book,
            nnodes = nnodes,
        ))
    end

    write_table_9_1_xlsx(output_xlsx, results)
    println("Excel output: ", output_xlsx)
finally
    gmsh.finalize()
end

println(to)
