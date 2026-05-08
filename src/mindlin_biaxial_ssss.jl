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

# T1 9.3：SSSS 矩形板双轴均匀受压。γ = Nᵧ/Nₓ。
table_biaxial = [
    (1.0, 0.0),
    (1.0, 0.25),
    (1.0, 0.50),
    (1.0, 1.0),
    (1.0, 2.0),
    (1.0, 4.0),
]

w(x,y,z) = 0.0

const to = TimerOutput()
const output_xlsx = normpath(joinpath(@__DIR__, "..", "results", "mindlin_biaxial_ssss.xlsx"))

function aspect_tag(r)
    return replace(@sprintf("%.2f", r), "."=>"p")
end

function mesh_file(r)
    return normpath(joinpath(@__DIR__, "..", "msh", "mindlin_biaxial_ssss_ab_$(aspect_tag(r)).msh"))
end

function exact_biaxial_k(r, γ; mmax=80, nmax=80)
    k_exact = Inf
    m_exact = 0
    n_exact = 0
    for m in 1:mmax, n in 1:nmax
        m² = (m/r)^2
        n² = n^2
        denominator = m² + γ*n²
        denominator <= 0.0 && continue
        k = (m² + n²)^2/denominator
        if k < k_exact
            k_exact = k
            m_exact = m
            n_exact = n
        end
    end
    return k_exact, m_exact, n_exact
end

exact_λcr(k_exact) = k_exact*π^2*Dᵇ/b^2
exact_σcr(k_exact) = exact_λcr(k_exact)/h

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
        "a/b", "gamma_Ny_over_Nx", "m", "n", "mesh", "lambda_cr",
        "k_num", "k_exact", "k_error", "sigma_cr_num", "sigma_cr_exact", "nodes",
    ]
    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "T1 9.3 biaxial")
        for (j, header) in enumerate(headers)
            sheet[excel_cell(1, j)] = header
        end
        for (i, row) in enumerate(rows)
            values = [
                row.r, row.γ, row.m_exact, row.n_exact, row.mesh, row.λcr,
                row.k_num, row.k_exact, row.k_err, row.σcr_num, row.σcr_exact, row.nnodes,
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

function solve_case(r, γ)
    σ₁₁(x,y,z) = 1.0
    σ₂₂(x,y,z) = γ
    σ₁₂(x,y,z) = 0.0

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
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
        @timeit to "calculate shape functions" set∇𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ=>elements
        𝑎ᵠʷ = ∫φwdΩ=>elements
        𝑎ᵠᵠ = [∫φφdΩ=>elements, ∫κκdΩ=>elements]
        @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)

        prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        𝑎ᴳʷʷ = ∫wwGdΩ2D=>elements
        @timeit to "assemble" 𝑎ᴳʷʷ(kgʷʷ)
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
        isempty(λ_positive) && error("no positive finite buckling eigenvalue found for a/b = $r, γ = $γ")
        λcr = first(λ_positive)
        k_num = λcr*b^2/(π^2*Dᵇ)
        σcr_num = λcr/h
    end

    return λcr, k_num, σcr_num, nʷ
end

gmsh.initialize()
try
    println("T1 9.3: SSSS rectangular plate under biaxial in-plane compression")
    println("E = 30e6 psi, ν = 0.3, h/b = 0.01")
    results = []
    @printf("%7s %8s %3s %3s %18s %12s %14s %14s %12s %14s %14s %8s\n",
        "a/b", "γ", "m", "n", "msh", "λcr", "k_num", "k_exact", "k_err",
        "σcr_num", "σcr_exact", "nodes")
    println("-"^145)

    for (r, γ) in table_biaxial
        k_exact, m_exact, n_exact = exact_biaxial_k(r, γ)
        σcr_exact = exact_σcr(k_exact)
        λcr, k_num, σcr_num, nnodes = solve_case(r, γ)
        k_err = abs(k_num-k_exact)/abs(k_exact)
        @printf("%7.2f %8.3f %3d %3d %18s %12.6e %14.10f %14.10f %12.4e %14.6f %14.6f %8d\n",
            r, γ, m_exact, n_exact, basename(mesh_file(r)), λcr, k_num, k_exact,
            k_err, σcr_num, σcr_exact, nnodes)
        push!(results, (
            r = r, γ = γ, m_exact = m_exact, n_exact = n_exact,
            mesh = basename(mesh_file(r)), λcr = λcr, k_num = k_num,
            k_exact = k_exact, k_err = k_err, σcr_num = σcr_num,
            σcr_exact = σcr_exact, nnodes = nnodes,
        ))
    end

    write_xlsx(output_xlsx, results)
    println("Excel output: ", output_xlsx)
finally
    gmsh.finalize()
end

println(to)
