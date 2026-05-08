using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫αwwdΓ

using LinearAlgebra
using Printf
using TimerOutputs
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 1e-1
a = 1.0
b = 1.0
δ0 = 1.0
δ = 1.0
γ = 1.0/δ0
Dᵇ = E*h^3/12/(1-ν^2)
k_exact_biaxial = 2.0
λcr_exact_biaxial = k_exact_biaxial*π^2*Dᵇ/b^2
msh_file = "msh/mindlin_ssss.msh"

w(x,y,z) = 0.0
σ₁₁_biaxial(x,y,z) = 1.0
σ₂₂_biaxial(x,y,z) = γ
σ₁₂_biaxial(x,y,z) = 0.0
σ₁₁_shear(x,y,z) = 1.0
σ₂₂_shear(x,y,z) = 0.0
σ₁₂_shear(x,y,z) = 1.0/δ

const to = TimerOutput()

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

function solve_λcr(kᵠᵠ, kᵠʷ, kʷʷ, kgʷʷ)
    Keff = kʷʷ - kᵠʷ'*(kᵠᵠ\kᵠʷ)
    λ_positive = positive_finite_eigenvalues(Keff, kgʷʷ)
    isempty(λ_positive) && error("no positive finite buckling eigenvalue found")
    return first(λ_positive)
end

gmsh.initialize()
@timeit to "open msh file" gmsh.open(msh_file)
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
kgʷʷ_biaxial = zeros(nʷ,nʷ)
kgʷʷ_shear = zeros(nʷ,nʷ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ, ∫wwGdΩ2D" begin
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

    prescribe!(elements, :σ₁₁=>σ₁₁_biaxial, :σ₂₂=>σ₂₂_biaxial, :σ₁₂=>σ₁₂_biaxial)
    𝑎ᴳʷʷ_biaxial = ∫wwGdΩ2D=>elements
    @timeit to "assemble" 𝑎ᴳʷʷ_biaxial(kgʷʷ_biaxial)

    prescribe!(elements, :σ₁₁=>σ₁₁_shear, :σ₂₂=>σ₂₂_shear, :σ₁₂=>σ₁₂_shear)
    𝑎ᴳʷʷ_shear = ∫wwGdΩ2D=>elements
    @timeit to "assemble" 𝑎ᴳʷʷ_shear(kgʷʷ_shear)
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

@timeit to "solve combined buckling eigenvalues" begin
    K = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
    KG_biaxial = [
        zeros(2*nᵠ,2*nᵠ) zeros(2*nᵠ,nʷ)
        zeros(nʷ,2*nᵠ) kgʷʷ_biaxial
    ]
    KG_shear = [
        zeros(2*nᵠ,2*nᵠ) zeros(2*nᵠ,nʷ)
        zeros(nʷ,2*nᵠ) kgʷʷ_shear
    ]
    λcr_biaxial = solve_λcr(kᵠᵠ, kᵠʷ, kʷʷ, kgʷʷ_biaxial)
    λcr_shear = solve_λcr(kᵠᵠ, kᵠʷ, kʷʷ, kgʷʷ_shear)
    k_num_biaxial = λcr_biaxial*b^2/(π^2*Dᵇ)
    k_num_shear = λcr_shear*b^2/(π^2*Dᵇ)
    rel_error_biaxial = abs(k_num_biaxial - k_exact_biaxial)/abs(k_exact_biaxial)
end

gmsh.finalize()

println(to)
println("case 1: SSSS square plate, combined biaxial compression")
println("source: T1 9.3")
println("msh_file: ", msh_file)
@printf("δ0 = σ₁₁/σ₂₂: %.6f\n", δ0)
@printf("γ = σ₂₂/σ₁₁: %.6f\n", γ)
println("λcr_biaxial: ", λcr_biaxial)
println("λcr_exact_biaxial: ", λcr_exact_biaxial)
println("k_num_biaxial: ", k_num_biaxial)
println("k_exact_biaxial: ", k_exact_biaxial)
println("rel_error_biaxial: ", rel_error_biaxial)
println()
println("case 2: SSSS square plate, combined uniaxial compression and pure shear")
println("source: T2 9.7 curve/energy-method reference")
@printf("δ = σ₁₁/σ₁₂: %.6f\n", δ)
println("λcr_shear_combined: ", λcr_shear)
println("k_num_shear_combined: ", k_num_shear)
