using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫αwwdΓ

using LinearAlgebra
using Printf
using TimerOutputs
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 1e-2
a = 1.0
b = 1.0
Dᵇ = E*h^3/12/(1-ν^2)
k_ref = 9.340
λcr_ref = k_ref*π^2*Dᵇ/b^2
msh_file = "msh/mindlin_ssss.msh"

w(x,y,z) = 0.0
σ₁₁(x,y,z) = 0.0
σ₂₂(x,y,z) = 0.0
σ₁₂(x,y,z) = 1.0

const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open(msh_file)
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
kgʷʷ = zeros(nʷ,nʷ)

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

    prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
    𝑎ᴳʷʷ = ∫wwGdΩ2D=>elements
    @timeit to "assemble" 𝑎ᴳʷʷ(kgʷʷ)
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
        zeros(2*nᵠ,2*nᵠ) zeros(2*nᵠ,nʷ)
        zeros(nʷ,2*nᵠ) kgʷʷ
    ]
    Keff = kʷʷ - kᵠʷ'*(kᵠᵠ\kᵠʷ)
    λ = eigvals(Keff, kgʷʷ)
    λ_positive = sort!(collect(real(λᵢ) for λᵢ in λ if isfinite(real(λᵢ)) && isfinite(imag(λᵢ)) && abs(imag(λᵢ)) < 1.0e-7 && real(λᵢ) > 0.0))
    isempty(λ_positive) && error("no positive finite buckling eigenvalue found")
    λcr = first(λ_positive)
    k_num = λcr*b^2/(π^2*Dᵇ)
    rel_error = abs(k_num - k_ref)/abs(k_ref)
end

gmsh.finalize()

println(to)
println("case: SSSS square plate, pure in-plane shear")
println("source: T2 9.7, pure shear table value for a/b = 1")
println("msh_file: ", msh_file)
@printf("a/b: %.6f\n", a/b)
@printf("h/b: %.6f\n", h/b)
println("λcr: ", λcr)
println("λcr_ref: ", λcr_ref)
println("k_num: ", k_num)
println("k_ref: ", k_ref)
println("rel_error: ", rel_error)
