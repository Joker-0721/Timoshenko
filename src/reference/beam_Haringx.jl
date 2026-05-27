using ApproxOperator
using LinearAlgebra
using TimerOutputs
import Gmsh: gmsh
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Timoshenko:∫wwdΩ, ∫φwdΩ, ∫κκdΩ, ∫φφdΩ, ∫φwGHdΩ, ∫φφGHdΩ, ∫wφGHdΩ, ∫wqdΩ

E = 1.0e8
ν = 0.3
κ = 5.0 / 6.0
h = 1.0e-3
b = 1.0
α = 1.0e8 * E

G = E / (2.0 * (1.0 + ν))
A = b * h
I = b * h^3 / 12.0
EI = E * I
kGA = κ * G * A

const to = TimerOutput()
Pcr = NaN

gmsh.initialize()
try
    @timeit to "open msh file" gmsh.open("msh/beam.msh")
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)

    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(nᵠ, nᵠ)
    kᵠʷ = zeros(nᵠ, nʷ)
    kgwᵠ = zeros(nʷ, nᵠ)
    kgᵠw = zeros(nᵠ, nʷ)
    kgᵠᵠ = zeros(nᵠ, nᵠ)
    fʷ = zeros(nʷ)

    @timeit to "assemble material stiffness" begin
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :EI => EI, :kGA => kGA)
        @timeit to "calculate shape functions" set∇𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ => elements
        𝑎ᵠʷ = ∫φwdΩ => elements
        𝑎ᵠᵠ = [
            ∫φφdΩ => elements,
            ∫κκdΩ => elements,
        ]
        𝑎gʷᵠ = ∫wφGHdΩ => elements
        𝑎gᵠw = ∫φwGHdΩ => elements
        𝑎gᵠᵠ = ∫φφGHdΩ => elements
        𝑓ʷ = ∫wqdΩ => elements
        𝑓ᵠ = ∫φmdΩ => elements

        @timeit to "assemble" begin
            𝑎ʷʷ(kʷʷ)
            𝑎ᵠʷ(kᵠʷ)
            𝑎ᵠᵠ(kᵠᵠ)
            𝑎gʷᵠ(kgwᵠ)
            𝑎gᵠw(kgᵠw)
            𝑎gᵠᵠ(kgᵠᵠ)
            𝑓ʷ(fʷ)
            𝑓ᵠ(fᵠ)
        end
    end

    @timeit to "assemble boundary penalty on w" begin
        elements_1 = getElements(nodes, entities["Γ¹"])
        elements_2 = getElements(nodes, entities["Γ²"])
        prescribe!(elements_1, :α => α, :g => (x, y, z) -> 0.0)
        prescribe!(elements_2, :α => α, :g => (x, y, z) -> 0.0)
        set𝝭!(elements_1)
        set𝝭!(elements_2)
        Γ = elements_1 ∪ elements_2
        (ApproxOperator.Timoshenko.∫αwwdΓ => Γ)(kʷʷ, fʷ)
    end

    @timeit to "solve buckling EVP (Haringx)" begin
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        KG = [-kgᵠᵠ kgᵠw; kgwᵠ zeros(nʷ, nʷ)]
        λ = eigvals(K, KG)
        λ = real.(λ[isfinite.(λ) .& (abs.(imag.(λ)) .< 1.0e-8)])
        λ = sort(λ[λ .> 0.0])
        if isempty(λ)
            error("No positive finite Haringx buckling eigenvalue found.")
        end
        global Pcr = first(λ)
    end
finally
    gmsh.finalize()
end

println(to)
println("Haringx P_cr = ", Pcr)
