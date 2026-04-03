using ApproxOperator
using LinearAlgebra
using TimerOutputs
import Gmsh: gmsh
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Timoshenko:∫wwdΩ, ∫φwdΩ, ∫κκdΩ, ∫φφdΩ, ∫wwGdΩ


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
    kgʷʷ = zeros(nʷ, nʷ)
    fʷ = zeros(nʷ)

    @timeit to "assemble material stiffness" begin
        elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :EI => EI, :kGA => kGA)
        set∇𝝭!(elements)
        (ApproxOperator.Timoshenko.∫wwdΩ => elements)(kʷʷ)
        (ApproxOperator.Timoshenko.∫φwdΩ => elements)(kᵠʷ)
        (ApproxOperator.Timoshenko.∫κκdΩ => elements)(kᵠᵠ)
        (ApproxOperator.Timoshenko.∫φφdΩ => elements)(kᵠᵠ)
        (ApproxOperator.Timoshenko.∫wwGdΩ => elements)(kgʷʷ)
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

    @timeit to "solve buckling EVP (Engesser)" begin
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        KG = [zeros(nᵠ, nᵠ) zeros(nᵠ, nʷ); zeros(nʷ, nᵠ) kgʷʷ]
        λ = eigvals(K, KG)
        λ = real.(λ[isfinite.(λ) .& (abs.(imag.(λ)) .< 1.0e-8)])
        λ = sort(λ[λ .> 0.0])
        if isempty(λ)
            error("No positive finite Engesser buckling eigenvalue found.")
        end
        global Pcr = first(λ)
    end
finally
    gmsh.finalize()
end

println(to)
println("Engesser P_cr = ", Pcr)
