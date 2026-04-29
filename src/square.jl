using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Timoshenko:
    ∫κκdΩ,
    ∫wwdΩ,
    ∫φφdΩ,
    ∫φwdΩ,
    L₂,
    L₂φ,
    w_exact_cf,
    φ_exact_cf

using TimerOutputs
import Gmsh: gmsh

# square.jl:
# formula-based field comparison using Timoshenko beam exact solution

E = 1.0e8
ν = 0.3
κ = 5.0 / 6.0
h = 1.0e-3
L = 1.0
q₀ = 1.0
α = 1.0e8 * E

G = E / (2.0 * (1.0 + ν))
A = h
I = h^3 / 12.0
EI = E * I
kGA = κ * G * A

w(x, y, z) = w_exact_cf(x, E, I, κ, G, A, L, q₀)
φ(x, y, z) = φ_exact_cf(x, E, I, L, q₀)
q(x, y, z) = q₀
m(x, y, z) = 0.0
M(x, y, z) = 0.0
V(x, y, z) = 0.0
g(x, y, z) = 0.0
g₁(x, y, z) = 0.0

const to = TimerOutput()
L₂_w = NaN
L₂_φ = NaN

gmsh.initialize()
try
    @timeit to "open msh file" gmsh.open("msh/beam_seg2_16.msh")
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(nᵠ, nᵠ)
    kᵠʷ = zeros(nᵠ, nʷ)
    fʷ = zeros(nʷ)
    fᵠ = zeros(nᵠ)

    @timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ" begin
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
        prescribe!(elements, :EI => EI, :kGA => kGA, :q => q, :m => m)
        @timeit to "calculate shape functions" set∇𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ => elements
        𝑎ᵠʷ = ∫φwdΩ => elements
        𝑎ᵠᵠ = [
            ∫φφdΩ => elements,
            ∫κκdΩ => elements,
        ]
        𝑓ʷ = ∫wqdΩ => elements
        𝑓ᵠ = ∫φmdΩ => elements
        @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
        @timeit to "assemble" 𝑓ʷ(fʷ)
        @timeit to "assemble" 𝑓ᵠ(fᵠ)
    end

    @timeit to "calculate ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ" begin
        @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
        @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
        prescribe!(elements_1, :α => α, :g => g, :g₁ => g₁, :M => M, :V => V)
        prescribe!(elements_2, :α => α, :g => g, :g₁ => g₁, :M => M, :V => V)
        @timeit to "calculate shape functions" set𝝭!(elements_1)
        @timeit to "calculate shape functions" set𝝭!(elements_2)
        𝑎ʷ = ∫αwwdΓ => elements_1
        𝑎ᵠ = ∫αφφdΓ => elements_1
        𝑓V = ∫wVdΓ => elements_2
        𝑓M = ∫φMdΓ => elements_2
        @timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)
        @timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
        @timeit to "assemble" 𝑓V(fʷ)
        @timeit to "assemble" 𝑓M(fᵠ)
    end

    @timeit to "solve" d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]
    push!(nodes, :d => d[nᵠ + 1:end], :φ => d[1:nᵠ])

    @timeit to "calculate formula-based field error" begin
        @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 10)
        prescribe!(elements, :u => w, :φ => φ)
        @timeit to "calculate shape functions" set𝝭!(elements)
        global L₂_w = L₂(elements)
        global L₂_φ = L₂φ(elements)
    end
finally
    gmsh.finalize()
end

println(to)
println("Formula-based exact field from Timoshenko_Report.tex")
println("L₂ field error of w: ", L₂_w)
println("L₂ field error of φ: ", L₂_φ)
