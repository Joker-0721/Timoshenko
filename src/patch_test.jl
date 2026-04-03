using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Timoshenko:
    ∫κκdΩ,
    ∫wwdΩ,
    ∫φφdΩ,
    ∫φwdΩ,
    w_exact_ss,
    φ_exact_ss

using TimerOutputs
import Gmsh: gmsh

# patch_test.jl:
# direct benchmark comparison using Timoshenko beam exact numbers

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

w(x, y, z) = w_exact_ss(x, E, I, κ, G, A, L, q₀)
φ(x, y, z) = φ_exact_ss(x, E, I, L, q₀)
q(x, y, z) = q₀
m(x, y, z) = 0.0
M(x, y, z) = 0.0
V(x, y, z) = 0.0
g(x, y, z) = 0.0

# Direct exact numbers from Timoshenko_Report.tex (SS, h/L = 0.001)
x_w_ref = L / 2.0
x_φ_ref = L / 4.0
w_exact_value = 1.5625039
φ_exact_value = 3.4375

const to = TimerOutput()
L₂_w = NaN
L₂_φ = NaN
w_fem = NaN
φ_fem = NaN

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

    @timeit to "calculate ∫αwwdΓ, ∫wVdΓ, ∫φMdΓ" begin
        @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
        @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
        prescribe!(elements_1, :α => α, :g => g, :M => M, :V => V)
        prescribe!(elements_2, :α => α, :g => g, :M => M, :V => V)
        @timeit to "calculate shape functions" set𝝭!(elements_1)
        @timeit to "calculate shape functions" set𝝭!(elements_2)
        Γ = elements_1 ∪ elements_2
        𝑎ʷ = ∫αwwdΓ => Γ
        𝑓V = ∫wVdΓ => Γ
        𝑓M = ∫φMdΓ => Γ
        @timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)
        @timeit to "assemble" 𝑓V(fʷ)
        @timeit to "assemble" 𝑓M(fᵠ)
    end

    @timeit to "solve" d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]
    push!(nodes, :d => d[nᵠ + 1:end], :φ => d[1:nᵠ])

    @timeit to "calculate direct benchmark error" begin
        i_w = argmin(abs.([node.x - x_w_ref for node in nodes]))
        i_φ = argmin(abs.([node.x - x_φ_ref for node in nodes]))
        global w_fem = nodes[i_w].d
        global φ_fem = nodes[i_φ].φ
        global L₂_w = abs(w_fem - w_exact_value) / abs(w_exact_value)
        global L₂_φ = abs(φ_fem - φ_exact_value) / abs(φ_exact_value)
    end
finally
    gmsh.finalize()
end

println(to)
println("Direct exact benchmark comparison from Timoshenko_Report.tex")
println("w_exact(mid) = ", w_exact_value)
println("w_fem(mid)   = ", w_fem)
println("φ_exact(quarter) = ", φ_exact_value)
println("φ_fem(quarter)   = ", φ_fem)
println("relative error of w: ", L₂_w)
println("relative error of φ: ", L₂_φ)
