using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫φφdΩ, ∫wwdΩ, ∫∇w∇wdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh
include("cal_area_support.jl")

E = 10.92e6
ν = 0.3
h = 1e-0
Dᵇ = E/12/(1-ν^2)
Dˢ = 5/6*E/h^2/(2*(1+ν))
integrationOrder = 2

const to = TimerOutput()

gmsh.initialize()
# ──────────────────────────────────────────────────────────
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})
# ndiv_φ = 4, nʷ = 11, nˢ = 21
# ndiv_φ = 4, nʷ = 11, nˢ = 21
ndiv_φ = 32
ndiv_w = ndiv_φ-1
ndiv_q = ndiv_φ
nʷ = 1034
# nˢ = 21
# ─── Deflection W ─────────────────────────────────────────
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_$nʷ.msh")
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_w.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp_w = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
@timeit to "get entities" entities = getPhysicalGroups()
elements_support = getElements(nodes_w, entities["Ω"], 1)
nʷ = length(nodes_w)
s_w, var_A = cal_area_support(elements_support)
s₁ = 1.5*s_w*ones(nʷ)
s₂ = 1.5*s_w*ones(nʷ)
s₃ = 1.5*s_w*ones(nʷ)
push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
# ─── Rotation Φ ───────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
nᵠ = length(nodes_φ)
@timeit to "get entities" entities = getPhysicalGroups()
elements_support = getElements(nodes_φ, entities["Ω"], 1)
s_φ, var_A = cal_area_support(elements_support)
s₁ = 1.5*s_φ*ones(nᵠ)
s₂ = 1.5*s_φ*ones(nᵠ)
s₃ = 1.5*s_φ*ones(nᵠ)
push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
# ─── Shear ────────────────────────────────────────────────
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_$nˢ.msh")
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_q.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
elements_support = getElements(nodes, entities["Ω"], 1)
s_q, var_A = cal_area_support(elements_support)

nˢ = length(nodes)
kˢˢ = zeros(2 * nˢ, 2 * nˢ)
kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
kˢʷ = zeros(2 * nˢ, nʷ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
kʷʷ = zeros(nʷ, nʷ)
fᵠ = zeros(2 * nᵠ)
fʷ = zeros(nʷ)
fˢ = zeros(2 * nˢ)

@timeit to "calculate ∫QQdΩ ∫∇QwdΩ" begin
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)

    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set𝝭!(elements_w)

    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

    @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)
 
    𝑎ˢˢ = ∫QQdΩ=>elements_q
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements_q,elements_w),
        ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
    ]
    𝑎ʷʷ = ∫wwdΩ=>elements_w
    𝑓ʷ = ∫wqdΩ=>elements_w
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
end

nₑ = length(elements_q)

@timeit to "calculate ∫MMdΩ ∫MφdΩ" begin
    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)

    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

    𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
end

αʷ = 1e8
@timeit to "calculate ∫QwdΓ" begin
    @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
    prescribe!(elements_w_1, :α=>αʷ)
    prescribe!(elements_w_2, :α=>αʷ)
    prescribe!(elements_w_3, :α=>αʷ)
    prescribe!(elements_w_4, :α=>αʷ)
    @timeit to "calculate shape functions" set𝝭!(elements_q_1)
    @timeit to "calculate shape functions" set𝝭!(elements_q_2)
    @timeit to "calculate shape functions" set𝝭!(elements_q_3)
    @timeit to "calculate shape functions" set𝝭!(elements_q_4)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
    𝑎 = ∫QwdΓ => (elements_q_1 ∪ elements_q_2 ∪ elements_q_3 ∪ elements_q_4, elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4)
    𝑎ʷ = ∫αwwdΓ => elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4
    @timeit to "assemble" 𝑎(kˢʷ)
end

println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ, nˢ = $nˢ")

k̃ᵠᵠ = - kˢᵠ'*(kˢˢ\kˢᵠ)
k̃ʷʷ = - kˢʷ'*(kˢˢ\kˢʷ)
k̃ᵠʷ = - kˢᵠ'*(kˢˢ\kˢʷ)
# ─── Eigen Test For βʷ ────────────────────────────────────
print("nʷ≤⌊[nˢ]⌋-1:         ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
βʷ² = eigvals(kˢʷ*(kʷʷ\kˢʷ')*s_q^(-2))
βʷ² = real.(βʷ²)
βʷ²⁺ = βʷ²[βʷ² .≥ 1e5*eps()]
βʷ⁺ = βʷ²⁺.^0.5
nʷ⁺ = length(βʷ⁺)
βʷ⁺ = min(βʷ⁺...)
println("βʷ⁺ = $βʷ⁺, nʷ⁺ = $nʷ⁺")

# ─── Eigen Test For βᵞ ────────────────────────────────────
print("2nˢ≤nʷ+2nᵠ-min(nʷ,nᶠ):")
nᶠ = floor(0.5*((1+8*nᵠ)^0.5-3))
nᶜ = ceil(0.5*((1+8*nᵠ)^0.5-3))
n = 0.5*(nᶠ+2)*(nᶠ+3)
n_diff_f = 0.5*(nʷ+2nᵠ-min(nʷ,n))-nˢ
n_diff_f≥0.0 ? println("✓:$n_diff_f") : println("×:$n_diff_f")
println("nʷ = $nʷ, n = $n")
print("2nˢ≤nʷ+2nᵠ-min(nʷ,nᶜ):")
n = 0.5*(nᶜ+2)*(nᶜ+3)
n_diff_c = 0.5*(nʷ+2nᵠ-min(nʷ,n))-nˢ
n_diff_c≥0.0 ? println("✓:$n_diff_c") : println("×:$n_diff_c")
println("nʷ = $nʷ, n = $n")

# βᵞ² = eigvals([k̃ᵠᵠ k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ]/Dˢ*s_φ^(-2))
βᵞ² = eigvals([k̃ᵠᵠ k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ]/Dˢ*(1/ndiv_φ)^(-2))
βᵞ² = real.(βᵞ²)
βᵞ²⁺ = βᵞ²[βᵞ² .≥ 1e5*eps()]
βᵞ⁺ = βᵞ²⁺.^0.5
nᵞ⁺ = length(βᵞ⁺)
βᵞ⁺ = min(βᵞ⁺...)
println("βᵞ⁺ = $βᵞ⁺, nᵞ⁺ = $nᵞ⁺")

logβʷ⁺ = log10(βʷ⁺)
logβᵞ⁺ = log10(βᵞ⁺)
println("$logβʷ⁺, $logβᵞ⁺")