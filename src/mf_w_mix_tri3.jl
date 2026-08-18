using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ

using TimerOutputs, LinearAlgebra, WriteVTK
import Gmsh: gmsh
include("cal_area_support.jl")

E = 1.0
ν = 0.3
ρ = 1.0
h = 1e-2
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
σ₁₁ = 1e0
σ₂₂ = 0.0
σ₁₂ = 0.0
a = 1.0

const to = TimerOutput()

integrationOrder = 2

type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})

ndiv_φ = 16
ndiv_w = 14
ndiv_q = 16

gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/patchtest_tri3_$ndiv_w.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp_w = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
elements_support = getElements(nodes_w, entities["Ω"], 1)
s_w, var_A = cal_area_support(elements_support)
nʷ = length(nodes_w)
s₁ = 1.5*s_w*ones(nʷ)
s₂ = 1.5*s_w*ones(nʷ)
s₃ = 1.5*s_w*ones(nʷ)
push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
@timeit to "get entities" entities_φ = getPhysicalGroups()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
nᵠ = length(nodes_φ)
elements_support = getElements(nodes_φ, entities_φ["Ω"], 1)
s_φ, var_A = cal_area_support(elements_support)
s₁ = 1.5*s_φ*ones(nᵠ)
s₂ = 1.5*s_φ*ones(nᵠ)
s₃ = 1.5*s_φ*ones(nᵠ)
push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_q.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

nˢ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kˢˢ = zeros(2*nˢ,2*nˢ)
kᴳʷʷ = zeros(nʷ,nʷ)
kᴳᵠᵠ = zeros(2*nᵠ,2*nᵠ)
mʷʷ = zeros(nʷ,nʷ)
mᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
kˢʷ = zeros(2*nˢ,nʷ)
kˢᵠ = zeros(2*nˢ,2*nᵠ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"],integrationOrder)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)

    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set𝝭!(elements_w)

    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)

    @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
    prescribe!(elements_m, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_m)

    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

    @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)

    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

    @timeit to "get elements" elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
    @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

    𝑎ˢˢ = ∫QQdΩ=>elements_q
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements_q,elements_w),
        ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
    ]
    𝑎ᵐᵐ = ∫MMdΩ=>elements_m
    𝑎ᵐᵠ = [
        ∫∇MφdΩ=>(elements_m,elements_φ),
        ∫MφdΓ=>(elements_m_Γ,elements_φ_Γ),
    ]
    𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements_w
    𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements_φ
    𝑎ᵐʷʷ = ∫ρwwdΩ=>elements_w
    𝑎ᵐᵠᵠ = ∫ρφφdΩ=>elements_φ
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ᴳʷʷ(kᴳʷʷ)
    @timeit to "assemble" 𝑎ᴳᵠᵠ(kᴳᵠᵠ)
    @timeit to "assemble" 𝑎ᵐʷʷ(mʷʷ)
    @timeit to "assemble" 𝑎ᵐᵠᵠ(mᵠᵠ)
end

@timeit to "calculate  ∫QwdΓ" begin
    @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
    prescribe!(elements_w_1, :α=>αʷ, :g=>w)
    prescribe!(elements_w_2, :α=>αʷ, :g=>w)
    prescribe!(elements_w_3, :α=>αʷ, :g=>w)
    prescribe!(elements_w_4, :α=>αʷ, :g=>w)
    @timeit to "calculate shape functions" set𝝭!(elements_q_1)
    @timeit to "calculate shape functions" set𝝭!(elements_q_2)
    @timeit to "calculate shape functions" set𝝭!(elements_q_3)
    @timeit to "calculate shape functions" set𝝭!(elements_q_4)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
    𝑎 = ∫QwdΓ => (elements_q_1 ∪ elements_q_2 ∪ elements_q_3 ∪ elements_q_4, elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4)
    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
end

@timeit to "calculate ∫MφdΓ" begin
    @timeit to "get elements" elements_m_1 = getElements(entities["Γ¹"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_2 = getElements(entities["Γ²"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_3 = getElements(entities["Γ³"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_4 = getElements(entities["Γ⁴"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_φ_1 = getElements(nodes_φ, entities["Γ¹"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_2 = getElements(nodes_φ, entities["Γ²"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_3 = getElements(nodes_φ, entities["Γ³"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_4 = getElements(nodes_φ, entities["Γ⁴"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    prescribe!(elements_φ_1, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_2, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_3, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_4, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
    𝑎 = ∫MφdΓ => (elements_m_1 ∪ elements_m_2 ∪ elements_m_3 ∪ elements_m_4, elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4)
    𝑎ᵅ = ∫αφφdΓ => elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4
    # @timeit to "assemble" 𝑎(kᵐᵠ,fᵐ)
    # @timeit to "assemble" 𝑎ᵅ(kᵅᵠᵠ, fᵅᵠ)
end

gmsh.finalize()

kᵠᵠ .+= - kˢᵠ'*(kˢˢ\kˢᵠ) - kᵐᵠ'*(kᵐᵐ\kᵐᵠ)
kᵠʷ .+= - kˢᵠ'*(kˢˢ\kˢʷ)
kʷʷ .+= - kˢʷ'*(kˢˢ\kˢʷ)

k = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
kᴳ = [kᴳᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) kᴳʷʷ]
m = [mᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) mʷʷ]

λ,v = eigen(k,kᴳ)
# λ,v = eigen(k,m)

index = findfirst(real.(λ).>1e-8)

d₁ = real.(v[2nᵠ+1:2nᵠ+nʷ,index])
d₂ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+1])
d₃ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+2])
d₄ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+3])
d₅ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+4])
d₆ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+5])
d₇ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+6])
d₈ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+7])
d₉ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+8])
d₁₀ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+9])
d₁₁ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+10])
d₁₂ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+11])
push!(nodes,
    :d₁=>d₁,
    :d₂=>d₂,
    :d₃=>d₃,
    :d₄=>d₄,
    :d₅=>d₅,
    :d₆=>d₆,
    :d₇=>d₇,
    :d₈=>d₈,
    :d₉=>d₉,
    :d₁₀=>d₁₀,
    :d₁₁=>d₁₁,
    :d₁₂=>d₁₂,
)


xs = [node.x for node in nodes]'
ys = [node.y for node in nodes]'
zs = [node.z for node in nodes]'
points = [xs; ys; zs]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]

vtk_grid("./vtk/fem.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    vtk["v₁"] = [node.d₁ for node in nodes]
    vtk["v₂"] = [node.d₂ for node in nodes]
    vtk["v₃"] = [node.d₃ for node in nodes]
    vtk["v₄"] = [node.d₄ for node in nodes]
    vtk["v₅"] = [node.d₅ for node in nodes]
    vtk["v₆"] = [node.d₆ for node in nodes]
    vtk["v₇"] = [node.d₇ for node in nodes]
    vtk["v₈"] = [node.d₈ for node in nodes]
    vtk["v₉"] = [node.d₉ for node in nodes]
    vtk["v₁₀"] = [node.d₁₀ for node in nodes]
    vtk["v₁₁"] = [node.d₁₁ for node in nodes]
    vtk["v₁₂"] = [node.d₁₂ for node in nodes]
end

# (λ.*ρ/Dˢ).^0.5
println(λ[index]*a^2/(π^2*Dᵇ)*h)