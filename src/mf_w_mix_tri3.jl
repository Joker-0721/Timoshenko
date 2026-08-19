using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ

using TimerOutputs, LinearAlgebra, WriteVTK
import Gmsh: gmsh
include("cal_area_support.jl")

E = 1.0
ν = 0.3
ρ = 1.0
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
σ₁₁ = 1e0
σ₂₂ = 0.0
σ₁₂ = 0.0
a = 1.0
αʷ = 0e6
αᵠ = 0e3

const to = TimerOutput()

integrationOrder = 2

type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})

ndiv = 10
ndiv_φ = ndiv
ndiv_w = ndiv-1
ndiv_q = ndiv

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
@timeit to "get entities" entities = getPhysicalGroups()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
nᵠ = length(nodes_φ)
elements_support = getElements(nodes_φ, entities["Ω"], 1)
s_φ, var_A = cal_area_support(elements_support)
s₁ = 1.5*s_φ*ones(nᵠ)
s₂ = 1.5*s_φ*ones(nᵠ)
s₃ = 1.5*s_φ*ones(nᵠ)
push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_q.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
nˢ = length(nodes)

nₑ = length(elements_support)
nᵐ = nₑ*ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[],𝑿ₛ[]))

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
kᵐᵐ = zeros(3*nᵐ,3*nᵐ)
kᵐᵠ = zeros(3*nᵐ,2*nᵠ)
kᵐʷ = zeros(3*nᵐ,nʷ)
kˢᵐ = zeros(2*nˢ,3*nᵐ)
fˢ = zeros(2*nˢ)
fᵐ = zeros(3*nᵐ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"],integrationOrder)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)

    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)

    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)

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
    𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    𝑎ᵐᵐ = ∫MMdΩ=>elements_m
    𝑎ᵐᵠ = [
        ∫∇MφdΩ=>(elements_m,elements_φ),
        ∫MφdΓ=>(elements_m_Γ,elements_φ_Γ),
    ]
    𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements_w
    𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements_φ
    𝑎ᵐʷʷ = ∫ρwwdΩ=>elements_w
    𝑎ᵐᵠᵠ = ∫ρφφdΩ=>elements_φ
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
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
    prescribe!(elements_w_1, :α=>αʷ, :g=>0.0)
    prescribe!(elements_w_2, :α=>αʷ, :g=>0.0)
    prescribe!(elements_w_3, :α=>αʷ, :g=>0.0)
    prescribe!(elements_w_4, :α=>αʷ, :g=>0.0)
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
    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ʷ(kʷʷ)
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
    prescribe!(elements_φ_1, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_2, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_3, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_4, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
    𝑎 = ∫MφdΓ => (elements_m_1 ∪ elements_m_2 ∪ elements_m_3 ∪ elements_m_4, elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4)
    𝑎ᵅ = ∫αφφdΓ => elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4
    # @timeit to "assemble" 𝑎(kᵐᵠ,fᵐ)
    # @timeit to "assemble" 𝑎ᵅ(kᵠᵠ)
end

gmsh.finalize()

kᵠᵠ .+= - kˢᵠ'*(kˢˢ\kˢᵠ) - kᵐᵠ'*(kᵐᵐ\kᵐᵠ)
kᵠʷ .+= - kˢᵠ'*(kˢˢ\kˢʷ)
kʷʷ .+= - kˢʷ'*(kˢˢ\kˢʷ)

k = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
kᴳ = [kᴳᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) kᴳʷʷ]
m = [mᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) mʷʷ]

λ,v = eigen(k,kᴳ)
# # λ,v = eigen(k,m)

index = findfirst(real.(λ).>1e-8)

n_index = 12
d = zeros(nˢ,n_index)
𝗠 = zeros(21)
for (i,xᵢ) in enumerate(nodes)
    x = xᵢ.x
    y = xᵢ.y
    indices = sp_w(x,y,0.0)
    ni = length(indices)
    𝓒 = [nodes_w[i] for i in indices]
    data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
    ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    𝓖 = [ξ]
    a = eval(type_w)(𝓒,𝓖)
    set𝝭!(a)
    for j in 1:n_index
        u = 0.0
        N = ξ[:𝝭]
        for (k,xₖ) in enumerate(𝓒)
            I = xₖ.𝐼
            u += N[k]*v[2nᵠ+I,index+j-1]
        end
        d[i,j] = u
    end
end

push!(nodes,
    :d₁=>d[:,1],
    :d₂=>d[:,2],
    :d₃=>d[:,3],
    :d₄=>d[:,4],
    :d₅=>d[:,5],
    :d₆=>d[:,6],
    :d₇=>d[:,7],
    :d₈=>d[:,8],
    :d₉=>d[:,9],
    :d₁₀=>d[:,10],
    :d₁₁=>d[:,11],
    :d₁₂=>d[:,12],
)



xs = [node.x for node in nodes]'
ys = [node.y for node in nodes]'
zs = [node.z for node in nodes]'
points = [xs; ys; zs]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_q]

vtk_grid("./vtk/mf_w_φ.vtu", points, cells;
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