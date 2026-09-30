using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ

using TimerOutputs, LinearAlgebra, WriteVTK, DelimitedFiles
import Gmsh: gmsh
include("cal_area_support.jl")

E = 1.0
ν = 0.3
ρ = 1.0
h = 1e-2
G = E/(2*(1+ν))
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
σ₁₁ = 1e0
σ₂₂ = 0.0
σ₁₂ = 0.0
a = 1.0
αʷ = 0e6
αᵠ = 0e3

const to = TimerOutput()
open("../date/new_vib/FSCS_mf_w_mix_tri3.csv", "w") do io
write(io, "ndiv,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40\n")

ndivs = [5,10,20]
for ndiv in ndivs

integrationOrder = 4

type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :tri3
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})

ndiv_φ = ndiv
ndiv_w = ndiv-1
ndiv_q = ndiv

gmsh.initialize()
@timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv_w.msh")
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

@timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
nᵠ = length(nodes_φ)

@timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv_q.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
nˢ = length(nodes)

# nₑ 必須取自 M 場所在的網格（= 當前 gmsh 模型 = q/φ 網格），
# 不能再沿用 w 網格的 elements_support，否則 ndiv_w ≠ ndiv_q 時 kᵐᵐ 尺寸不足
elements_support = getElements(nodes, entities["Ω"], 1)
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

    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], integrationOrder)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)

    @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
    prescribe!(elements_m, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_m)

    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

    @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)

    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], integrationOrder,  normal=true)
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
    𝑎 = ∫QwdΓ => (
            elements_q_2 ∪ 
            elements_q_3 ∪ 
            elements_q_4,
            elements_w_2 ∪ 
            elements_w_3 ∪ 
            elements_w_4)

    𝑎ʷ = ∫αwwdΓ =>
    elements_w_2 ∪ elements_w_3 ∪ elements_w_4

    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ʷ(kʷʷ)
end

@timeit to "calculate ∫MφdΓ" begin
    @timeit to "get elements" elements_m_1 = getElements(entities["Γ¹"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_2 = getElements(entities["Γ²"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_3 = getElements(entities["Γ³"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_4 = getElements(entities["Γ⁴"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_φ_1 = getElements(nodes_φ, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_φ_2 = getElements(nodes_φ, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_φ_3 = getElements(nodes_φ, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_φ_4 = getElements(nodes_φ, entities["Γ⁴"], integrationOrder, normal=true)
    prescribe!(elements_φ_1, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_2, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_3, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_4, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
    𝑎 = ∫MφdΓ => (
        elements_m_3,
        elements_φ_3)

    𝑎ᵅ = ∫αφφdΓ =>
    elements_φ_3

    @timeit to "assemble" 𝑎(kᵐᵠ,fᵐ)
    @timeit to "assemble" 𝑎ᵅ(kᵠᵠ)
end

gmsh.finalize()

kᵠᵠ .+= - kˢᵠ'*(kˢˢ\kˢᵠ) - kᵐᵠ'*(kᵐᵐ\kᵐᵠ)
kᵠʷ .+= - kˢᵠ'*(kˢˢ\kˢʷ)
kʷʷ .+= - kˢʷ'*(kˢˢ\kˢʷ)

k = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
kᴳ = [kᴳᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) kᴳʷʷ]
m = [mᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) mʷʷ]

# λ,v = eigen(k,kᴳ)
λ,v = eigen(k,m)

index = findfirst(real.(λ).>1e-8)

n_index = 40
index = clamp(index, 1, max(1, length(λ) - n_index))
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





# (λ.*ρ/Dˢ).^0.5
# println(λ[index]*a^2/(π^2*Dᵇ)*h)

# println(λ[index:end])
κ = real.((λ[index:index+40-1]).^0.5 .* a .* (ρ/G)^0.5)
println("FSCS n=$ndiv_w/$ndiv Ω = ", round.(κ, digits=6))
write(io, "$ndiv_w," * join(string.(κ), ",") * "\n")


end
end
