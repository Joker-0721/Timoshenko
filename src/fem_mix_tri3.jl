using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ

using TimerOutputs, LinearAlgebra, WriteVTK
import Gmsh: gmsh

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
gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/patchtest_tri3_10.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
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
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"],integrationOrder)
    @timeit to "get elements" elements_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements_Γ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements
    𝑎ˢᵠ = ∫QφdΩ=>elements
    𝑎ˢˢ = ∫QQdΩ=>elements
    𝑎ˢʷ = [
        ∫∇QwdΩ=>elements,
        ∫QwdΓ=>elements_Γ,
    ]
    𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements
    𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements
    𝑎ᵐʷʷ = ∫ρwwdΩ=>elements
    𝑎ᵐᵠᵠ = ∫ρφφdΩ=>elements
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ᴳʷʷ(kᴳʷʷ)
    @timeit to "assemble" 𝑎ᴳᵠᵠ(kᴳᵠᵠ)
    @timeit to "assemble" 𝑎ᵐʷʷ(mʷʷ)
    @timeit to "assemble" 𝑎ᵐᵠᵠ(mᵠᵠ)
end

@timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(elements_1, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ʷ = ∫αwwdΓ=>elements_2∪elements_3∪elements_4
    # 𝑎ᵠ = ∫αφφdΓ=>elements_3
    @timeit to "assemble" 𝑎ʷ(kʷʷ)
    # @timeit to "assemble" 𝑎ᵠ(kᵠᵠ)
    # @timeit to "assemble" 𝑎ʷ(mʷʷ)
    # @timeit to "assemble" 𝑎ᵠ(mᵠᵠ)
    # @timeit to "assemble" 𝑎ʷ(kᴳʷʷ)
end

gmsh.finalize()

kᵠᵠ .+= - kˢᵠ'*(kˢˢ\kˢᵠ)
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