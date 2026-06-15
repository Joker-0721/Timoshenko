using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using LinearAlgebra
using TimerOutputs, WriteVTK, XLSX
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 1e-0
Dᵇ = E * h^3 / 12 / (1 - ν^2)
Dˢ = 5 / 6 * E * h / (2 * (1 + ν))

w(x, y, z) = 1.0 + x + y
w₁(x, y, z) = 1.0
w₂(x, y, z) = 1.0
w₁₁(x, y, z) = 0.0
w₂₂(x, y, z) = 0.0
φ₁(x, y, z) = 1.0 + x + y
φ₂(x, y, z) = 1.0 + x + y
φ₁₁(x, y, z) = 1.0
φ₁₂(x, y, z) = 1.0
φ₂₁(x, y, z) = 1.0
φ₂₂(x, y, z) = 1.0
φ₁₁₁(x, y, z) = 0.0
φ₁₁₂(x, y, z) = 0.0
φ₂₂₁(x, y, z) = 0.0
φ₂₂₂(x, y, z) = 0.0
φ₁₂₁(x, y, z) = 0.0
φ₁₂₂(x, y, z) = 0.0

# r = 1
# w(x,y,z) = (x+y)^r
# w₁(x,y,z) = r*(x+y)^abs(r-1)
# w₂(x,y,z) = r*(x+y)^abs(r-1)
# w₁₁(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
# w₂₂(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
# φ₁(x,y,z) = r*(x+y)^abs(r-1)
# φ₂(x,y,z) = r*(x+y)^abs(r-1)
# φ₁₁(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₁₂(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₂₁(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₂₂(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₁₁₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₁₁₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₂₂₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₂₂₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₁₂₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₁₂₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)

M₁₁(x, y, z) = -Dᵇ * (φ₁₁(x, y, z) + ν * φ₂₂(x, y, z))
M₁₂(x, y, z) = -Dᵇ * (1 - ν) * 0.5 * (φ₁₂(x, y, z) + φ₂₁(x, y, z))
M₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁(x, y, z) + φ₂₂(x, y, z))
M₁₁₁(x, y, z) = -Dᵇ * (φ₁₁₁(x, y, z) + ν * φ₂₂₁(x, y, z))
M₁₂₂(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₂(x, y, z)
M₁₂₁(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₁(x, y, z)
M₂₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁₂(x, y, z) + φ₂₂₂(x, y, z))

Q₁(x, y, z) = Dˢ * (w₁(x, y, z) - φ₁(x, y, z))
Q₂(x, y, z) = Dˢ * (w₂(x, y, z) - φ₂(x, y, z))
Q₁₁(x, y, z) = Dˢ * (w₁₁(x, y, z) - φ₁₁(x, y, z))
Q₂₂(x, y, z) = Dˢ * (w₂₂(x, y, z) - φ₂₂(x, y, z))
q(x, y, z) = -Q₁₁(x, y, z) - Q₂₂(x, y, z)
m₁(x, y, z) = M₁₁₁(x, y, z) + M₁₂₂(x, y, z) - Q₁(x, y, z)
m₂(x, y, z) = M₁₂₁(x, y, z) + M₂₂₂(x, y, z) - Q₂(x, y, z)


const to = TimerOutput()

gmsh.initialize()
# @timeit to "open msh file" gmsh.open("msh/patchtest_3.msh")
# @timeit to "get nodes" nodes_s = get𝑿ᵢ()

integrationOrder = 2
# ──────────────────────────────────────────────────────────
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})
ndiv_φ = 16
ndiv_w = 14
ndiv = 16
sʷ = 1.5
sᵠ = 1.5
# XLSX.openxlsx("xls/patchtest.xlsx", mode="w") do xf
# for ndiv = 4:16
# ndiv_w = ndiv
# row = ndiv
# ─── Deflection W ─────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_high_un_tri3_$ndiv_w.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp_w = RegularGrid(xʷ, yʷ, zʷ, n=3, γ=5)
nʷ = length(nodes_w)
s = 1 / ndiv_w
s₁ = sʷ * s * ones(nʷ)
s₂ = sʷ * s * ones(nʷ)
s₃ = sʷ * s * ones(nʷ)
push!(nodes_w, :s₁ => s₁, :s₂ => s₂, :s₃ => s₃)
# ─── Rotation Φ ───────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_high_un_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ, yᵠ, zᵠ, n=3, γ=5)
nᵠ = length(nodes_φ)
s = 1 / ndiv_φ
s₁ = sᵠ * s * ones(nᵠ)
s₂ = sᵠ * s * ones(nᵠ)
s₃ = sᵠ * s * ones(nᵠ)
push!(nodes_φ, :s₁ => s₁, :s₂ => s₂, :s₃ => s₃)
# ─── Shear ────────────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_high_un_tri3_$ndiv.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

nˢ = length(nodes)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kʷʷ = zeros(nʷ, nʷ)
kˢˢ = zeros(2 * nˢ, 2 * nˢ)
kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
kˢʷ = zeros(2 * nˢ, nʷ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
fᵠ = zeros(2 * nᵠ)
fʷ = zeros(nʷ)
fˢ = zeros(2 * nˢ)
kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵅʷʷ = zeros(nʷ, nʷ)
fᵅᵠ = zeros(2 * nᵠ)
fᵅʷ = zeros(nʷ)

@timeit to "calculate ∫QQdΩ ∫∇QwdΩ" begin
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(elements_q, :E => E, :ν => ν, :h => h, :m₁ => m₁, :m₂ => m₂, :q => q)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)

    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements_w, :E => E, :ν => ν, :h => h, :q => q)
    @timeit to "calculate shape functions" set𝝭!(elements_w)

    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

    @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)

    𝑎ˢˢ = ∫QQdΩ => elements_q
    𝑎ˢʷ = [
        ∫∇QwdΩ => (elements_q, elements_w),
        ∫QwdΓ => (elements_q_Γ, elements_w_Γ),
    ]
    𝑓ʷ = ∫wqdΩ => elements_w
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
end

nₑ = length(elements_q)
nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))
nᵐ = nₑ * nᵖ
kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ)
kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)
kᵐʷ = zeros(3 * nᵐ, nʷ)
kˢᵐ = zeros(2 * nˢ, 3 * nᵐ)
fᵐ = zeros(3 * nᵐ)

@timeit to "calculate ∫MMdΩ ∫MφdΩ" begin
    @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
    prescribe!(elements_m, :E => E, :ν => ν, :h => h, :m₁ => m₁, :m₂ => m₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_m)

    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :m₁ => m₁, :m₂ => m₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)

    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

    @timeit to "get elements" elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
    @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

    𝑎ᵠᵠ = ∫κκdΩ => elements_φ
    𝑎ᵐᵐ = ∫MMdΩ => elements_m
    𝑎ᵐᵠ = [
        ∫∇MφdΩ => (elements_m, elements_φ),
        ∫MφdΓ => (elements_m_Γ, elements_φ_Γ),
    ]
    𝑎ˢᵠ = ∫QφdΩ => (elements_q, elements_φ)
    𝑓ᵠ = ∫φmdΩ => elements_φ
    @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

αʷ = 0e1;αᵠ = 0e1;
@timeit to "calculate ∫QwdΓ" begin
    @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
    prescribe!(elements_w_1, :α => αʷ, :g => w)
    prescribe!(elements_w_2, :α => αʷ, :g => w)
    prescribe!(elements_w_3, :α => αʷ, :g => w)
    prescribe!(elements_w_4, :α => αʷ, :g => w)
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
    @timeit to "assemble" 𝑎(kˢʷ, fˢ)
    @timeit to "assemble" 𝑎ʷ(kᵅʷʷ, fᵅʷ)
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
    prescribe!(elements_φ_1, :α => αᵠ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_φ_2, :α => αᵠ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_φ_3, :α => αᵠ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_φ_4, :α => αᵠ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
    𝑎 = ∫MφdΓ => (elements_m_1 ∪ elements_m_2 ∪ elements_m_3 ∪ elements_m_4, elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4)
    𝑎ᵅ = ∫αφφdΓ => elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4
    @timeit to "assemble" 𝑎(kᵐᵠ, fᵐ)
    @timeit to "assemble" 𝑎ᵅ(kᵅᵠᵠ, fᵅᵠ)
end

# ──────────────────────────────────────────────────────────
# dᵠ = zeros(2*nᵠ)
# dʷ = zeros(nʷ)
# dˢ = zeros(2*nˢ)
# dᵐ = zeros(3*nᵐ)
# for node in nodes_φ
#     x = node.x
#     y = node.y
#     z = node.z
#     dᵠ[2*node.𝐼-1] = φ₁(x,y,z)
#     dᵠ[2*node.𝐼]   = φ₂(x,y,z)
# end

# for node in nodes_w
#     x = node.x
#     y = node.y
#     z = node.z
#     dʷ[node.𝐼] = w(x,y,z)
# end

# for node in nodes
#     x = node.x
#     y = node.y
#     z = node.z
#     dˢ[2*node.𝐼-1] = Q₁(x,y,z)
#     dˢ[2*node.𝐼]   = Q₂(x,y,z)
# end
# for i in 1:nₑ
#     dᵐ[9*i-8] = M₁₁(0.0,0.0,0.0)
#     dᵐ[9*i-7] = M₂₂(0.0,0.0,0.0)
#     dᵐ[9*i-6] = M₁₂(0.0,0.0,0.0)
#     dᵐ[9*i-5] = M₁₁(1.0,0.0,0.0)-M₁₁(0.0,0.0,0.0)
#     dᵐ[9*i-4] = M₂₂(1.0,0.0,0.0)-M₂₂(0.0,0.0,0.0)
#     dᵐ[9*i-3] = M₁₂(1.0,0.0,0.0)-M₁₂(0.0,0.0,0.0)
#     dᵐ[9*i-2] = M₁₁(0.0,1.0,0.0)-M₁₁(0.0,0.0,0.0)
#     dᵐ[9*i-1] = M₂₂(0.0,1.0,0.0)-M₂₂(0.0,0.0,0.0)
#     dᵐ[9*i]   = M₁₂(0.0,1.0,0.0)-M₁₂(0.0,0.0,0.0)
# end
# println(kᵠᵠ*dᵠ+kᵐᵠ'*dᵐ+kˢᵠ'*dˢ - fᵠ)
# println(Q₁(1.0,0.0,0.0)-M₁₁₁(1.0,0.0,0.0)- M₁₂₂(1.0,0.0,0.0))
# println(Q₂(1.0,0.0,0.0)-M₁₂₁(1.0,0.0,0.0)- M₂₂₂(1.0,0.0,0.0))
# println(kˢᵠ'*dˢ+kᵐᵠ'*dᵐ)
# println(fᵠ)
# println(kˢᵠ'*dˢ+kᵐᵠ'*dᵐ - fᵠ)
# println(kᵠᵠ*dᵠ+kᵠʷ*dʷ+kˢᵠ'*dˢ+kᵐᵠ'*dᵐ - fᵠ)
# println(kᵠʷ'*dᵠ+kʷʷ*dʷ+kˢʷ'*dˢ+kᵐʷ'*dᵐ - fʷ)
# println(kˢᵠ*dᵠ+kˢʷ*dʷ+kˢˢ*dˢ+kˢᵐ*dᵐ - fˢ)
# println(kᵐᵠ*dᵠ+kᵐʷ*dʷ+kˢᵐ'*dˢ+kᵐᵐ*dᵐ - fᵐ)
# println(norm(kᵠᵠ*dᵠ+kᵠʷ*dʷ+kˢᵠ'*dˢ+kᵐᵠ'*dᵐ - fᵠ))
# println(norm(kᵠʷ'*dᵠ+kʷʷ*dʷ+kˢʷ'*dˢ+kᵐʷ'*dᵐ - fʷ))
# println(norm(kˢᵠ*dᵠ+kˢʷ*dʷ+kˢˢ*dˢ+kˢᵐ*dᵐ - fˢ))
# println(norm(kᵐᵠ*dᵠ+kᵐʷ*dʷ+kˢᵐ'*dˢ+kᵐᵐ*dᵐ - fᵐ))

# println(kˢˢ*dˢ - (fˢ-kˢᵠ*dᵠ-kˢʷ*dʷ))
# println(dˢ - kˢˢ\(fˢ-kˢᵠ*dᵠ-kˢʷ*dʷ))
# println(kᵐᵐ*dᵐ - (fᵐ-kᵐᵠ*dᵠ-kᵐʷ*dʷ))
# println(dᵐ - kᵐᵐ\(fᵐ-kᵐᵠ*dᵠ-kᵐʷ*dʷ))
# println(kᵠᵠ*dᵠ+kᵠʷ*dʷ+kˢᵠ'*dˢ+kᵐᵠ'*dᵐ - fᵠ)
# println(kˢᵠ'*dˢ+kᵐᵠ'*dᵐ - fᵠ)
# println(kˢᵠ'*(kˢˢ\(fˢ-kˢᵠ*dᵠ-kˢʷ*dʷ))+kᵐᵠ'*(kᵐᵐ\(fᵐ-kᵐᵠ*dᵠ-kᵐʷ*dʷ)) - fᵠ)
# println((kˢᵠ'*(kˢˢ\fˢ)-kˢᵠ'*(kˢˢ\kˢᵠ)*dᵠ-kˢᵠ'*(kˢˢ\kˢʷ)*dʷ)+(kᵐᵠ'*(kᵐᵐ\fᵐ)-kᵐᵠ'*(kᵐᵐ\kᵐᵠ)*dᵠ-kᵐᵠ'*(kᵐᵐ\kᵐʷ)*dʷ) - fᵠ)
# println(-(kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ))*dᵠ-(kˢᵠ'*(kˢˢ\kˢʷ)+kᵐᵠ'*(kᵐᵐ\kᵐʷ))*dʷ - (fᵠ-kˢᵠ'*(kˢˢ\fˢ)-kᵐᵠ'*(kᵐᵐ\fᵐ)))
# println(kˢʷ'*dˢ+kᵐʷ'*dᵐ-fʷ)
# println(kˢʷ'*(kˢˢ\(fˢ-kˢᵠ*dᵠ-kˢʷ*dʷ))+kᵐʷ'*(kᵐᵐ\(fᵐ-kᵐᵠ*dᵠ-kᵐʷ*dʷ))-fʷ)
# println((kˢʷ'*(kˢˢ\fˢ)-kˢʷ'*(kˢˢ\kˢᵠ)*dᵠ-kˢʷ'*(kˢˢ\kˢʷ)*dʷ)+(kᵐʷ'*(kᵐᵐ\fᵐ)-kᵐʷ'*(kᵐᵐ\kᵐᵠ)*dᵠ-kᵐʷ'*(kᵐᵐ\kᵐʷ)*dʷ)-fʷ)
# println(-(kˢʷ'*(kˢˢ\kˢᵠ)+kᵐʷ'*(kᵐᵐ\kᵐᵠ))*dᵠ-(kˢʷ'*(kˢˢ\kˢʷ)+kᵐʷ'*(kᵐᵐ\kᵐʷ))*dʷ-(fʷ-kˢʷ'*(kˢˢ\fˢ)-kᵐʷ'*(kᵐᵐ\fᵐ)))
k = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ kˢᵠ'*(kˢˢ\kˢʷ)+kᵐᵠ'*(kᵐᵐ\kᵐʷ); kˢʷ'*(kˢˢ\kˢᵠ)+kᵐʷ'*(kᵐᵐ\kᵐᵠ) kˢʷ'*(kˢˢ\kˢʷ)+kᵐʷ'*(kᵐᵐ\kᵐʷ)-kᵅʷʷ]
f = [fᵠ - kˢᵠ' * (kˢˢ \ fˢ) - kᵐᵠ' * (kᵐᵐ \ fᵐ) + fᵅᵠ; fʷ - kˢʷ' * (kˢˢ \ fˢ) - kᵐʷ' * (kᵐᵐ \ fᵐ) + fᵅʷ]
# println(kᵅʷʷ*dʷ-fᵅʷ)
# println(kᵅᵠᵠ*dᵠ-fᵅᵠ)
# println(norm(k*[dᵠ;dʷ]-f))
# ──────────────────────────────────────────────────────────
# @timeit to "solve" d = [kᵠᵠ+kᵅᵠᵠ kᵠʷ kˢᵠ' kᵐᵠ';kᵠʷ' kʷʷ+kᵅʷʷ kˢʷ' kᵐʷ';kˢᵠ kˢʷ kˢˢ kˢᵐ;kᵐᵠ kᵐʷ kˢᵐ' kᵐᵐ]\[fᵠ+fᵅᵠ;fʷ+fᵅʷ;fˢ;fᵐ]
# push!(nodes_φ,:d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ])
# push!(nodes_w,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])
# push!(nodes,:q₁=>d[2*nᵠ+nʷ+1:2:2*nᵠ+nʷ+2*nˢ], :q₂=>d[2*nᵠ+nʷ+2:2:2*nᵠ+nʷ+2*nˢ])
# push!(nodes,:m₁₁=>d[2*nᵠ+nʷ+2*nˢ+1:3:end],:m₂₂=>d[2*nᵠ+nʷ+2*nˢ+2:3:end],:m₁₂=>d[2*nᵠ+nʷ+2*nˢ+3:3:end])

## ──────────────────────────────────────────────────────────
d = k \ f
dᵠ = d[1:2*nᵠ]
dʷ = d[2*nᵠ+1:2*nᵠ+nʷ]
dˢ = kˢˢ \ (fˢ - kˢᵠ * dᵠ - kˢʷ * dʷ)
dᵐ = kᵐᵐ \ (fᵐ - kᵐᵠ * dᵠ - kᵐʷ * dʷ)

push!(nodes_φ, :d₁ => dᵠ[1:2:2*nᵠ], :d₂ => dᵠ[2:2:2*nᵠ])
push!(nodes_w, :d => dʷ)
push!(nodes, :q₁ => dˢ[1:2:end], :q₂ => dˢ[2:2:end])
push!(nodes, :m₁₁ => dᵐ[1:3:end], :m₂₂ => dᵐ[2:3:end], :m₁₂ => dᵐ[3:3:end])

@timeit to "calculate error" begin
    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), 10, sp_φ)
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), 10, sp_w)
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], 10)
    @timeit to "get elements" elements_m = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :φ₁ => φ₁, :φ₂ => φ₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)
    prescribe!(elements_w, :E => E, :ν => ν, :h => h, :w => w)
    @timeit to "calculate shape functions" set𝝭!(elements_w)
    prescribe!(elements_q, :E => E, :ν => ν, :h => h, :Q₁ => Q₁, :Q₂ => Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements_q)
    prescribe!(elements_m, :E => E, :ν => ν, :h => h, :M₁₁ => Q₁, :Q₂ => Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements_q)
end

# points = zeros(3, nᵛ)
# for node in nodes_q
#     I = node.𝐼
#     points[1,I] = node.x
#     points[2,I] = node.y
#     points[3,I] = node.z
# end
# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [node.𝐼 for node in elm.𝓒]) for elm in elements]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [node.𝐼 for node in elm.𝓒]) for elm in elements_q]
# vtk_grid("vtk/square.vtu", points, cells) do vtk
#     vtk["Q₁"] = [node.q₁ for node in nodes_q]
#     vtk["Q₂"] = [node.q₂ for node in nodes_q]
#     vtk["Q̄₁"] = [Q₁(node.x,node.y,node.z) for node in nodes_q]
#     vtk["Q̄₂"] = [Q₂(node.x,node.y,node.z) for node in nodes_q]
# end

# println(to)

@timeit to "calculate error" begin
    L₂_w = L₂w(elements_w)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_q)
end
println("h=$h,αʷ=$αʷ,αᵠ=$αᵠ,nᵠ=$nᵠ,nʷ=$nʷ,nˢ=$nˢ,nᵐ=$nᵐ")
println("L₂ error of w: ", log10(L₂_w))
println("L₂ error of φ: ", log10(L₂_φ))
println("L₂ error of Q: ", log10(L₂_Q))

# ──────────────────────────────────────────────────────────
#     sheet = xf[1]
#     XLSX.rename!(sheet, "new_sheet")
#     sheet["A1"] = "type w"
#     sheet["B1"] = "nʷ"
#     sheet["C1"] = "type φ"
#     sheet["D1"] = "nᵠ"
#     sheet["E1"] = "type Q"
#     sheet["F1"] = "nˢ"
#     sheet["G1"] = "type M"
#     sheet["H1"] = "nᵐ"
#     sheet["I1"] = "L₂w"
#     sheet["J1"] = "L₂φ"
#     sheet["K1"] = "L₂Q"
#     sheet["A$row"] = "$type_w"
#     sheet["B$row"] = nʷ
#     sheet["C$row"] = "$type_φ"
#     sheet["D$row"] = nᵠ
#     sheet["E$row"] = "$type_Q"
#     sheet["F$row"] = nˢ
#     sheet["G$row"] = "$type_M"
#     sheet["H$row"] = nᵐ
#     sheet["I$row"] = log10(L₂_w)
#     sheet["J$row"] = log10(L₂_φ)
#     sheet["K$row"] = log10(L₂_Q)
# end
# end
gmsh.finalize()

