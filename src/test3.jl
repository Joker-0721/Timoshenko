using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ, ∫wVdΓ, ∫φMdΓ
using WriteVTK
using LinearAlgebra
using TimerOutputs
import Gmsh: gmsh

# -----------------------------
# case switch (reference solution from 圆板精確解.md)
# -----------------------------

    E = 10.92
    ν = 0.3
    h = 1e0
    R = 5.0
    fz = 1.0
    α = 1e8 * E 
    Dᵇ = E*h^3/12/(1-ν^2)
    Dˢ = 5/6*E*h/(2*(1+ν))
    
 

# 制造解-----------------------------
# r = 1
# w(x, y, z) = (x + y)^r
# w₁(x, y, z) = r * (x + y)^abs(r - 1)
# w₂(x, y, z) = r * (x + y)^abs(r - 1)
# w₁₁(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# w₂₂(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₁(x, y, z) = r * (x + y)^abs(r - 1)
# φ₂(x, y, z) = r * (x + y)^abs(r - 1)
# φ₁₁(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₁₂(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₂₁(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₂₂(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₁₁₁(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₁₁₂(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₂₂₁(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₂₂₂(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₁₂₁(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₁₂₂(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)

# M₁₁(x, y, z) = -Dᵇ * (φ₁₁(x, y, z) + ν * φ₂₂(x, y, z))
# M₁₂(x, y, z) = -Dᵇ * (1 - ν) * 0.5 * (φ₁₂(x, y, z) + φ₂₁(x, y, z))
# M₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁(x, y, z) + φ₂₂(x, y, z))
# M₁₁₁(x, y, z) = -Dᵇ * (φ₁₁₁(x, y, z) + ν * φ₂₂₁(x, y, z))
# M₁₂₂(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₂(x, y, z)
# M₁₂₁(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₁(x, y, z)
# M₂₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁₂(x, y, z) + φ₂₂₂(x, y, z))

# Q₁(x, y, z) = Dˢ * (w₁(x, y, z) - φ₁(x, y, z))
# Q₂(x, y, z) = Dˢ * (w₂(x, y, z) - φ₂(x, y, z))
# Q₁₁(x, y, z) = Dˢ * (w₁₁(x, y, z) - φ₁₁(x, y, z))
# Q₂₂(x, y, z) = Dˢ * (w₂₂(x, y, z) - φ₂₂(x, y, z))
# q(x, y, z) = -Q₁₁(x, y, z) - Q₂₂(x, y, z)
# m₁(x, y, z) = M₁₁₁(x, y, z) + M₁₂₂(x, y, z) - Q₁(x, y, z)
# m₂(x, y, z) = M₁₂₁(x, y, z) + M₂₂₂(x, y, z) - Q₂(x, y, z)

# 解析解------------------------------------------------------------------------------

#  ξ(x, y) = sqrt(x^2 + y^2) / R
# w(x, y, z) = fz * R^4 / (64 * Dᵇ) * (1 - ξ(x, y)^2) * ((1 - ξ(x, y)^2) + (8 * (h / R)^2) / (3 * k * (1 - ν)))
# Mr(x, y, z) = fz * R^2 / 16 * (1 + ν) * (1 - ((3 + ν) / (1 + ν)) * ξ(x, y)^2)
# Mθ(x, y, z) = fz * R^2 / 16 * (1 + ν) * (1 - ((1 + 3ν) / (1 + ν)) * ξ(x, y)^2)
# Pi_int() = (fz^2 * R^6 * ν) / (384 * Dᵇ) * (1 + (4 * (h / R)^2) / (k * (1 - ν)))
# Tr(x, y, z) = -fz * sqrt(x^2 + y^2) / 2

# 补充解析解------------------------------------------------------------------------------

r(x, y) = sqrt(x^2 + y^2)
ξ(x, y) = r(x, y) / R

w(x, y, z) = fz * R^4 / (64 * Dᵇ) * (1 - ξ(x, y)^2) * ((1 - ξ(x, y)^2) + 8 * h^2 / (3 * (5/6) * R^2 * (1 - ν)))


φᵣ(x, y, z) = fz * r(x, y) * (r(x, y)^2 - R^2) / (16 * Dᵇ)
φ₁(x,y,z) = φᵣ(x,y,z)*x/r(x,y)
φ₂(x,y,z) = φᵣ(x,y,z)*y/r(x,y)

κᵣᵣ(x, y, z) = -fz * R^2 * (3 * ξ(x, y)^2 - 1) / (16 * Dᵇ)
κᵩᵩ(x, y, z) = -fz * R^2 * (ξ(x, y)^2 - 1) / (16 * Dᵇ)

γᵣ(x, y, z) = -fz * r(x, y) / Dˢ

Mᵣᵣ(x, y, z) = -fz * R^2 / 16 * ((3 + ν) * ξ(x, y)^2 - (1 + ν))
Mᵩᵩ(x, y, z) = -fz * R^2 / 16 * ((1 + 3ν) * ξ(x, y)^2 - (1 + ν))

Qᵣ(x, y, z) = -fz * r(x, y) / 2
# --------------------------------------------------------------------------------


# 自然
# （示意图边界验收时用零边界；解析解对比时也依然是 clamped 零边界）
# w0(x, y, z) = 0.0
# φ10(x, y, z) = 0.0
# φ20(x, y, z) = 0.0

# V(x, y, z) = 0.0
# M₁(x, y, z) = 0.0
# M₂(x, y, z) = 0.0



const to = TimerOutput()
n = 15
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("./msh/circular.msh")
# @timeit to "open msh file" gmsh.open("./msh/circular_tri3_9.msh")
@timeit to "open msh file" gmsh.open("./msh/circular_tri3_$n.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)

kʷʷ = zeros(nʷ, nʷ)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
fʷ = zeros(nʷ)
fᵠ = zeros(2 * nᵠ)

@timeit to "assemble domain" begin
    elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :E => E, :ν => ν, :h => h, :q => fz)
    set∇𝝭!(elements)

    𝑎ʷʷ = ∫wwdΩ => elements
    𝑎ᵠʷ = ∫φwdΩ => elements
    𝑎ᵠᵠ = [
        ∫φφdΩ => elements,
        ∫κκdΩ => elements,
    ]
    𝑓ʷ = ∫wqdΩ => elements
    # 𝑓ᵠ = ∫φmdΩ => elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)

    global elements_domain = elements 

end



# println("after domain: ‖fʷ‖₂=", norm(fʷ), "  ‖fᵠ‖₂=", norm(fᵠ))

@timeit to "Γᵉ" begin
    elements = getElements(nodes, entities["Γᵉ"])
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g => w,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 1.0,
        :n₁₂ => 0.0,
        :n₂₂ => 1.0,
    )
    𝑎ʷ = ∫αwwdΓ => elements
    𝑎ᵠ = ∫αφφdΓ => elements
    @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
 
end

@timeit to "Γˡ" begin
    elements = getElements(nodes, entities["Γˡ"],normal=true)
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 1.0,
        :n₁₂ => 0.0,
        :n₂₂ => 0.0,
    )
    𝑎ᵠ = ∫αφφdΓ => elements
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end

@timeit to "Γᵇ" begin
    elements = getElements(nodes, entities["Γᵇ"],normal=true)
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 0.0,
        :n₁₂ => 0.0,
        :n₂₂ => 1.0,
    )
    𝑎ᵠ = ∫αφφdΓ => elements
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end

# println("after bnd: ‖fʷ‖₂=", norm(fʷ), "  ‖fᵠ‖₂=", norm(fᵠ))

@timeit to "solve" d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]

# println("check d norms: ‖d‖₂=", norm(d), "  ‖f‖₂=", norm([fᵠ; fʷ]))
# println("check rhs norms: ‖fʷ‖₂=", norm(fʷ), "  ‖fᵠ‖₂=", norm(fᵠ))

push!(nodes, :d => d[2*nᵠ+1:end], :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])

# println("check node fields:")
# println("  node1: d=", nodes[1].d, " d₁=", nodes[1].d₁, " d₂=", nodes[1].d₂)
# println("  node2: d=", nodes[2].d, " d₁=", nodes[2].d₁, " d₂=", nodes[2].d₂)


@timeit to "calculate error" begin
    elements_err = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_err, :E => E, :ν => ν, :h => h, :u => w, :φ₁ => φ₁, :φ₂ => φ₂)
    set𝝭!(elements_err)
    global L₂_w = L₂w(elements_err)
    global L₂_φ = L₂φ(elements_err)
end


gmsh.finalize()

println(to)

println("α penalty: ", α)
println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)

# 图--------------------------------------------------------------------------------

# 坐标------------------------------------------------------------------------------
# nₚ = length(nodes)
# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d*4
#     # points[3,i] = us[i]*4
# end

nₚ = length(nodes)
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = node.d/15
    # points[3,i] = us[i]*4
end

# 二维------------------------------------------------------------------------------
# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
# vtk_grid("./vtk/hmd_2d/Tri3_d_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
#     # vtk["精确解"] = us
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# 三维------------------------------------------------------------------------------

# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/circular_tri3_"*string(ndiv), points, cells) do vtk
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]


# 三维误差-------------------------------------------------------------------------------------------------
# vtk_grid("./vtk/circular_tri3_" * string(ndiv) * ".vtu", points, cells;
#          ascii=true, append=false, compress=false) do vtk


#     vtk["L2_w", WriteVTK.VTKFieldData()] = [L₂_w]
#     vtk["L2_phi", WriteVTK.VTKFieldData()] = [L₂_φ]
# end
# -------------------------------------------------------------------------------------------------

# 三维变形------------------------------------------------------------------------------
# vtk_grid("./vtk/circular_Clamped_tri3_" * string(ndiv) * ".vtu", points, cells;
#          ascii=true, append=false, compress=false) do vtk
vtk_grid("./vtk/circular_Clamped_tri3_$n.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    # 挠度 w
    vtk["w"] = [node.d for node in nodes]
    # 转角 φ₁
    vtk["phi_1"] = [node.d₁ for node in nodes]
    # 转角 φ₂  
    vtk["phi_2"] = [node.d₂ for node in nodes]
end
# -------------------------------------------------------------------------------------------------