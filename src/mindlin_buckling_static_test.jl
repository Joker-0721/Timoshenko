using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ

using LinearAlgebra
using Printf
using TimerOutputs
using WriteVTK
import Gmsh: gmsh

# Static consistency test for the buckling operator:
#
#     (K - λtest KG) d = f
#
# with a manufactured nodal field d_exact and f := (K - λtest KG)d_exact.
# This checks the assembled buckling blocks, DOF ordering, boundary penalty,
# and VTK export without changing the buckling formula used in mindlin.jl.

E = 200e9
ν = 0.3
h = 1.0
a = 1.0
b = 1.0
Dᵇ = E*h^3/12/(1 - ν^2)
k_exact = 4.0
λ_exact = k_exact*π^2*Dᵇ/b^2
λtest = 0.1*λ_exact
α = 1e8*E

σ₁₁(x, y, z) = 1.0
σ₂₂(x, y, z) = 0.0
σ₁₂(x, y, z) = 0.0

w_exact(x, y, z) = sin(π*x/a)*sin(π*y/b)
φ₁_exact(x, y, z) = -(π/a)*cos(π*x/a)*sin(π*y/b)
φ₂_exact(x, y, z) = -(π/b)*sin(π*x/a)*cos(π*y/b)
w_boundary(x, y, z) = 0.0

const to = TimerOutput()

function exact_vector(nodes)
    nʷ = length(nodes)
    d = zeros(3*nʷ)
    for node in nodes
        I = node.𝐼
        d[2*I - 1] = φ₁_exact(node.x, node.y, node.z)
        d[2*I] = φ₂_exact(node.x, node.y, node.z)
        d[2*nʷ + I] = w_exact(node.x, node.y, node.z)
    end
    return d
end

function relative_error(x, xref)
    denom = norm(xref)
    return denom == 0.0 ? norm(x - xref) : norm(x - xref)/denom
end

gmsh.initialize()
try
    integrationOrder = 2
    @timeit to "open msh file" gmsh.open("./msh/patchtest_quad4_16.msh")
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
    kᵠʷ = zeros(2*nᵠ, nʷ)
    kᴳʷʷ = zeros(nʷ, nʷ)
    kᴳᵠᵠ = zeros(2*nᵠ, 2*nᵠ)

    @timeit to "assemble domain" begin
        elements = getElements(nodes, entities["Ω"], integrationOrder)
        prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
        set∇𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ=>elements
        𝑎ᵠʷ = ∫φwdΩ=>elements
        𝑎ᵠᵠ = [
            ∫φφdΩ=>elements,
            ∫κκdΩ=>elements,
        ]
        𝑎ʷʷ(kʷʷ)
        𝑎ᵠʷ(kᵠʷ)
        𝑎ᵠᵠ(kᵠᵠ)

        prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements
        𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements
        𝑎ᴳʷʷ(kᴳʷʷ)
        𝑎ᴳᵠᵠ(kᴳᵠᵠ)

        global elements_domain = elements
    end

    @timeit to "assemble w boundary penalty" begin
        elements_1 = getElements(nodes, entities["Γ¹"], integrationOrder)
        elements_2 = getElements(nodes, entities["Γ²"], integrationOrder)
        elements_3 = getElements(nodes, entities["Γ³"], integrationOrder)
        elements_4 = getElements(nodes, entities["Γ⁴"], integrationOrder)
        prescribe!(elements_1, :α=>α, :g=>w_boundary)
        prescribe!(elements_2, :α=>α, :g=>w_boundary)
        prescribe!(elements_3, :α=>α, :g=>w_boundary)
        prescribe!(elements_4, :α=>α, :g=>w_boundary)
        set𝝭!(elements_1)
        set𝝭!(elements_2)
        set𝝭!(elements_3)
        set𝝭!(elements_4)
        𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
        fᵅ = zeros(nʷ)
        𝑎ʷ(kʷʷ, fᵅ)
    end

    @timeit to "solve static buckling operator test" begin
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        KG = [
            kᴳᵠᵠ zeros(2*nᵠ, nʷ)
            zeros(nʷ, 2*nᵠ) kᴳʷʷ
        ]
        A = K - λtest*KG
        d_exact = exact_vector(nodes)
        rhs = A*d_exact
        d = A\rhs

        global rel_d = relative_error(d, d_exact)
        global rel_w = relative_error(d[2*nᵠ+1:end], d_exact[2*nᵠ+1:end])
        global rel_φ = relative_error(d[1:2*nᵠ], d_exact[1:2*nᵠ])
        global solved_d = d
    end

    push!(nodes,
        :d => solved_d[2*nᵠ+1:end],
        :d₁ => solved_d[1:2:2*nᵠ],
        :d₂ => solved_d[2:2:2*nᵠ],
    )

    nₚ = length(nodes)
    points = zeros(3, nₚ)
    for (i, node) in enumerate(nodes)
        points[1, i] = node.x
        points[2, i] = node.y
        points[3, i] = node.d/5
    end
    cells = [MeshCell(VTKCellTypes.VTK_QUAD, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]
    mkpath("./vtk")
    vtk_grid("./vtk/mindlin_buckling_static_test.vtu", points, cells;
             ascii=true, append=false, compress=false) do vtk
        vtk["w"] = [node.d for node in nodes]
        vtk["phi_1"] = [node.d₁ for node in nodes]
        vtk["phi_2"] = [node.d₂ for node in nodes]
        vtk["relative_error_d", WriteVTK.VTKFieldData()] = [rel_d]
        vtk["relative_error_w", WriteVTK.VTKFieldData()] = [rel_w]
        vtk["relative_error_phi", WriteVTK.VTKFieldData()] = [rel_φ]
        vtk["lambda_test", WriteVTK.VTKFieldData()] = [λtest]
    end
finally
    gmsh.finalize()
end

println(to)
@printf("λtest = %.6e\n", λtest)
@printf("relative error of d   = %.6e\n", rel_d)
@printf("relative error of w   = %.6e\n", rel_w)
@printf("relative error of phi = %.6e\n", rel_φ)
println("VTK output: ./vtk/mindlin_buckling_static_test.vtu")
