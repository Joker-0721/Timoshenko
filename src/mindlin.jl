using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ

using LinearAlgebra
using TimerOutputs
using WriteVTK
import Gmsh: gmsh

E = 200e9
ν = 0.3
h = 1e-0
a = 1.0
b = 1.0
Dᵇ = E*h^3/12/(1-ν^2)
k_exact = 4.0
λcr_exact = k_exact*π^2*Dᵇ/b^2

σ₁₁ = 1.0
σ₂₂ = 0.0
σ₁₂ = 0.0

const to = TimerOutput()

gmsh.initialize()
integrationOrder = 3
integrationOrder_shear = 2
@timeit to "open msh file" gmsh.open("./msh/patchtest_quad4_4.msh")
# @timeit to "open msh file" gmsh.open("./msh/bui_2011_square_17x17.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
kᴳʷʷ = zeros(nʷ,nʷ)
kᴳᵠᵠ = zeros(2*nᵠ,2*nᵠ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"],integrationOrder)
    @timeit to "get shear elements" elements_s = getElements(nodes, entities["Ω"],integrationOrder_shear)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_s, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shear shape functions" set∇𝝭!(elements_s)
    𝑎ʷʷ = ∫wwdΩ=>elements_s
    𝑎ᵠʷ = ∫φwdΩ=>elements_s
    𝑎ᵠᵠ = [
        ∫φφdΩ=>elements_s,
        ∫κκdΩ=>elements,
    ]
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)

    prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
    𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements
    𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements
    @timeit to "assemble" 𝑎ᴳʷʷ(kᴳʷʷ)
    @timeit to "assemble" 𝑎ᴳᵠᵠ(kᴳᵠᵠ)

    global elements_domain = elements
    global elements_shear = elements_s
end

@timeit to "calculate ∫αwwdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(elements_1, :α=>1e8*E, :g=>0.0)
    prescribe!(elements_2, :α=>1e8*E, :g=>0.0)
    prescribe!(elements_3, :α=>1e8*E, :g=>0.0)
    prescribe!(elements_4, :α=>1e8*E, :g=>0.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ʷ(kʷʷ)
    # @timeit to "assemble" 𝑎ʷ(kᴳʷʷ)
end

@timeit to "solve buckling eigenvalue" begin
    K = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
    Kᴳ = [
        kᴳᵠᵠ zeros(2*nᵠ,nʷ)
        zeros(nʷ,2*nᵠ) kᴳʷʷ
    ]
    Keff = kʷʷ - kᵠʷ'*(kᵠᵠ\kᵠʷ)
    # λ = eigvals(Keff, kᴳʷʷ)
    F = eigen(K, Kᴳ)
    λ = F.values
    V = F.vectors

    mode_ids = sort!(
        collect(i for i in eachindex(λ)
            if isfinite(real(λ[i])) &&
               isfinite(imag(λ[i])) &&
               abs(imag(λ[i])) < 1.0e-7 &&
               real(λ[i]) > 0.0),
        by = i -> real(λ[i]),
    )
    isempty(mode_ids) && error("no positive finite buckling eigenvalue found")
    sort!(mode_ids, by = i -> abs(real(λ[i])*b^2/(π^2*Dᵇ) - k_exact))

    println(λ)
    λcr = real(λ[first(mode_ids)])
    k_num = λcr*b^2/(π^2*Dᵇ)
    rel_error = abs(k_num - k_exact)/abs(k_exact)
end
 
gmsh.finalize()

println(to)

println("λcr: ", λcr)
println("λcr_exact: ", λcr_exact)
println("k_num: ", k_num)
println("k_exact: ", k_exact)
println("rel_error: ", rel_error)


# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]
cells = [MeshCell(VTKCellTypes.VTK_QUAD, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]

nₚ = length(nodes)
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = 0.0
end

vtk_grid("./vtk/mindlin_Q4int_modes.vtu", points, cells;
# vtk_grid("./vtk/mindlin_T3_modes.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    vtk_mode_ids = collect(eachindex(λ))
    vtk["lambda_real", WriteVTK.VTKFieldData()] = real.(λ)
    vtk["lambda_imag", WriteVTK.VTKFieldData()] = imag.(λ)
    vtk["lambda_isfinite", WriteVTK.VTKFieldData()] = Float64.(isfinite.(real.(λ)) .& isfinite.(imag.(λ)))
    vtk["lambda_isinf", WriteVTK.VTKFieldData()] = Float64.(isinf.(real.(λ)) .| isinf.(imag.(λ)))
    vtk["lambda_isnan", WriteVTK.VTKFieldData()] = Float64.(isnan.(real.(λ)) .| isnan.(imag.(λ)))

    for (mode_rank, mode_id) in enumerate(vtk_mode_ids)
        dm_raw = real.(V[:, mode_id])
        dm = [isfinite(x) ? x : 0.0 for x in dm_raw]
        w = dm[2*nᵠ+1:end]
        phi_1 = dm[1:2:2*nᵠ]
        phi_2 = dm[2:2:2*nᵠ]
        dm_vtk = zeros(3,nₚ)
        dm_vtk[3,:] .= w

        # 挠度 w
        vtk["w$(mode_rank)"] = w
        # 转角 φ₁
        # vtk["phi_1_$(mode_rank)"] = phi_1
        # 转角 φ₂
        # vtk["phi_2_$(mode_rank)"] = phi_2
    end
end
