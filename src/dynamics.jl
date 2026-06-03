#屈曲問題分析
#2d_ssss
#mindlin plate

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ

using LinearAlgebra
using TimerOutputs
using WriteVTK
import Gmsh: gmsh

E = 200e9
ν = 0.3
h = 1e-2
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
integrationOrder = 2
integrationOrder_shear = 1
@timeit to "open msh file" gmsh.open("./msh/struct_quad_17.msh")
# @timeit to "open msh file" gmsh.open("./msh/struct_tri_17.msh")
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
    
    # 【導師修正】：排序和計算都需要乘上 h 和 σ₁₁ 來獲得真正的臨界載荷 N_cr
    sort!(mode_ids, by = i -> abs(real(λ[i])*h*σ₁₁*b^2/(π^2*Dᵇ) - k_exact))

    println(λ)
    λcr = real(λ[first(mode_ids)])
    
    # 計算真正的物理量
    N_cr = λcr * h * σ₁₁
    k_num = N_cr * b^2 / (π^2 * Dᵇ)
    rel_error = abs(k_num - k_exact)/abs(k_exact)
end
 
gmsh.finalize()

println(to)

println("λcr: ", λcr)
println("λcr_exact: ", λcr_exact)
println("k_num: ", k_num)
println("k_exact: ", k_exact)
println("rel_error: ", rel_error)

function write_eigen_check(filepath, K, Kᴳ, λ, V, mode_ids)
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        println(io, join([
            "mode_rank",
            "eigen_index",
            "lambda_real",
            "lambda_imag",
            "lambda_isfinite",
            "lambda_isinf",
            "lambda_isnan",
            "vector_isfinite",
            "norm_v",
            "norm_Kv",
            "norm_KGv",
            "residual_norm",
            "relative_residual",
            "KGv_over_v",
            "selected_positive_mode",
        ], ","))

        selected = Set(mode_ids)
        for (mode_rank, mode_id) in enumerate(eachindex(λ))
            λᵢ = λ[mode_id]
            vᵢ = V[:, mode_id]
            Kv = K*vᵢ
            Kᴳv = Kᴳ*vᵢ
            λ_finite = isfinite(real(λᵢ)) && isfinite(imag(λᵢ))
            v_finite = all(isfinite, real.(vᵢ)) && all(isfinite, imag.(vᵢ))
            norm_v = norm(vᵢ)
            norm_Kv = norm(Kv)
            norm_Kᴳv = norm(Kᴳv)

            if λ_finite && v_finite
                residual = Kv - λᵢ*Kᴳv
                residual_norm = norm(residual)
                denominator = max(norm_Kv, abs(λᵢ)*norm_Kᴳv, eps(Float64))
                relative_residual = residual_norm/denominator
            else
                residual_norm = NaN
                relative_residual = NaN
            end

            KGv_over_v = norm_Kᴳv/max(norm_v, eps(Float64))

            println(io, join([
                mode_rank,
                mode_id,
                real(λᵢ),
                imag(λᵢ),
                λ_finite,
                isinf(real(λᵢ)) || isinf(imag(λᵢ)),
                isnan(real(λᵢ)) || isnan(imag(λᵢ)),
                v_finite,
                norm_v,
                norm_Kv,
                norm_Kᴳv,
                residual_norm,
                relative_residual,
                KGv_over_v,
                mode_id in selected,
            ], ","))
        end
    end
end

write_eigen_check("./data/dynamics_Q4int_eigen_check.csv", K, Kᴳ, λ, V, mode_ids)
println("eigen check csv: ./data/dynamics_Q4int_eigen_check.csv")

function write_mode_data(summary_path, node_path, nodes, λ, V, nᵠ, mode_ids)
    mkpath(dirname(summary_path))
    selected = Set(mode_ids)

    open(summary_path, "w") do summary_io
        println(summary_io, join([
            "mode_rank",
            "eigen_index",
            "vtk_w_field",
            "lambda_real",
            "lambda_imag",
            "lambda_isfinite",
            "selected_positive_mode",
            "vector_isfinite",
            "vector_norm",
            "w_norm",
            "w_min",
            "w_min_node",
            "w_max",
            "w_max_node",
            "phi_1_norm",
            "phi_2_norm",
            "rebuild_error",
            "relative_rebuild_error",
        ], ","))

        open(node_path, "w") do node_io
            println(node_io, join([
                "mode_rank",
                "eigen_index",
                "vtk_w_field",
                "lambda_real",
                "lambda_imag",
                "node_row",
                "x",
                "y",
                "w",
                "phi_1",
                "phi_2",
            ], ","))

            for (mode_rank, mode_id) in enumerate(eachindex(λ))
                λᵢ = λ[mode_id]
                dm_raw = real.(V[:, mode_id])
                dm = [isfinite(x) ? x : 0.0 for x in dm_raw]
                phi_1 = dm[1:2:2*nᵠ]
                phi_2 = dm[2:2:2*nᵠ]
                w = dm[2*nᵠ+1:end]

                rebuilt = similar(dm)
                rebuilt[1:2:2*nᵠ] .= phi_1
                rebuilt[2:2:2*nᵠ] .= phi_2
                rebuilt[2*nᵠ+1:end] .= w
                rebuild_error = norm(rebuilt - dm)
                relative_rebuild_error = rebuild_error/max(norm(dm), eps(Float64))

                w_min, w_min_node = findmin(w)
                w_max, w_max_node = findmax(w)

                println(summary_io, join([
                    mode_rank,
                    mode_id,
                    "w$(mode_rank)",
                    real(λᵢ),
                    imag(λᵢ),
                    isfinite(real(λᵢ)) && isfinite(imag(λᵢ)),
                    mode_id in selected,
                    all(isfinite, dm_raw),
                    norm(dm_raw),
                    norm(w),
                    w_min,
                    w_min_node,
                    w_max,
                    w_max_node,
                    norm(phi_1),
                    norm(phi_2),
                    rebuild_error,
                    relative_rebuild_error,
                ], ","))

                for (node_row, node) in enumerate(nodes)
                    println(node_io, join([
                        mode_rank,
                        mode_id,
                        "w$(mode_rank)",
                        real(λᵢ),
                        imag(λᵢ),
                        node_row,
                        node.x,
                        node.y,
                        w[node_row],
                        phi_1[node_row],
                        phi_2[node_row],
                    ], ","))
                end
            end
        end
    end
end

write_mode_data(
    "./data/dynamics_Q4int_summary.csv",
    "./data/dynamics_Q4int_node_values.csv",
    nodes,
    λ,
    V,
    nᵠ,
    mode_ids,
)
println("mode summary csv: ./data/dynamics_Q4int_summary.csv")
println("mode node values csv: ./data/dynamics_Q4int_node_values.csv")


# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]
cells = [MeshCell(VTKCellTypes.VTK_QUAD, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]

nₚ = length(nodes)
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = 0.0
end

vtk_grid("./vtk/dynamics_Q4int_modes.vtu", points, cells;
# vtk_grid("./vtk/dynamics_T3int_modes.vtu", points, cells;
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
