# 屈曲問題分析 —— 完美對齊 Dawe & Roufaeil (1982) 算例
# 2d_ssss (Hard S2 邊界)
# mindlin plate

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ, ∫αφφdΓ

using LinearAlgebra
using TimerOutputs
using WriteVTK
import Gmsh: gmsh

E = 200e9
ν = 0.3
h = 1e-2  # 【論文參數研究點】：可在此切換 1e-2 (0.01), 5e-2 (0.05), 1e-1 (0.10)
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
@timeit to "open msh file" gmsh.open("./msh/st_q_17.msh")
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

@timeit to "calculate boundary conditions (Hard SSSS - S2)" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    
    # 賦予邊界大剛度懲罰
# Γ¹: 下邊界 (y = 0), 外法線 n = (0, -1) ➔ n₁=0, n₂=-1 ➔ n₁₁=1, n₁₂=0, n₂₂=0
    prescribe!(elements_1, :α=>1e8*E, :g=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>0.0)
    
    # Γ²: 右邊界 (x = 1), 外法線 n = (1, 0)  ➔ n₁=1, n₂=0  ➔ n₁₁=0, n₁₂=0, n₂₂=1
    prescribe!(elements_2, :α=>1e8*E, :g=>0.0, :n₁₁=>0.0, :n₁₂=>0.0, :n₂₂=>1.0)
    
    # Γ³: 上邊界 (y = 1), 外法線 n = (0, 1)  ➔ n₁=0, n₂=1  ➔ n₁₁=1, n₁₂=0, n₂₂=0
    prescribe!(elements_3, :α=>1e8*E, :g=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>0.0)
    
    # Γ⁴: 左邊界 (x = 0), 外法線 n = (-1, 0) ➔ n₁=-1, n₂=0 ➔ n₁₁=0, n₁₂=0, n₂₂=1
    prescribe!(elements_4, :α=>1e8*E, :g=>0.0, :n₁₁=>0.0, :n₁₂=>0.0, :n₂₂=>1.0)
    
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    
    # 1. 限制面外撓度 w = 0
    𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble w BC" 𝑎ʷ(kʷʷ)
    
    # 【導師修改一】：調用 thick_plate.jl 的轉角邊界函數，限制切向轉角以實現 Hard SSSS (S2)
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble phi BC" 𝑎ᵠ(kᵠᵠ)
end

@timeit to "solve buckling eigenvalue" begin
    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    
    # 【導師修改二】：拋棄 kᴳᵠᵠ，強制左上角為零矩陣，只保留真正的面外屈曲幾何剛度 kᴳʷʷ
    Kᴳ = [
        zeros(2*nᵠ, 2*nᵠ)    zeros(2*nᵠ, nʷ)
        zeros(nʷ, 2*nᵠ)      kᴳʷʷ
    ]
    
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
    
    # 根據物理公式正確排序
    sort!(mode_ids, by = i -> abs(real(λ[i])*h*σ₁₁*b^2/(π^2*Dᵇ) - k_exact))

    λcr = real(λ[first(mode_ids)])
    
    # 計算真正的物理量
    N_cr = λcr * h * σ₁₁
    k_num = N_cr * b^2 / (π^2 * Dᵇ)
    rel_error = abs(k_num - k_exact)/abs(k_exact)
end
 
gmsh.finalize()

println(to)

println("λcr (應力放大係數): ", λcr)
println("k_num (數值屈曲係數): ", k_num)
println("k_exact (薄板解析解): ", k_exact)
println("rel_error (相對於薄板極限的偏差): ", rel_error)

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

write_eigen_check("./data/buckling/buckling_Q4int_eigen_check.csv", K, Kᴳ, λ, V, mode_ids)
println("eigen check csv: ./data/buckling/buckling_Q4int_eigen_check.csv")

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
    "./data/buckling/buckling_Q4int_summary.csv",
    "./data/buckling/buckling_Q4int_node_values.csv",
    nodes,
    λ,
    V,
    nᵠ,
    mode_ids,
)
println("mode summary csv: ./data/buckling/buckling_Q4int_summary.csv")
println("mode node values csv: ./data/buckling/buckling_Q4int_node_values.csv")


# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]
# 建立四邊形網格連接關係
cells = [MeshCell(VTKCellTypes.VTK_QUAD, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]

nₚ = length(nodes)
points = zeros(3, nₚ)
for (i, node) in enumerate(nodes)
    points[1, i] = node.x
    points[2, i] = node.y
    points[3, i] = 0.0
end

# 確保 vtk 資料夾存在
mkpath("./vtk")

# 只挑選前 10 個最關鍵的有效正實數屈曲模態（避免 Inf 垃圾解污染選單）
max_output_modes = min(10, length(mode_ids))
selected_vtk_modes = mode_ids[1:max_output_modes]

vtk_grid("./vtk/buckling_Q4int_split.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    # 寫入全局特徵值資料 (Field Data)
    lambdas_real_box = [real(λ[m_id]) for m_id in selected_vtk_modes]
    vtk["Selected_Lambda_Real", WriteVTK.VTKFieldData()] = lambdas_real_box

    # 循環解開每個選定模態的自由度
    for (mode_rank, mode_id) in enumerate(selected_vtk_modes)
        dm_raw = real.(V[:, mode_id])
        dm = [isfinite(x) ? x : 0.0 for x in dm_raw]
        
        # 【核心拆解】：根據矩陣自由度編排，精確切出獨立的力學物理場
        phi_1 = dm[1:2:2*nᵠ]       # x 方向轉角
        phi_2 = dm[2:2:2*nᵠ]       # y 方向轉角
        w     = dm[2*nᵠ+1:end]     # 面外橫向撓度
        
        # 1. 寫入面外撓度 w (純量場) —— 用來做 Warp By Scalar 隆起泡泡
        vtk["Mode_$(mode_rank)_w", VTKPointData()] = w
        
        # 2. 寫入轉角 φ₁ (純量場) —— 觀察面內截面繞 y 軸的旋轉應變
        vtk["Mode_$(mode_rank)_phi1", VTKPointData()] = phi_1
        
        # 3. 寫入轉角 φ₂ (純量場) —— 觀察面內截面繞 x 軸的旋轉應變
        vtk["Mode_$(mode_rank)_phi2", VTKPointData()] = phi_2
        
        # 4. 【高級技巧】：將 phi1 和 phi2 組裝成一個 3D 向量場 (z分量設為0)
        # 這樣你可以在 ParaView 裡直接呼叫 Glyph (箭頭濾鏡) 觀察截面旋轉的流向！
        phi_vec = zeros(3, nₚ)
        phi_vec[1, :] .= phi_1
        phi_vec[2, :] .= phi_2
        vtk["Mode_$(mode_rank)_phi_vector", VTKPointData()] = phi_vec
    end
end

println("【成功】已生成包含獨立 w, phi1, phi2 與轉角向量場的 VTU 檔案：./vtk/buckling_Q4int.vtu")
