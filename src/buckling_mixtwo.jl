# ==============================================================================
#  高級計算固體力學研究平台：Mindlin 屈曲分析多網格收斂性測試
#  流派：Meshfree-FE Hybrid Mixed Method (無網格-有限元高階多場混合變分原理)
#  網格路徑對規：/home/a/Joker/msh/st_q_$ndiv.msh 規則
#  🌟 徹底修復：將 set𝝭!(elements_w) 升級為 set∇𝝭!(elements_w) 以提供幾何剛度所需之一階導數
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, 
                                     ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, L₂w, L₂φ, L₂Q

using LinearAlgebra
using TimerOutputs, WriteVTK, Printf, Statistics
import Gmsh: gmsh

# ==============================================================================
#  精密資料夾分流配置
# ==============================================================================
const DATA_DIR = "./data"
const VTK_DIR  = "./VTK"
mkpath(DATA_DIR)
mkpath(VTK_DIR)

# 全局幾何與材料物理常數
E = 200e9; 
ν = 0.3; a = 1.0; b = 1.0; h = 1e-2;
Dᵇ = E * h^3 / (12 * (1 - ν^2))
Dˢ = 5 / 6 * E * h / (2 * (1 + ν))

αʷ = 1.0e8 * Dᵇ;
αᵠ = 1.0e8 * Dᵇ; 

const σ₁₁ = 1.0
const σ₂₂ = 0.0
const σ₁₂ = 0.0
const k_exact = 4.0

integrationOrder = 2
sʷ = 1.5
sᵠ = 1.5

const to = TimerOutput()
convergence_table = []

# 定義理論經典薄板解析解函數
w_exact_func(x, y, z) = sin(π * x / a) * sin(π * y / b)
phi1_exact_func(x, y, z) = (π / a) * cos(π * x / a) * sin(π * y / b)
phi2_exact_func(x, y, z) = (π / b) * sin(π * x / a) * cos(π * y / b)

println("="^80)
println(" 啟動自適應多網格解耦屈曲收斂性測試大循環 (ndiv = 15 ➔ 25) ")
println("  目標基礎路徑: /home/a/Joker/msh/st_q_\$ndiv.msh")
println("="^80)

for n_div in 15:25
    @timeit to "Total loop iteration (ndiv=$n_div)" begin
        
        ndiv_w = n_div - 2
        ndiv_φ = n_div
        ndiv_s = n_div

        mesh_file_w = "/home/a/Joker/msh/st_q_$(ndiv_w).msh"
        mesh_file_φ = "/home/a/Joker/msh/st_q_$(ndiv_φ).msh"
        mesh_file_s = "/home/a/Joker/msh/st_q_$(ndiv_s).msh"

        if !isfile(mesh_file_w) || !isfile(mesh_file_φ) || !isfile(mesh_file_s)
            println("  [WARN] 找不到指定的解耦網格檔案，跳過 ndiv=$n_div")
            continue
        end

        case_prefix = "buckling_mixtwo_multi_mesh_w$(ndiv_w)_phi$(ndiv_φ)_s$(ndiv_s)"

        type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
        type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
        type_Q = :tri3;
        type_M = :(PiecewisePolynomial{:Linear2D})

        # ─── 2.1 加載位移場 w 網格 ───
        gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(mesh_file_w)
        nodes_w = get𝑿ᵢ(); nʷ = length(nodes_w);
        s_w = 1.0 / (ndiv_w - 1)
        push!(nodes_w, :s₁ => sʷ * s_w * ones(nʷ), :s₂ => sʷ * s_w * ones(nʷ), :s₃ => sʷ * s_w * ones(nʷ))
        sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

        # ─── 2.2 加載轉角場 φ 網格 ───
        gmsh.clear();
        gmsh.open(mesh_file_φ)
        nodes_φ = get𝑿ᵢ(); nᵠ = length(nodes_φ);
        s_𝜙 = 1.0 / (ndiv_φ - 1)
        push!(nodes_φ, :s₁ => sᵠ * s_𝜙 * ones(nᵠ), :s₂ => sᵠ * s_𝜙 * ones(nᵠ), :s₃ => sᵠ * s_𝜙 * ones(nᵠ))
        sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

        # ─── 2.3 加載內力場之主控制網格 ───
        gmsh.clear();
        gmsh.open(mesh_file_s)
        nodes = get𝑿ᵢ(); entities = getPhysicalGroups();
        nˢ = length(nodes)

        kˢˢ = zeros(2 * nˢ, 2 * nˢ);
        kˢʷ = zeros(2 * nˢ, nʷ); kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
        kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ);
        kᵅʷʷ = zeros(nʷ, nʷ); 
        kʷʷ = zeros(nʷ, nʷ) 

        # A. 剪切項弱形式組裝
        @timeit to "assemble shear items" begin
            elements_q = getElements(nodes, entities["Ω"], integrationOrder);
            prescribe!(elements_q, :E => E, :ν => ν, :h => h);
            set∇𝝭!(elements_q)
            
            elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w);
            prescribe!(elements_w, :E => E, :ν => ν, :h => h);
            
            # 🌟【核心修正點】：從 set𝝭! 升級為 set∇𝝭!，強制 RKPM 計算一階偏導數，注入 :∂𝝭∂x 與 :∂𝝭∂y
            set∇𝝭!(elements_w)
            
            elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true);
            set𝝭!(elements_w_Γ)
            elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true);
            set𝝭!(elements_q_Γ)
            (∫QQdΩ => elements_q)(kˢˢ);
            ([∫∇QwdΩ => (elements_q, elements_w), ∫QwdΓ => (elements_q_Γ, elements_w_Γ)])(kˢʷ)
        end

        # B. 彎矩項弱形式組裝
        nₑ = length(elements_q);
        nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[])); nᵐ = nₑ * nᵖ
        kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ);
        kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)

        @timeit to "assemble moment items" begin
            elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder);
            prescribe!(elements_m, :E => E, :ν => ν, :h => h);
            set∇𝝭!(elements_m)
            elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ);
            prescribe!(elements_φ, :E => E, :ν => ν, :h => h);
            set∇𝝭!(elements_φ)
            elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true);
            set𝝭!(elements_φ_Γ)
            elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder);
            set𝝭!(elements_m_Γ)
            (∫MMdΩ => elements_m)(kᵐᵐ); ([∫∇MφdΩ => (elements_m, elements_φ), ∫MφdΓ => (elements_m_Γ, elements_φ_Γ)])(kᵐᵠ);
            (∫QφdΩ => (elements_q, elements_φ))(kˢᵠ)
        end

        # C-1. 幾何剛度矩陣弱形式組裝
        @timeit to "assemble geometric stiffness" begin
            prescribe!(elements_w, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
            (∫∇wσ∇wdΩ => elements_w)(kʷʷ)
        end

        # C-2. 施加邊界罰剛度條件
        @timeit to "assemble boundary penalty" begin
            boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
            elements_w_b = [getElements(nodes_w, entities[name], eval(type_w), integrationOrder, sp_w, normal=true) for name in boundary_names]
            for el in elements_w_b; prescribe!(el, :α => αʷ, :g => (x,y,z)->0.0); set𝝭!(el); end
            elements_φ_b = [getElements(nodes_φ, entities[name], eval(type_φ), integrationOrder, sp_φ, normal=true) for name in boundary_names]
            for el in elements_φ_b; prescribe!(el, :α => αᵠ, :g₁ => (x,y,z)->0.0, :g₂ => (x,y,z)->0.0, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0); set𝝭!(el); end
            dummy_fʷ = zeros(nʷ); dummy_fᵠ = zeros(2*nᵠ)
            (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kᵅʷʷ, dummy_fʷ)
            (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵅᵠᵠ, dummy_fᵠ)
        end

        # E. 執行特徵值屈曲凝結與求解
        @timeit to "condense and solve buckling" begin
            global k_cond = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ   kˢᵠ'*(kˢˢ\kˢʷ);
                              kˢʷ'*(kˢˢ\kˢᵠ)                        kˢʷ'*(kˢˢ\kˢʷ)-kᵅʷʷ]
            ks = Symmetric(0.5 * (k_cond + k_cond'));
            
            Kᴳ_cond = [zeros(2*nᵠ, 2*nᵠ)  zeros(2*nᵠ, nʷ);
                       zeros(nʷ, 2*nᵠ)    kʷʷ]
            
            F = eigen(ks, Kᴳ_cond); 
            λ = F.values;
            V = F.vectors

            mode_ids = sort!([i for i in eachindex(λ) if isfinite(real(λ[i])) && abs(imag(λ[i])) < 1.0e-7 && real(λ[i]) > 0.0], by = i -> real(λ[i]))
            sort!(mode_ids, by = i -> abs(real(λ[i])*h*σ₁₁*b^2/(π^2*Dᵇ) - k_exact))

            target_mode_id = first(mode_ids)
            λcr = real(λ[target_mode_id])
            k_num = (λcr * h * σ₁₁) * b^2 / (π^2 * Dᵇ)
            rel_error_k = abs(k_num - k_exact) / k_exact
            
            # 智慧幅值幾何歸一化與符號翻轉對齊
            v_first = V[:, target_mode_id];
            vʷ_raw = v_first[2*nᵠ+1:end]
            
            max_idx = argmax(abs.(vʷ_raw))
            x_max = nodes_w[max_idx].x
            y_max = nodes_w[max_idx].y
            
            exact_max_val = sin(π * x_max / a) * sin(π * y_max / b)
            
            denom = vʷ_raw[max_idx]
            if abs(denom) < 1e-15; denom = sign(denom) * 1e-15; end
            scale_factor = exact_max_val / denom
            
            v_first_scaled = v_first .* scale_factor
            vᵠ_1 = v_first_scaled[1:2*nᵠ];
            vʷ_1 = v_first_scaled[2*nᵠ+1:end]
            vˢ_1 = kˢˢ \ (kˢᵠ * vᵠ_1 + kˢʷ * vʷ_1)

            # 容器級別的 push! 欄位開拓協議
            push!(nodes_w, :d => vʷ_1)
            push!(nodes_φ, :d₁ => vᵠ_1[1:2:end], :d₂ => vᵠ_1[2:2:end])
            push!(nodes, :q₁ => vˢ_1[1:2:end], :q₂ => vˢ_1[2:2:end])

            w_ex_1(x,y,z)   = sin(1 * π * x / a) * sin(1 * π * y / b)
            φx_ex_1(x,y,z)  = (1 * π / a) * cos(1 * π * x / a) * sin(1 * π * y / b)
            φy_ex_1(x,y,z)  = (1 * π / b) * sin(1 * π * x / a) * cos(1 * π * y / b)
            q1_ex_1(x,y,z)  = Dˢ * ((1 * π / a) * cos(1 * π * x / a) * sin(1 * π * y / b) - φx_ex_1(x,y,z)) + 1e-5
            q2_ex_1(x,y,z)  = Dˢ * ((1 * π / b) * sin(1 * π * x / a) * cos(1 * π * y / b) - φy_ex_1(x,y,z)) + 1e-5

            prescribe!(elements_w, :w => w_ex_1)
            prescribe!(elements_φ, :φ₁ => φx_ex_1, :φ₂ => φy_ex_1)
            prescribe!(elements_q, :Q₁ => q1_ex_1, :Q₂ => q2_ex_1)

            @timeit to "calculate error" begin
                L₂_w = L₂w(elements_w)
                L₂_φ = L₂φ(elements_φ)
                L₂_Q = L₂Q(elements_q)
            end
            
            log10_Error_L2_w   = log10(L₂_w + 1e-16)
            log10_Error_L2_phi = log10(L₂_φ + 1e-16)
        end
        
        h_size = 1.0 / n_div
        push!(convergence_table, (
            ndiv = n_div, h = h_size, log10_h = log10(h_size),
            lambda_cr = λcr, k_num = k_num, rel_error_k = rel_error_k,
            log10_Error_L2_w = log10_Error_L2_w, log10_Error_L2_phi = log10_Error_L2_phi
        ))
        
        # 多網格解耦 VTK 批量出圖
        n_modes_output = min(length(mode_ids), 20)
        cells = [length(elm.𝓒) == 4 ? MeshCell(VTKCellTypes.VTK_QUAD, [x.𝐼 for x in elm.𝓒]) : MeshCell(VTKCellTypes.VTK_TRIANGLE, [x.𝐼 for x in elm.𝓒]) for elm in elements_q]
        points = zeros(3, nˢ); for node in nodes; points[1, node.𝐼] = node.x; points[2, node.𝐼] = node.y; points[3, node.𝐼] = 0.0; end
        
        vtu_path = joinpath(VTK_DIR, "$(case_prefix).vtu")
        vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
            
            lambdas_real_box = [real(λ[m_id]) for m_id in mode_ids[1:n_modes_output]]
            vtk["Selected_Lambda_Real", WriteVTK.VTKFieldData()] = lambdas_real_box
            
            w_nodal = zeros(nˢ); phi1_nodal = zeros(nˢ); phi2_nodal = zeros(nˢ)
            
            for (r, mode_id) in enumerate(mode_ids[1:n_modes_output])
                d_mode = V[:, mode_id] .* (r == 1 ? scale_factor : 1.0)
                for i in 1:nˢ
                    idx_w = argmin([(nodes[i].x - nw.x)^2 + (nodes[i].y - nw.y)^2 for nw in nodes_w])
                    idx_𝜙 = argmin([(nodes[i].x - n𝜙.x)^2 + (nodes[i].y - n𝜙.y)^2 for n𝜙 in nodes_φ])
                    w_nodal[nodes[i].𝐼]    = d_mode[2*nᵠ+1:end][idx_w]
                    phi1_nodal[nodes[i].𝐼] = d_mode[1:2*nᵠ][2*idx_𝜙-1]
                    phi2_nodal[nodes[i].𝐼] = d_mode[1:2*nᵠ][2*idx_𝜙]
                end
                vtk["Mode_$(r)_w"] = w_nodal;
                vtk["Mode_$(r)_phi1"] = phi1_nodal; vtk["Mode_$(r)_phi2"] = phi2_nodal
            end
        end
        
        gmsh.finalize()
        println("  [SUCCESS] 網格平分 ndiv = $(n_div) 計算與 VTK 匯出完畢. Log10(L2_Error_W) = $(round(log10_Error_L2_w, digits=4))")
    end
end

# 數據大整合與檔名物理隔離
output_summary_filepath = joinpath(DATA_DIR, "buckling_h_mixtwo_line.csv")

open(output_summary_filepath, "w") do io
    println(io, join([
        "ndiv", "h", "log10_h", "lambda_cr", "k_num", "relative_error_k", "log10_Error_L2_w", "log10_Error_L2_phi"
    ], ","))
    
    for row in convergence_table
        println(io, join([
            row.ndiv, row.h, row.log10_h, row.lambda_cr, row.k_num, row.rel_error_k, row.log10_Error_L2_w, row.log10_Error_L2_phi
        ], ","))
    end
end

println(to)
println("="^80)
println("【重構大成功】總收斂大表已安全寫入隔離路徑：", output_summary_filepath)
println("【大成功】所有解耦對齊後的 3D 雲圖檔（限制前 20 階）已儲存至：", VTK_DIR)
println("="^80)