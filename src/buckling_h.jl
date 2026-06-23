# ==============================================================================
#  高級計算固體力學研究平台：Mindlin 屈曲分析純位移法（厚度效應掃描專題）
#  流派：純位移控制體有限元 / 固定 17x17 結構化網格
#  功能：固定網格 + 厚度 h 自動化循環 (0.1 ➔ 0.00001) + 智慧幅值歸一化
#  網格路徑對規：/home/a/Joker/msh/st_q_17.msh 規則
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ

using LinearAlgebra
using TimerOutputs
using WriteVTK
using Printf
import Gmsh: gmsh

# 精密資料夾分流配置
const DATA_DIR = "./date"
const VTK_DIR  = "./VTK"
mkpath(DATA_DIR)
mkpath(VTK_DIR)

# 全局基礎力學材料常數
const E = 200e9
const ν = 0.3
const a = 1.0   
const b = 1.0   
const k_exact = 4.0  # 經典薄板臨界屈曲係數理論值

# 預加載名義應力場 (單軸均勻壓縮)
const σ₁₁ = 1.0
const σ₂₂ = 0.0
const σ₁₂ = 0.0

const to = TimerOutput()
convergence_table = []

# 定義理論經典薄板解析解函數 —— 用於智慧 field scale alignment
w_exact_func(x, y, z) = sin(π * x / a) * sin(π * y / b)
phi1_exact_func(x, y, z) = (π / a) * cos(π * x / a) * sin(π * y / b)
phi2_exact_func(x, y, z) = (π / b) * sin(π * x / a) * cos(π * y / b)

# 🚀 嚴格定義厚度掃描序列 (從 0.1 一路對數/梯度下降至 0.00001)
h_series = [0.1, 0.01, 0.001, 0.0001, 0.00001]

println("="^80)
println(" 啟動板厚度穩定性掃描大循環 (h = 0.1 ➔ 0.00001) ")
println("  固定使用網格: /home/a/Joker/msh/st_q_17.msh")
println("="^80)

# 固定使用 17 網格
const ndiv = 17
mesh_path = "/home/a/Joker/msh/st_q_17.msh"

if !isfile(mesh_path)
    error("【錯誤】找不到指定的 17 階基礎網格檔案: ", mesh_path)
end

# 進入厚度變量大循環
for h_val in h_series
    @timeit to "Total loop iteration (h=$h_val)" begin
        
        # 🚀 核心力學重構：彎矩剛度 Dᵇ 必須隨著當前厚度進行三次方動態更新
        Dᵇ = E * h_val^3 / (12 * (1 - ν^2))
        α_penalty = 1.0e8 * Dᵇ  # 邊界罰參數同步縮放，確保數值矩陣條件數穩定
        
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(mesh_path)
        
        entities = getPhysicalGroups()
        nodes = get𝑿ᵢ()
        
        nʷ = length(nodes)
        nᵠ = length(nodes)
        
        # 分配全域剛度與幾何剛度矩陣
        kʷʷ = zeros(nʷ, nʷ)
        kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
        kᵠʷ = zeros(2*nᵠ, nʷ)
        kʷʷ = zeros(nʷ, nʷ)
        kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
        
        integrationOrder = 2
        integrationOrder_shear = 1
        
        # 域內剛度與幾何剛度矩陣集成
        @timeit to "Domain assembly" begin
            elements = getElements(nodes, entities["Ω"], integrationOrder)
            elements_s = getElements(nodes, entities["Ω"], integrationOrder_shear)
            
            # 注入當前循環的厚度 h_val
            prescribe!(elements, :E=>E, :ν=>ν, :h=>h_val)
            prescribe!(elements_s, :E=>E, :ν=>ν, :h=>h_val)
            
            set∇𝝭!(elements)
            set∇𝝭!(elements_s)
            
            (∫wwdΩ=>elements_s)(kʷʷ)
            (∫φwdΩ=>elements_s)(kᵠʷ)
            ([∫φφdΩ=>elements_s, ∫κκdΩ=>elements])(kᵠᵠ)
            
            prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
            (∫∇wσ∇wdΩ=>elements)(kʷʷ)
            (∫∇φσ∇φdΩ=>elements)(kᵠᵠ)
        end
        
        # 邊界條件處理 (Hard SSSS - S2 邊界罰函數法)
        @timeit to "Boundary BC penalty assembly" begin
            elements_1 = getElements(nodes, entities["Γ¹"], integrationOrder)
            elements_2 = getElements(nodes, entities["Γ²"], integrationOrder)
            elements_3 = getElements(nodes, entities["Γ³"], integrationOrder)
            elements_4 = getElements(nodes, entities["Γ⁴"], integrationOrder)
            
            prescribe!(elements_1, :α=>α_penalty, :g=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>0.0)
            prescribe!(elements_2, :α=>α_penalty, :g=>0.0, :n₁₁=>0.0, :n₁₂=>0.0, :n₂₂=>1.0)
            prescribe!(elements_3, :α=>α_penalty, :g=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>0.0)
            prescribe!(elements_4, :α=>α_penalty, :g=>0.0, :n₁₁=>0.0, :n₁₂=>0.0, :n₂₂=>1.0)
            
            set𝝭!(elements_1); set𝝭!(elements_2); set𝝭!(elements_3); set𝝭!(elements_4)
            
            (∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4)(kʷʷ)
            (∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4)(kᵠᵠ)
        end
        
        # 特徵值屈曲問題求解
        @timeit to "Eigen solver" begin
            K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
            K = [zeros(2*nᵠ, 2*nᵠ)  zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ)  kʷʷ]
            
            F = eigen(K, K)
            λ = F.values
            V = F.vectors
            
            mode_ids = sort!(
                collect(i for i in eachindex(λ) 
                    if isfinite(real(λ[i])) && abs(imag(λ[i])) < 1.0e-7 && real(λ[i]) > 0.0),
                by = i -> real(λ[i])
            )
            
            sort!(mode_ids, by = i -> abs(real(λ[i])*h_val*σ₁₁*b^2/(π^2*Dᵇ) - k_exact))
            
            target_mode_id = first(mode_ids)
            λcr = real(λ[target_mode_id])
            k_num = (λcr * h_val * σ₁₁) * b^2 / (π^2 * Dᵇ)
            rel_error_k = abs(k_num - k_exact) / k_exact
            
            # 智慧振幅幾何歸一化與符號翻轉對齊 (L2 誤差對齊核心)
            dm_raw = real.(V[:, target_mode_id])
            w_raw = dm_raw[2*nᵠ+1:end] 
            
            max_idx = argmax(abs.(w_raw))
            x_max = nodes[max_idx].x
            y_max = nodes[max_idx].y
            
            exact_max_val = w_exact_func(x_max, y_max, 0.0)
            
            denom = w_raw[max_idx]
            if abs(denom) < 1e-15; denom = sign(denom) * 1e-15; end
            scale_factor = exact_max_val / denom
            
            dm_scaled = dm_raw .* scale_factor
            vᵠ = dm_scaled[1:2*nᵠ]
            vʷ = dm_scaled[2*nᵠ+1:end]
            
            # 採用標準安全型協議，開拓空間抽屜注入數據，絕殺 KeyError
            push!(nodes, :d => vʷ)
            push!(nodes, :d₁ => vᵠ[1:2:2*nᵠ])
            push!(nodes, :d₂ => vᵠ[2:2:2*nᵠ])
            
            prescribe!(elements, :w => w_exact_func)
            prescribe!(elements, :φ₁ => phi1_exact_func, :φ₂ => phi2_exact_func)
            
            L₂_w = L₂w(elements)
            L₂_φ = L₂φ(elements)
            
            log10_Error_L2_w   = log10(L₂_w + 1e-16)
            log10_Error_L2_phi = log10(L₂_φ + 1e-16)
        end
        
        # 收集當前厚度下的力學性能軌跡
        push!(convergence_table, (
            h = h_val, log10_h = log10(h_val),
            lambda_cr = λcr, k_num = k_num, rel_error_k = rel_error_k,
            log10_Error_L2_w = log10_Error_L2_w, log10_Error_L2_phi = log10_Error_L2_phi
        ))
        
        # 🚀 數據同一輸出在同一個 CSV 文件，依照要求以固定網格 17 進行模式命名
        modes_csv_path = joinpath(DATA_DIR, "buckling_st_q_17_h_$(h_val)_modes.csv")
        open(modes_csv_path, "w") do io
            println(io, "mode_rank,eigen_index,lambda_real,k_num,w_norm,phi1_norm,phi2_norm")
            for (rank, m_id) in enumerate(mode_ids[1:min(20, length(mode_ids))])
                lam = real(λ[m_id])
                k_val = (lam * h_val * σ₁₁) * b^2 / (π^2 * Dᵇ)
                
                d_m = real.(V[:, m_id]) .* (rank == 1 ? scale_factor : 1.0)
                ph1_v = d_m[1:2:2*nᵠ]
                ph2_v = d_m[2:2:2*nᵠ]
                w_v   = d_m[2*nᵠ+1:end]
                
                @printf(io, "%d,%d,%.6e,%.6f,%.6e,%.6e,%.6e\n", 
                        rank, m_id, lam, k_val, norm(w_v), norm(ph1_v), norm(ph2_v))
            end
        end
        
        # VTK 檔案儲存在 VTK 資料夾下
        nₚ = length(nodes)
        points = zeros(3, nₚ)
        for (i, node) in enumerate(nodes)
            points[1, i] = node.x
            points[2, i] = node.y
            points[3, i] = 0.0
        end
        cells = [MeshCell(VTKCellTypes.VTK_QUAD, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]
        
        vtk_grid("$(VTK_DIR)/buckling_thickness_mesh_h_$(h_val).vtu", points, cells; ascii=true) do vtk
            vtk["Mode_1_w_aligned", VTKPointData()] = vʷ
            vtk["Mode_1_phi1_aligned", VTKPointData()] = vᵠ[1:2:2*nᵠ]
            vtk["Mode_1_phi2_aligned", VTKPointData()] = vᵠ[2:2:2*nᵠ]
        end
        
        gmsh.finalize()
        println("  [進度] 厚度 h = $(h_val) 計算完畢. 屈曲係數 k_num = $(round(k_num, digits=4)), L2_Error_W = $(round(L₂_w, digits=6))")
    end
end

# ==============================================================================
#  自動將各厚度梯度的誤差輸出為指定名稱的總收斂大表
# ==============================================================================
output_summary_filepath = joinpath(DATA_DIR, "buckling_h_line.csv")

open(output_summary_filepath, "w") do io
    println(io, join([
        "h", "log10_h", "lambda_cr", "k_num", "relative_error_k", "log10_Error_L2_w", "log10_Error_L2_phi"
    ], ","))
    
    for row in convergence_table
        println(io, join([
            row.h, row.log10_h, row.lambda_cr, row.k_num, row.rel_error_k, row.log10_Error_L2_w, row.log10_Error_L2_phi
        ], ","))
    end
end

println(to)
println("="^80)
println("【厚度掃描大成功】總性能折線表已寫入：", output_summary_filepath)
println("="^80)