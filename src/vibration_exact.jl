# ==============================================================================
#  2D Mindlin Plate Exact Spectral Analysis - Master All-In-One Loop
#  流派：純淨版理論解析解全譜對齊數據庫 (Exact Master Base) - 修正版
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
# 🌟 核心修正：補上原本遺漏的 GmshImport 函數引進，解決 UndefVarError 錯誤
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements

using LinearAlgebra, WriteVTK, Printf
import Gmsh: gmsh

# ==================== 1. 物理參數設定 ====================
a = 1.0; E = 1.0e6; ν = 0.3; b = 1.0; h_over_b = 0.001; h = h_over_b * b; ρ = 1.0
Dᵇ = E * h^3 / (12 * (1 - ν^2))

data_dir = normpath(joinpath(@__DIR__, "..", "date"))
if !ispath(data_dir); mkpath(data_dir); end

function exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
    ν_ex = 0.3; E_ex = Dᵇ * 12 * (1 - ν_ex^2) / h^3
    Dˢ_ex = (5/6) * E_ex * h / (2 * (1 + ν_ex))
    α_m = m * π / a; β_n = n * π / b; λ² = α_m^2 + β_n^2
    K_ex = zeros(3,3)
    K_ex[1,1] = Dˢ_ex * λ²; K_ex[1,2] = Dˢ_ex * α_m; K_ex[1,3] = Dˢ_ex * β_n
    K_ex[2,1] = Dˢ_ex * α_m; K_ex[2,2] = Dᵇ * α_m^2 + Dᵇ * ((1-ν_ex)/2) * β_n^2 + Dˢ_ex; K_ex[2,3] = Dᵇ * ((1+ν_ex)/2) * α_m * β_n
    K_ex[3,1] = Dˢ_ex * β_n; K_ex[3,2] = K_ex[2,3]; K_ex[3,3] = Dᵇ * β_n^2 + Dᵇ * ((1-ν_ex)/2) * α_m^2 + Dˢ_ex
    M_ex = zeros(3,3)
    M_ex[1,1] = ρ * h; M_ex[2,2] = ρ * h^3 / 12; M_ex[3,3] = ρ * h^3 / 12
    vals = eigvals(K_ex, M_ex)
    ω_exact = sqrt(minimum(real.(vals)))
    w_h_exact = ω_exact * a^2 * sqrt(ρ * h / Dᵇ)
    return w_h_exact
end

function generate_mn_pairs(max_modes)
    pairs = Tuple{Int,Int}[]
    max_mn = ceil(Int, sqrt(max_modes)) + 2
    for total in 2:(2*max_mn)
        for m_harm in 1:max_mn
            n_harm = total - m_harm
            if n_harm >= 1 && n_harm <= max_mn
                push!(pairs, (m_harm, n_harm))
            end
        end
    end
    sort!(pairs, by = p -> p[1]^2 + p[2]^2)
    return pairs[1:min(max_modes, length(pairs))]
end

ndiv_series = 9:25
master_results = []

# ==================== 2. 自動化網格大迴圈 ====================
println("🚀 開始生成理論解析解 (Exact) 跨網格對齊數據庫...")

for n_div in ndiv_series
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()
    
    mesh_file = normpath(joinpath(@__DIR__, "..", "msh", "st_q_$(n_div).msh"))
    if !isfile(mesh_file)
        gmsh.finalize()
        continue
    end
    
    gmsh.open(mesh_file)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    n_pts = length(nodes)
    
    # 依據當前網格規模對齊的理論模態組
    mn_pairs = generate_mn_pairs(n_pts)
    n_vtu_output = min(length(mn_pairs), 20) # 鎖定前 20 階

    for mode_idx in 1:n_vtu_output
        m_val, n_val = mn_pairs[mode_idx]
        w_h_ex = exact_omega_sq(m_val, n_val, a, b, Dᵇ, ρ, h)
        
        # 打包儲存
        push!(master_results, (n_div, mode_idx, m_val, n_val, w_h_ex))
    end
    
    gmsh.finalize()
end

# ==================== 3. 匯出單一匯總 CSV ====================
csv_master_path = joinpath(data_dir, "vibration_exact_master_all_mesh.csv")
open(csv_master_path, "w") do io
    println(io, join(["n_div", "mode_rank", "m", "n", "w_h_exact"], ","))
    for data in master_results
        n_div, rank, m, n, w_exact = data
        println(io, join([n_div, rank, m, n, @sprintf("%.4f", w_exact)], ","))
    end
end
println("📊 成功導出 Exact 總表至: $csv_master_path")