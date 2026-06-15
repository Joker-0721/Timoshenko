# ==============================================================================
#  2D Mindlin Plate Free Vibration Analysis - Fully Simply Supported (SSSS)
#  流派：Meshfree-FE Hybrid Mixed Method (無網格-有限元高階多場混合變分原理)
#  特性：
#    1. CSV 摘要報告：全資料完整填入 (包含所有特徵模態之 w, phi1, phi2 能量拆解)
#    2. VTK 可視化雲圖：保持精簡，僅導出前 20 階低頻核心振態
#    3. 自動掃描：自動呼叫雷達探測器，動態尋找最優高斯階數與無網格半徑
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, 
                                     ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ, L₂Q,
                                     ∫ρhwwdΩ, ∫ρIφφdΩ  

using LinearAlgebra
using TimerOutputs, WriteVTK, Printf
import Gmsh: gmsh

# ==================== 1. 材料與幾何參數設定 ====================
E = 1.0e6
ν = 0.3
a = 1.0
b = 1.0
h_over_b = 0.001          # 厚薄比 (0.001 代表超薄板極限)
h = h_over_b * b
ρ = 1.0

# 計算板彎曲與剪切剛度
Dᵇ = E * h^3 / (12 * (1 - ν^2))
Dˢ = 5 / 6 * E * h / (2 * (1 + ν))

# 四邊簡支 SSSS 罰參數設定 (放開轉角約束)
αʷ = 1.0e8 * Dᵇ     
αᵠ = 0.0            

eigen_imag_tol = 1.0e-7
omega_sq_tol = 1.0e-12

# 要分析的網格檔案
mesh_files = [
    normpath(joinpath(@__DIR__, "..", "msh", "st_q_17.msh")),
]

output_dir = normpath(joinpath(@__DIR__, "..", "vtk"))
mkpath(output_dir)

# 相對殘差評估核心函數
function relative_residual(K, M, ω²ᵢ, vᵢ)
    Kv = K * vᵢ
    Mv = M * vᵢ
    residual = Kv - ω²ᵢ * Mv
    denominator = max(norm(Kv), abs(ω²ᵢ) * norm(Mv), eps(Float64))
    return norm(residual) / denominator
end

# 薄板 Navier 頻率理論解析解公式
function exact_omega_sq(m, n, a, b, D, ρ, h)
    Ω_exact = π^2 * (m^2 + n^2)
    ω_exact = Ω_exact * sqrt(D / (ρ * h)) / a^2
    return ω_exact, Ω_exact
end

# VTK 單元構建輔助函數
function vtk_cell(elm)
    node_ids = [xᵢ.𝐼 for xᵢ in elm.𝓒]
    if length(node_ids) == 4
        return MeshCell(VTKCellTypes.VTK_QUAD, node_ids)
    else
        return MeshCell(VTKCellTypes.VTK_TRIANGLE, node_ids)
    end
end

const to = TimerOutput()

using Statistics # 引入用於計算平均值的套件

"""
    auto_find_optimal_parameters(...)
    透過二維網格掃描，自動尋找讓前3 階動力殘差最小、且矩陣最穩定的 (integrationOrder, s) 黃金組合
"""
function auto_find_optimal_parameters(mesh_file, entities, nodes, nˢ, Dᵇ, ρ, h, E, ν, type_w, type_φ, type_Q, type_M, s_size, αʷ, αᵠ)
    println("\n🔍 [雷達掃描啟動] 開始全自動探測最優數值積分階數與無網格半徑...")
    
    gauss_orders = [2, 3, 4]             
    s_radius_range = 1.2:0.1:2.0         
    
    best_order = 2
    best_s = 1.5
    min_avg_residual = Inf
    best_condition = Inf

    println(@sprintf("  %-18s %-8s %-15s %-12s %-10s", "探測坐標(Order, s)", "狀態", "平均相對殘差", "剛度條件數", "評級"))
    println("  " * "-" ^ 70) # 🌟 修正點：* 70 變更為 ^ 70

    for order in gauss_orders
        for s_try in s_radius_range
            try
                gmsh.clear()
                gmsh.open(mesh_file)
                nodes_w = get𝑿ᵢ()
                nʷ = length(nodes_w)
                push!(nodes_w, :s₁ => s_try * s_size * ones(nʷ), :s₂ => s_try * s_size * ones(nʷ), :s₃ => s_try * s_size * ones(nʷ))
                sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

                gmsh.clear()
                gmsh.open(mesh_file)
                nodes_φ = get𝑿ᵢ()
                nᵠ = length(nodes_φ)
                push!(nodes_φ, :s₁ => s_try * s_size * ones(nᵠ), :s₂ => s_try * s_size * ones(nᵠ), :s₃ => s_try * s_size * ones(nᵠ))
                sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

                kˢˢ = zeros(2 * nˢ, 2 * nˢ); kˢʷ = zeros(2 * nˢ, nʷ); kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
                kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ); kᵅʷʷ = zeros(nʷ, nʷ)
                mʷʷ = zeros(nʷ, nʷ); mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

                elements_q = getElements(nodes, entities["Ω"], order)
                prescribe!(elements_q, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_q)
                elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), order, sp_w)
                prescribe!(elements_w, :E => E, :ν => ν, :h => h); set𝝭!(elements_w)
                elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), order, sp_w, normal=true); set𝝭!(elements_w_Γ)
                elements_q_Γ = getElements(nodes, entities["Γ"], order, normal=true); set𝝭!(elements_q_Γ)
                (∫QQdΩ => elements_q)(kˢˢ)
                ([∫∇QwdΩ => (elements_q, elements_w), ∫QwdΓ => (elements_q_Γ, elements_w_Γ)])(kˢʷ)

                elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), order)
                prescribe!(elements_m, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_m)
                elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), order, sp_φ)
                prescribe!(elements_φ, :E => E, :ν => ν, :h => h); set∇𝝭!(elements_φ)
                elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), order, sp_φ, normal=true); set𝝭!(elements_φ_Γ)
                elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), order); set𝝭!(elements_m_Γ)

                nₑ = length(elements_q); nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[])); nᵐ = nₑ * nᵖ
                kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ); kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)
                (∫MMdΩ => elements_m)(kᵐᵐ); ([∫∇MφdΩ => (elements_m, elements_φ), ∫MφdΓ => (elements_m_Γ, elements_φ_Γ)])(kᵐᵠ)
                (∫QφdΩ => (elements_q, elements_φ))(kˢᵠ)

                boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
                elements_w_b = [getElements(nodes_w, entities[name], eval(type_w), order, sp_w, normal=true) for name in boundary_names]
                for el in elements_w_b; prescribe!(el, :α => αʷ, :g => (x,y,z)->0.0); set𝝭!(el); end
                elements_φ_b = [getElements(nodes_φ, entities[name], eval(type_φ), order, sp_φ, normal=true) for name in boundary_names]
                for el in elements_φ_b; prescribe!(el, :α => αᵠ, :g₁ => (x,y,z)->0.0, :g₂ => (x,y,z)->0.0, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0); set𝝭!(el); end
                
                dummy_fʷ = zeros(nʷ); dummy_fᵠ = zeros(2*nᵠ)
                (∫αwwdΓ => elements_w_b[1] ∪ elements_w_b[2] ∪ elements_w_b[3] ∪ elements_w_b[4])(kᵅʷʷ, dummy_fʷ)
                (∫αφφdΓ => elements_φ_b[1] ∪ elements_φ_b[2] ∪ elements_φ_b[3] ∪ elements_φ_b[4])(kᵅᵠᵠ, dummy_fᵠ)

                prescribe!(elements_w, :ρ => ρ, :h => h); (∫ρhwwdΩ => elements_w)(mʷʷ)
                prescribe!(elements_φ, :ρ => ρ, :h => h); (∫ρIφφdΩ => elements_φ)(mᵠᵠ)
                M_tmp = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]

                K_cond = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ   kˢᵠ'*(kˢˢ\kˢʷ); kˢʷ'*(kˢˢ\kˢᵠ)  kˢʷ'*(kˢˢ\kˢʷ)-kᵅʷʷ]
                ks_s = Symmetric(0.5 * (K_cond + K_cond'))
                ms_s = Symmetric(0.5 * (M_tmp + M_tmp'))
                
                cond_val = cond(ks_s)
                if cond_val > 1.0e12
                    @sprintf("  (%d, %.1f)          %-8s %-15s %-12.2e %-10s", order, s_try, "代數病態", "N/A", cond_val, "❌ 淘汰") |> println
                    continue
                end

                F = eigen(ks_s, ms_s)
                mode_filter = [i for i in eachindex(F.values) if real(F.values[i]) > 1e-12 && abs(imag(F.values[i])) < 1e-7]
                if length(mode_filter) < 3
                    @sprintf("  (%d, %.1f)          %-8s %-15s %-12.2e %-10s", order, s_try, "秩虧噪聲", "N/A", cond_val, "❌ 淘汰") |> println
                    continue
                end
                
                avg_res = mean([relative_residual(ks_s, ms_s, real(F.values[mid]), F.vectors[:, mid]) for mid in mode_filter[1:3]])

                rating = "   優良"
                if avg_res < min_avg_residual
                    min_avg_residual = avg_res
                    best_order = order
                    best_s = s_try
                    best_condition = cond_val
                    rating = "🌟 目前最優"
                end
                @sprintf("  (%d, %.1f)          %-8s %-15.4e %-12.2e %-10s", order, s_try, "計算成功", avg_res, cond_val, rating) |> println
            catch e
                @sprintf("  (%d, %.1f)          %-8s %-15s %-12s %-10s", order, s_try, "矩陣奇異", "N/A", "N/A", "❌ 崩潰") |> println
                continue
            end
        end
    end
    println("  " * "-" ^ 70) # 🌟 修正點：* 70 變更為 ^ 70
    println("🎯 [雷達探測成功] 最佳數值參數組合已鎖定：")
    println("  -> 最佳高斯積分階數 (integrationOrder) = ", best_order)
    println("  -> 最佳無網格影響域半徑 (s) = ", best_s)
    println("  -> 此時系統最優平均動力相對殘差 = ", min_avg_residual)
    println("="^70 * "\n")
    
    return best_order, best_s
end

# ==================== 2. 主程式網格循環求解 ====================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

for mesh_file in mesh_files
    mesh_name = splitext(basename(mesh_file))[1]
    println("\n" * "="^70)
    println("混合離散動力分析：$(mesh_name) 網格 (SSSS 全簡支狀態)")
    println("="^70)

    case_prefix = "vibration_mix_$(mesh_name)"
    
    type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_Q = :tri3
    type_M = :(PiecewisePolynomial{:Linear2D})

    m_match = match(r"(\d+)", mesh_name)
    ndiv = m_match === nothing ? 16 : parse(Int, m_match.captures[1]) - 1
    s_size = 1.0 / ndiv

    gmsh.clear()
    gmsh.open(mesh_file)
    nodes = get𝑿ᵢ()
    entities = getPhysicalGroups()
    nˢ = length(nodes)

    # 🌟【核心修正通道】使用 Base.invokelatest 呼叫探測雷達，百分之百消除 Julia 1.12 的世界年齡警告與報錯
    optimal_order, optimal_s = Base.invokelatest(
        auto_find_optimal_parameters,
        mesh_file, entities, nodes, nˢ, Dᵇ, ρ, h, E, ν, 
        type_w, type_φ, type_Q, type_M, s_size, αʷ, αᵠ
    )

    integrationOrder = optimal_order
    sʷ = optimal_s
    sᵠ = optimal_s

    gmsh.clear()
    gmsh.open(mesh_file)
    nodes_w = get𝑿ᵢ()
    nʷ = length(nodes_w)
    push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ), :s₂ => sʷ * s_size * ones(nʷ), :s₃ => sʷ * s_size * ones(nʷ))
    sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

    gmsh.clear()
    gmsh.open(mesh_file)
    nodes_φ = get𝑿ᵢ()
    nᵠ = length(nodes_φ)
    push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ), :s₂ => sᵠ * s_size * ones(nᵠ), :s₃ => sᵠ * s_size * ones(nᵠ))
    sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

    kˢˢ = zeros(2 * nˢ, 2 * nˢ)
    kˢʷ = zeros(2 * nˢ, nʷ)
    kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
    kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
    kᵅʷʷ = zeros(nʷ, nʷ)

    # A. 剪切項弱形式組裝
    @timeit to "assemble shear items" begin
        elements_q = getElements(nodes, entities["Ω"], integrationOrder)
        prescribe!(elements_q, :E => E, :ν => ν, :h => h)
        set∇𝝭!(elements_q)

        elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
        prescribe!(elements_w, :E => E, :ν => ν, :h => h)
        set𝝭!(elements_w)

        elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
        set𝝭!(elements_w_Γ)

        elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
        set𝝭!(elements_q_Γ)

        𝑎ˢˢ = ∫QQdΩ => elements_q
        𝑎ˢʷ = [∫∇QwdΩ => (elements_q, elements_w), ∫QwdΓ => (elements_q_Γ, elements_w_Γ)]
        
        𝑎ˢˢ(kˢˢ)
        𝑎ˢʷ(kˢʷ)
    end

    # B. 彎矩項弱形式組裝
    nₑ = length(elements_q)
    nᵖ = ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))
    nᵐ = nₑ * nᵖ
    kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ)
    kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)

    @timeit to "assemble moment items" begin
        elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
        prescribe!(elements_m, :E => E, :ν => ν, :h => h)
        set∇𝝭!(elements_m)

        elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
        prescribe!(elements_φ, :E => E, :ν => ν, :h => h)
        set∇𝝭!(elements_φ)

        elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        set𝝭!(elements_φ_Γ)

        elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
        set𝝭!(elements_m_Γ)

        𝑎ᵐᵐ = ∫MMdΩ => elements_m
        𝑎ᵐᵠ = [∫∇MφdΩ => (elements_m, elements_φ), ∫MφdΓ => (elements_m_Γ, elements_φ_Γ)]
        𝑎ˢᵠ = ∫QφdΩ => (elements_q, elements_φ)

        𝑎ᵐᵐ(kᵐᵐ)
        𝑎ᵐᵠ(kᵐᵠ)
        𝑎ˢᵠ(kˢᵠ)
    end

    # C. 施加四邊邊界自適應罰剛度
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

    # 2.3 組裝無網格動力一致質量矩陣
    mʷʷ = zeros(nʷ, nʷ)
    mᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)

    @timeit to "assemble mass matrix" begin
        prescribe!(elements_w, :ρ => ρ, :h => h)
        (∫ρhwwdΩ => elements_w)(mʷʷ)

        prescribe!(elements_φ, :ρ => ρ, :h => h) 
        (∫ρIφφdΩ => elements_φ)(mᵠᵠ)

        global M_total = [mᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) mʷʷ]
    end

    # 2.4 執行靜態凝聚與廣義特徵值求解
    @timeit to "static condensation" begin
        global k = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ   kˢᵠ'*(kˢˢ\kˢʷ); kˢʷ'*(kˢˢ\kˢᵠ)                        kˢʷ'*(kˢˢ\kˢʷ)-kᵅʷʷ]
    end

    @timeit to "solve eigenvalue" begin
        ks = Symmetric(0.5 * (k + k'))
        ms = Symmetric(0.5 * (M_total + M_total'))
        F = eigen(ks, ms)
        ω² = F.values
        V = F.vectors
    end

    # 篩選合法物理振動模態
    mode_ids = sort!(
        collect(i for i in eachindex(ω²) if real(ω²[i]) > omega_sq_tol && abs(imag(ω²[i])) < eigen_imag_tol && isfinite(real(ω²[i]))),
        by = i -> real(ω²[i]),
    )
    
    n_modes_output = min(length(mode_ids), 20)

    # 計算相對殘差
    global residuals = [
        relative_residual(ks, ms, real(ω²[mid]), V[:, mid])
        for mid in mode_ids
    ]

    println("成功解出特徵值！有效物理模態總數: $(length(mode_ids))")

    # ──────────────────────────────────────────────────────────
    # 2.5 全域量化核心：遍歷所有模態進行【能量拆解】與【理論 Navier 配對】
    # ──────────────────────────────────────────────────────────
    modal_energies = []
    mac_matches = []
    
    println("正在執行全頻譜場變量能量定量拆解...")
    for r in 1:length(mode_ids)
        mode_id = mode_ids[r]
        d_mode = V[:, mode_id]
        
        vᵠ = d_mode[1:2*nᵠ]
        vʷ = d_mode[2*nᵠ+1:end]
        
        E_kinetic_w = vʷ' * mʷʷ * vʷ
        
        vᵠ₁ = vᵠ[1:2:end]
        vᵠ₂ = vᵠ[2:2:end]
        
        mᵠ₁ᵠ₁ = mᵠᵠ[1:2:end, 1:2:end]
        mᵠ₂ᵠ₂ = mᵠᵠ[2:2:end, 2:2:end]
        
        E_kinetic_𝜙1 = vᵠ₁' * mᵠ₁ᵠ₁ * vᵠ₁
        E_kinetic_𝜙2 = vᵠ₂' * mᵠ₂ᵠ₂ * vᵠ₂
        
        E_kinetic_total = E_kinetic_w + E_kinetic_𝜙1 + E_kinetic_𝜙2
        
        W_pct = (E_kinetic_w / E_kinetic_total) * 100
        Phi1_pct = (E_kinetic_𝜙1 / E_kinetic_total) * 100
        Phi2_pct = (E_kinetic_𝜙2 / E_kinetic_total) * 100
        push!(modal_energies, (W_pct, Phi1_pct, Phi2_pct))
        
        ω²ᵢ = real(ω²[mode_id])
        ωᵢ = sqrt(ω²ᵢ)
        Ω_FEM = ωᵢ * a^2 * sqrt(ρ * h / Dᵇ)
        
        best_err = Inf
        best_mn = (0,0)
        for m_harm in 1:15 
            for n_harm in 1:15
                _, Ω_ex = exact_omega_sq(m_harm, n_harm, a, b, Dᵇ, ρ, h)
                err = abs(Ω_FEM - Ω_ex) / Ω_ex * 100
                if err < best_err
                    best_err = err
                    best_mn = (m_harm, n_harm)
                end
            end
        end
        push!(mac_matches, (best_mn[1], best_mn[2]))
    end

    # ──────────────────────────────────────────────────────────
    # 2.6 圖形場 VTU 導出 (保持僅輸出前 20 階，維持輕量化)
    # ──────────────────────────────────────────────────────────
    elements_Ω = getElements(nodes, entities["Ω"], 1)
    cells = [vtk_cell(elm) for elm in elements_Ω]
    points = zeros(3, nˢ)
    for node in nodes
        points[1, node.𝐼] = node.x
        points[2, node.𝐼] = node.y
        points[3, node.𝐼] = 0.0
    end

    vtu_path = joinpath(output_dir, "$(case_prefix).vtu")

    vtk_grid(vtu_path, points, cells; ascii=true, append=false, compress=false) do vtk
        ω²_selected = [real(ω²[mid]) for mid in mode_ids[1:n_modes_output]]
        Ω_selected = [sqrt(w2) * a^2 * sqrt(ρ * h / Dᵇ) for w2 in ω²_selected]
        
        vtk["Omega_dimensionless", WriteVTK.VTKFieldData()] = Ω_selected
        vtk["color_range_min", WriteVTK.VTKFieldData()] = [-1.0]
        vtk["color_range_max", WriteVTK.VTKFieldData()] = [1.0]

        for (mode_rank, mode_id) in enumerate(mode_ids[1:n_modes_output])
            d_mode = V[:, mode_id]
            dᵠ_mode = d_mode[1:2*nᵠ]
            dʷ_mode = d_mode[2*nᵠ+1:end]

            w_nodal = zeros(nˢ)
            phi1_nodal = zeros(nˢ)
            phi2_nodal = zeros(nˢ)
            
            for i in 1:nˢ
                w_nodal[i]    = dʷ_mode[i]
                phi1_nodal[i] = dᵠ_mode[2*i-1]
                phi2_nodal[i] = dᵠ_mode[2*i]
            end

            vtk["Mode_$(mode_rank)_w"]    = w_nodal
            vtk["Mode_$(mode_rank)_phi1"] = phi1_nodal
            vtk["Mode_$(mode_rank)_phi2"] = phi2_nodal
        end
    end
    println("  [VTK] 前 $n_modes_output 階低頻幾何圖形場成功導出: ", vtu_path)

    # ──────────────────────────────────────────────────────────
    # 2.7 CSV 全數據摘要導出 (包含全部模態與 w, phi1, phi2 量化數值欄位)
    # ──────────────────────────────────────────────────────────
    csv_summary_path = joinpath(output_dir, "$(case_prefix)_modes.csv")
    
    open(csv_summary_path, "w") do io
        println(io, join([
            "mode_rank", "eigen_index", "omega_sq_real", "omega_real",
            "Omega_FEM", "Omega_exact", "matched_m", "matched_n",
            "W_energy_pct", "Phi1_energy_pct", "Phi2_energy_pct", "relative_residual"
        ], ","))
        
        for r in 1:length(mode_ids)
            mode_id = mode_ids[r]
            ω²_real = real(ω²[mode_id])
            ω_real  = sqrt(ω²_real)
            Ω_FEM   = ω_real * a^2 * sqrt(ρ * h / Dᵇ)
            res_val = residuals[r]
            
            m, n = mac_matches[r]
            _, Ω_exact = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h) # 🌟 修正點：補齊為標準 7 參數
            W_pct, Phi1_pct, Phi2_pct = modal_energies[r]

            println(io, join([
                r, mode_id,
                @sprintf("%.6e", ω²_real), @sprintf("%.6e", ω_real),
                @sprintf("%.4f", Ω_FEM), @sprintf("%.4f", Ω_exact),
                m, n,
                @sprintf("%.4f", W_pct), @sprintf("%.4f", Phi1_pct), @sprintf("%.4f", Phi2_pct),
                @sprintf("%.2e", res_val)
            ], ","))
        end
    end
    println("  [CSV] 全頻譜全資料數據摘要（含變量能量解析）成功導出: ", csv_summary_path)

    # 終端精簡列印前 10 階，便於即時預覽
    println("\n結果摘要 - 前 10 階能量控制權分流預覽:")
    println("  " * @sprintf("%-6s %-6s %-10s %-10s %-10s %-12s %-12s", "mode", "(m,n)", "Ω_FEM", "Ω_exact", "W_eng(%)", "φ1_eng(%)", "φ2_eng(%)"))
    for r in 1:min(n_modes_output, 10)
        m, n = mac_matches[r]
        ω_real = sqrt(real(ω²[mode_ids[r]]))
        Ω_FEM = ω_real * a^2 * sqrt(ρ * h / Dᵇ)
        _, Ω_exact = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
        W_pct, Phi1_pct, Phi2_pct = modal_energies[r]
        println("  " * @sprintf("%-6d (%-d,%-d)   %-10.4f %-10.4f %-10.2f %-12.2f %-12.2f", r, m, n, Ω_FEM, Ω_exact, W_pct, Phi1_pct, Phi2_pct))
    end

end

gmsh.finalize()
println("\n" * "="^70)
println("全譜數據滿載分析與圖形分流成功結束！")
println("="^70)
println(to)