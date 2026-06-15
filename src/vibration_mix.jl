# ==============================================================================
#  2D Mindlin Plate Free Vibration Analysis - Fully Simply Supported (SSSS)
#  流派：Meshfree-FE Hybrid Mixed Method (無網格-有限元高階多場混合變分原理)
#  特性：
#    1. CSV 摘要報告：全資料完整填入 (包含所有特徵模態之 w, phi1, phi2 能量拆解)
#    2. VTK 可視化雲圖：保持精簡，僅導出前 20 階低頻核心振態
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

# ==================== 2. 主程式網格循環求解 ====================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

for mesh_file in mesh_files
    mesh_name = splitext(basename(mesh_file))[1]
    println("\n" * "="^70)
    println("混合離散動力分析：$(mesh_name) 網格 (SSSS 全簡支狀態)")
    println("="^70)

    case_prefix = "vibration_mix_$(mesh_name)"
    integrationOrder = 2

    # 2.1 物理場離散類型與空間半徑設定 (嚴格遵從 Inf-Sup 穩定性規範)
    type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
    type_Q = :tri3
    type_M = :(PiecewisePolynomial{:Linear2D})
    
    sʷ = 1.5
    sᵠ = 1.5

    m_match = match(r"(\d+)", mesh_name)
    ndiv = m_match === nothing ? 16 : parse(Int, m_match.captures[1]) - 1
    s_size = 1.0 / ndiv

    # 建立位移場 w 的粒子支援域
    gmsh.clear()
    gmsh.open(mesh_file)
    nodes_w = get𝑿ᵢ()
    nʷ = length(nodes_w)
    push!(nodes_w, :s₁ => sʷ * s_size * ones(nʷ), :s₂ => sʷ * s_size * ones(nʷ), :s₃ => sʷ * s_size * ones(nʷ))
    sp_w = RegularGrid(nodes_w.x, nodes_w.y, nodes_w.z, n=3, γ=5)

    # 建立轉角場 φ 的粒子支援域
    gmsh.clear()
    gmsh.open(mesh_file)
    nodes_φ = get𝑿ᵢ()
    nᵠ = length(nodes_φ)
    push!(nodes_φ, :s₁ => sᵠ * s_size * ones(nᵠ), :s₂ => sᵠ * s_size * ones(nᵠ), :s₃ => sᵠ * s_size * ones(nᵠ))
    sp_φ = RegularGrid(nodes_φ.x, nodes_φ.y, nodes_φ.z, n=3, γ=5)

    gmsh.clear()
    gmsh.open(mesh_file)
    nodes = get𝑿ᵢ()
    entities = getPhysicalGroups()
    nˢ = length(nodes)

    # 剛度分塊矩陣初始化
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

    # 2.3 組裝無網格動力一致質量矩陣 (保留完整分塊以利後續能量拆解)
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
        global k = -[kˢᵠ'*(kˢˢ\kˢᵠ)+kᵐᵠ'*(kᵐᵐ\kᵐᵠ)-kᵅᵠᵠ   kˢᵠ'*(kˢˢ\kˢʷ); 
                     kˢʷ'*(kˢˢ\kˢᵠ)                        kˢʷ'*(kˢˢ\kˢʷ)-kᵅʷʷ]
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
    
    # 圖形導出上限保持為前 20 階
    n_modes_output = min(length(mode_ids), 20)

    # 計算相對殘差
    global residuals = [
        relative_residual(ks, ms, real(ω²[mid]), V[:, mid])
        for mid in mode_ids
    ]

    println("成功解出特徵值！有效物理模態總數: $(length(mode_ids))")

    # ──────────────────────────────────────────────────────────
    # 2.5 🌟 全域量化核心：遍歷所有模態進行【能量拆解】與【理論 Navier 配對】
    # ──────────────────────────────────────────────────────────
    modal_energies = []
    mac_matches = []
    
    println("正在執行全頻譜場變量能量定量拆解...")
    for r in 1:length(mode_ids)
        mode_id = mode_ids[r]
        d_mode = V[:, mode_id]
        
        # 提取分塊特徵向量
        vᵠ = d_mode[1:2*nᵠ]
        vʷ = d_mode[2*nᵠ+1:end]
        
        # A. 位移場 w 攜帶的平移動能
        E_kinetic_w = vʷ' * mʷʷ * vʷ
        
        # B. 轉角場 φ₁ (x方向) 與 φ₂ (y方向) 攜帶的旋轉動能細分
        vᵠ₁ = vᵠ[1:2:end]
        vᵠ₂ = vᵠ[2:2:end]
        mᵠ₁ᵠ₁ = mᵠᵠ[1:2:end, 1:2:end]
        mᵠ₂ᵠ₂ = mᵠᵠ[2:2:end, 2:2:end]
        
        E_kinetic_𝜙1 = vᵠ₁' * mᵠ₁ᵠ₁ * vᵠ₁
        E_kinetic_𝜙2 = vᵠ₂' * mᵠ₂ᵠ₂ * vᵠ₂
        
        E_kinetic_total = E_kinetic_w + E_kinetic_𝜙1 + E_kinetic_𝜙2
        
        # 計算各自對當前特徵值的控制份額 (%)
        W_pct = (E_kinetic_w / E_kinetic_total) * 100
        Phi1_pct = (E_kinetic_𝜙1 / E_kinetic_total) * 100
        Phi2_pct = (E_kinetic_𝜙2 / E_kinetic_total) * 100
        push!(modal_energies, (W_pct, Phi1_pct, Phi2_pct))
        
        # 進行 Navier 理論解自動配對
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
    # 2.7 🌟 CSV 全數據摘要導出 (包含全部模態與 w, phi1, phi2 量化數值欄位)
    # ──────────────────────────────────────────────────────────
    csv_summary_path = joinpath(output_dir, "$(case_prefix)_modes.csv")
    
    open(csv_summary_path, "w") do io
        # 寫入包含能量佔比的豐富 CSV 標頭
        println(io, join([
            "mode_rank", "eigen_index", "omega_sq_real", "omega_real", 
            "Omega_FEM", "Omega_exact", "matched_m", "matched_n", 
            "W_energy_pct", "Phi1_energy_pct", "Phi2_energy_pct", "relative_residual"
        ], ","))
        
        # 🌟 滿載循環：直接寫入解出的所有合法特徵對
        for r in 1:length(mode_ids)
            mode_id = mode_ids[r]
            ω²_real = real(ω²[mode_id])
            ω_real  = sqrt(ω²_real)
            Ω_FEM   = ω_real * a^2 * sqrt(ρ * h / Dᵇ)
            res_val = residuals[r]
            
            # 讀取對應的配對與能量數值
            m, n = mac_matches[r]
            _, Ω_exact = exact_omega_sq(m, n, a, b, Dᵇ, ρ, h)
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