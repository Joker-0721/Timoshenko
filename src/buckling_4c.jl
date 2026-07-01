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
using Printf
import Gmsh: gmsh

# 全局數據庫路徑規範
const DATA_DIR = "./data"
const UNIFIED_CSV = joinpath(DATA_DIR, "buckling_FEM_cccc.csv")   # 輸出檔名可改為 CCCC
mkpath(DATA_DIR)

# 全局力學材料常數
const E = 200e9
const ν = 0.3
const h = 1e-2  
const a = 1.0   
const b = 1.0   
const Dᵇ = E * h^3 / (12 * (1 - ν^2))

const α_base = 1.0e8 * Dᵇ 
const σ₁₁ = 1.0
const σ₂₂ = 0.0
const σ₁₂ = 0.0

# 四邊固支（CCCC）第一模態擬合參考值 (Reddy, §7.5)
const k_exact_modes = [10.31, 23.92, 23.92, 39.57, 50.80, 50.80]

const integrationOrder = 2
const integrationOrder_shear = 1
const to = TimerOutput()

println("="^80)
println(" 執行 Pure FEM 模組：四邊固支 (CCCC)，網格加密大循環 (ndiv = 9 ➔ 25) ")
println("="^80)

for n_div in 9:25
    @timeit to "Pure_FEM Loop (ndiv=$n_div)" begin
        
        mesh_path = "/home/a/Joker/msh/st_q_$(n_div).msh"
        if !isfile(mesh_path)
            println("  [WARN] 找不到指定的網格檔案，跳過此輪: ", mesh_path)
            continue
        end
        
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(mesh_path)
        
        entities = getPhysicalGroups()
        nodes = get𝑿ᵢ()
        
        nʷ = length(nodes)
        nᵠ = length(nodes)
        
        # ❶ 建立所有元素
        elements   = getElements(nodes, entities["Ω"], integrationOrder)
        elements_s = getElements(nodes, entities["Ω"], integrationOrder_shear)
        
        boundary_names = ["Γ¹", "Γ²", "Γ³", "Γ⁴"]
        el_b_w_arrays   = [getElements(nodes, entities[nm], integrationOrder) for nm in boundary_names]
        el_b_phi_arrays = [getElements(nodes, entities[nm], integrationOrder) for nm in boundary_names]
        
        # ❷ 賦予物理常數
        prescribe!(elements,   :E=>E, :ν=>ν, :h=>h)
        prescribe!(elements_s, :E=>E, :ν=>ν, :h=>h)
        for el_b_w in el_b_w_arrays
            prescribe!(el_b_w, :α=>α_base, :g=>0.0)
        end
        for el_b_phi in el_b_phi_arrays
            prescribe!(el_b_phi, :α=>α_base, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
        end
        
        # ❸ 設定形函數導數
        set∇𝝭!(elements)
        set∇𝝭!(elements_s)
        for el_b_w in el_b_w_arrays
            set𝝭!(el_b_w)
        end
        for el_b_phi in el_b_phi_arrays
            set𝝭!(el_b_phi)
        end
        
        # ❹ 初始化剛度矩陣
        kʷʷ = zeros(nʷ, nʷ)
        kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
        kᵠʷ = zeros(2*nᵠ, nʷ)
        kᴳʷʷ = zeros(nʷ, nʷ)
        kᴳᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
        
        # ❺ 組裝彈性剛度（區域積分）
        (∫wwdΩ=>elements_s)(kʷʷ)
        (∫φwdΩ=>elements_s)(kᵠʷ)
        ([∫φφdΩ=>elements_s, ∫κκdΩ=>elements])(kᵠᵠ)
        
        # ❻ 組裝幾何剛度
        prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        (∫∇wσ∇wdΩ=>elements)(kᴳʷʷ)
        (∫∇φσ∇φdΩ=>elements)(kᴳᵠᵠ)
        
        # ❼ 組裝邊界罰函數（固支：w=0 且 φx=0, φy=0）
        for el_b_w in el_b_w_arrays
            (∫αwwdΓ=>el_b_w)(kʷʷ)
        end
        for el_b_phi in el_b_phi_arrays
            (∫αφφdΓ=>el_b_phi)(kᵠᵠ)    # 強制 φx=0, φy=0
        end
        
        # ❽ 組裝整體矩陣並求解
        K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
        Kᴳ = [kᴳᵠᵠ zeros(2*nᵠ, nʷ); zeros(nʷ, 2*nᵠ) kᴳʷʷ]
        
        F = eigen(K, Kᴳ)
        λ = F.values
        
        mode_ids = sort!(
            collect(i for i in eachindex(λ) 
                if isfinite(real(λ[i])) && abs(imag(λ[i])) < 1.0e-7 && real(λ[i]) > 0.0),
            by = i -> real(λ[i])
        )
        
        h_size = 1.0 / n_div
        log10_h = log10(h_size)
        
        # ❾ 輸出結果至 CSV
        file_exists = isfile(UNIFIED_CSV)
        open(UNIFIED_CSV, "a") do io
            if !file_exists
                println(io, "method,ndiv,h,log10_h,mode_rank,lambda_cr,k_num,k_exact,error_y")
            end
            for rank in 1:min(6, length(mode_ids))
                m_id = mode_ids[rank]
                lam = real(λ[m_id])
                k_num = (lam * h * σ₁₁) * b^2 / (π^2 * Dᵇ)
                k_ex = k_exact_modes[rank]
                error_y = (k_num / k_ex) - 1.0
                @printf(io, "Pure_FEM_CCCC,%d,%.6e,%.6f,%d,%.6e,%.6f,%.6f,%.6e\n", 
                        n_div, h_size, log10_h, rank, lam, k_num, k_ex, error_y)
            end
        end
        
        gmsh.finalize()
        println("  [Pure_FEM_CCCC] ndiv = $(n_div) 前 6 階數據導出成功。")
    end
end