using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ


using TimerOutputs, LinearAlgebra, WriteVTK, DelimitedFiles
import Gmsh: gmsh
include("cal_area_support.jl")


E = 1.0
ν = 0.3
ρ = 1.0
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
σ₁₁ = 1e0
σ₂₂ = 0.0
σ₁₂ = 0.0
a = 1.0
αʷ = 0e6
αᵠ = 0e3

# ============================================================
#  邊界條件設定
#  w_edges     : 參與 w 邊界積分（∫QwdΓ, ∫αwwdΓ）的邊號
#  assemble_w  : 是否執行 w 邊界組裝（SSSS 為 false）
#  φ_edges     : 參與 φ 邊界積分（∫MφdΓ, ∫αφφdΓ）的邊號
#  assemble_φ  : 是否執行 φ 邊界組裝（SSSS、FSSS 為 false）
# ============================================================
bcs = [
    (name="CCCC", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[1,2,3,4], assemble_φ=true),
    (name="SSSS", w_edges=[1,2,3,4], assemble_w=true, φ_edges=[1,2,3,4], assemble_φ=false),
    (name="CSCS", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[1,3],     assemble_φ=true),
    (name="SCSC", w_edges=[1,2,3,4], assemble_w=true,  φ_edges=[2,4],     assemble_φ=true),
    (name="FSCS", w_edges=[2,3,4],   assemble_w=true,  φ_edges=[3],       assemble_φ=true),
    (name="FSSS", w_edges=[2,3,4],   assemble_w=true,  φ_edges=[3],       assemble_φ=false),
]

to = TimerOutput()

# 為每個邊界條件開啟獨立 CSV 檔
ios = Dict{String,IOStream}()
for bc in bcs
    mkpath("../date/mf/vibration/plate")
    ios[bc.name] = open("../date/mf/vibration/plate/vib_mf_w_φ_$(bc.name)2_0.001.csv", "w")
    write(ios[bc.name], "nʷ,nᵠ,nˢ,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀\n")
end

integrationOrder = 2

type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})

ndivs = 20:30
for ndiv in ndivs
    ndiv_φ = ndiv
    ndiv_w = ndiv-2
    ndiv_q = ndiv

    gmsh.initialize()
    @timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv_w.msh")
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes_w = get𝑿ᵢ()
    xʷ = nodes_w.x
    yʷ = nodes_w.y
    zʷ = nodes_w.z
    sp_w = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
    elements_support = getElements(nodes_w, entities["Ω"], 1)
    s_w, var_A = cal_area_support(elements_support)
    nʷ = length(nodes_w)
    s₁ = 1.5*s_w*ones(nʷ)
    s₂ = 1.5*s_w*ones(nʷ)
    s₃ = 1.5*s_w*ones(nʷ)
    push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

    @timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv_φ.msh")
    @timeit to "get nodes" nodes_φ = get𝑿ᵢ()
    @timeit to "get entities" entities = getPhysicalGroups()
    xᵠ = nodes_φ.x
    yᵠ = nodes_φ.y
    zᵠ = nodes_φ.z
    sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
    nᵠ = length(nodes_φ)
    elements_support = getElements(nodes_φ, entities["Ω"], 1)
    s_φ, var_A = cal_area_support(elements_support)
    s₁ = 1.5*s_φ*ones(nᵠ)
    s₂ = 1.5*s_φ*ones(nᵠ)
    s₃ = 1.5*s_φ*ones(nᵠ)
    push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

    @timeit to "open msh file" gmsh.open("../msh/patchtest_tri3_$ndiv_q.msh")
    @timeit to "get nodes" nodes = get𝑿ᵢ()
    @timeit to "get entities" entities = getPhysicalGroups()
    nˢ = length(nodes)

    nₑ = length(elements_support)
    nᵐ = nₑ*ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[],𝑿ₛ[]))

    kʷʷ = zeros(nʷ,nʷ)
    kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kˢˢ = zeros(2*nˢ,2*nˢ)
    kᴳʷʷ = zeros(nʷ,nʷ)
    kᴳᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    mʷʷ = zeros(nʷ,nʷ)
    mᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kᵠʷ = zeros(2*nᵠ,nʷ)
    kˢʷ = zeros(2*nˢ,nʷ)
    kˢᵠ = zeros(2*nˢ,2*nᵠ)
    kᵐᵐ = zeros(3*nᵐ,3*nᵐ)
    kᵐᵠ = zeros(3*nᵐ,2*nᵠ)
    kᵐʷ = zeros(3*nᵐ,nʷ)
    kˢᵐ = zeros(2*nˢ,3*nᵐ)
    fˢ = zeros(2*nˢ)
    fᵐ = zeros(3*nᵐ)

    @timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
        @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"],integrationOrder)
        prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
        @timeit to "calculate shape functions" set∇𝝭!(elements_q)

        @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
        prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
        @timeit to "calculate shape functions" set∇𝝭!(elements_w)

        @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
        prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
        @timeit to "calculate shape functions" set∇𝝭!(elements_φ)

        @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
        prescribe!(elements_m, :E=>E, :ν=>ν, :h=>h)
        @timeit to "calculate shape functions" set∇𝝭!(elements_m)

        @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

        @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
        @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)

        @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

        @timeit to "get elements" elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
        @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

        𝑎ˢˢ = ∫QQdΩ=>elements_q
        𝑎ˢʷ = [
            ∫∇QwdΩ=>(elements_q,elements_w),
            ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
        ]
        𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
        𝑎ᵐᵐ = ∫MMdΩ=>elements_m
        𝑎ᵐᵠ = [
            ∫∇MφdΩ=>(elements_m,elements_φ),
            ∫MφdΓ=>(elements_m_Γ,elements_φ_Γ),
        ]
        𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements_w
        𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements_φ
        𝑎ᵐʷʷ = ∫ρwwdΩ=>elements_w
        𝑎ᵐᵠᵠ = ∫ρφφdΩ=>elements_φ
        @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
        @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
        @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
        @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
        @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
        @timeit to "assemble" 𝑎ᴳʷʷ(kᴳʷʷ)
        @timeit to "assemble" 𝑎ᴳᵠᵠ(kᴳᵠᵠ)
        @timeit to "assemble" 𝑎ᵐʷʷ(mʷʷ)
        @timeit to "assemble" 𝑎ᵐᵠᵠ(mᵠᵠ)
    end

    @timeit to "calculate  ∫QwdΓ" begin
        @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
        @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
        @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
        @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
        @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
        @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
        prescribe!(elements_w_1, :α=>αʷ, :g=>0.0)
        prescribe!(elements_w_2, :α=>αʷ, :g=>0.0)
        prescribe!(elements_w_3, :α=>αʷ, :g=>0.0)
        prescribe!(elements_w_4, :α=>αʷ, :g=>0.0)
        @timeit to "calculate shape functions" set𝝭!(elements_q_1)
        @timeit to "calculate shape functions" set𝝭!(elements_q_2)
        @timeit to "calculate shape functions" set𝝭!(elements_q_3)
        @timeit to "calculate shape functions" set𝝭!(elements_q_4)
        @timeit to "calculate shape functions" set𝝭!(elements_w_1)
        @timeit to "calculate shape functions" set𝝭!(elements_w_2)
        @timeit to "calculate shape functions" set𝝭!(elements_w_3)
        @timeit to "calculate shape functions" set𝝭!(elements_w_4)
        # w 邊界組裝移至下方 BC 迴圈中執行
    end

    @timeit to "calculate ∫MφdΓ" begin
        @timeit to "get elements" elements_m_1 = getElements(entities["Γ¹"], entities["Γ"], elements_m_Γ)
        @timeit to "get elements" elements_m_2 = getElements(entities["Γ²"], entities["Γ"], elements_m_Γ)
        @timeit to "get elements" elements_m_3 = getElements(entities["Γ³"], entities["Γ"], elements_m_Γ)
        @timeit to "get elements" elements_m_4 = getElements(entities["Γ⁴"], entities["Γ"], elements_m_Γ)
        @timeit to "get elements" elements_φ_1 = getElements(nodes_φ, entities["Γ¹"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        @timeit to "get elements" elements_φ_2 = getElements(nodes_φ, entities["Γ²"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        @timeit to "get elements" elements_φ_3 = getElements(nodes_φ, entities["Γ³"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        @timeit to "get elements" elements_φ_4 = getElements(nodes_φ, entities["Γ⁴"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        prescribe!(elements_φ_1, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
        prescribe!(elements_φ_2, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
        prescribe!(elements_φ_3, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
        prescribe!(elements_φ_4, :α=>αᵠ, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
        @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
        @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
        @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
        @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
        # φ 邊界組裝移至下方 BC 迴圈中執行
    end

    gmsh.finalize()

    # ----------------------------------------------------------
    #  儲存域積分（含全邊界通用項）後的矩陣快照
    #  每個 BC 開始前從此狀態還原，避免互相污染
    # ----------------------------------------------------------
    kˢʷ₀ = copy(kˢʷ)
    kᵐᵠ₀ = copy(kᵐᵠ)

    # 各邊元素陣列，便於依 BC 動態組合
    elements_q_edges = [elements_q_1, elements_q_2, elements_q_3, elements_q_4]
    elements_w_edges = [elements_w_1, elements_w_2, elements_w_3, elements_w_4]
    elements_m_edges = [elements_m_1, elements_m_2, elements_m_3, elements_m_4]
    elements_φ_edges = [elements_φ_1, elements_φ_2, elements_φ_3, elements_φ_4]

    # VTK 格點與 cells（與 BC 無關，只算一次）
    xs = [node.x for node in nodes]'
    ys = [node.y for node in nodes]'
    zs = [node.z for node in nodes]'
    points = [xs; ys; zs]
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_q]

    # ==============================================================
    #  邊界條件迴圈
    # ==============================================================
    for bc in bcs
        @timeit to "BC $(bc.name)" begin

            # 還原 BC 相關矩陣至域積分後狀態
            kˢʷ_bc = copy(kˢʷ₀)
            kᵐᵠ_bc = copy(kᵐᵠ₀)
            kʷʷ_bc = zeros(nʷ, nʷ)
            kᵠᵠ_bc = zeros(2*nᵠ, 2*nᵠ)
            kᵠʷ_bc = zeros(2*nᵠ, nʷ)
            fˢ_bc = zeros(2*nˢ)
            fᵐ_bc = zeros(3*nᵐ)

            # ---------- w 邊界條件 ----------
            q_w = reduce(∪, [elements_q_edges[i] for i in bc.w_edges])
            w_w = reduce(∪, [elements_w_edges[i] for i in bc.w_edges])
            𝑎_w  = ∫QwdΓ   => (q_w, w_w)
            𝑎ʷ_w = ∫αwwdΓ => w_w
            if bc.assemble_w
                @timeit to "assemble" 𝑎_w(kˢʷ_bc, fˢ_bc)
                @timeit to "assemble" 𝑎ʷ_w(kʷʷ_bc)
            end

            # ---------- φ 邊界條件 ----------
            m_φ = reduce(∪, [elements_m_edges[i] for i in bc.φ_edges])
            φ_φ = reduce(∪, [elements_φ_edges[i] for i in bc.φ_edges])
            𝑎_φ   = ∫MφdΓ   => (m_φ, φ_φ)
            𝑎ᵅ_φ  = ∫αφφdΓ => φ_φ
            if bc.assemble_φ
                @timeit to "assemble" 𝑎_φ(kᵐᵠ_bc, fᵐ_bc)
                @timeit to "assemble" 𝑎ᵅ_φ(kᵠᵠ_bc)
            end

            # ---------- Schur 補餘 ----------
            kᵠᵠ_bc .+= - kˢᵠ'*(kˢˢ\kˢᵠ) - kᵐᵠ_bc'*(kᵐᵐ\kᵐᵠ_bc)
            kᵠʷ_bc .+= - kˢᵠ'*(kˢˢ\kˢʷ_bc)
            kʷʷ_bc .+= - kˢʷ_bc'*(kˢˢ\kˢʷ_bc)

            # ---------- 組裝整體矩陣 ----------
            k_mat  = [kᵠᵠ_bc kᵠʷ_bc; kᵠʷ_bc' kʷʷ_bc]
            kᴳ_mat = [kᴳᵠᵠ zeros(2*nᵠ,nʷ); zeros(nʷ,2*nᵠ) kᴳʷʷ]
            m_mat  = [mᵠᵠ  zeros(2*nᵠ,nʷ); zeros(nʷ,2*nᵠ) mʷʷ]

            # ---------- 特徵值求解 ----------
            # λ, v = eigen(k_mat, kᴳ_mat)
            λ, v = eigen(k_mat, m_mat)
            index = findfirst(real.(λ).>1e-8)

            # ---------- 計算位移場 d（用於 VTK 輸出） ----------
            n_index = 12
            d = zeros(nˢ, n_index)
            𝗠 = zeros(21)
            for (i, xᵢ) in enumerate(nodes)
                x = xᵢ.x
                y = xᵢ.y
                indices = sp_w(x, y, 0.0)
                ni = length(indices)
                𝓒 = [nodes_w[i] for i in indices]
                data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
                ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
                𝓖 = [ξ]
                a_elm = eval(type_w)(𝓒, 𝓖)
                set𝝭!(a_elm)
                for j in 1:n_index
                    u = 0.0
                    N = ξ[:𝝭]
                    for (k_idx, xₖ) in enumerate(𝓒)
                        I = xₖ.𝐼
                        u += N[k_idx]*v[2*nᵠ+I, index+j-1]
                    end
                    d[i,j] = u
                end
            end

            push!(nodes,
                :d₁=>d[:,1],  :d₂=>d[:,2],  :d₃=>d[:,3],  :d₄=>d[:,4],
                :d₅=>d[:,5],  :d₆=>d[:,6],  :d₇=>d[:,7],  :d₈=>d[:,8],
                :d₉=>d[:,9],  :d₁₀=>d[:,10],:d₁₁=>d[:,11],:d₁₂=>d[:,12],
            )

            # ---------- VTK 輸出 ----------
            mkpath("./vtk/mf/vibration/$(bc.name)")
            vtk_grid("./vtk/mf/vibration/$(bc.name)/vib_mf_w_φ_$(bc.name)2_0.001.vtu", points, cells;
                     ascii=false, append=false, compress=false) do vtk
                vtk["v₁"]  = [node.d₁  for node in nodes]
                vtk["v₂"]  = [node.d₂  for node in nodes]
                vtk["v₃"]  = [node.d₃  for node in nodes]
                vtk["v₄"]  = [node.d₄  for node in nodes]
                vtk["v₅"]  = [node.d₅  for node in nodes]
                vtk["v₆"]  = [node.d₆  for node in nodes]
                vtk["v₇"]  = [node.d₇  for node in nodes]
                vtk["v₈"]  = [node.d₈  for node in nodes]
                vtk["v₉"]  = [node.d₉  for node in nodes]
                vtk["v₁₀"] = [node.d₁₀ for node in nodes]
                vtk["v₁₁"] = [node.d₁₁ for node in nodes]
                vtk["v₁₂"] = [node.d₁₂ for node in nodes]
            end

            # ---------- 特徵值結果輸出 ----------
            # k = (λ[index]*a^2/(π^2*Dᵇ)*h)
            # println("$(bc.name): ", k)
            # write(ios[bc.name], "$ndiv_w,$ndiv_φ,$ndiv_q,$k\n")
            # flush(ios[bc.name])

            n_modes = 10
            # 由 index 開始取前 10 個非零特徵值，轉成無量綱頻率參數 k
            k_vals = [sqrt(real(λ[index+j-1])) * a^2 * sqrt(ρ*h/Dᵇ) for j in 1:n_modes]
            println("$(bc.name): ", k_vals)
            write(ios[bc.name], "$ndiv_w,$ndiv_φ,$ndiv_q,$(join(k_vals, ","))\n")
            flush(ios[bc.name])


            
        end
    end
end

# 關閉所有 CSV 檔
for bc in bcs
    close(ios[bc.name])
end
