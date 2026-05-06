using ApproxOperator
using LinearAlgebra
import Gmsh: gmsh
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Timoshenko: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ

# Mindlin-Reissner 方板材料參數。
E = 1.0e8
ν = 0.3
h = 1.0e-3

# 幾何剛度 KG 使用的參考面內壓應力場。
# 目前算例採雙向單位壓縮：σ₁₁ = σ₂₂ = 1，σ₁₂ = 0。
σ₁₁ref(x, y, z) = 1.0
σ₂₂ref(x, y, z) = 1.0
σ₁₂ref(x, y, z) = 0.0

# 優先沿用 beam_Engesser.jl 的相對路徑寫法；若從 repo 根目錄執行，
# 則退回到腳本所在位置推得的 patchtest.msh 路徑。
const MESH_FILE = isfile("msh/patchtest.msh") ? "msh/patchtest.msh" : joinpath(@__DIR__, "..", "msh", "patchtest.msh")

# patchtest.msh 的四條邊界 physical group。
# 四邊簡支在這裡只限制橫向位移 w，轉角 φ₁、φ₂ 保持自由。
const BOUNDARY_NAMES = ("Γ¹", "Γ²", "Γ³", "Γ⁴")

Pcr = NaN

# 求解 Kϕ = λ KGϕ，並只保留正的、有限的、實數特徵值。
function get_positive_finite_eigenvalues(k, kg)
    λ = eigvals(k, kg)
    λᵖ = Float64[]
    for λᵢ in λ
        if isfinite(real(λᵢ)) && isfinite(imag(λᵢ)) && abs(imag(λᵢ)) < 1.0e-8 && real(λᵢ) > 0.0
            push!(λᵖ, real(λᵢ))
        end
    end
    sort(λᵖ)
end

gmsh.initialize()
try
    # 1. 讀取 mesh 與 physical group。
    println("[Mindlin square] open msh file")
    gmsh.open(MESH_FILE)

    println("[Mindlin square] get entities")
    entities = getPhysicalGroups()
    for name in ("Ω", BOUNDARY_NAMES...)
        haskey(entities, name) || error("missing physical group: $name")
    end

    println("[Mindlin square] get nodes")
    nodes = get𝑿ᵢ()
    nʷ = length(nodes)
    nᵠ = length(nodes)

    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
    kᵠʷ = zeros(2*nᵠ, nʷ)
    kgʷʷ = zeros(nʷ, nʷ)

    # 2. 用 Mindlin patch test 的分塊寫法組裝材料剛度。
    # kʷʷ：剪切項中的 w-w 區塊。
    # kᵠʷ：剪切項中的 φ-w 耦合區塊。
    # kᵠᵠ：彎曲項與剪切項中的 φ-φ 區塊。
    println("[Mindlin square] assemble material stiffness")
    elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    set∇𝝭!(elements)
    𝑎ʷʷ = ∫wwdΩ => elements
    𝑎ᵠʷ = ∫φwdΩ => elements
    𝑎ᵠᵠ = [
        ∫φφdΩ => elements,
        ∫κκdΩ => elements,
    ]
    𝑎ʷʷ(kʷʷ)
    𝑎ᵠʷ(kᵠʷ)
    𝑎ᵠᵠ(kᵠᵠ)

    # 3. 組裝幾何剛度 KG，只需要 w-w 區塊。
    # ∫wwGdΩ2D 的原始實作寫入 [w, φ₁, φ₂] 排列的大矩陣，
    # 因此先組臨時矩陣，再抽出 w-w 區塊供分塊廣義特徵值問題使用。
    println("[Mindlin square] assemble geometric stiffness")
    prescribe!(elements, :σ₁₁=>σ₁₁ref, :σ₂₂=>σ₂₂ref, :σ₁₂=>σ₁₂ref)
    kg = zeros(3*nʷ, 3*nʷ)
    (∫wwGdΩ2D => elements)(kg)
    kgʷʷ .= kg[1:3:end, 1:3:end]

    # 4. 施加四邊簡支邊界條件：四條邊上的 w = 0。
    # 轉角自由度 φ₁、φ₂ 不約束。
    println("[Mindlin square] apply simply supported boundary")
    elements_1 = getElements(nodes, entities["Γ¹"])
    elements_2 = getElements(nodes, entities["Γ²"])
    elements_3 = getElements(nodes, entities["Γ³"])
    elements_4 = getElements(nodes, entities["Γ⁴"])
    elements_Γ = elements_1 ∪ elements_2 ∪ elements_3 ∪ elements_4
    fixed_nodes = sort(unique([node.𝐼 for element in elements_Γ for node in element.𝓒]))

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    KG = [
        zeros(2*nᵠ, 2*nᵠ) zeros(2*nᵠ, nʷ)
        zeros(nʷ, 2*nᵠ) kgʷʷ
    ]

    ndofs = 2*nᵠ + nʷ
    fixed_dofs = 2*nᵠ .+ fixed_nodes
    free_dofs = setdiff(1:ndofs, fixed_dofs)

    # 5. 解約束後的廣義特徵值問題 Kff v = λ KGff v。
    println("[Mindlin square] solve eigenvalue problem")
    λᵖ = get_positive_finite_eigenvalues(K[free_dofs, free_dofs], KG[free_dofs, free_dofs])
    isempty(λᵖ) && error("no positive finite eigenvalue found")
    global Pcr = first(λᵖ)

    println("nodes = ", nʷ)
    println("total dofs = ", ndofs)
    println("fixed w dofs = ", length(fixed_dofs))
    println("positive eigenvalues = ", λᵖ[1:min(5, length(λᵖ))])
finally
    gmsh.finalize()
end

println("Mindlin square λ_cr = ", Pcr)
