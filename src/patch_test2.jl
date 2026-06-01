# Mindlin 板 Patch Test（線性精確性檢驗）
# 使用 Gmsh 建立 2×2 Q4 網格，正確傳入 Pair 給 getElements

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using LinearAlgebra
using ApproxOperator
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φwdΩ, ∫φφdΩ
import ApproxOperator.GmshImport: get𝑿ᵢ, getElements
import Gmsh: gmsh

function run_patch_test()
    # ========== 參數設定 ==========
    E = 1.0e8
    ν = 0.3
    h = 0.1
    a = 1.0
    b = 1.0
    α = 1.0e6 * E

    w_exact(x, y) = x + y
    φx_exact(x, y) = 1.0
    φy_exact(x, y) = 1.0

    # ========== 建立網格 ==========
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()
    model = gmsh.model
    model.add("patch_test")

    model.occ.addRectangle(0, 0, 0, a, b)
    model.occ.synchronize()

    # 獲取所有面與線，並設定為結構化網格
    surfaces = model.getEntities(2)
    surf_tags = Int64[tag for (dim, tag) in surfaces]   # 確保為 Vector{Int64}
    if isempty(surf_tags)
        error("找不到任何曲面實體")
    end

    # 對所有線邊界設定分段，保證 2x2 網格 (需要 3 個節點)
    curves = model.getEntities(1)
    for (dim, tag) in curves
        model.mesh.setTransfiniteCurve(tag, 3)
    end

    # 對所有面設定跨限幾何及合成四邊形單元
    for tag in surf_tags
        model.mesh.setTransfiniteSurface(tag)
        model.mesh.setRecombine(2, tag)
    end

    model.mesh.generate(2)

    # ========== 讀取節點與單元 ==========
    nodes = get𝑿ᵢ()
    # 建構符合 getElements 所需的 Pair
    entity_pair = 2 => surf_tags   # 等價於 Pair(2, surf_tags)
    # 積分階先用 3（彎曲）
    elements = getElements(nodes, entity_pair, 3)

    n_nodes = length(nodes)
    nʷ = n_nodes
    nᵠ = n_nodes

    Kww = zeros(nʷ, nʷ)
    Kφφ = zeros(2nᵠ, 2nᵠ)
    Kφw = zeros(2nᵠ, nʷ)
    Fw  = zeros(nʷ)
    Fφ  = zeros(2nᵠ)

    prescribe!(elements, :E => E, :ν => ν, :h => h)

    # 彎曲部分 (高階積分)
    set∇𝝭!(elements)
    𝑎_κκ = ∫κκdΩ => elements
    𝑎_κκ(Kφφ)

    # 剪切部分 (低階積分)
    set∇𝝭!(elements)
    𝑎_φφ = ∫φφdΩ => elements
    𝑎_φφ(Kφφ)
    𝑎_φw = ∫φwdΩ => elements
    𝑎_φw(Kφw)

    # ========== 罰函數邊界條件 ==========
    for node in nodes
        if node.x ≈ 0.0 || node.x ≈ a || node.y ≈ 0.0 || node.y ≈ b
            row_w = node.𝐼
            Kww[row_w, row_w] += α
            Fw[row_w] += α * w_exact(node.x, node.y)

            row_φx = 2*(node.𝐼-1)+1
            Kφφ[row_φx, row_φx] += α
            Fφ[row_φx] += α * φx_exact(node.x, node.y)

            row_φy = 2*(node.𝐼-1)+2
            Kφφ[row_φy, row_φy] += α
            Fφ[row_φy] += α * φy_exact(node.x, node.y)
        end
    end

    # ========== 求解 ==========
    K = [Kφφ Kφw; Kφw' Kww]
    F = vcat(Fφ, Fw)
    U = K \ F

    # ========== 檢查中心節點 ==========
    center = findfirst(node -> node.x ≈ 0.5 && node.y ≈ 0.5, nodes)
    if center === nothing
        error("找不到中心節點 (0.5,0.5)")
    end
    idx = nodes[center].𝐼

    w_num   = U[2nᵠ + idx]
    φx_num  = U[2*(idx-1)+1]
    φy_num  = U[2*(idx-1)+2]

    w_ex   = w_exact(0.5, 0.5)
    φx_ex  = φx_exact(0.5, 0.5)
    φy_ex  = φy_exact(0.5, 0.5)

    err_w   = abs(w_num - w_ex)
    err_φx  = abs(φx_num - φx_ex)
    err_φy  = abs(φy_num - φy_ex)
    tol = 1e-7

    println("========== Mindlin 板 Patch Test 結果 ==========")
    println("網格: 2×2 Q4 單元 (3×3 節點)")
    println("精確解: w = x + y,  φx = -1, φy = -1")
    println("中心節點 (0.5,0.5) 比較:")
    println("  w_numerical  = $w_num")
    println("  w_exact      = $w_ex   → 誤差 = $err_w")
    println("  φx_numerical = $φx_num")
    println("  φx_exact     = $φx_ex  → 誤差 = $err_φx")
    println("  φy_numerical = $φy_num")
    println("  φy_exact     = $φy_ex  → 誤差 = $err_φy")
    println("----------------------------------------------")
    if err_w < tol && err_φx < tol && err_φy < tol
        println("✓ Patch test 通過！公式在線性條件下具有精確性。")
    else
        println("✗ Patch test 失敗！誤差過大，請檢查積分階數或罰因數。")
    end
    println("===============================================")

    gmsh.finalize()
end

run_patch_test()