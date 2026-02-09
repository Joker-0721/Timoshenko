using LinearAlgebra

### 兩構件 + 中間內鉸 + 兩端固支：Timoshenko 梁線性屈曲（廣義特徵值）
# - 兩構件：左段長 L1、右段長 L2
# - 內鉸：中間接點 w 連續、ϕ 允許左右不連續（複製轉角 DOF：ϕ2L、ϕ2R）
# - 兩端固支：w=0 且 ϕ=0（用 DOF 消去強制）
#
# 特徵值問題（以 P 為特徵值）：
#   K * d = P * KG * d
# 其中 KG 只作用在 w 子空間（幾何剛度/屈曲矩陣）。

# -----------------------------
# Step 1: 材料/截面/幾何
# -----------------------------
E = 2.1e0      # 楊氏模數 (Pa)（此例沿用原模板的尺度）
ν = 0.3        # 卜松比
L = 1.0        # 總長 (m)
L1 = L / 2     # 左構件長度 (m)
L2 = L / 2     # 右構件長度 (m)

h = 1e-1       # 梁高 (m)
b = 0.01       # 梁寬 (m)
I = b * h^3 / 12
A = b * h
κ = 5 / 6
G = E / (2 * (1 + ν))
EI = E * I
κGA = κ * A * G

# 每一構件的元素數（可調；ne1=ne2=1 就是最精簡的 2 元素模型）
ne1 = 25
ne2 = 25

# -----------------------------
# Step 2: 幾何節點與元素連接
# -----------------------------
# 幾何節點：左段 [0,L1] 與右段 [L1,L]
# w 的節點在鉸點處共用（不重複）
left_nodes = collect(range(0.0, L1, length=ne1 + 1))
right_nodes = collect(range(L1, L1 + L2, length=ne2 + 1))

# 合併（去掉 right_nodes 的第一個點，避免重複鉸點座標）
nodes = vcat(left_nodes, right_nodes[2:end])

n_nodes = length(nodes)                 # 幾何節點數
hinge_node = ne1 + 1                    # 鉸點（幾何節點索引），座標 x = L1
nₑ = ne1 + ne2                          # 總元素數

# -----------------------------
# Step 3: DOF 設計（內鉸複製 ϕ）
# -----------------------------
# w：每幾何節點 1 個 DOF => nʷ = n_nodes
# ϕ：除了鉸點要分成左/右兩個 DOF，其餘節點 1 個 DOF
#    => nᵠ = n_nodes + 1
nʷ = n_nodes
nᵠ = n_nodes + 1

# 幾何節點 -> w DOF（直接 1..n_nodes）
w_dof(node::Int) = node

# 幾何節點 -> ϕ DOF：鉸點特例，其餘節點在鉸點右側需要 +1 位移（因為多了一個 ϕ DOF）
function ϕ_dof(node::Int, side::Symbol)
    if node < hinge_node
        return node
    elseif node == hinge_node
        return side === :left ? hinge_node : hinge_node + 1
    else
        return node + 1
    end
end

# 元素 e（1-based）對應的幾何節點（i,i+1）
# 元素 1..(hinge_node-1) 在左構件；元素 hinge_node..(n_nodes-1) 在右構件

# -----------------------------
# Step 4: 初始化全域分塊矩陣（命名對齊 cantilever 範本）
# -----------------------------
kʷʷ = zeros(nʷ, nʷ)
kʷᵠ = zeros(nʷ, nᵠ)
kᵠᵠ = zeros(nᵠ, nᵠ)
kᴳ = zeros(nʷ, nʷ)

# -----------------------------
# Step 5: Gauss 積分（2 點）
# -----------------------------
ξ = [-1 / √3, 1 / √3]
w = [1.0, 1.0]

# 參考元素上的線性形狀函數與導數（對 ξ）
N1(ξ) = 0.5 * (1 - ξ)
N2(ξ) = 0.5 * (1 + ξ)
dN1_dξ = -0.5
dN2_dξ = 0.5

# -----------------------------
# Step 6: 逐元素組裝
# -----------------------------
for e in 1:nₑ
    i = e
    j = e + 1

    x1 = nodes[i]
    x2 = nodes[j]
    le = x2 - x1
    J = le / 2

    # 這個元素位於鉸點左側或右側？（決定鉸點的 ϕ 用 left 還是 right）
    element_side = (j <= hinge_node) ? :left : :right

    w_idx = (w_dof(i), w_dof(j))
    ϕ_idx = (
        ϕ_dof(i, element_side),
        ϕ_dof(j, element_side),
    )

    for (ξᵢ, wᵢ) in zip(ξ, w)
        # shape functions
        N = (N1(ξᵢ), N2(ξᵢ))

        # dN/dx = dN/dξ * dξ/dx = dN/dξ * (2/le)
        dNdx = (dN1_dξ * (2 / le), dN2_dξ * (2 / le))

        # ---- Kww: 剪力項中的 w_x * w_x ----
        for a in 1:2, b in 1:2
            kʷʷ[w_idx[a], w_idx[b]] += κGA * dNdx[a] * dNdx[b] * J * wᵢ
        end

        # ---- Kwϕ: 剪力項中的 - w_x * ϕ ----
        for a in 1:2, b in 1:2
            kʷᵠ[w_idx[a], ϕ_idx[b]] += -κGA * dNdx[a] * N[b] * J * wᵢ
        end

        # ---- Kϕϕ: 彎曲 EI ϕ_x ϕ_x + 剪力 κGA ϕ ϕ ----
        for a in 1:2, b in 1:2
            kᵠᵠ[ϕ_idx[a], ϕ_idx[b]] += (EI * dNdx[a] * dNdx[b] + κGA * N[a] * N[b]) * J * wᵢ
        end

        # ---- KG: 幾何剛度（屈曲矩陣）∫ (w_x)^T (w_x) dx ----
        for a in 1:2, b in 1:2
            kᴳ[w_idx[a], w_idx[b]] += dNdx[a] * dNdx[b] * J * wᵢ
        end
    end
end

# -----------------------------
# Step 7: 組成廣義特徵值問題
#   [Kww  Kwϕ; Kwϕ'  Kϕϕ] d = P [KG 0; 0 0] d
# -----------------------------
K = [kʷʷ kʷᵠ; kʷᵠ' kᵠᵠ]
Kᴳ = [kᴳ zeros(nʷ, nᵠ); zeros(nᵠ, nʷ) zeros(nᵠ, nᵠ)]

# -----------------------------
# Step 8: 兩端固支（w=0, ϕ=0）用 DOF 消去
# -----------------------------
ndof = nʷ + nᵠ

# 全域 DOF 排列：前 nʷ 個是 w，後 nᵠ 個是 ϕ
w_global(node::Int) = w_dof(node)
ϕ_global(node::Int, side::Symbol) = nʷ + ϕ_dof(node, side)

# 固支端：節點 1 與節點 n_nodes
constrained = Int[]
append!(constrained, [w_global(1), w_global(n_nodes)])
append!(constrained, [ϕ_global(1, :left), ϕ_global(n_nodes, :right)])

keep = trues(ndof)
for c in constrained
    keep[c] = false
end

free = findall(keep)
Kred = K[free, free]
Kᴳred = Kᴳ[free, free]

# -----------------------------
# Step 9: 解特徵值並取最小正實根
# -----------------------------
try
    λ = eigvals(Kred, Kᴳred)
    λr = real.(λ)

    # 過濾：有限、正值
    candidates = filter(x -> isfinite(x) && x > 0.0, λr)
    if isempty(candidates)
        error("找不到有限且為正的屈曲特徵值；可能是邊界條件或矩陣組裝有問題。")
    end

    Pcr = minimum(candidates)
    println("兩端固支 + 中間內鉸：第一屈曲荷重數值解 Pcr (N)：", Pcr)

    # -----------------------------
    # Step 10: 解析/近似對照（A 選項）
    # -----------------------------
    # 對於「兩端固支 + 中點內鉸」的 Euler-Bernoulli 參考，第一模態常可用
    # pinned-pinned 係數作量級近似：Pcr_Euler ≈ π^2 EI / L^2。
    Pcr_Euler = π^2 * E * I / L^2
    Pcr_Timoshenko = Pcr_Euler / (1 + (Pcr_Euler / (κ * A * G)))
    println("參考（Euler 近似）Pcr_Euler (N)：", Pcr_Euler)
    println("參考（Timoshenko 修正）Pcr_Timoshenko (N)：", Pcr_Timoshenko)

    # 補充：沒有內鉸的固定-固定 Euler buckling（僅供量級對照）
    Pcr_Euler_fixed_fixed = 4 * π^2 * E * I / L^2
    println("補充（無內鉸）固定-固定 Euler 屈曲荷重 (N)：", Pcr_Euler_fixed_fixed)
catch e
    println("特徵值計算失敗：", e)
end
