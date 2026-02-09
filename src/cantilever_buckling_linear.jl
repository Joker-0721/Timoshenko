
using LinearAlgebra

### 步驟 1：定義梁的材料/幾何參數與離散
E = 2.1e0  # 楊氏模數 (Pa)
ν = 0.3     # 卜松比
L = 1.0     # 梁長 (m)
h = 1e-1    # 梁高 (m)
b = 0.01    # 梁寬 (m)
I = b * h^3 / 12 # 截面二次矩 (m^4)
A = b * h      # 截面面積 (m^2)
κ = 5/6     # 鐵摩辛柯剪力修正係數
G = E / (2 * (1 + ν)) # 剪力模數 (Pa)
EI = E*I
κGA = κ * A * G

nₑ = 50 # 元素數
nₚ = nₑ + 1

# 產生節點座標
nodes = range(0, L, length=nₚ)
# 定義自由度數（此例每個節點有 w 與 ϕ 各一個 DOF）
nʷ = nₚ
nᵠ = nₚ

### 步驟 2：初始化全域矩陣（使用稠密矩陣）
kʷʷ = zeros(nʷ, nʷ)
kʷᵠ = zeros(nʷ, nᵠ)
kᵠᵠ = zeros(nᵠ, nᵠ)
kᴳ = zeros(nʷ, nʷ)

### 步驟 3：設定 2 點 Gauss 積分點與權重
ξ = [-1/√3, 1/√3]
w = [1.0, 1.0]

### 步驟 4：逐元素計算並組裝全域矩陣
for i in 1:nₑ
    le = nodes[i+1] - nodes[i]
    for (ξᵢ, wᵢ) in zip(ξ, w)
        N₁ = 0.5 * (1 - ξᵢ)
        N₂ = 0.5 * (1 + ξᵢ)
        B₁ = -le/2
        B₂ = le/2

        # 座標映射的 Jacobian
        J = le / 2

        # 元素剛度矩陣的各分量（組裝到全域矩陣）
        kʷʷ[i,i] += κGA*B₁*B₁*J*wᵢ
        kʷʷ[i,i+1] += κGA*B₁*B₂*J*wᵢ
        kʷʷ[i+1,i] += κGA*B₂*B₁*J*wᵢ
        kʷʷ[i+1,i+1] += κGA*B₂*B₂*J*wᵢ

        kʷᵠ[i,i] -= κGA*B₁*N₁*J*wᵢ
        kʷᵠ[i,i+1] -= κGA*B₁*N₂*J*wᵢ
        kʷᵠ[i+1,i] -= κGA*B₂*N₁*J*wᵢ
        kʷᵠ[i+1,i+1] -= κGA*B₂*N₂*J*wᵢ

        kᵠᵠ[i,i] += (EI*B₁*B₁ + κGA*N₁*N₁)*J*wᵢ
        kᵠᵠ[i,i+1] += (EI*B₁*B₂ + κGA*N₁*N₂)*J*wᵢ
        kᵠᵠ[i+1,i] += (EI*B₂*B₁ + κGA*N₂*N₁)*J*wᵢ
        kᵠᵠ[i+1,i+1] += (EI*B₂*B₂ + κGA*N₂*N₂)*J*wᵢ

        # 元素幾何剛度矩陣（屈曲矩陣）
        kᴳ[i,i] += B₁*B₁*J*wᵢ
        kᴳ[i,i+1] += B₁*B₂*J*wᵢ
        kᴳ[i+1,i] += B₂*B₁*J*wᵢ
        kᴳ[i+1,i+1] += B₂*B₂*J*wᵢ
    end
end

α = 1e8
kᵅʷ = zeros(nʷ, nʷ)
kᵅᵠ = zeros(nᵠ, nᵠ)
kᵅʷ[1,1] += α*κGA
kᵅᵠ[1,1] += α*(EI+κGA)
# kᵅʷ[nʷ,nʷ] += α*κGA
# kᵅᵠ[nᵠ,nᵠ] += α*(EI+κGA)


### 步驟 6：用 LinearAlgebra.eigvals 解廣義特徵值問題
# 解 K*v = λ*B*v（此處 B 與幾何剛度矩陣相關）
try
    # eigenvalues = eigvals(k, -kᴳ)
    # eigenvectors = eigvecs(k, -kᴳ)
    eigenvalues = eigvals([kʷʷ+kᵅʷ kʷᵠ;kʷᵠ' kᵠᵠ+kᵅᵠ], [kᴳ+kᵅʷ zeros(nʷ, nᵠ); zeros(nᵠ, nʷ) kᵅᵠ])
    # eigenvalues = eigvals([-kᴳ zeros(nʷ, nᵠ); zeros(nᵠ, nʷ) zeros(nᵠ, nᵠ)],[kʷʷ kʷᵠ;kʷᵠ' kᵠᵠ])
    
    # 取最小的正實特徵值（對應第一屈曲荷重）
    # println("Computed eigenvalues: ", eigenvalues)
    # println(eigenvectors)
    buckling_loads = real(eigenvalues)
    first_buckling_load = minimum(filter(x -> x > 0.0, buckling_loads))
    
    println("Numerical solution for the first buckling load (N): ", first_buckling_load)
catch e
    println("Eigenvalue computation failed: ", e)
end

### Step 7: Calculate and display analytical solution
# Analytical solution for a cantilever Timoshenko beam
Pcr_Euler = (π^2 * E * I) / (4 * L^2)
Pcr_Timoshenko = Pcr_Euler / (1 + (Pcr_Euler / (κ * A * G)))

println("Analytical Euler buckling load (N): ", Pcr_Euler)
println("Analytical Timoshenko buckling load (N): ", Pcr_Timoshenko)
