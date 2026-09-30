# ============================================================================
#  plate_buckling_exact.jl — 方板單軸屈曲「標準解」（公式 → 數值求根 → VTK）
#
#  出圖鏈（使用者 2026-09-24 指令，見 AGENTS.md B 節）：
#      公式給 w(x,y) 每階模態 → 數值求根 → 寫 .vtu（WriteVTK）
#      → Makie 讀 VTK 出單張 PNG（RdBu、contour 12、無色條/邊框/標題）
#      → LaTeX 逐張 \includegraphics
#
#  板 a = b = 1，單軸壓縮沿 x：λ = N_x b²/D，無量綱係數 k = λ/π²
#  控制方程：D∇⁴w + N_x w,xx = 0
#
#  邊界字串次序與論文一致：底(y=0) 右(x=1) 頂(y=1) 左(x=0)
#      SSSS  Navier 雙三角閉式      k = (m²+n²)²/m²
#      CSCS  底/頂固支、左/右簡支    Levy（x 向 sin 展開，Y 兩端固支）
#      FSCS  底自由、頂固支、左/右簡支 Levy（Y 自由–固支）
#      FSSS  底自由、其餘簡支        Levy（Y 自由–簡支）
#      SCSC  左/右固支、底/頂簡支    Levy（y 向 sin 展開，X 兩端固支）
#      CCCC  四邊固支（無簡支對邊）  Rayleigh–Ritz（固支樑函數基）
#
#  Levy 基（對 s 解析，無 s→0 退化偽根）：
#      Y(y) = A·C(s₁,y) + B·S(s₁,y) + C·C(s₂,y) + D·S(s₂,y)
#      C(s,y) = cos(√s·y)，S(s,y) = sin(√s·y)/√s（s→0 取極限 y）
#      端部邊界：固支 Y=Y'=0；簡支 Y=Y''=0；
#               自由 M∝Y''−νμ²Y=0、V∝Y'''−(2−ν)μ²Y'=0
#
#  執行：julia plate_buckling_exact.jl        （vtu 寫到 <repo>/vtk/std/buckling/）
# ============================================================================
using LinearAlgebra, Printf, WriteVTK

const NU      = 0.3
const NP      = 21                        # 每向節點數（80×80 單元）
const NMODE   = 12
const XS      = range(0.0, 1.0, length = NP)
const VFIELDS = ["v₁","v₂","v₃","v₄","v₅","v₆",
                 "v₇","v₈","v₉","v₁₀","v₁₁","v₁₂"]
const REF = Dict("SSSS" => 4.0000, "CCCC" => 10.0700, "CSCS" => 7.6900,
                 "SCSC" => 6.7500, "FSCS" => 1.7000,  "FSSS" => 1.4400)

# ---------------------------------------------------------------------------
#  通用 Levy 基函數
# ---------------------------------------------------------------------------
cv(s, y) = s == 0 ? one(ComplexF64) : cos(sqrt(complex(s)) * y)
function sv(s, y)
    if abs(s) * y * y < 1e-12
        return complex(y) * (1 - s*y^2/6 + s^2*y^4/120)
    end
    return sin(sqrt(complex(s)) * y) / sqrt(complex(s))
end

# 第 k 階導數（k = 0,1,2,3）
Ck_s(s, y, k) = k == 0 ? cv(s, y) : k == 1 ? -s * sv(s, y) :
                k == 2 ? -s * cv(s, y) : s^2 * sv(s, y)
Sk_s(s, y, k) = k == 0 ? sv(s, y) : k == 1 ? cv(s, y) :
                k == 2 ? -s * sv(s, y) : -s * cv(s, y)
# 對 s 的偏導（僅用於 s₁≈s₂ 的重根情形）
dC_s(s, y, k) = k == 0 ? -(y/2)*sv(s, y) : k == 1 ? -(y*cv(s, y) + sv(s, y))/2 :
                k == 2 ? -cv(s, y) + (s*y/2)*sv(s, y) : (s*y*cv(s, y) + 3s*sv(s, y))/2
dS_s(s, y, k) = k == 0 ? (y*cv(s, y) - sv(s, y))/(2s) : k == 1 ? -(y/2)*sv(s, y) :
                k == 2 ? -(y*cv(s, y) + sv(s, y))/2 : -cv(s, y) + (s*y/2)*sv(s, y)

# 4 個基值 [C(s₁), S(s₁), ΔC/Δs, ΔS/Δs]
# 用差商取代直接的第二對基以消除 s₁=s₂ 時嘅退化偽零（如 SCSC 嘅 Λ=4ν²）
function basis_row(s1, s2, y, k)
    Δ = s2 - s1
    v = Vector{ComplexF64}(undef, 4)
    v[1] = Ck_s(s1, y, k); v[2] = Sk_s(s1, y, k)
    if abs(Δ) > 1e-7 * max(1.0, abs(s1))
        v[3] = (Ck_s(s2, y, k) - v[1]) / Δ
        v[4] = (Sk_s(s2, y, k) - v[2]) / Δ
    else
        v[3] = dC_s(s1, y, k); v[4] = dS_s(s1, y, k)
    end
    return v
end

# 4×4 邊界行列式矩陣（行序：w(0)、w(1)、第二條件(0)、第二條件(1)）
function levy_matrix(s1, s2, mueff, bc0, bc1; ν::Float64 = NU)
    M = zeros(ComplexF64, 4, 4)
    for (i, (kind, y, ord)) in enumerate(((bc0, 0.0, 0), (bc1, 1.0, 0),
                                          (bc0, 0.0, 1), (bc1, 1.0, 1)))
        local row
        if kind == :C
            row = basis_row(s1, s2, y, ord == 0 ? 0 : 1)
        elseif kind == :S
            row = basis_row(s1, s2, y, ord == 0 ? 0 : 2)
        elseif kind == :F
            r0 = basis_row(s1, s2, y, 0); r1 = basis_row(s1, s2, y, 1)
            r2 = basis_row(s1, s2, y, 2); r3 = basis_row(s1, s2, y, 3)
            row = ord == 0 ? (r2 .- ν * mueff^2 .* r0) :
                             (r3 .- (2 - ν) * mueff^2 .* r1)
        else
            error("未知邊界 $kind")
        end
        sc = maximum(abs, row)
        M[i, :] = sc > 0 ? row ./ sc : row
    end
    return M
end

function levy_det(get_s, Lam, mueff, bc0, bc1; ν::Float64 = NU)
    s1, s2 = get_s(Lam)
    return real(det(levy_matrix(s1, s2, mueff, bc0, bc1; ν = ν)))
end

# 掃描 + 二分求根（Λ ∈ (0, kmax·π²]）
function scan_roots(get_s, mueff, bc0, bc1; kmax::Float64 = 60.0,
                    dk::Float64 = 0.1, ν::Float64 = NU)
    roots = Float64[]
    Lhi = kmax * π^2
    dΛ  = dk * π^2
    Λ   = 1e-6
    dprev = levy_det(get_s, Λ, mueff, bc0, bc1; ν = ν)
    while Λ < Lhi
        Λn = Λ + dΛ
        dn = levy_det(get_s, Λn, mueff, bc0, bc1; ν = ν)
        if isfinite(dprev) && isfinite(dn) && sign(dn) != sign(dprev)
            lo, hi = Λ, Λn
            dlo = dprev
            for _ in 1:80
                mid = (lo + hi) / 2
                dm = levy_det(get_s, mid, mueff, bc0, bc1; ν = ν)
                if sign(dm) == sign(dlo)
                    lo = mid; dlo = dm
                else
                    hi = mid
                end
            end
            push!(roots, (lo + hi) / 2)
        end
        Λ = Λn; dprev = dn
    end
    return roots
end

# type = :x（左/右簡支 → w = Σ Y(y) sin(mπx)）
# type = :y（底/頂簡支 → w = Σ X(x) sin(nπy)）
# 註：基函數 C(s,y)=cos(√s·y) 對應 Y'' = −s·Y，故 s 由下列二次方程解出
#     type :x  Y''''−2μ²Y''+(μ⁴−Λμ²)Y=0  →  s²+2μ²s+(μ⁴−Λμ²)=0
#     type :y  X''''+(Λ−2ν²)X''+ν⁴X=0    →  s²−(Λ−2ν²)s+ν⁴=0
function get_s_x(μ)
    return Lam -> begin
        r = μ * sqrt(Lam)
        return (-μ^2 + r, -μ^2 - r)
    end
end
function get_s_y(ν_)
    return Lam -> begin
        p2 = Lam - 2ν_^2
        d  = sqrt(complex(p2^2 - 4ν_^4))
        ((p2 + d) / 2, (p2 - d) / 2)
    end
end

function levy_modes(typ::Symbol, bc0, bc1; wrange = 1:8, kmax = 60.0, dk = 0.1)
    out = Tuple{Float64,Int}[]          # (k, 波數)
    for w in wrange
        mueff = w * π
        get_s = typ == :x ? get_s_x(mueff) : get_s_y(mueff)
        for r in scan_roots(get_s, mueff, bc0, bc1; kmax = kmax, dk = dk)
            push!(out, (r / π^2, w))
        end
    end
    sort!(out, by = first)
    return out
end

# 由根取 Levy 剖面（長度 NP 的實向量）
function levy_profile(typ::Symbol, k::Float64, w::Int, bc0, bc1)
    μ = w * π
    get_s = typ == :x ? get_s_x(μ) : get_s_y(μ)
    s1, s2 = get_s(k * π^2)
    M = levy_matrix(s1, s2, μ, bc0, bc1)
    V = svd(M).V
    c = V[:, end]
    return real([sum(c[j] * basis_row(s1, s2, y, 0)[j] for j in 1:4) for y in XS])
end

function grid_from_profile(typ::Symbol, w::Int, prof::Vector{Float64})
    W = zeros(NP, NP)
    for ix in 1:NP, iy in 1:NP
        W[ix, iy] = typ == :x ? prof[iy] * sin(w * π * XS[ix]) :
                                prof[ix] * sin(w * π * XS[iy])
    end
    return W
end

function normalize!(W)
    s = maximum(abs, W)
    s > 0 && (W ./= s)
    if W[cld(NP, 2), cld(NP, 2)] < 0
        W .= -W
    end
    return W
end

# ---------------------------------------------------------------------------
#  SSSS：Navier 閉式
# ---------------------------------------------------------------------------
function ssss_modes()
    lst = sort(vec([(((m^2 + n^2)^2 / m^2), m, n) for m in 1:8, n in 1:8]), by = first)
    out = Tuple{Float64,String,Matrix{Float64}}[]
    for (k, m, n) in lst[1:NMODE]
        W = [sin(m * π * XS[ix]) * sin(n * π * XS[iy]) for ix in 1:NP, iy in 1:NP]
        push!(out, (k, "$m,$n", normalize!(W)))
    end
    return out
end

# ---------------------------------------------------------------------------
#  CCCC：Rayleigh–Ritz（固支樑函數基）
# ---------------------------------------------------------------------------
function gauss_legendre(n::Int)
    β = [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* (E.vectors[1, :]).^2      # 節點、權重
end
const GLX, GLW = gauss_legendre(200)

# 固支–固支樑精確特徵值 β_n：cos β cosh β = 1
function cc_beta(n::Int)
    β = (n + 0.5) * π
    for _ in 1:200
        g  = cos(β) * cosh(β) - 1
        gp = -sin(β) * cosh(β) + cos(β) * sinh(β)
        d  = g / gp
        β -= d
        abs(d) < 1e-15 * max(1.0, abs(β)) && break
    end
    return β
end
const CCB = [cc_beta(n) for n in 1:16]

function cc_coefs(β::Float64)
    den = sinh(β) - sin(β)
    oms = (-exp(-β) + cos(β) - sin(β)) / den      # 1 − σ
    ops = ( exp(β)  - sin(β) - cos(β)) / den      # 1 + σ
    σ   = (cosh(β) - cos(β)) / den
    return oms, ops, σ
end
# 以下三式用 e^{±βt} 形式改寫，避免 cosh/sinh 大數相減（高階 β 才穩定）
function beam_cc(n::Int, t::Float64)
    β = CCB[n]; oms, ops, σ = cc_coefs(β)
    p1 = (oms * exp(β * t) + ops * exp(-β * t)) / 2
    return p1 - (cos(β * t) - σ * sin(β * t))
end
function beam_cc_d(n::Int, t::Float64)
    β = CCB[n]; oms, ops, σ = cc_coefs(β)
    return β * ((oms * exp(β * t) - ops * exp(-β * t)) / 2 + sin(β * t) + σ * cos(β * t))
end
function beam_cc_dd(n::Int, t::Float64)
    β = CCB[n]; oms, ops, σ = cc_coefs(β)
    p1 = (oms * exp(β * t) + ops * exp(-β * t)) / 2
    return β^2 * (p1 + cos(β * t) - σ * sin(β * t))
end
beam_ss(n, t)    = sin(n * π * t)
beam_ss_d(n, t)  = n * π * cos(n * π * t)
beam_ss_dd(n, t) = -(n * π)^2 * sin(n * π * t)

# 固支基（條件數良好）：ψ_n(t) = sin(πt)·sin(nπt)，滿足 ψ=ψ'=0 於 t=0,1
beam_cl(n, t)    = sin(π * t) * sin(n * π * t)
beam_cl_d(n, t)  = π * (cos(π * t) * sin(n * π * t) + n * sin(π * t) * cos(n * π * t))
beam_cl_dd(n, t) = π^2 * (2n * cos(π * t) * cos(n * π * t) -
                          (1 + n^2) * sin(π * t) * sin(n * π * t))

function ritz_plate(kind::Symbol, N::Int)
    f  = kind == :C ? beam_cl    : beam_ss
    d1 = kind == :C ? beam_cl_d  : beam_ss_d
    d2 = kind == :C ? beam_cl_dd : beam_ss_dd
    A0 = zeros(N, N); A1 = zeros(N, N); A2 = zeros(N, N)
    for p in 1:N, q in 1:N
        s0 = 0.0; s1 = 0.0; s2 = 0.0
        for (t, w) in zip(GLX, GLW)
            s0 += w * f(p, t) * f(q, t)
            s1 += w * d1(p, t) * d1(q, t)
            s2 += w * d2(p, t) * d2(q, t)
        end
        A0[p, q] = s0; A1[p, q] = s1; A2[p, q] = s2
    end
    # 各基函數歸一化（改善高階樑函數的病態），廣義特徵值不變
    d = sqrt.(diag(A0))
    A0 ./= d * d'; A1 ./= d * d'; A2 ./= d * d'
    K = kron(A2, A0) + 2 * kron(A1, A1) + kron(A0, A2)
    G = kron(A1, A0)
    E = eigen(Symmetric(K), Symmetric(G))
    return E.values ./ π^2, E.vectors
end

function cccc_modes(N::Int = 14)
    ks, V = ritz_plate(:C, N)
    out = Tuple{Float64,Matrix{Float64}}[]
    for j in 1:NMODE
        C = reshape(V[:, j], N, N)
        W = zeros(NP, NP)
        for ix in 1:NP, iy in 1:NP
            s = 0.0
            for p in 1:N, q in 1:N
                s += C[p, q] * beam_cl(p, XS[ix]) * beam_cl(q, XS[iy])
            end
            W[ix, iy] = s
        end
        push!(out, (ks[j], normalize!(W)))
    end
    return out
end

# ---------------------------------------------------------------------------
#  VTU 輸出（與既有 makie.jl 讀法相容：inline base64、三角單元、欄位 v₁…v₁₂）
# ---------------------------------------------------------------------------
function write_vtu(path::String, modes::Vector{Matrix{Float64}})
    mkpath(dirname(path))
    pts = zeros(2, NP * NP)
    for ix in 1:NP, iy in 1:NP
        p = (ix - 1) * NP + iy
        pts[1, p] = XS[ix] - 0.5      # 移到 [-0.5, 0.5]
        pts[2, p] = XS[iy] - 0.5
    end
    cells = Vector{MeshCell}()
    for ix in 1:(NP - 1), iy in 1:(NP - 1)
        p1 = (ix - 1) * NP + iy
        p2 = ix * NP + iy
        p3 = ix * NP + iy + 1
        p4 = (ix - 1) * NP + iy + 1
        push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, [p1, p2, p3]))
        push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, [p1, p3, p4]))
    end
    vtk_grid(path, pts, cells; ascii = false, append = false, compress = false) do vtk
        for j in 1:NMODE
            # 注意：點序為 p=(ix-1)*NP+iy，而 modes[j] 係 W[ix,iy]（第一維 = x）。
            # vec(W) 會將值按 (iy 慢、ix 快) 排，同點序相反 → 必須 permutedims 才對位。
            vtk[VFIELDS[j]] = vec(permutedims(modes[j]))
        end
    end
    return path
end

# ---------------------------------------------------------------------------
#  主流程
# ---------------------------------------------------------------------------
function report(bc::String, ks::Vector{Float64}, wave::Vector{String})
    @printf("%-5s 前 12 階 k = ", bc)
    for v in ks
        @printf("%9.4f", v)
    end
    println()
    @printf("      (波數)     ")
    for v in wave
        @printf("%9s", v)
    end
    println()
    e = REF[bc]
    rel = abs(ks[1] - e) / e * 100
    @printf("      最低階 %.4f   經典參考值 %.4f   相對差 %.4f%%   %s\n\n",
            ks[1], e, rel, rel < 0.05 ? "OK" : "**待核**")
end

function main()
    println("="^100)
    println("方板單軸屈曲標準解（公式法）：ν = $(NU)，a = b = 1，k = N_x b²/(π²D)")
    println("參考值來源：Timoshenko & Gere (1961) / Bui et al. (2011)（論文 CH4 表 CH4_compare_table 所用）")
    println("="^100)

    L = Dict{String,Any}()
    L["SSSS"] = (:closed, ssss_modes())
    for (bc, b0, b1) in (("CSCS", :C, :C), ("FSCS", :F, :C), ("FSSS", :F, :S))
        lst = levy_modes(:x, b0, b1)
        # 統一標度：SVD 零空間向量係任意標度，必須歸一（max|w| = 1）才與 SSSS/CCCC 同色標可比
        ms = [(k, w, normalize!(grid_from_profile(:x, w, levy_profile(:x, k, w, b0, b1))))
              for (k, w) in lst[1:NMODE]]
        L[bc] = (:levy_x, ms)
    end
    lst = levy_modes(:y, :C, :C)
    L["SCSC"] = (:levy_y, [(k, w, normalize!(grid_from_profile(:y, w, levy_profile(:y, k, w, :C, :C))))
                           for (k, w) in lst[1:NMODE]])
    L["CCCC"] = (:ritz, [(k, 0, m) for (k, m) in cccc_modes()])

    # 自檢：SS 樑基 Ritz 應重現 SSSS 閉式 4.0000 / 6.2500 / 11.1111
    ks_ss, _ = ritz_plate(:S, 10)
    @printf("自檢 SS 樑基 Ritz：k = %.4f %.4f %.4f %.4f  （閉式 4.0000 6.2500 11.1111 16.0000）\n\n",
            ks_ss[1], ks_ss[2], ks_ss[3], ks_ss[4])

    for bc in ("CCCC", "SSSS", "CSCS", "SCSC", "FSCS", "FSSS")
        tag, lst = L[bc]
        ks = [Float64(x[1]) for x in lst]
        wv = tag == :ritz ? fill("—", NMODE) :
             tag == :closed ? [x[2] for x in lst] : [string(Int(x[2])) for x in lst]
        report(bc, ks[1:NMODE], wv)
        modes = [x[3] for x in lst[1:NMODE]]
        out = joinpath(@__DIR__, "..", "vtk", "std", "buckling_n20", bc, "std_$(bc).vtu")
        write_vtu(out, modes)
        println("      → VTU: $out\n")
    end
    println("="^100)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
