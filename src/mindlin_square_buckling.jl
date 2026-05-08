using ApproxOperator
using LinearAlgebra
using Printf

import Gmsh: gmsh
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ

# Mindlin-Reissner 板屈曲 benchmark 腳本。
#
# 解析解只使用 Timoshenko/Gere《Theory of Elastic Stability》中文版掃描本
# 第 9 章中的公式：
# - T1(1).pdf 第 1 頁：板屈曲控制方程 (9-1)。
# - T1(1).pdf 第 6 頁：四邊簡支矩形板受壓公式 (9-7) 與表 9-1。
#
# 數值模型：
# - 用 Gmsh 生成矩形板網格，並使用 ASCII physical groups：
#   left/right/top/bottom/domain。
# - 用 ApproxOperator.MindlinPlate 組裝與原始腳本一致的 Mindlin 分塊矩陣：
#     K = [Kphiphi  Kphiw
#          Kwphi    Kww]
#     KG 只作用在橫向位移 w 上。
# - 因為 KG 在轉角自由度區塊為奇異矩陣，求解屈曲特徵值前先靜力凝聚轉角自由度。

const DEFAULT_E = 1.0e8
const DEFAULT_NU = 0.3
const DEFAULT_B = 1.0
const DEFAULT_MESH_N = 12
const EXACT_SEARCH_MODES = 60

struct BucklingCase
    name::String
    load::Symbol
    boundary::String
    aspect::Float64
    gamma::Float64
    mesh_n::Int
end

struct LoadState
    σ11::Float64
    σ22::Float64
    σ12::Float64
end

function bending_rigidity(E::Float64, ν::Float64, h::Float64)
    return E*h^3/(12.0*(1.0 - ν^2))
end

function load_state(case::BucklingCase)
    if case.load == :uniaxial
        return LoadState(1.0, 0.0, 0.0)
    elseif case.load == :biaxial
        return LoadState(1.0, case.gamma, 0.0)
    elseif case.load == :shear
        return LoadState(0.0, 0.0, 1.0)
    else
        error("unsupported load case: $(case.load)")
    end
end

function exact_uniaxial_ssss(aspect::Float64; max_mode::Int=EXACT_SEARCH_MODES)
    return minimum((m/aspect + aspect/m)^2 for m in 1:max_mode)
end

function exact_biaxial_ssss(aspect::Float64, gamma::Float64; max_mode::Int=EXACT_SEARCH_MODES)
    values = Float64[]
    for m in 1:max_mode, n in 1:max_mode
        mx = m^2/aspect^2
        ny = n^2
        denom = mx + gamma*ny
        denom > 0.0 || continue
        push!(values, (mx + ny)^2/denom)
    end
    isempty(values) && error("no exact biaxial value found")
    return minimum(values)
end

function exact_k(case::BucklingCase)
    if case.boundary != "SSSS"
        return missing
    elseif case.load == :uniaxial
        return exact_uniaxial_ssss(case.aspect)
    elseif case.load == :biaxial
        return exact_biaxial_ssss(case.aspect, case.gamma)
    else
        return missing
    end
end

function generate_rect_mesh!(case::BucklingCase)
    a = case.aspect*DEFAULT_B
    b = DEFAULT_B
    nx = max(round(Int, case.mesh_n*case.aspect), 1) + 1
    ny = max(case.mesh_n, 1) + 1

    gmsh.clear()
    gmsh.model.add("mindlin_rect")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(a, 0.0, 0.0)
    p3 = gmsh.model.geo.addPoint(a, b, 0.0)
    p4 = gmsh.model.geo.addPoint(0.0, b, 0.0)

    bottom = gmsh.model.geo.addLine(p1, p2)
    right = gmsh.model.geo.addLine(p2, p3)
    top = gmsh.model.geo.addLine(p3, p4)
    left = gmsh.model.geo.addLine(p4, p1)
    loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(left, ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(right, ny)
    gmsh.model.geo.mesh.setTransfiniteSurface(surface)
    gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()

    left_group = gmsh.model.addPhysicalGroup(1, [left])
    gmsh.model.setPhysicalName(1, left_group, "left")
    right_group = gmsh.model.addPhysicalGroup(1, [right])
    gmsh.model.setPhysicalName(1, right_group, "right")
    top_group = gmsh.model.addPhysicalGroup(1, [top])
    gmsh.model.setPhysicalName(1, top_group, "top")
    bottom_group = gmsh.model.addPhysicalGroup(1, [bottom])
    gmsh.model.setPhysicalName(1, bottom_group, "bottom")
    domain_group = gmsh.model.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, domain_group, "domain")

    gmsh.model.mesh.generate(2)
end

function constrained_nodes(nodes, entities, boundary::String)
    length(boundary) == 4 || error("boundary must have four letters: left/right/top/bottom")
    side_names = ("left", "right", "top", "bottom")
    fixed_w = Set{Int}()
    fixed_phi = Set{Int}()

    for (bc, side) in zip(collect(boundary), side_names)
        bc == 'F' && continue
        elements = getElements(nodes, entities[side])
        side_nodes = collect(node.𝐼 for element in elements for node in element.𝓒)
        if bc == 'S'
            union!(fixed_w, side_nodes)
        elseif bc == 'C'
            union!(fixed_w, side_nodes)
            union!(fixed_phi, side_nodes)
        else
            error("unsupported boundary code '$bc' in $boundary")
        end
    end

    return sort!(collect(fixed_w)), sort!(collect(fixed_phi))
end

function positive_finite_eigenvalues(k::AbstractMatrix, kg::AbstractMatrix)
    λ = eigvals(k, kg)
    values = Float64[]
    for λᵢ in λ
        if isfinite(real(λᵢ)) && isfinite(imag(λᵢ)) && abs(imag(λᵢ)) < 1.0e-7 && real(λᵢ) > 0.0
            push!(values, real(λᵢ))
        end
    end
    return sort!(values)
end

function assemble_mindlin_case(case::BucklingCase, h_over_b::Float64; E::Float64=DEFAULT_E, ν::Float64=DEFAULT_NU)
    generate_rect_mesh!(case)

    h = h_over_b*DEFAULT_B
    state = load_state(case)
    entities = getPhysicalGroups()
    for name in ("domain", "left", "right", "top", "bottom")
        haskey(entities, name) || error("missing physical group: $name")
    end

    nodes = get𝑿ᵢ()
    nʷ = length(nodes)
    nᵠ = length(nodes)

    kʷʷ = zeros(nʷ, nʷ)
    kᵠᵠ = zeros(2*nᵠ, 2*nᵠ)
    kᵠʷ = zeros(2*nᵠ, nʷ)
    kgʷʷ = zeros(nʷ, nʷ)

    elements = getElements(nodes, entities["domain"])
    prescribe!(elements, :E => E, :ν => ν, :h => h)
    set∇𝝭!(elements)
    (∫wwdΩ => elements)(kʷʷ)
    (∫φwdΩ => elements)(kᵠʷ)
    ([∫φφdΩ => elements, ∫κκdΩ => elements])(kᵠᵠ)

    σ₁₁ref(x, y, z) = state.σ11
    σ₂₂ref(x, y, z) = state.σ22
    σ₁₂ref(x, y, z) = state.σ12
    prescribe!(elements, :σ₁₁ => σ₁₁ref, :σ₂₂ => σ₂₂ref, :σ₁₂ => σ₁₂ref)
    (∫wwGdΩ2D => elements)(kgʷʷ)

    K = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ]
    KG = [
        zeros(2*nᵠ, 2*nᵠ) zeros(2*nᵠ, nʷ)
        zeros(nʷ, 2*nᵠ) kgʷʷ
    ]

    fixed_w_nodes, fixed_phi_nodes = constrained_nodes(nodes, entities, case.boundary)
    fixed_phi_dofs = reduce(vcat, ([2*i - 1, 2*i] for i in fixed_phi_nodes), init=Int[])
    fixed_w_dofs = 2*nᵠ .+ fixed_w_nodes
    fixed_dofs = sort!(unique!(vcat(fixed_phi_dofs, fixed_w_dofs)))
    free_dofs = setdiff(1:size(K, 1), fixed_dofs)

    return K, KG, free_dofs, h, length(nodes), length(fixed_dofs), 2*nᵠ
end

function solve_buckling(case::BucklingCase, h_over_b::Float64; E::Float64=DEFAULT_E, ν::Float64=DEFAULT_NU)
    K, KG, free_dofs, h, nnodes, nfixed, nphi_dofs = assemble_mindlin_case(case, h_over_b; E=E, ν=ν)
    free_phi = [dof for dof in free_dofs if dof <= nphi_dofs]
    free_w = [dof for dof in free_dofs if dof > nphi_dofs]
    isempty(free_w) && error("all transverse displacement dofs are constrained for $(case.name)")

    # KG 只作用在 w 上。先靜力凝聚轉角自由度，避免直接求解含有零 KG
    # 轉角區塊的奇異廣義特徵值問題。
    Kpp = K[free_phi, free_phi]
    Kpw = K[free_phi, free_w]
    Kwp = K[free_w, free_phi]
    Kww = K[free_w, free_w]
    KGww = KG[free_w, free_w]
    Keff = isempty(free_phi) ? Kww : Kww - Kwp*(Kpp\Kpw)

    λ = positive_finite_eigenvalues(Keff, KGww)
    isempty(λ) && error("no positive finite eigenvalue found for $(case.name), h/b=$h_over_b")

    λcr = first(λ)
    D = bending_rigidity(E, ν, h)
    k_num = λcr*DEFAULT_B^2/(π^2*D)
    return (; λcr, k_num, nnodes, nfixed)
end

function benchmark_cases()
    cases = BucklingCase[]

    # T1(1).pdf 第 6 頁，公式 (9-7) 與表 9-1：
    # 四邊簡支矩形板承受均勻單軸壓縮。
    for aspect in (1.0, 1.5, 2.0, 3.0, 4.0)
        push!(cases, BucklingCase("book-9-1 uniaxial SSSS a/b=$(aspect)", :uniaxial, "SSSS", aspect, 0.0, DEFAULT_MESH_N))
    end

    # T1(1).pdf 第 1 頁，方程 (9-1)，其中 Ny = gamma*Nx 且 Nxy = 0。
    for aspect in (1.0, 1.5, 2.0, 3.0, 4.0)
        push!(cases, BucklingCase("book-9-1 biaxial SSSS gamma=1 a/b=$(aspect)", :biaxial, "SSSS", aspect, 1.0, DEFAULT_MESH_N))
    end

    # 保留這個項目是為了讓輸出範圍明確。三份掃描本中尚未轉錄出
    # 四邊簡支矩形板純剪切屈曲的明確係數。
    push!(cases, BucklingCase("book scan pure shear rectangular", :shear, "SSSS", 1.0, 0.0, DEFAULT_MESH_N))

    return cases
end

function print_header()
    println("Mindlin-Reissner buckling vs. Timoshenko/Gere thin-plate formulas")
    println("Exact source: Chinese Theory of Elastic Stability scans, chapter 9.")
    println("Mesh source: Gmsh generated rectangle with physical groups left/right/top/bottom/domain.")
    println("Normalization: k = lambda*b^2/(pi^2*D), D = E*h^3/(12*(1-nu^2)).")
    println()
    @printf("%-46s %7s %8s %8s %12s %12s %12s %8s %8s\n",
        "case", "h/b", "a/b", "bc", "k_exact", "k_num", "rel_error", "nodes", "fixed")
    println("-"^136)
end

function run_benchmarks()
    h_over_b_values = (0.001, 0.01, 0.05, 0.1)
    case_limit = parse(Int, get(ENV, "MINDLIN_CASE_LIMIT", "0"))
    thickness_limit = parse(Int, get(ENV, "MINDLIN_THICKNESS_LIMIT", "0"))

    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        cases = benchmark_cases()
        case_limit > 0 && (cases = cases[1:min(case_limit, length(cases))])
        thickness_limit > 0 && (h_over_b_values = h_over_b_values[1:min(thickness_limit, length(h_over_b_values))])

        print_header()
        for case in cases
            k_exact = exact_k(case)
            if ismissing(k_exact)
                @printf("%-46s %7s %8.3f %8s %12s %12s %12s %8s %8s\n",
                    case.name, "-", case.aspect, case.boundary,
                    "skipped", "skipped", "no book formula/table", "-", "-")
                continue
            end
            for h_over_b in h_over_b_values
                result = solve_buckling(case, h_over_b)
                rel_error = abs(result.k_num - k_exact)/abs(k_exact)
                @printf("%-46s %7.4f %8.3f %8s %12.6f %12.6f %12.4e %8d %8d\n",
                    case.name, h_over_b, case.aspect, case.boundary,
                    k_exact, result.k_num, rel_error, result.nnodes, result.nfixed)
            end
        end
    finally
        gmsh.finalize()
    end
end

run_benchmarks()
