using ApproxOperator
import ApproxOperator: Tri3toTri6!, Seg2toSeg3!
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫wwGdΩ2D, ∫κκdΩBui, ∫wwdΩ, ∫ψψdΩBui, ∫ψwdΩBui, ∫ψxψxGdΩ2DBui, ∫ψyψyGdΩ2DBui, ∫αwwdΓ

using LinearAlgebra
using Printf
using TimerOutputs
using XLSX
import Gmsh: gmsh

E = 30.0e6
ν = 0.3
b = 1.0
h_over_b = 0.01
h = h_over_b*b
Dᵇ = E*h^3/12/(1-ν^2)

# T1 表 9-1：均勻受壓簡支矩形板，E = 30e6 psi, ν = 0.3, h/b = 0.01。
# 後兩欄是書上列印值，計算時另用解析公式得到更多位數。
table_9_1 = [
    (0.20, 27.00, 73200.0),
    (0.30, 13.20, 35800.0),
    (0.40,  8.41, 22800.0),
    (0.50,  6.25, 16900.0),
    (0.60,  5.14, 13900.0),
    (0.70,  4.53, 12300.0),
    (0.80,  4.20, 11400.0),
    (0.90,  4.04, 11000.0),
    (1.00,  4.00, 10800.0),
    (1.10,  4.04, 11000.0),
    (1.20,  4.13, 11200.0),
    (1.30,  4.28, 11600.0),
    (1.40,  4.47, 12100.0),
    (1.41,  4.49, 12200.0),
]

function exact_uniaxial_k(r; mmax=80)
    k_exact = Inf
    m_exact = 0
    for m in 1:mmax
        k = (m/r + r/m)^2
        if k < k_exact
            k_exact = k
            m_exact = m
        end
    end
    return k_exact, m_exact
end

exact_ξcr(k_exact) = k_exact*π^2*Dᵇ/b^2
exact_σcr(k_exact) = exact_ξcr(k_exact)/h

w(x,y,z) = 0.0
σ₁₁(x,y,z) = 1.0
σ₂₂(x,y,z) = 0.0
σ₁₂(x,y,z) = 0.0

const to = TimerOutput()
const output_xlsx = normpath(joinpath(@__DIR__, "..", "results", "mindlin_uniaxial_ssss_table_9_1_k_nullspace_diagnostic.xlsx"))
const eigen_imag_tol = 1.0e-7
const spectrum_first_positive_count = 40
const spectrum_closest_exact_count = 10
const spectrum_selection_mode = "closest_exact"
const k_nullspace_count = 40
const k_nullspace_near_zero_tol = 1.0e-10

function aspect_tag(r)
    return replace(@sprintf("%.2f", r), "."=>"p")
end

function mesh_file(r)
    return normpath(joinpath(@__DIR__, "..", "msh", "mindlin_uniaxial_ssss_ab_$(aspect_tag(r)).msh"))
end

function excel_column_name(j)
    name = ""
    while j > 0
        j, r = divrem(j - 1, 26)
        name = string(Char('A' + r), name)
    end
    return name
end

function excel_cell(row, col)
    return string(excel_column_name(col), row)
end

function write_sheet_rows!(sheet, headers, rows)
    for (j, header) in enumerate(headers)
        sheet[excel_cell(1, j)] = header
    end
    for (i, row) in enumerate(rows)
        for (j, value) in enumerate(row)
            sheet[excel_cell(i + 1, j)] = value
        end
    end
end

function table_9_1_values(row)
    return [
        row.r,
        row.m_exact,
        row.mesh,
        row.kg_strategy,
        row.selection_mode,
        row.mode_index,
        row.rank_by_xi,
        row.xi_cr,
        row.xi_exact,
        row.xi_ratio,
        row.log10_xi_ratio,
        row.w_participation,
        row.k_num,
        row.k_exact,
        row.k_book,
        row.k_err,
        row.σcr_num,
        row.σcr_exact,
        row.σcr_book,
        row.nnodes,
    ]
end

function spectrum_values(row)
    return [
        row.r,
        row.mesh,
        row.kg_strategy,
        row.spectrum_group,
        row.group_rank,
        row.mode_index,
        row.rank_by_xi,
        row.xi,
        row.xi_exact,
        row.xi_ratio,
        row.log10_xi_ratio,
        row.distance_to_exact,
        row.w_participation,
        row.k,
        row.k_exact,
    ]
end

function block_values(row)
    return [
        row.r,
        row.mesh,
        row.kg_strategy,
        row.mode_index,
        row.rank_by_xi,
        row.xi_selected,
        row.xi_exact,
        row.xi_ratio,
        row.norm_kww,
        row.norm_kqq,
        row.norm_kqw,
        row.norm_kgww,
        row.active_norm_kgqq,
        row.assembled_norm_kgqq,
        row.ratio_kqq_kww,
        row.ratio_kqw_geom_mean,
        row.active_ratio_kgqq_kgww,
        row.assembled_ratio_kgqq_kgww,
        row.energy_kww,
        row.energy_kqq,
        row.energy_kqw,
        row.energy_K_total,
        row.energy_kgww,
        row.active_energy_kgqq,
        row.assembled_energy_kgqq,
        row.active_energy_KG_total,
        row.assembled_energy_KG_total,
        row.rayleigh_xi_active,
        row.w_participation,
    ]
end

function kqq_nullspace_values(row)
    return [
        row.r,
        row.mesh,
        row.mode_rank,
        row.lambda_kqq,
        row.abs_lambda_kqq,
        row.lambda_ratio_to_max,
        row.is_near_zero,
    ]
end

function k_nullspace_values(row)
    return [
        row.r,
        row.mesh,
        row.mode_rank,
        row.lambda_K,
        row.abs_lambda_K,
        row.lambda_ratio_to_max,
        row.w_participation,
        row.psi_participation,
        row.is_near_zero,
    ]
end

function selected_mode_energy_values(row)
    return [
        row.r,
        row.mesh,
        row.kg_strategy,
        row.mode_index,
        row.rank_by_xi,
        row.psi_norm_ratio,
        row.w_norm_ratio,
        row.energy_kqq,
        row.energy_kqw,
        row.energy_kww,
        row.energy_K_total,
        row.abs_energy_kqq_ratio,
        row.abs_energy_kqw_ratio,
        row.abs_energy_kww_ratio,
    ]
end

safe_ratio(numerator, denominator) = denominator == 0.0 ? NaN : numerator/denominator

function k_nullspace_diagnostics(r, kᵠᵠ, kᵠʷ, kʷʷ, nᵠ, nʷ)
    mesh_name = basename(mesh_file(r))

    E_kqq = eigen(Symmetric(kᵠᵠ))
    λ_kqq = real.(E_kqq.values)
    λ_kqq_max = maximum(abs.(λ_kqq))
    λ_kqq_tol = k_nullspace_near_zero_tol*λ_kqq_max
    kqq_order = sortperm(abs.(λ_kqq))
    kqq_near_zero_count = count(abs(λ) <= λ_kqq_tol for λ in λ_kqq)
    kqq_rows = NamedTuple[]
    for (mode_rank, idx) in enumerate(first(kqq_order, min(length(kqq_order), k_nullspace_count)))
        λ = λ_kqq[idx]
        push!(kqq_rows, (
            r = r,
            mesh = mesh_name,
            mode_rank = mode_rank,
            lambda_kqq = λ,
            abs_lambda_kqq = abs(λ),
            lambda_ratio_to_max = safe_ratio(abs(λ), λ_kqq_max),
            is_near_zero = abs(λ) <= λ_kqq_tol,
        ))
    end

    K = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
    E_K = eigen(Symmetric(K))
    λ_K = real.(E_K.values)
    λ_K_max = maximum(abs.(λ_K))
    λ_K_tol = k_nullspace_near_zero_tol*λ_K_max
    k_order = sortperm(abs.(λ_K))
    k_near_zero_count = count(abs(λ) <= λ_K_tol for λ in λ_K)
    k_rows = NamedTuple[]
    for (mode_rank, idx) in enumerate(first(k_order, min(length(k_order), k_nullspace_count)))
        λ = λ_K[idx]
        mode = @view E_K.vectors[:, idx]
        mode_norm = sum(abs2, mode)
        psi_participation = mode_norm <= 0.0 ? NaN : sum(abs2, @view mode[1:2*nᵠ])/mode_norm
        w_participation = mode_norm <= 0.0 ? NaN : sum(abs2, @view mode[2*nᵠ+1:2*nᵠ+nʷ])/mode_norm
        push!(k_rows, (
            r = r,
            mesh = mesh_name,
            mode_rank = mode_rank,
            lambda_K = λ,
            abs_lambda_K = abs(λ),
            lambda_ratio_to_max = safe_ratio(abs(λ), λ_K_max),
            w_participation = w_participation,
            psi_participation = psi_participation,
            is_near_zero = abs(λ) <= λ_K_tol,
        ))
    end

    return (
        kqq_rows = kqq_rows,
        k_rows = k_rows,
        kqq_near_zero_count = kqq_near_zero_count,
        k_near_zero_count = k_near_zero_count,
    )
end

function selected_mode_energy_diagnostics(
    r, kg_strategy, selected_mode, F, nᵠ, nʷ,
    kʷʷ, kᵠᵠ, kᵠʷ,
)
    v = @view F.vectors[:, selected_mode.mode_index]
    vᵠ = @view v[1:2*nᵠ]
    vʷ = @view v[2*nᵠ+1:2*nᵠ+nʷ]

    mode_norm = sum(abs2, v)
    psi_norm_ratio = mode_norm <= 0.0 ? NaN : sum(abs2, vᵠ)/mode_norm
    w_norm_ratio = mode_norm <= 0.0 ? NaN : sum(abs2, vʷ)/mode_norm

    energy_kqq = real(dot(vᵠ, kᵠᵠ*vᵠ))
    energy_kqw = 2.0*real(dot(vᵠ, kᵠʷ*vʷ))
    energy_kww = real(dot(vʷ, kʷʷ*vʷ))
    energy_K_total = energy_kqq + energy_kqw + energy_kww
    abs_energy_total = abs(energy_kqq) + abs(energy_kqw) + abs(energy_kww)

    return (
        r = r,
        mesh = basename(mesh_file(r)),
        kg_strategy = kg_strategy,
        mode_index = selected_mode.mode_index,
        rank_by_xi = selected_mode.rank_by_xi,
        psi_norm_ratio = psi_norm_ratio,
        w_norm_ratio = w_norm_ratio,
        energy_kqq = energy_kqq,
        energy_kqw = energy_kqw,
        energy_kww = energy_kww,
        energy_K_total = energy_K_total,
        abs_energy_kqq_ratio = safe_ratio(abs(energy_kqq), abs_energy_total),
        abs_energy_kqw_ratio = safe_ratio(abs(energy_kqw), abs_energy_total),
        abs_energy_kww_ratio = safe_ratio(abs(energy_kww), abs_energy_total),
    )
end

function solve_spectrum(K, KG, r, k_exact, kg_strategy, nᵠ)
    F = eigen(K, KG)
    ξ_exact = exact_ξcr(k_exact)
    candidates_raw = NamedTuple[]
    for (i, ξᵢ) in pairs(F.values)
        ξᵣ = real(ξᵢ)
        ξᵢₘ = imag(ξᵢ)
        if !(isfinite(ξᵣ) && isfinite(ξᵢₘ) && abs(ξᵢₘ) < eigen_imag_tol && ξᵣ > 0.0)
            continue
        end
        mode = @view F.vectors[:, i]
        mode_norm = sum(abs2, mode)
        if !(isfinite(mode_norm) && mode_norm > 0.0)
            continue
        end
        w_participation = sum(abs2, @view mode[2*nᵠ+1:end])/mode_norm
        k = ξᵣ*b^2/(π^2*Dᵇ)
        ξ_ratio = ξᵣ/ξ_exact
        log10_ξ_ratio = log10(abs(ξ_ratio))
        distance_to_exact = abs(log10_ξ_ratio)
        push!(candidates_raw, (
            mode_index = i,
            xi = ξᵣ,
            k = k,
            xi_ratio = ξ_ratio,
            log10_xi_ratio = log10_ξ_ratio,
            distance_to_exact = distance_to_exact,
            w_participation = w_participation,
        ))
    end
    sort!(candidates_raw, by = c -> c.xi)
    isempty(candidates_raw) && error("no positive finite buckling eigenvalue found for a/b = $r, kg_strategy = $kg_strategy")
    candidates = NamedTuple[]
    for (rank_by_xi, candidate) in enumerate(candidates_raw)
        push!(candidates, (
            mode_index = candidate.mode_index,
            rank_by_xi = rank_by_xi,
            xi = candidate.xi,
            k = candidate.k,
            xi_ratio = candidate.xi_ratio,
            log10_xi_ratio = candidate.log10_xi_ratio,
            distance_to_exact = candidate.distance_to_exact,
            w_participation = candidate.w_participation,
        ))
    end

    closest_candidates = sort(candidates, by = c -> c.distance_to_exact)
    selected_mode = first(closest_candidates)
    first_positive_modes = first(candidates, min(length(candidates), spectrum_first_positive_count))
    closest_exact_modes = first(closest_candidates, min(length(closest_candidates), spectrum_closest_exact_count))
    mesh_name = basename(mesh_file(r))
    spectrum_rows = NamedTuple[]
    for (group_rank, mode) in enumerate(first_positive_modes)
        push!(spectrum_rows, (
            r = r,
            mesh = mesh_name,
            kg_strategy = kg_strategy,
            spectrum_group = "first_positive",
            group_rank = group_rank,
            mode_index = mode.mode_index,
            rank_by_xi = mode.rank_by_xi,
            xi = mode.xi,
            xi_exact = ξ_exact,
            xi_ratio = mode.xi_ratio,
            log10_xi_ratio = mode.log10_xi_ratio,
            distance_to_exact = mode.distance_to_exact,
            w_participation = mode.w_participation,
            k = mode.k,
            k_exact = k_exact,
        ))
    end
    for (group_rank, mode) in enumerate(closest_exact_modes)
        push!(spectrum_rows, (
            r = r,
            mesh = mesh_name,
            kg_strategy = kg_strategy,
            spectrum_group = "closest_exact",
            group_rank = group_rank,
            mode_index = mode.mode_index,
            rank_by_xi = mode.rank_by_xi,
            xi = mode.xi,
            xi_exact = ξ_exact,
            xi_ratio = mode.xi_ratio,
            log10_xi_ratio = mode.log10_xi_ratio,
            distance_to_exact = mode.distance_to_exact,
            w_participation = mode.w_participation,
            k = mode.k,
            k_exact = k_exact,
        ))
    end

    return (
        F = F,
        selected_mode = selected_mode,
        spectrum_rows = spectrum_rows,
        xi_cr = selected_mode.xi,
        k_num = selected_mode.k,
        σcr_num = selected_mode.xi/h,
        w_participation = selected_mode.w_participation,
    )
end

function matrix_block_diagnostics(
    r, k_exact, kg_strategy, selected_mode, F, nᵠ, nʷ,
    kʷʷ, kᵠᵠ, kᵠʷ, kgʷʷ, active_kgᵠᵠ, assembled_kgᵠᵠ,
)
    v = @view F.vectors[:, selected_mode.mode_index]
    vᵠ = @view v[1:2*nᵠ]
    vʷ = @view v[2*nᵠ+1:2*nᵠ+nʷ]

    norm_kww = norm(kʷʷ)
    norm_kqq = norm(kᵠᵠ)
    norm_kqw = norm(kᵠʷ)
    norm_kgww = norm(kgʷʷ)
    active_norm_kgqq = norm(active_kgᵠᵠ)
    assembled_norm_kgqq = norm(assembled_kgᵠᵠ)

    energy_kww = real(dot(vʷ, kʷʷ*vʷ))
    energy_kqq = real(dot(vᵠ, kᵠᵠ*vᵠ))
    energy_kqw = 2.0*real(dot(vᵠ, kᵠʷ*vʷ))
    energy_K_total = energy_kww + energy_kqq + energy_kqw
    energy_kgww = real(dot(vʷ, kgʷʷ*vʷ))
    active_energy_kgqq = real(dot(vᵠ, active_kgᵠᵠ*vᵠ))
    assembled_energy_kgqq = real(dot(vᵠ, assembled_kgᵠᵠ*vᵠ))
    active_energy_KG_total = energy_kgww + active_energy_kgqq
    assembled_energy_KG_total = energy_kgww + assembled_energy_kgqq

    return (
        r = r,
        mesh = basename(mesh_file(r)),
        kg_strategy = kg_strategy,
        mode_index = selected_mode.mode_index,
        rank_by_xi = selected_mode.rank_by_xi,
        xi_selected = selected_mode.xi,
        xi_exact = exact_ξcr(k_exact),
        xi_ratio = selected_mode.xi_ratio,
        norm_kww = norm_kww,
        norm_kqq = norm_kqq,
        norm_kqw = norm_kqw,
        norm_kgww = norm_kgww,
        active_norm_kgqq = active_norm_kgqq,
        assembled_norm_kgqq = assembled_norm_kgqq,
        ratio_kqq_kww = safe_ratio(norm_kqq, norm_kww),
        ratio_kqw_geom_mean = norm_kqq*norm_kww <= 0.0 ? NaN : norm_kqw/sqrt(norm_kqq*norm_kww),
        active_ratio_kgqq_kgww = safe_ratio(active_norm_kgqq, norm_kgww),
        assembled_ratio_kgqq_kgww = safe_ratio(assembled_norm_kgqq, norm_kgww),
        energy_kww = energy_kww,
        energy_kqq = energy_kqq,
        energy_kqw = energy_kqw,
        energy_K_total = energy_K_total,
        energy_kgww = energy_kgww,
        active_energy_kgqq = active_energy_kgqq,
        assembled_energy_kgqq = assembled_energy_kgqq,
        active_energy_KG_total = active_energy_KG_total,
        assembled_energy_KG_total = assembled_energy_KG_total,
        rayleigh_xi_active = safe_ratio(energy_K_total, active_energy_KG_total),
        w_participation = selected_mode.w_participation,
    )
end

function write_table_9_1_xlsx(
    filepath, rows, spectrum_rows, block_rows,
    kqq_nullspace_rows, k_nullspace_rows, selected_mode_energy_rows,
)
    mkpath(dirname(filepath))
    table_headers = [
        "a/b", "m", "mesh", "kg_strategy", "selection_mode", "mode_index", "rank_by_xi",
        "xi_cr", "xi_exact", "xi_ratio", "log10_xi_ratio", "w_participation",
        "k_num", "k_exact", "k_book", "k_error", "sigma_cr_num",
        "sigma_cr_exact", "sigma_cr_book", "nodes",
    ]
    spectrum_headers = [
        "a/b", "mesh", "kg_strategy", "spectrum_group", "group_rank", "mode_index", "rank_by_xi",
        "xi", "xi_exact", "xi_ratio", "log10_xi_ratio", "distance_to_exact",
        "w_participation", "k", "k_exact",
    ]
    block_headers = [
        "a/b", "mesh", "kg_strategy", "mode_index", "rank_by_xi", "xi_selected", "xi_exact",
        "xi_ratio", "norm_kww", "norm_kqq", "norm_kqw", "norm_kgww",
        "active_norm_kgqq", "assembled_norm_kgqq", "ratio_kqq_kww",
        "ratio_kqw_geom_mean", "active_ratio_kgqq_kgww",
        "assembled_ratio_kgqq_kgww", "energy_kww", "energy_kqq",
        "energy_kqw", "energy_K_total", "energy_kgww",
        "active_energy_kgqq", "assembled_energy_kgqq",
        "active_energy_KG_total", "assembled_energy_KG_total",
        "rayleigh_xi_active", "w_participation",
    ]
    kqq_nullspace_headers = [
        "a/b", "mesh", "mode_rank", "lambda_kqq", "abs_lambda_kqq",
        "lambda_ratio_to_max", "is_near_zero",
    ]
    k_nullspace_headers = [
        "a/b", "mesh", "mode_rank", "lambda_K", "abs_lambda_K",
        "lambda_ratio_to_max", "w_participation", "psi_participation",
        "is_near_zero",
    ]
    selected_mode_energy_headers = [
        "a/b", "mesh", "kg_strategy", "mode_index", "rank_by_xi",
        "psi_norm_ratio", "w_norm_ratio", "energy_kqq", "energy_kqw",
        "energy_kww", "energy_K_total", "abs_energy_kqq_ratio",
        "abs_energy_kqw_ratio", "abs_energy_kww_ratio",
    ]

    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Table 9-1")
        write_sheet_rows!(sheet, table_headers, table_9_1_values.(rows))

        spectrum_sheet = XLSX.addsheet!(xf, "Spectrum")
        write_sheet_rows!(spectrum_sheet, spectrum_headers, spectrum_values.(spectrum_rows))

        block_sheet = XLSX.addsheet!(xf, "MatrixBlocks")
        write_sheet_rows!(block_sheet, block_headers, block_values.(block_rows))

        kqq_sheet = XLSX.addsheet!(xf, "KqqNullspace")
        write_sheet_rows!(kqq_sheet, kqq_nullspace_headers, kqq_nullspace_values.(kqq_nullspace_rows))

        k_sheet = XLSX.addsheet!(xf, "KNullspace")
        write_sheet_rows!(k_sheet, k_nullspace_headers, k_nullspace_values.(k_nullspace_rows))

        selected_energy_sheet = XLSX.addsheet!(xf, "SelectedModeEnergy")
        write_sheet_rows!(
            selected_energy_sheet,
            selected_mode_energy_headers,
            selected_mode_energy_values.(selected_mode_energy_rows),
        )
    end
end

function solve_table_9_1_case(r, k_exact)
    gmsh.clear()
    @timeit to "open msh file" gmsh.open(mesh_file(r))
    @timeit to "get entities" entities = getPhysicalGroups()
    @timeit to "get nodes" nodes = get𝑿ᵢ()
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 3)
    @timeit to "convert Ω elements to Tri6" begin
        elements, midside_nodes = Tri3toTri6!(nodes, elements)
    end
    @timeit to "get boundary elements" begin
        elements_1 = getElements(nodes, entities["Γ¹"], 3)
        elements_2 = getElements(nodes, entities["Γ²"], 3)
        elements_3 = getElements(nodes, entities["Γ³"], 3)
        elements_4 = getElements(nodes, entities["Γ⁴"], 3)
    end
    @timeit to "convert Γ elements to Seg3" begin
        elements_1 = Seg2toSeg3!(nodes, elements_1, midside_nodes)
        elements_2 = Seg2toSeg3!(nodes, elements_2, midside_nodes)
        elements_3 = Seg2toSeg3!(nodes, elements_3, midside_nodes)
        elements_4 = Seg2toSeg3!(nodes, elements_4, midside_nodes)
    end

    nʷ = length(nodes)
    nᵠ = length(nodes)
    kʷʷ = zeros(nʷ,nʷ)
    kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kᵠʷ = zeros(2*nᵠ,nʷ)
    kgᵠᵠ = zeros(2*nᵠ,2*nᵠ)
    kgʷʷ = zeros(nʷ,nʷ)

    @timeit to "calculate ∫κκdΩBui, ∫wwdΩ, ∫ψψdΩBui, ∫ψwdΩBui, ∫wwGdΩ2D, ∫ψψGdΩ2DBui" begin
        prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
        @timeit to "calculate shape functions" set∇²𝝭!(elements)
        𝑎ʷʷ = ∫wwdΩ=>elements
        𝑎ᵠʷ = ∫ψwdΩBui=>elements
        𝑎ᵠᵠ = [
            ∫ψψdΩBui=>elements,
            ∫κκdΩBui=>elements,
        ]
        @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
        @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
        @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)

        prescribe!(elements, :σ₁₁=>σ₁₁, :σ₂₂=>σ₂₂, :σ₁₂=>σ₁₂)
        𝑎ᴳʷʷ = ∫wwGdΩ2D=>elements
        𝑎ᴳᵠˣᵠˣ = ∫ψxψxGdΩ2DBui=>elements
        𝑎ᴳᵠʸᵠʸ = ∫ψyψyGdΩ2DBui=>elements
        @timeit to "assemble" 𝑎ᴳʷʷ(kgʷʷ)
        @timeit to "assemble" 𝑎ᴳᵠˣᵠˣ(view(kgᵠᵠ, 1:2:2*nᵠ, 1:2:2*nᵠ))
        @timeit to "assemble" 𝑎ᴳᵠʸᵠʸ(view(kgᵠᵠ, 2:2:2*nᵠ, 2:2:2*nᵠ))
    end

    @timeit to "calculate ∫αwwdΓ" begin
        prescribe!(elements_1, :α=>1e8*E, :g=>w)
        prescribe!(elements_2, :α=>1e8*E, :g=>w)
        prescribe!(elements_3, :α=>1e8*E, :g=>w)
        prescribe!(elements_4, :α=>1e8*E, :g=>w)
        @timeit to "calculate shape functions" set𝝭!(elements_1)
        @timeit to "calculate shape functions" set𝝭!(elements_2)
        @timeit to "calculate shape functions" set𝝭!(elements_3)
        @timeit to "calculate shape functions" set𝝭!(elements_4)
        𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
        fᵅ = zeros(nʷ)
        @timeit to "assemble" 𝑎ʷ(kʷʷ,fᵅ)
    end

    strategy_results = NamedTuple[]
    k_diag = nothing
    @timeit to "solve buckling eigenvalue" begin
        K = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
        k_diag = k_nullspace_diagnostics(r, kᵠᵠ, kᵠʷ, kʷʷ, nᵠ, nʷ)
        zero_kgᵠᵠ = zeros(2*nᵠ,2*nᵠ)
        for (kg_strategy, active_kgᵠᵠ) in (
            ("full_bui", kgᵠᵠ),
            ("kgww_only", zero_kgᵠᵠ),
        )
            KG = [
                active_kgᵠᵠ zeros(2*nᵠ,nʷ)
                zeros(nʷ,2*nᵠ) kgʷʷ
            ]
            solved = solve_spectrum(K, KG, r, k_exact, kg_strategy, nᵠ)
            block_diag = matrix_block_diagnostics(
                r, k_exact, kg_strategy, solved.selected_mode, solved.F, nᵠ, nʷ,
                kʷʷ, kᵠᵠ, kᵠʷ, kgʷʷ, active_kgᵠᵠ, kgᵠᵠ,
            )
            selected_energy = selected_mode_energy_diagnostics(
                r, kg_strategy, solved.selected_mode, solved.F, nᵠ, nʷ,
                kʷʷ, kᵠᵠ, kᵠʷ,
            )
            push!(strategy_results, (
                kg_strategy = kg_strategy,
                xi_cr = solved.xi_cr,
                k_num = solved.k_num,
                σcr_num = solved.σcr_num,
                w_participation = solved.w_participation,
                selected_mode = solved.selected_mode,
                spectrum_rows = solved.spectrum_rows,
                block_diag = block_diag,
                selected_energy = selected_energy,
            ))
        end
    end

    return nʷ, strategy_results, k_diag
end

gmsh.initialize()
try
    println("T1 表 9-1: SSSS rectangular plate under uniaxial in-plane compression")
    println("E = 30e6 psi, ν = 0.3, h/b = 0.01")
    println("Required mesh names: msh/mindlin_uniaxial_ssss_ab_<a_over_b>.msh, e.g. ab_1p00.")
    results = []
    spectrum_results = []
    block_results = []
    kqq_nullspace_results = []
    k_nullspace_results = []
    selected_mode_energy_results = []
    @printf("%7s %3s %18s %11s %14s %10s %10s %12s %12s %12s %14s %14s %14s %12s %8s\n",
        "a/b", "m", "msh", "kg_strategy", "selection", "mode", "rank", "ξcr", "ξ_exact",
        "ξ/ξex", "w_part", "k_num", "k_exact", "k_err", "nodes")
    println("-"^188)

    for (r, k_book, σcr_book) in table_9_1
        k_exact, m_exact = exact_uniaxial_k(r)
        σcr_exact = exact_σcr(k_exact)
        ξ_exact = exact_ξcr(k_exact)
        nnodes, strategy_rows, k_diag = solve_table_9_1_case(r, k_exact)
        append!(kqq_nullspace_results, k_diag.kqq_rows)
        append!(k_nullspace_results, k_diag.k_rows)
        if isapprox(r, 1.0; atol=1.0e-12)
            println("K-only nullspace diagnostics a/b=1.0:")
            @printf("  Kqq near-zero count = %d\n", k_diag.kqq_near_zero_count)
            @printf("  K near-zero count   = %d\n", k_diag.k_near_zero_count)
        end
        for strategy_row in strategy_rows
            selected_mode = strategy_row.selected_mode
            block_diag = strategy_row.block_diag
            selected_energy = strategy_row.selected_energy
            append!(spectrum_results, strategy_row.spectrum_rows)
            push!(block_results, block_diag)
            push!(selected_mode_energy_results, selected_energy)
            k_err = abs(strategy_row.k_num-k_exact)/abs(k_exact)
            @printf("%7.2f %3d %18s %11s %14s %10d %10d %12.6e %12.6e %12.4e %14.6e %14.10f %14.10f %12.4e %8d\n",
                r, m_exact, basename(mesh_file(r)), strategy_row.kg_strategy, spectrum_selection_mode,
                selected_mode.mode_index, selected_mode.rank_by_xi, strategy_row.xi_cr, ξ_exact,
                selected_mode.xi_ratio, strategy_row.w_participation, strategy_row.k_num,
                k_exact, k_err, nnodes)
            if isapprox(r, 1.0; atol=1.0e-12)
                println("Matrix block diagnostics a/b=1.0, kg_strategy=$(strategy_row.kg_strategy):")
                @printf("  norm_kww=%12.6e norm_kqq=%12.6e norm_kqw=%12.6e norm_kgww=%12.6e active_norm_kgqq=%12.6e assembled_norm_kgqq=%12.6e\n",
                    block_diag.norm_kww, block_diag.norm_kqq, block_diag.norm_kqw,
                    block_diag.norm_kgww, block_diag.active_norm_kgqq,
                    block_diag.assembled_norm_kgqq)
                @printf("  energy_K_total=%12.6e active_energy_KG_total=%12.6e assembled_energy_KG_total=%12.6e rayleigh_xi_active=%12.6e\n",
                    block_diag.energy_K_total, block_diag.active_energy_KG_total,
                    block_diag.assembled_energy_KG_total, block_diag.rayleigh_xi_active)
                @printf("  selected %-9s psi_participation=%12.6e w_participation=%12.6e\n",
                    strategy_row.kg_strategy, selected_energy.psi_norm_ratio,
                    selected_energy.w_norm_ratio)
            end
            push!(results, (
                r = r,
                m_exact = m_exact,
                mesh = basename(mesh_file(r)),
                kg_strategy = strategy_row.kg_strategy,
                selection_mode = spectrum_selection_mode,
                mode_index = selected_mode.mode_index,
                rank_by_xi = selected_mode.rank_by_xi,
                xi_cr = strategy_row.xi_cr,
                xi_exact = ξ_exact,
                xi_ratio = selected_mode.xi_ratio,
                log10_xi_ratio = selected_mode.log10_xi_ratio,
                w_participation = strategy_row.w_participation,
                k_num = strategy_row.k_num,
                k_exact = k_exact,
                k_book = k_book,
                k_err = k_err,
                σcr_num = strategy_row.σcr_num,
                σcr_exact = σcr_exact,
                σcr_book = σcr_book,
                nnodes = nnodes,
            ))
        end
    end

    write_table_9_1_xlsx(
        output_xlsx, results, spectrum_results, block_results,
        kqq_nullspace_results, k_nullspace_results, selected_mode_energy_results,
    )
    println("Excel output: ", output_xlsx)
finally
    gmsh.finalize()
end

println(to)
