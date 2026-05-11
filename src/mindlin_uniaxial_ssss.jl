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
const output_xlsx = normpath(joinpath(@__DIR__, "..", "results", "mindlin_uniaxial_ssss_table_9_1.xlsx"))
const eigen_imag_tol = 1.0e-7
const spectrum_first_positive_count = 40
const spectrum_closest_exact_count = 10
const spectrum_selection_mode = "closest_exact"

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
        row.mode_index,
        row.rank_by_xi,
        row.xi_selected,
        row.xi_exact,
        row.xi_ratio,
        row.norm_kww,
        row.norm_kqq,
        row.norm_kqw,
        row.norm_kgww,
        row.norm_kgqq,
        row.ratio_kqq_kww,
        row.ratio_kqw_geom_mean,
        row.ratio_kgqq_kgww,
        row.energy_kww,
        row.energy_kqq,
        row.energy_kqw,
        row.energy_K_total,
        row.energy_kgww,
        row.energy_kgqq,
        row.energy_KG_total,
        row.rayleigh_xi,
        row.w_participation,
    ]
end

safe_ratio(numerator, denominator) = denominator == 0.0 ? NaN : numerator/denominator

function matrix_block_diagnostics(
    r, k_exact, selected_mode, F, nᵠ, nʷ,
    kʷʷ, kᵠᵠ, kᵠʷ, kgʷʷ, kgᵠᵠ,
)
    v = @view F.vectors[:, selected_mode.mode_index]
    vᵠ = @view v[1:2*nᵠ]
    vʷ = @view v[2*nᵠ+1:2*nᵠ+nʷ]

    norm_kww = norm(kʷʷ)
    norm_kqq = norm(kᵠᵠ)
    norm_kqw = norm(kᵠʷ)
    norm_kgww = norm(kgʷʷ)
    norm_kgqq = norm(kgᵠᵠ)

    energy_kww = real(dot(vʷ, kʷʷ*vʷ))
    energy_kqq = real(dot(vᵠ, kᵠᵠ*vᵠ))
    energy_kqw = 2.0*real(dot(vᵠ, kᵠʷ*vʷ))
    energy_K_total = energy_kww + energy_kqq + energy_kqw
    energy_kgww = real(dot(vʷ, kgʷʷ*vʷ))
    energy_kgqq = real(dot(vᵠ, kgᵠᵠ*vᵠ))
    energy_KG_total = energy_kgww + energy_kgqq

    return (
        r = r,
        mesh = basename(mesh_file(r)),
        mode_index = selected_mode.mode_index,
        rank_by_xi = selected_mode.rank_by_xi,
        xi_selected = selected_mode.xi,
        xi_exact = exact_ξcr(k_exact),
        xi_ratio = selected_mode.xi_ratio,
        norm_kww = norm_kww,
        norm_kqq = norm_kqq,
        norm_kqw = norm_kqw,
        norm_kgww = norm_kgww,
        norm_kgqq = norm_kgqq,
        ratio_kqq_kww = safe_ratio(norm_kqq, norm_kww),
        ratio_kqw_geom_mean = norm_kqq*norm_kww <= 0.0 ? NaN : norm_kqw/sqrt(norm_kqq*norm_kww),
        ratio_kgqq_kgww = safe_ratio(norm_kgqq, norm_kgww),
        energy_kww = energy_kww,
        energy_kqq = energy_kqq,
        energy_kqw = energy_kqw,
        energy_K_total = energy_K_total,
        energy_kgww = energy_kgww,
        energy_kgqq = energy_kgqq,
        energy_KG_total = energy_KG_total,
        rayleigh_xi = safe_ratio(energy_K_total, energy_KG_total),
        w_participation = selected_mode.w_participation,
    )
end

function write_table_9_1_xlsx(filepath, rows, spectrum_rows, block_rows)
    mkpath(dirname(filepath))
    table_headers = [
        "a/b", "m", "mesh", "selection_mode", "mode_index", "rank_by_xi",
        "xi_cr", "xi_exact", "xi_ratio", "log10_xi_ratio", "w_participation",
        "k_num", "k_exact", "k_book", "k_error", "sigma_cr_num",
        "sigma_cr_exact", "sigma_cr_book", "nodes",
    ]
    spectrum_headers = [
        "a/b", "mesh", "spectrum_group", "group_rank", "mode_index", "rank_by_xi",
        "xi", "xi_exact", "xi_ratio", "log10_xi_ratio", "distance_to_exact",
        "w_participation", "k", "k_exact",
    ]
    block_headers = [
        "a/b", "mesh", "mode_index", "rank_by_xi", "xi_selected", "xi_exact",
        "xi_ratio", "norm_kww", "norm_kqq", "norm_kqw", "norm_kgww",
        "norm_kgqq", "ratio_kqq_kww", "ratio_kqw_geom_mean",
        "ratio_kgqq_kgww", "energy_kww", "energy_kqq", "energy_kqw",
        "energy_K_total", "energy_kgww", "energy_kgqq", "energy_KG_total",
        "rayleigh_xi", "w_participation",
    ]

    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Table 9-1")
        write_sheet_rows!(sheet, table_headers, table_9_1_values.(rows))

        spectrum_sheet = XLSX.addsheet!(xf, "Spectrum")
        write_sheet_rows!(spectrum_sheet, spectrum_headers, spectrum_values.(spectrum_rows))

        block_sheet = XLSX.addsheet!(xf, "MatrixBlocks")
        write_sheet_rows!(block_sheet, block_headers, block_values.(block_rows))
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

    @timeit to "solve buckling eigenvalue" begin
        K = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
        KG = [
            kgᵠᵠ zeros(2*nᵠ,nʷ)
            zeros(nʷ,2*nᵠ) kgʷʷ
        ]
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
        isempty(candidates_raw) && error("no positive finite buckling eigenvalue found for a/b = $r")
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

        ξcr = selected_mode.xi
        w_participation = selected_mode.w_participation
        k_num = selected_mode.k
        σcr_num = ξcr/h
        block_diag = matrix_block_diagnostics(
            r, k_exact, selected_mode, F, nᵠ, nʷ,
            kʷʷ, kᵠᵠ, kᵠʷ, kgʷʷ, kgᵠᵠ,
        )
    end

    return ξcr, k_num, σcr_num, nʷ, w_participation, selected_mode, spectrum_rows, block_diag
end

gmsh.initialize()
try
    println("T1 表 9-1: SSSS rectangular plate under uniaxial in-plane compression")
    println("E = 30e6 psi, ν = 0.3, h/b = 0.01")
    println("Required mesh names: msh/mindlin_uniaxial_ssss_ab_<a_over_b>.msh, e.g. ab_1p00.")
    results = []
    spectrum_results = []
    block_results = []
    @printf("%7s %3s %18s %14s %10s %10s %12s %12s %12s %14s %14s %14s %12s %8s\n",
        "a/b", "m", "msh", "selection", "mode", "rank", "ξcr", "ξ_exact",
        "ξ/ξex", "w_part", "k_num", "k_exact", "k_err", "nodes")
    println("-"^172)

    for (r, k_book, σcr_book) in table_9_1
        k_exact, m_exact = exact_uniaxial_k(r)
        σcr_exact = exact_σcr(k_exact)
        ξ_exact = exact_ξcr(k_exact)
        ξcr, k_num, σcr_num, nnodes, w_participation, selected_mode, spectrum_rows, block_diag = solve_table_9_1_case(r, k_exact)
        append!(spectrum_results, spectrum_rows)
        push!(block_results, block_diag)
        k_err = abs(k_num-k_exact)/abs(k_exact)
        @printf("%7.2f %3d %18s %14s %10d %10d %12.6e %12.6e %12.4e %14.6e %14.10f %14.10f %12.4e %8d\n",
            r, m_exact, basename(mesh_file(r)), spectrum_selection_mode, selected_mode.mode_index,
            selected_mode.rank_by_xi, ξcr, ξ_exact, selected_mode.xi_ratio, w_participation,
            k_num, k_exact, k_err, nnodes)
        if isapprox(r, 1.0; atol=1.0e-12)
            println("Matrix block diagnostics a/b=1.0:")
            @printf("  norm_kww=%12.6e norm_kqq=%12.6e norm_kqw=%12.6e norm_kgww=%12.6e norm_kgqq=%12.6e\n",
                block_diag.norm_kww, block_diag.norm_kqq, block_diag.norm_kqw,
                block_diag.norm_kgww, block_diag.norm_kgqq)
            @printf("  energy_K_total=%12.6e energy_KG_total=%12.6e rayleigh_xi=%12.6e\n",
                block_diag.energy_K_total, block_diag.energy_KG_total, block_diag.rayleigh_xi)
        end
        push!(results, (
            r = r,
            m_exact = m_exact,
            mesh = basename(mesh_file(r)),
            selection_mode = spectrum_selection_mode,
            mode_index = selected_mode.mode_index,
            rank_by_xi = selected_mode.rank_by_xi,
            xi_cr = ξcr,
            xi_exact = ξ_exact,
            xi_ratio = selected_mode.xi_ratio,
            log10_xi_ratio = selected_mode.log10_xi_ratio,
            w_participation = w_participation,
            k_num = k_num,
            k_exact = k_exact,
            k_book = k_book,
            k_err = k_err,
            σcr_num = σcr_num,
            σcr_exact = σcr_exact,
            σcr_book = σcr_book,
            nnodes = nnodes,
        ))
    end

    write_table_9_1_xlsx(output_xlsx, results, spectrum_results, block_results)
    println("Excel output: ", output_xlsx)
finally
    gmsh.finalize()
end

println(to)
