using LinearAlgebra
using Printf

const ROOT = dirname(dirname(@__FILE__))
const DATA_DIR = joinpath(ROOT, "data")

function element_matrices(E, nu, L, h, kappa, n_elem; rho=1.0)
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h^3 / 12.0
    ndof = 2 * (n_elem + 1)
    K = zeros(Float64, ndof, ndof)
    M = zeros(Float64, ndof, ndof)
    KG = zeros(Float64, ndof, ndof)
    le = L / n_elem

    EI = E * I
    kGA = kappa * G * A

    Kb = zeros(Float64, 4, 4)
    Kb[2, 2] = 1.0
    Kb[2, 4] = -1.0
    Kb[4, 2] = -1.0
    Kb[4, 4] = 1.0
    Kb .*= EI / le

    Bs = [-1.0 / le, -0.5, 1.0 / le, -0.5]
    Ks = (kGA * le) .* (Bs * Bs')
    Kloc = Kb + Ks

    Mloc = zeros(Float64, 4, 4)
    mt = rho * A * le
    mr = rho * I * le
    Mloc[1, 1] = mt / 3.0
    Mloc[3, 3] = mt / 3.0
    Mloc[1, 3] = mt / 6.0
    Mloc[3, 1] = mt / 6.0
    Mloc[2, 2] = mr / 3.0
    Mloc[4, 4] = mr / 3.0
    Mloc[2, 4] = mr / 6.0
    Mloc[4, 2] = mr / 6.0

    KGloc = zeros(Float64, 4, 4)
    KGloc[1, 1] = 1.0 / le
    KGloc[1, 3] = -1.0 / le
    KGloc[3, 1] = -1.0 / le
    KGloc[3, 3] = 1.0 / le

    for e in 1:n_elem
        dofs = [2e - 1, 2e, 2e + 1, 2e + 2]
        K[dofs, dofs] .+= Kloc
        M[dofs, dofs] .+= Mloc
        KG[dofs, dofs] .+= KGloc
    end

    return K, M, KG, A, I
end

function free_dofs(n_elem, bc)
    ndof = 2 * (n_elem + 1)
    fixed = if bc == "cantilever"
        [1, 2]
    elseif bc == "fixed-fixed"
        [1, 2, ndof - 1, ndof]
    elseif bc in ("simply-supported", "pinned-pinned")
        [1, ndof - 1]
    else
        error("Unknown boundary condition: $bc")
    end

    mask = trues(ndof)
    mask[fixed] .= false
    return findall(mask)
end

output_scale(E, I, A, L, rho) = L^2 * sqrt(rho * A / (E * I))

function solve_vibration_eigenpairs(E, nu, L, h, kappa, n_elem, bc;
                                    n_modes=15, rho=1.0, output="beta2")
    K, M, _, A, I = element_matrices(E, nu, L, h, kappa, n_elem; rho=rho)
    dofs = free_dofs(n_elem, bc)
    Kff = K[dofs, dofs]
    Mff = M[dofs, dofs]

    F = eigen(Symmetric(Kff), Symmetric(Mff))
    vals = real.(F.values)
    vecs = real.(F.vectors)
    keep = findall(v -> isfinite(v) && v > 1.0e-12, vals)
    order = keep[sortperm(vals[keep])]
    order = order[1:min(n_modes, length(order))]

    lambdas = vals[order]
    omegas = sqrt.(lambdas)
    scale = output_scale(E, I, A, L, rho)
    outputs = omegas .* scale
    if output == "beta"
        outputs = sqrt.(outputs)
    elseif output != "beta2"
        error("Unknown vibration output: $output")
    end

    full_modes = Vector{Vector{Float64}}()
    residuals = Float64[]
    for mode_id in order
        v = vecs[:, mode_id]
        lambda = vals[mode_id]
        push!(residuals, norm(Kff * v - lambda * Mff * v) /
                         max(norm(Kff * v), abs(lambda) * norm(Mff * v), eps(Float64)))

        full = zeros(Float64, 2 * (n_elem + 1))
        full[dofs] .= v
        push!(full_modes, full)
    end

    return outputs, lambdas, full_modes, residuals
end

function normalized_w_modes(full_modes, L)
    n_nodes = div(length(first(full_modes)), 2)
    x = collect(range(0.0, L; length=n_nodes))
    modes = []
    for full in full_modes
        w = copy(full[1:2:end])
        max_abs = maximum(abs.(w))
        if max_abs > 0.0
            w ./= max_abs
        end
        if w[cld(length(w), 2)] < 0.0
            w .*= -1.0
        end
        push!(modes, (x=x, w=w))
    end
    return modes
end

function write_vibration_csvs(prefix, bc, outputs, lambdas, modes, residuals)
    mkpath(DATA_DIR)
    eigen_path = joinpath(DATA_DIR, "$(prefix)_eigenvalues.csv")
    mode_path = joinpath(DATA_DIR, "$(prefix)_mode_shapes.csv")

    open(eigen_path, "w") do io
        println(io, "bc,mode,lambda_omega_sq,omega,beta_output,relative_residual")
        for i in eachindex(outputs)
            omega = sqrt(lambdas[i])
            @printf(io, "%s,%d,%.17g,%.17g,%.17g,%.17g\n",
                    bc, i, lambdas[i], omega, outputs[i], residuals[i])
        end
    end

    open(mode_path, "w") do io
        println(io, "bc,mode,node,x,w_normalized")
        for (mode_id, mode) in enumerate(modes)
            for node in eachindex(mode.x)
                @printf(io, "%s,%d,%d,%.17g,%.17g\n",
                        bc, mode_id, node, mode.x[node], mode.w[node])
            end
        end
    end

    return eigen_path, mode_path
end

function run_example()
    E = 1.0
    nu = 0.3
    L = 1.0
    kappa = 5.0 / 6.0
    h = 0.01 * L
    n_elem = 40
    rho = 1.0

    for bc in ("fixed-fixed", "simply-supported")
        outputs, lambdas, full_modes, residuals = solve_vibration_eigenpairs(
            E, nu, L, h, kappa, n_elem, bc;
            n_modes=4, rho=rho, output="beta",
        )
        modes = normalized_w_modes(full_modes, L)
        prefix = "julia_ch10_vibration_$(replace(bc, "-" => "_"))"
        eigen_path, mode_path = write_vibration_csvs(prefix, bc, outputs, lambdas, modes, residuals)

        println("Boundary condition: $bc")
        for i in eachindex(outputs)
            @printf("  mode %d: lambda=%.8e omega=%.8e beta=%.8e residual=%.3e\n",
                    i, lambdas[i], sqrt(lambdas[i]), outputs[i], residuals[i])
        end
        println("  wrote: $eigen_path")
        println("  wrote: $mode_path")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
