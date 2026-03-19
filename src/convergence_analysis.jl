"""
Convergence Analysis for Timoshenko Beam
Elements: 20, 30, 40, 50, 60, 70, 80, 90, 100
With Analytical Solution Comparison
"""

using LinearAlgebra

# Parameters (shared across all analyses)
E = 2.1e11  # Pa
nu = 0.3
G = E / (2 * (1 + nu))
A = 0.02    # m²
I = 10e-5   # m⁴
kappa = 5/6
L = 10.0    # m
rho = 1.0
q = 1000.0  # N/m

# Analytical solutions
EI = E * I
kappaGA = kappa * G * A

# Static bending analytical solutions
w_analytical_simply = 5 * q * L^4 / (384 * EI) + q * L^2 / (8 * kappaGA)
w_analytical_clamped = q * L^4 / (8 * EI) + q * L^2 / (2 * kappaGA)

# Buckling analytical solutions
P_cr_euler_pinned = π^2 * EI / L^2
P_cr_euler_cantilever = π^2 * EI / (4 * L^2)
P_analytical_pinned = P_cr_euler_pinned / (1 + P_cr_euler_pinned / (kappa * A * G))
P_analytical_cantilever = P_cr_euler_cantilever / (1 + P_cr_euler_cantilever / (kappa * A * G))

# Vibration analytical solutions (Euler-Bernoulli)
beta = [1.8751, 4.6941, 7.8548, 10.9955]
omega_analytical = beta.^2 .* sqrt.(EI / (rho * A)) ./ L^2

element_counts = [20, 30, 40, 50, 60, 70, 80, 90, 100]

println("=" ^ 90)
println("CONVERGENCE ANALYSIS - TIMOSHENKO BEAM WITH ANALYTICAL SOLUTION COMPARISON")
println("=" ^ 90)

# =====================
# 1. STATIC BENDING
# =====================
println("\n" * "=" ^ 90)
println("1. STATIC BENDING ANALYSIS")
println("=" ^ 90)
println("Analytical Solution:")
println("  Simply-supported: w = ", w_analytical_simply, " m")
println("  Clamped:          w = ", w_analytical_clamped, " m")
println("\n%6s | %15s | %15s | %8s | %15s | %15s | %8s", 
        "Elements", "FEM (m)", "Analytical", "Err(%)", "FEM (m)", "Analytical", "Err(%)")
println("-" ^ 90)

for n_elements in element_counts
    n_nodes = n_elements + 1
    ndof = 2 * n_nodes
    nodes = range(0, L, length=n_nodes)
    
    K = zeros(ndof, ndof)
    F = zeros(ndof)
    
    for i in 1:n_elements
        le = nodes[i+1] - nodes[i]
        
        K_b = zeros(4, 4)
        K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
        K_b .*= EI / le
        
        B_s = [-1/le, -1/2, 1/le, -1/2]
        K_s = (kappaGA * le) * (B_s * B_s')
        
        K_local = K_b + K_s
        
        dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
        K[dofs, dofs] .+= K_local
        F[dofs] .+= [q*le/2, 0.0, q*le/2, 0.0]
    end
    
    # Simply-supported
    free_dofs_s = setdiff(1:ndof, [1, ndof-1])
    U1 = zeros(ndof)
    U1[free_dofs_s] = K[free_dofs_s, free_dofs_s] \ F[free_dofs_s]
    w_simply = maximum(abs.(U1[1:2:end]))
    error_simply = abs(w_simply - w_analytical_simply) / w_analytical_simply * 100
    
    # Clamped
    free_dofs_c = setdiff(1:ndof, [1, 2])
    U2 = zeros(ndof)
    U2[free_dofs_c] = K[free_dofs_c, free_dofs_c] \ F[free_dofs_c]
    w_clamped = maximum(abs.(U2[1:2:end]))
    error_clamped = abs(w_clamped - w_analytical_clamped) / w_analytical_clamped * 100
    
    println(" %6d | %15.8f | %15.8f | %7.4f%% | %15.8f | %15.8f | %7.4f%%", 
            n_elements, w_simply, w_analytical_simply, error_simply, 
            w_clamped, w_analytical_clamped, error_clamped)
end

# =====================
# 2. BUCKLING
# =====================
println("\n" * "=" ^ 90)
println("2. BUCKLING ANALYSIS")
println("=" ^ 90)
println("Analytical Solution:")
println("  Pinned-Pinned: P_cr = ", P_analytical_pinned, " N")
println("  Cantilever:    P_cr = ", P_analytical_cantilever, " N")
println("\n%6s | %15s | %15s | %8s | %15s | %15s | %8s", 
        "Elements", "FEM (N)", "Analytical", "Err(%)", "FEM (N)", "Analytical", "Err(%)")
println("-" ^ 90)

for n_elements in element_counts
    n_nodes = n_elements + 1
    ndof = 2 * n_nodes
    
    K = zeros(ndof, ndof)
    K_G = zeros(ndof, ndof)
    
    for i in 1:n_elements
        element_L = L / n_elements
        
        K_b = zeros(4, 4)
        K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
        K_b .*= EI / element_L
        
        B_s = [-1/element_L, -1/2, 1/element_L, -1/2]
        K_s = (kappa * G * A * element_L) * (B_s * B_s')
        
        K_local = K_b + K_s
        
        KG_local = zeros(4, 4)
        KG_local[1,1] = 1/element_L
        KG_local[1,3] = -1/element_L
        KG_local[3,1] = -1/element_L
        KG_local[3,3] = 1/element_L
        
        dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
        K[dofs, dofs] .+= K_local
        K_G[dofs, dofs] .+= KG_local
    end
    
    # Pinned-Pinned
    free_dofs_pinned = setdiff(1:ndof, [1, ndof-1])
    eigenvalues_pinned = eigvals(K[free_dofs_pinned, free_dofs_pinned], K_G[free_dofs_pinned, free_dofs_pinned])
    P_cr_pinned = minimum(filter(x -> x > 0 && isfinite(x), real(eigenvalues_pinned)))
    error_pinned = abs(P_cr_pinned - P_analytical_pinned) / P_analytical_pinned * 100
    
    # Cantilever
    free_dofs_cantilever = 3:ndof
    eigenvalues_cant = eigvals(K[free_dofs_cantilever, free_dofs_cantilever], K_G[free_dofs_cantilever, free_dofs_cantilever])
    P_cr_cant = minimum(filter(x -> x > 0 && isfinite(x), real(eigenvalues_cant)))
    error_cantilever = abs(P_cr_cant - P_analytical_cantilever) / P_analytical_cantilever * 100
    
    println(" %6d | %15.2f | %15.2f | %7.4f%% | %15.2f | %15.2f | %7.4f%%", 
            n_elements, P_cr_pinned, P_analytical_pinned, error_pinned,
            P_cr_cant, P_analytical_cantilever, error_cantilever)
end

# =====================
# 3. VIBRATIONS (Mode 1-4)
# =====================
println("\n" * "=" ^ 90)
println("3. VIBRATION ANALYSIS (Mode 1-4)")
println("=" ^ 90)
println("Analytical Solution (Euler-Bernoulli):")
for m in 1:4
    println("  Mode ", m, ": ω = ", omega_analytical[m], " rad/s")
end
println("\n%6s | %12s | %12s | %8s | %12s | %12s | %8s | %12s | %12s | %8s", 
        "Elements", "FEM(M1)", "Anal(M1)", "Err(%)", "FEM(M2)", "Anal(M2)", "Err(%)", "FEM(M3)", "Anal(M3)", "Err(%)")
println("-" ^ 90)

for n_elements in element_counts
    n_nodes = n_elements + 1
    ndof = 2 * n_nodes
    
    K = zeros(ndof, ndof)
    M = zeros(ndof, ndof)
    
    for i in 1:n_elements
        element_L = L / n_elements
        
        K_b = zeros(4, 4)
        K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
        K_b .*= EI / element_L
        
        B_s = [-1/element_L, -1/2, 1/element_L, -1/2]
        K_s = (kappa * G * A * element_L) * (B_s * B_s')
        
        K_local = K_b + K_s
        
        M_local = zeros(4, 4)
        M_local[1,1] = M_local[3,3] = rho * A * element_L / 3
        M_local[1,3] = M_local[3,1] = rho * A * element_L / 6
        M_local[2,2] = M_local[4,4] = rho * I * element_L / 3
        M_local[2,4] = M_local[4,2] = rho * I * element_L / 6
        
        dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
        K[dofs, dofs] .+= K_local
        M[dofs, dofs] .+= M_local
    end
    
    # Cantilever
    free_dofs = 3:ndof
    eigenvalues = eigvals(K[free_dofs, free_dofs], M[free_dofs, free_dofs])
    ω² = real(eigenvalues)
    ω_positive = sqrt.(sort(filter(x -> x > 0 && isfinite(x), ω²)))
    
    error1 = abs(ω_positive[1] - omega_analytical[1]) / omega_analytical[1] * 100
    error2 = abs(ω_positive[2] - omega_analytical[2]) / omega_analytical[2] * 100
    error3 = abs(ω_positive[3] - omega_analytical[3]) / omega_analytical[3] * 100
    
    println(" %6d | %12.4f | %12.4f | %7.4f%% | %12.4f | %12.4f | %7.4f%% | %12.4f | %12.4f | %7.4f%%", 
            n_elements, ω_positive[1], omega_analytical[1], error1,
            ω_positive[2], omega_analytical[2], error2,
            ω_positive[3], omega_analytical[3], error3)
end

# Mode 4 table
println("\n%6s | %12s | %12s | %8s", "Elements", "FEM(M4)", "Anal(M4)", "Err(%)")
println("-" ^ 60)
for n_elements in element_counts
    n_nodes = n_elements + 1
    ndof = 2 * n_nodes
    
    K = zeros(ndof, ndof)
    M = zeros(ndof, ndof)
    
    for i in 1:n_elements
        element_L = L / n_elements
        
        K_b = zeros(4, 4)
        K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
        K_b .*= EI / element_L
        
        B_s = [-1/element_L, -1/2, 1/element_L, -1/2]
        K_s = (kappa * G * A * element_L) * (B_s * B_s')
        
        K_local = K_b + K_s
        
        M_local = zeros(4, 4)
        M_local[1,1] = M_local[3,3] = rho * A * element_L / 3
        M_local[1,3] = M_local[3,1] = rho * A * element_L / 6
        M_local[2,2] = M_local[4,4] = rho * I * element_L / 3
        M_local[2,4] = M_local[4,2] = rho * I * element_L / 6
        
        dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
        K[dofs, dofs] .+= K_local
        M[dofs, dofs] .+= M_local
    end
    
    free_dofs = 3:ndof
    eigenvalues = eigvals(K[free_dofs, free_dofs], M[free_dofs, free_dofs])
    ω² = real(eigenvalues)
    ω_positive = sqrt.(sort(filter(x -> x > 0 && isfinite(x), ω²)))
    
    error4 = abs(ω_positive[4] - omega_analytical[4]) / omega_analytical[4] * 100
    
    println(" %6d | %12.4f | %12.4f | %7.4f%%", 
            n_elements, ω_positive[4], omega_analytical[4], error4)
end

println("\n" * "=" ^ 90)
println("CONVERGENCE ANALYSIS COMPLETED")
println("=" ^ 90)
