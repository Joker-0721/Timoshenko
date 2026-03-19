"""
Problem 16Buckling: Timoshenko Beam Buckling
From: MATLAB Codes for Finite Element Analysis (Ferreira 2009)

This problem solves buckling of Timoshenko beam.
- E = 210 GPa
- G = 84 GPa
- A = 0.02 m²
- I = 10e-5 m⁴
- κ = 5/6
- L = 10 m

Finds critical buckling load.
"""

using LinearAlgebra

# Parameters
E = 2.1e11  # Pa (210 GPa)
nu = 0.3
G = E / (2 * (1 + nu))
A = 0.02     # m²
I = 10e-5    # m⁴
kappa = 5/6  # Shear correction factor
L = 10.0     # m
n_elements = 40

println("=" ^ 60)
println("Problem 16Buckling: Timoshenko Beam Buckling")
println("=" ^ 60)
println("Parameters: E = $E Pa, G = $(round(G, digits=5)) Pa")
println("A = $A m², I = $I m⁴, κ = $(round(kappa, digits=4))")
println("Length: L = $L m")
println("Elements: $n_elements")
println("-" ^ 60)

# Number of nodes and DOFs
n_nodes = n_elements + 1
ndof = 2 * n_nodes  # w (transverse) and θ (rotation) at each node

# Initialize global matrices
K = zeros(ndof, ndof)      # Stiffness matrix
K_G = zeros(ndof, ndof)    # Geometric stiffness matrix

println("Assembling stiffness and geometric stiffness matrices...")
for i in 1:n_elements
    element_L = L / n_elements
    EI = E * I
    kGA = kappa * G * A
    
    # Exact 1-point integration stiffness matrix for Timoshenko beam
    K_b = zeros(4, 4)
    K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
    K_b .*= EI / element_L
    
    B_s = [-1/element_L, -1/2, 1/element_L, -1/2]
    K_s = (kGA * element_L) * (B_s * B_s')
    
    K_local = K_b + K_s
    
    # Geometric stiffness for Timoshenko beam (using linear w)
    KG_local = zeros(4, 4)
    KG_local[1,1] = 1/element_L
    KG_local[1,3] = -1/element_L
    KG_local[3,1] = -1/element_L
    KG_local[3,3] = 1/element_L
    
    dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
    K[dofs, dofs] .+= K_local
    K_G[dofs, dofs] .+= KG_local
end

# Analytical solutions for Euler buckling
P_cr_euler_pinned = π^2 * E * I / L^2
P_cr_euler_cantilever = π^2 * E * I / (4 * L^2)

# Timoshenko correction analytical
P_cr_timoshenko_pinned = P_cr_euler_pinned / (1 + P_cr_euler_pinned / (kappa * A * G))
P_cr_timoshenko_cantilever = P_cr_euler_cantilever / (1 + P_cr_euler_cantilever / (kappa * A * G))

# ------------------------------------------------------------------
# Boundary Condition 1: Pinned-Pinned (w=0 at both ends)
# ------------------------------------------------------------------
println("\n" * "=" ^ 60)
println("CASE 1: Pinned-Pinned Beam")
println("=" ^ 60)

free_dofs_pinned = setdiff(1:ndof, [1, ndof-1])

eigenvalues_pinned = eigvals(K[free_dofs_pinned, free_dofs_pinned], K_G[free_dofs_pinned, free_dofs_pinned])
λ_pinned = sort(filter(x -> x > 0 && isfinite(x), real(eigenvalues_pinned)))

println("Critical buckling loads (Pinned-Pinned, N):")
for i in 1:min(3, length(λ_pinned))
    println("  Mode $i: P_cr = $(round(λ_pinned[i], digits=2)) N")
end

println("  Analytical Timoshenko (Pinned-Pinned): P_cr = $(round(P_cr_timoshenko_pinned, digits=2)) N")

# ------------------------------------------------------------------
# Boundary Condition 2: Cantilever (w=0, θ=0 at x=0)
# ------------------------------------------------------------------
println("\n" * "=" ^ 60)
println("CASE 2: Cantilever Beam")
println("=" ^ 60)

free_dofs_cantilever = 3:ndof

eigenvalues_cant = eigvals(K[free_dofs_cantilever, free_dofs_cantilever], K_G[free_dofs_cantilever, free_dofs_cantilever])
λ_cantilever = sort(filter(x -> x > 0 && isfinite(x), real(eigenvalues_cant)))

println("Critical buckling loads (Cantilever, N):")
for i in 1:min(3, length(λ_cantilever))
    println("  Mode $i: P_cr = $(round(λ_cantilever[i], digits=2)) N")
end

println("  Analytical Timoshenko (Cantilever): P_cr = $(round(P_cr_timoshenko_cantilever, digits=2)) N")

println("\n" * "=" ^ 60)
println("Problem 16Buckling completed!")
println("=" ^ 60)
