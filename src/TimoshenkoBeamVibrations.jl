"""
Problem 16vibrations: Timoshenko Beam Free Vibrations
From: MATLAB Codes for Finite Element Analysis (Ferreira 2009)

This problem solves free vibration of Timoshenko beam.
- E = 210 GPa
- G = 84 GPa
- A = 0.02 m²
- I = 10e-5 m⁴
- κ = 5/6
- L = 10 m
- ρ = 1.0 (density, normalized for analytical comparison)

Finds natural frequencies and mode shapes.
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
rho = 1.0    # density (normalized to match original formulation)
n_elements = 40

println("=" ^ 60)
println("Problem 16vibrations: Timoshenko Beam Free Vibrations")
println("=" ^ 60)
println("Parameters: E = $E Pa, G = $(round(G, digits=5)) Pa")
println("A = $A m², I = $I m⁴, κ = $(round(kappa, digits=4)), ρ = $rho")
println("Length: L = $L m")
println("Elements: $n_elements")
println("-" ^ 60)

# Number of nodes and DOFs
n_nodes = n_elements + 1
ndof = 2 * n_nodes  # w (transverse) and θ (rotation) at each node

# Initialize global matrices
K = zeros(ndof, ndof)
M = zeros(ndof, ndof)

println("Assembling stiffness and mass matrices...")
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
    
    # Consistent mass matrix for Timoshenko beam (linear interpolations)
    M_local = zeros(4, 4)
    
    # Transverse mass
    M_local[1,1] = M_local[3,3] = rho * A * element_L / 3
    M_local[1,3] = M_local[3,1] = rho * A * element_L / 6
    
    # Rotational mass
    M_local[2,2] = M_local[4,4] = rho * I * element_L / 3
    M_local[2,4] = M_local[4,2] = rho * I * element_L / 6

    dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
    K[dofs, dofs] .+= K_local
    M[dofs, dofs] .+= M_local
end

# Case: Cantilever beam (clamped at x=0, free at x=L)
println("\n" * "=" ^ 60)
println("CASE: Cantilever Beam (clamped at x=0)")
println("=" ^ 60)

# Apply boundary conditions: w=0 and θ=0 at x=0 (node 1) -> DOFs 1, 2
free_dofs = 3:ndof

println("Solving eigenvalue problem...")
eigenvalues = eigvals(K[free_dofs, free_dofs], M[free_dofs, free_dofs])

# Filter positive eigenvalues and sort
ω² = real(eigenvalues)
ω_positive = sqrt.(sort(filter(x -> x > 0 && isfinite(x), ω²)))

println("\nNatural frequencies (rad/s):")
for i in 1:min(5, length(ω_positive))
    println("  Mode $i: ω = $(round(ω_positive[i], digits=4)) rad/s")
end

# Analytical solution for cantilever beam (Euler-Bernoulli for comparison)
println("\n" * "-" ^ 60)
println("Analytical solution (Euler-Bernoulli):")
println("-" ^ 60)

β = [1.8751, 4.6941, 7.8548, 10.9955]  # exact β values for cantilever
omega_analytical = β.^2 .* sqrt.(E*I / (rho*A)) ./ L^2

for i in 1:min(4, length(β))
    println("  Mode $i: ω = $(round(omega_analytical[i], digits=4)) rad/s")
end

println("\n" * "=" ^ 60)
println("Problem 16vibrations completed!")
println("=" ^ 60)
