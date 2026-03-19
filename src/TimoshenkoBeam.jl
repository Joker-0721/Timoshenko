"""
Problem 16: Timoshenko Beam in Bending
From: MATLAB Codes for Finite Element Analysis (Ferreira 2009)

This problem solves a Timoshenko beam under uniform load.
Using exact 1-point integration to prevent shear locking.
"""

using LinearAlgebra

# Parameters
E = 210e9   # Young's modulus (Pa)
ν = 0.3     # Poisson's ratio
L = 10.0    # Length (m)
A = 0.02    # Cross-section area (m^2)
I = 0.0001  # Moment of inertia (m^4)
G = E / (2 * (1 + ν))  # Shear modulus
EI = E * I
κ = 5/6     # Shear correction factor
κGA = κ * A * G
q = 1000.0  # Load (N/m)

n_elements = 40
n_nodes = n_elements + 1
ndof = 2 * n_nodes

# Generate node coordinates
nodes = range(0, L, length=n_nodes)

println("=" ^ 60)
println("Problem 16: Timoshenko Beam")
println("=" ^ 60)
println("Parameters: E = $E Pa, G = $(round(G, digits=2)) Pa")
println("A = $A m², I = $I m⁴, κ = $(round(κ, digits=4))")
println("Length: L = $L m, Load: q = $q N/m")
println("Elements: $n_elements")
println("-" ^ 60)

# Initialize global matrices
K = zeros(ndof, ndof)
F = zeros(ndof)

# Assemble global stiffness and force matrices
println("Assembling stiffness matrix...")
for i in 1:n_elements
    le = nodes[i+1] - nodes[i]
    
    # Bending contribution
    K_b = zeros(4, 4)
    K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
    K_b .*= EI / le
    
    # Shear contribution (1-point Gauss integrated)
    B_s = [-1/le, -1/2, 1/le, -1/2]
    K_s = (κGA * le) * (B_s * B_s')
    
    K_local = K_b + K_s
    
    # Dofs: [w_1, θ_1, w_2, θ_2]
    dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
    K[dofs, dofs] .+= K_local
    
    # Force vector for uniform load
    F[dofs] .+= [q*le/2, 0.0, q*le/2, 0.0]
end

# Case 1: Simply-supported beam
println("\n" * "=" ^ 60)
println("CASE 1: Simply-Supported Beam")
println("=" ^ 60)

# Apply BC: w=0 at both ends (simply-supported)
free_dofs_s = setdiff(1:ndof, [1, ndof-1])
U1 = zeros(ndof)
U1[free_dofs_s] = K[free_dofs_s, free_dofs_s] \ F[free_dofs_s]

w_disp = U1[1:2:end]
w_max_simply = maximum(abs.(w_disp))
println("\nMaximum displacement (simply-supported): w = $(round(w_max_simply, digits=8)) m")

# Case 2: Clamped beam (cantilever)
println("\n" * "=" ^ 60)
println("CASE 2: Clamped Beam")
println("=" ^ 60)

# Apply BC: w=0 and θ=0 at node 1 (clamped)
free_dofs_c = setdiff(1:ndof, [1, 2])
U2 = zeros(ndof)
U2[free_dofs_c] = K[free_dofs_c, free_dofs_c] \ F[free_dofs_c]

w_disp2 = U2[1:2:end]
w_max_clamped = maximum(abs.(w_disp2))
println("\nMaximum displacement (clamped): w = $(round(w_max_clamped, digits=8)) m")

# Compare with analytical solutions
println("\n" * "-" ^ 60)
println("Comparison with analytical solutions:")
println("-" ^ 60)

# Bernoulli solutions
w_analytical_simply = 5 * q * L^4 / (384 * EI)
w_analytical_clamped = q * L^4 / (8 * EI)

# Timoshenko correction for shear
w_shear_simply = q * L^2 / (8 * κGA)
w_shear_clamped = q * L^2 / (2 * κGA)

w_timoshenko_simply = w_analytical_simply + w_shear_simply
w_timoshenko_clamped = w_analytical_clamped + w_shear_clamped

println("  Simply-supported (Bernoulli): w = $(round(w_analytical_simply, digits=8)) m")
println("  Simply-supported (Timoshenko exact): w = $(round(w_timoshenko_simply, digits=8)) m")
println("  Simply-supported (FEM): w = $(round(w_max_simply, digits=8)) m")
println("")
println("  Clamped (Bernoulli): w = $(round(w_analytical_clamped, digits=8)) m")
println("  Clamped (Timoshenko exact): w = $(round(w_timoshenko_clamped, digits=8)) m")
println("  Clamped (FEM): w = $(round(w_max_clamped, digits=8)) m")

println("\n" * "=" ^ 60)
println("Problem 16 completed!")
println("=" ^ 60)
