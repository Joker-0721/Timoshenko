"""
Timoshenko Beam - Full Integration Version
Using 2-point Gauss integration for shear stiffness matrix.

This version demonstrates the shear locking problem that occurs
when using full integration in Timoshenko beam elements.
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

# Analytical solutions
w_analytical_simply = 5 * q * L^4 / (384 * EI) + q * L^2 / (8 * κGA)
w_analytical_clamped = q * L^4 / (8 * EI) + q * L^2 / (2 * κGA)

println("=" ^ 60)
println("Timoshenko Beam - FULL INTEGRATION (2-point Gauss)")
println("=" ^ 60)
println("Parameters: E = $E Pa, G = $(round(G, digits=2)) Pa")
println("A = $A m², I = $I m⁴, κ = $(round(κ, digits=4))")
println("Length: L = $L m, Load: q = $q N/m")
println("Analytical (Simply-supported): w = $(w_analytical_simply) m")
println("Analytical (Clamped): w = $(w_analytical_clamped) m")
println("-" ^ 60)

# Gauss points and weights for 2-point integration
gauss_pts = [-1/sqrt(3), 1/sqrt(3)]
gauss_wts = [1.0, 1.0]

# Test different element numbers for convergence analysis
element_counts = [20, 30, 40, 50, 100]

println("\n*** CONVERGENCE ANALYSIS - Simply-Supported Beam ***")
println("-" ^ 60)
println(rpad("Elements", 10), rpad("FEM (m)", 18), rpad("Analytical (m)", 18), "Error (%)")
println("-" ^ 60)

for n_elements in element_counts
    n_nodes = n_elements + 1
    ndof = 2 * n_nodes
    nodes = range(0, L, length=n_nodes)
    
    K = zeros(ndof, ndof)
    F = zeros(ndof)
    
    for i in 1:n_elements
        le = nodes[i+1] - nodes[i]
        
        # Bending contribution
        K_b = zeros(4, 4)
        K_b[2,2] = 1; K_b[2,4] = -1; K_b[4,2] = -1; K_b[4,4] = 1
        K_b .*= EI / le
        
        # Shear contribution (Full 2-point Gauss integration)
        K_s = zeros(4, 4)
        for gp in 1:2
            ξ = gauss_pts[gp]
            weight = gauss_wts[gp]
            N1 = (1 - ξ) / 2
            N2 = (1 + ξ) / 2
            dN1_dx = -1 / le
            dN2_dx = 1 / le
            B_s = [dN1_dx, N1, dN2_dx, N2]
            K_s .+= weight * (B_s * B_s')
        end
        K_s .*= κGA * le
        
        K_local = K_b + K_s
        dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
        K[dofs, dofs] .+= K_local
        F[dofs] .+= [q*le/2, 0.0, q*le/2, 0.0]
    end
    
    # Simply-supported BC: w=0 at both ends
    free_dofs = setdiff(1:ndof, [1, ndof-1])
    U = zeros(ndof)
    U[free_dofs] = K[free_dofs, free_dofs] \ F[free_dofs]
    
    w_max = maximum(abs.(U[1:2:end]))
    error = abs(w_max - w_analytical_simply) / w_analytical_simply * 100
    
    println(rpad("$n_elements", 10), rpad("$(round(w_max, digits=8))", 18), rpad("$(round(w_analytical_simply, digits=8))", 18), "$(round(error, digits=4))")
end

println("\n*** CONVERGENCE ANALYSIS - Clamped Beam ***")
println("-" ^ 60)
println(rpad("Elements", 10), rpad("FEM (m)", 18), rpad("Analytical (m)", 18), "Error (%)")
println("-" ^ 60)

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
        
        K_s = zeros(4, 4)
        for gp in 1:2
            ξ = gauss_pts[gp]
            weight = gauss_wts[gp]
            N1 = (1 - ξ) / 2
            N2 = (1 + ξ) / 2
            dN1_dx = -1 / le
            dN2_dx = 1 / le
            B_s = [dN1_dx, N1, dN2_dx, N2]
            K_s .+= weight * (B_s * B_s')
        end
        K_s .*= κGA * le
        
        K_local = K_b + K_s
        dofs = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)]
        K[dofs, dofs] .+= K_local
        F[dofs] .+= [q*le/2, 0.0, q*le/2, 0.0]
    end
    
    # Clamped BC: w=0 and θ=0 at node 1
    free_dofs = setdiff(1:ndof, [1, 2])
    U = zeros(ndof)
    U[free_dofs] = K[free_dofs, free_dofs] \ F[free_dofs]
    
    w_max = maximum(abs.(U[1:2:end]))
    error = abs(w_max - w_analytical_clamped) / w_analytical_clamped * 100
    
    println(rpad("$n_elements", 10), rpad("$(round(w_max, digits=8))", 18), rpad("$(round(w_analytical_clamped, digits=8))", 18), "$(round(error, digits=4))")
end

println("\n" * "=" ^ 60)
println("NOTE: Full integration shows SHEAR LOCKING")
println("      Error does NOT decrease with mesh refinement!")
println("=" ^ 60)
