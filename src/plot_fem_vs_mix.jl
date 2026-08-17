using GLMakie
using DelimitedFiles

BCs = ["CCCC", "CSCS", "SCSC", "SSSS", "FSCS", "FSSS"]
ns  = [7, 9, 11, 13, 15, 17]

# Bui data (from user)
bui = Dict(
    "CCCC" => Dict(7=>17.365, 9=>11.609, 11=>10.406, 13=>10.108, 15=>10.182, 17=>10.101),
    "CSCS" => Dict(7=>11.930, 9=>8.847,  11=>8.209,  13=>7.702,  15=>7.719,  17=>7.698),
    "SCSC" => Dict(7=>9.872,  9=>7.268,  11=>6.933,  13=>6.780,  15=>6.756,  17=>6.748),
    "SSSS" => Dict(7=>4.628,  9=>4.099,  11=>4.024,  13=>4.005,  15=>3.997,  17=>3.999),
    "FSCS" => Dict(7=>2.121,  9=>1.793,  11=>1.729,  13=>1.693,  15=>1.717,  17=>1.704),
    "FSSS" => Dict(7=>1.668,  9=>1.495,  11=>1.452,  13=>1.439,  15=>1.445,  17=>1.435),
)
exact = Dict(
    "CCCC" => 10.070, "CSCS" => 7.690, "SCSC" => 6.750,
    "SSSS" => 4.000,  "FSCS" => 1.700, "FSSS" => 1.440,
)

function readcsv(path)
    m = readdlm(path, ',')
    d = Dict{Int,Float64}()
    for i in 2:size(m,1)
        d[Int(m[i,1])] = Float64(m[i,2])
    end
    return d
end

for bc in BCs
    mix  = readcsv("D:/Joker/Timoshenko/date/mix_w_φ_$(bc).csv")
    roit = readcsv("D:/Joker/Timoshenko/date/mix_w_φ_$(bc)_roit.csv")
    fem  = readcsv("D:/Joker/Timoshenko/date/fem_$(bc).csv")
    ex   = exact[bc]

    # Normalized buckling load = k_num / k_exact
    r_mix  = [mix[n]  / ex for n in ns]
    r_roit = [roit[n] / ex for n in ns]
    r_fem  = [fem[n]  / ex for n in ns]
    r_bui  = [bui[bc][n] / ex for n in ns]

    # ---- Convergence: X = n (7x7 ... 17x17), Y = normalized load ----
    fig1 = Figure(size=(800, 600))
    ax1 = Axis(fig1[1,1],
               xlabel = "n (n×n nodes)",
               ylabel = L"k_{\mathrm{num}} / k_{\mathrm{exact}}",
               title = "Normalized buckling load $(bc)")
    scatterlines!(ax1, ns, r_mix,  label="mix",  color=:blue)
    scatterlines!(ax1, ns, r_roit, label="mix-roit", color=:orange)
    scatterlines!(ax1, ns, r_fem,  label="fem",  color=:red)
    scatterlines!(ax1, ns, r_bui,  label="Bui",  color=:green)
    hlines!(ax1, [1.0], color=:black, linestyle=:dash, label="exact (=1)")
    axislegend(ax1, position=:rt)
    save("D:/Joker/Timoshenko/Fig/fem_vs_mix/convergence_$(bc).png", fig1)
end

println("plots done")