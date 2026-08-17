using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, ∫∇wσ∇wdΩ, ∫∇φσ∇φdΩ, ∫ρwwdΩ, ∫ρφφdΩ

using TimerOutputs, LinearAlgebra, WriteVTK, DelimitedFiles
import Gmsh: gmsh
include("cal_area_support.jl")

E = 1.0
ν = 0.3
ρ = 1.0
h = 1e-2
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
σ₁₁ = 1e0
σ₂₂ = 0.0
σ₁₂ = 0.0
a = 1.0
ns = [7,9,11,13,15,17]
integrationOrder = 2

const to = TimerOutput()

open("mix_w_φ_CCCC_roit.csv", "w") do io
    write(io, "node_density,lambda_scaled\n")
    for i in ns
        gmsh.initialize()
        type_q = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
        type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
        type_w = :tri3
        
        ndiv_φ = i-1
        ndiv_w = i
        ndiv_q = i-1
        
        @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_q.msh")
        @timeit to "get nodes" nodes = get𝑿ᵢ()
        xʷ = nodes.x
        yʷ = nodes.y
        zʷ = nodes.z
        sp_q = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
        @timeit to "get entities" entities = getPhysicalGroups()
        elements_support = getElements(nodes, entities["Ω"], 1)
        nʷ = length(nodes)
        s_q, var_A = cal_area_support(elements_support)
        s₁ = 1.5*s_q*ones(nʷ)
        s₂ = 1.5*s_q*ones(nʷ)
        s₃ = 1.5*s_q*ones(nʷ)
        push!(nodes,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
        # ─── Rotation Φ ───────────────────────────────────────────
        @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_φ.msh")
        @timeit to "get nodes" nodes_φ = get𝑿ᵢ()
        xᵠ = nodes_φ.x
        yᵠ = nodes_φ.y
        zᵠ = nodes_φ.z
        sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
        nᵠ = length(nodes_φ)
        @timeit to "get entities" entities = getPhysicalGroups()
        elements_support = getElements(nodes_φ, entities["Ω"], 1)
        s_φ, var_A = cal_area_support(elements_support)
        s₁ = 1.5*s_φ*ones(nᵠ)
        s₂ = 1.5*s_φ*ones(nᵠ)
        s₃ = 1.5*s_φ*ones(nᵠ)
        push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
        # ─── Shear ────────────────────────────────────────────────
        @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_w.msh")
        @timeit to "get nodes" nodes_w = get𝑿ᵢ()
        @timeit to "get entities" entities = getPhysicalGroups()
        elements_support = getElements(nodes_w, entities["Ω"], 1)
        s_w, var_A = cal_area_support(elements_support)

        nʷ = length(nodes_w)
        nᵠ = length(nodes_φ)
        nˢ = length(nodes)
        kʷʷ = zeros(nʷ,nʷ)
        kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
        kˢˢ = zeros(2*nˢ,2*nˢ)
        kᴳʷʷ = zeros(nʷ,nʷ)
        kᴳᵠᵠ = zeros(2*nᵠ,2*nᵠ)
        mʷʷ = zeros(nʷ,nʷ)
        mᵠᵠ = zeros(2*nᵠ,2*nᵠ)
        kᵠʷ = zeros(2*nᵠ,nʷ)
        kˢʷ = zeros(2*nˢ,nʷ)
        kˢᵠ = zeros(2*nˢ,2*nᵠ)

        @timeit to "calculate ∫κκdΩ, ∫QφdΩ, ∫∇φσ∇φdΩ, ∫ρφφdΩ" begin
            @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], eval(type_q), integrationOrder, sp_q)
            prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
            @timeit to "calculate shape functions" set𝝭!(elements_q)

            @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
            prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
            @timeit to "calculate shape functions" set∇𝝭!(elements_φ)

            @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
            @timeit to "calculate shape functions" set∇𝝭!(elements_φ_Γ)

            𝑎ᵠᵠ = ∫κκdΩ=>elements_φ
            𝑎ˢᵠ = ∫QφdΩ=>(elements_φ,elements_q)
            𝑎ᴳᵠᵠ = ∫∇φσ∇φdΩ=>elements_φ
            𝑎ᵐᵠᵠ = ∫ρφφdΩ=>elements_φ

            @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
            @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
            @timeit to "assemble" 𝑎ᴳᵠᵠ(kᴳᵠᵠ)
            @timeit to "assemble" 𝑎ᵐᵠᵠ(mᵠᵠ)
        end

        @timeit to "calculate ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫∇wσ∇wdΩ, ∫ρwwdΩ" begin
            @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], eval(type_q), integrationOrder, sp_q)
            prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
            @timeit to "calculate shape functions" set∇𝝭!(elements_q)

            @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], integrationOrder)
            prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :ρ=>ρ, :σ₁₁=>σ₁₁,:σ₂₂=>σ₂₂,:σ₁₂=>σ₁₂)
            @timeit to "calculate shape functions" set∇𝝭!(elements_w)

            @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], eval(type_q), integrationOrder, sp_q, normal=true)
            @timeit to "calculate shape functions" set∇𝝭!(elements_q_Γ)

            @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], integrationOrder, normal=true)
            @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

            𝑎ˢˢ = ∫QQdΩ=>elements_q
            𝑎ˢʷ = [
            ∫∇QwdΩ=>(elements_q,elements_w),
            ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
            ]
            𝑎ᴳʷʷ = ∫∇wσ∇wdΩ=>elements_w
            𝑎ᵐʷʷ = ∫ρwwdΩ=>elements_w

            @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
            @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
            @timeit to "assemble" 𝑎ᴳʷʷ(kᴳʷʷ)
            @timeit to "assemble" 𝑎ᵐʷʷ(mʷʷ)
        end

        @timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
            @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
            @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
            @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
            @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
            prescribe!(elements_1, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
            prescribe!(elements_2, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
            prescribe!(elements_3, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
            prescribe!(elements_4, :α=>1e8*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
            @timeit to "calculate shape functions" set𝝭!(elements_1)
            @timeit to "calculate shape functions" set𝝭!(elements_2)
            @timeit to "calculate shape functions" set𝝭!(elements_3)
            @timeit to "calculate shape functions" set𝝭!(elements_4)
            𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
            𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
            @timeit to "assemble" 𝑎ʷ(kʷʷ)
            @timeit to "assemble" 𝑎ᵠ(kᵠᵠ)
            # @timeit to "assemble" 𝑎ʷ(mʷʷ)
            # @timeit to "assemble" 𝑎ᵠ(mᵠᵠ)
            # @timeit to "assemble" 𝑎ʷ(kᴳʷʷ)
        end

        gmsh.finalize()

        kᵠᵠ .+= - kˢᵠ'*(kˢˢ\kˢᵠ)
        kᵠʷ .+= - kˢᵠ'*(kˢˢ\kˢʷ)
        kʷʷ .+= - kˢʷ'*(kˢˢ\kˢʷ)

        k = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]
        kᴳ = [kᴳᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) kᴳʷʷ]
        m = [mᵠᵠ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) mʷʷ]

        λ,v = eigen(k,kᴳ)
        # λ,v = eigen(k,m)

        index = findfirst(real.(λ).>1e-8)

        d₁ = real.(v[2nᵠ+1:2nᵠ+nʷ,index])
        d₂ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+1])
        d₃ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+2])
        d₄ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+3])
        d₅ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+4])
        d₆ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+5])
        d₇ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+6])
        d₈ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+7])
        d₉ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+8])
        d₁₀ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+9])
        d₁₁ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+10])
        d₁₂ = real.(v[2nᵠ+1:2nᵠ+nʷ,index+11])
        push!(nodes_w,
            :d₁=>d₁,
            :d₂=>d₂,
            :d₃=>d₃,
            :d₄=>d₄,
            :d₅=>d₅,
            :d₆=>d₆,
            :d₇=>d₇,
            :d₈=>d₈,
            :d₉=>d₉,
            :d₁₀=>d₁₀,
            :d₁₁=>d₁₁,
            :d₁₂=>d₁₂,
        )

    # 二维------------------------------------------------------------------------------
        # xs = [node.x for node in nodes]'
        # ys = [node.y for node in nodes]'
        # zs = [node.z for node in nodes]'
        # points = [xs; ys; zs]
        # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_w]

        # vtk_grid("./vtk/fem/CCCC/fem_CCCC_$n.vtu", points, cells;
        #         ascii=true, append=false, compress=false) do vtk

        #     vtk["v₁"] = [node.d₁ for node in nodes]
        #     vtk["v₂"] = [node.d₂ for node in nodes]
        #     vtk["v₃"] = [node.d₃ for node in nodes]
        #     vtk["v₄"] = [node.d₄ for node in nodes]
        #     vtk["v₅"] = [node.d₅ for node in nodes]
        #     vtk["v₆"] = [node.d₆ for node in nodes]
        #     vtk["v₇"] = [node.d₇ for node in nodes]
        #     vtk["v₈"] = [node.d₈ for node in nodes]
        #     vtk["v₉"] = [node.d₉ for node in nodes]
        #     vtk["v₁₀"] = [node.d₁₀ for node in nodes]
        #     vtk["v₁₁"] = [node.d₁₁ for node in nodes]
        #     vtk["v₁₂"] = [node.d₁₂ for node in nodes]

    # 三维------------------------------------------------------------------------------
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_w]

        names = ["v₁", "v₂", "v₃", "v₄", "v₅", "v₆", "v₇", "v₈", "v₉", "v₁₀", "v₁₁", "v₁₂"]
        modes = [:d₁, :d₂, :d₃, :d₄, :d₅, :d₆, :d₇, :d₈, :d₉, :d₁₀, :d₁₁, :d₁₂]
        for m in 1:12
            s = modes[m]          # 得到 :d₁、:d₂、...、:d₁₂
            vals = Float64[getproperty(node, s) for node in nodes_w]
            scale = 1.0 / maximum(abs, vals)

            nₚ = length(nodes_w)
            points = zeros(3, nₚ)
            for (i, node) in enumerate(nodes_w)
                points[1,i] = node.x
                points[2,i] = node.y
                points[3,i] = getproperty(node, s) * scale
            end

            vtk_grid("./vtk/mix/CCCC/mix_$(i)_CCCC_roit_mode_$m.vtu", points, cells;
                    ascii=false, append=false, compress=false) do vtk
                vtk[names[m]] = vals
            end
        end

        # (λ.*ρ/Dˢ).^0.5
        println(λ[index]*a^2/(π^2*Dᵇ)*h)
        k = (λ[index]*a^2/(π^2*Dᵇ)*h)

        write(io, "$i,$k\n")

    end
end