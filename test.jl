####################
## TURNOVER TESTS ##
####################
# projectdir = "C:/Users/Arman Angaji/OneDrive - Universität zu Köln/Dokumente/Uni-Köln/Masterarbeit/Workspace/Julia_Master/MasterProject_Julia"
# using Pkg; Pkg.activate(projectdir)
# include("Turnover.jl")

using Test

using .Turnover
@testset "Applying turnover methods" begin

    using DataFrames

    b, d, μ = params = (b=1., d=0.8, μ=0.1)
    @testset "simulation" begin

        out = neutral_growth(100 ; b = b, d = d, μ = μ, return_obs=true, showprogress=false)
        neutral_growth!(out[:tumor], out[:obs], 120 ; params..., showprogress=false)
        neutral_growth!(deepcopy(out[:tumor]), 130 ; params..., showprogress=false)

        @test !(:obs in keys( neutral_growth(10 ; params..., return_obs=false, showprogress=false) ) )

        @test begin
            out_clones = neutral_growth_clones(100; params..., showprogress=false )
            t = DataFrame(out_clones.tumor)
            all( round.(mfreqs(t).*sum(t.n), digits=3) .== out_clones.clonesizes )
        end
    end

    out = neutral_growth(1000 ; params..., return_obs=true, showprogress=false)
    htumor = DataFrame(out[:tumor])

    @testset "estranged" begin

        estr_res = estranged(htumor) |> df -> (count(df.isestranged), nrow(df))
        estr_treeless_res = filter(h-> !iszero(h.n), htumor) |> estranged_treeless |> df -> (sum(df.isestranged), sum(df.isgreen))
        @test estr_res == round.(Int, estr_treeless_res)

        q = (d + b*(μ/2)^2) / (b*(1-(μ/2))^2)
        estranged_expected(htumor, q, out[:obs][:,1])

    end

    @testset "orphaned" begin
        orph_res = orphaned_green(htumor) |> df -> ( count( df.isorphaned ) , nrow(df) )
        orph_red_res = orphaned_red(htumor) |> df -> ( count( unique(df).isorphaned ) , length(unique(df.mutation)) )
        @test orph_res == orph_red_res

        orph_res = orphaned_green(htumor) |> df -> ( count( df.isorphaned ) , nrow(df) )
        orph_treeless_res = filter(h-> !iszero(h.n), htumor) |> orphaned_green_treeless |> df -> (sum(df.isorphaned), sum(df.isgreen))
        @test orph_res == round.(Int, orph_treeless_res)

        orph_res = orphaned_red(htumor) |> df -> ( count( df.isorphaned ) , nrow(df) )
        orph_treeless_res = filter(h-> !iszero(h.n), htumor) |> orphaned_red_treeless |> df -> (sum(df.isorphaned), sum(df.isgreen))
        @test orph_res == round.(Int, orph_treeless_res)

        orphaned_green_expected(htumor, d/b, out[:obs][:,2])

        htumor, clonesizes, nbirth, N, t = neutral_growth_clones(1000; params..., showprogress=false)
        @test isapprox( mean( (d/b) .^ vcat(nbirth[.!iszero.(clonesizes)]...) ),
                        mean( orphaned_red_expected(d/b, clonesizes, nbirth).isorphaned ) )

    end
end
