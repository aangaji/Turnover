####################
## TURNOVER TESTS ##
####################

using Test

using .Turnover
    @testset "Applying turnover methods" begin

        using DataFrames
        using TumorGrowth: DataFrame, clones_by_mutations
    
        b, d, μ = 1., 0.8, 0.1
        @testset "simulation" begin
            
            out = neutral_growth(100 ; b = b, d = d, μ = μ, return_obs=true, showprogress=false)
            neutral_growth!(out[:tumor], out[:obs], 120 ; b = b, d = d, μ = μ, showprogress=false)
            neutral_growth!(deepcopy(out[:tumor]), 130 ; b = b, d = d, μ = μ, showprogress=false)
            
            @test !(:obs in keys( neutral_growth(10 ; b = b, d = d, μ = μ, return_obs=false, showprogress=false) ) )
        end
    
        out = neutral_growth(1000 ; b = b, d = d, μ = μ, return_obs=true, showprogress=false)
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

        end
    end