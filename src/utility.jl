# function bisection(f, interval, iter; precision=1e-9)
#     xlow, xup = interval
#     sign(f(xlow)) == sign(f(xup)) && error("no zero in interval $interval")
#     any((isnan.(f.(interval)))) && error("f not defined on interval $interval")
#     for i=1:iter
#         xm = (xlow+xup)/2
#         abs(f(xm)) < precision && return (xm, i)
#         if sign(f(xm)) != sign(f(xlow))
#             xup = xm
#         else
#             xlow = xm
#         end
#     end
#     println(xup)
#     error("no convergence within $iter steps")
# end

function estim_N( tumor_slice )
    r = mean( minimum( pairwise(norm ∘ -, tumor_slice.position) + I*Inf, dims=2 ) )
    n = nrow(tumor_slice)
    R = sqrt( n * r^2 / (π/(2*√3)) )
    N = (R/r)^3 * π/(3*√2)
    return N
end

function get_turnover(tumorinfo; 
        useknown_N = false, useknown_T = false, f_min, 
        tumor_sample_col = :tumor,
        tumor_sample_func = df -> df, mut_freqs_func = TumorGrowth.mutation_freqs, subsample_func = df -> df)
    
    clade_turnover = Float64[]
    clone_turnover = Float64[]
    tumorinfo.__sample = similar(tumorinfo.tumor)
    
    @showprogress for sim in eachrow(tumorinfo)

        tumor = getproperty(sim, tumor_sample_col)
        tumor = tumor_sample_func( tumor ) 
        
        mutations, b, d= sim.mutations, sim.b, sim.d
        sim.__sample = tumor

        htypes = unique(tumor.mutations)
        survivors = unique!(vcat(htypes...))
        freqs = if useknown_N
            Dict( survivors .=> 1 ./ mutations.N_birth[survivors] )
        elseif useknown_T
            Dict( survivors .=> @. 1 / exp( (b-d)*mutations.t_birth[survivors] ) )
        else
            mut_freqs_func(tumor) |> seq -> Dict(seq.mutation .=> seq.frequency)
        end

        orphaned_tumor = DataFrame( mutations = unique( filter.(m-> freqs[m] > f_min, htypes) ) ) |> subsample_func
        estranged_tumor = DataFrame(mutations = filter( muts -> all(  freqs[m] > f_min for m in muts), htypes) ) |> subsample_func

        W_a = Turnover.orphaned_red_treeless(orphaned_tumor) |> df -> sum(df.isorphaned)/sum(df.isgreen)

        W_o = Turnover.estranged_treeless(estranged_tumor) |> df -> sum(df.isestranged)/sum(df.isgreen)

        push!(clade_turnover, W_a)
        push!(clone_turnover, W_o)

        sleep(0.01)
    end
    return (ds = tumorinfo.d, Wa = clade_turnover, Wo = clone_turnover)
end

function infer_params( tumorinfo; N, Wa, Wo, usecorrection=true, estimate_N = true, tumor_sample_col = :__sample,)
    dfits = []
    mufits = [] 

    @showprogress for i in 1:nrow(tumorinfo)
        b, mu = tumorinfo.b[i], tumorinfo.μ[i]
        W_a, W_o = Wa[i], Wo[i]
        
        if estimate_N
           N = estim_N(tumorinfo[!,tumor_sample_col][i]) 
        end
        
        corr = d -> usecorrection ? (1-d/b) : 1
        
        # d_solve = min(1., 2*log(N)*W_a)*b
        d_solve = missing
        try
            d_solve = fzero(x -> min(1., Turnover.W_orphaned(x/b; N=N*corr(x))) - W_a, 0.01, 0.99)
        catch e
        end
        
        mu_solve = missing
        try
#             mu_solve, n = bisection(x -> min(1., Turnover.W_estranged(d_solve; b=b, μ=x, T=log(N)/(b-d_solve)))- W_o, [0.01,0.99], 100)
            mu_solve = fzero(x -> min(1., 
                    Turnover.W_estranged(d_solve; b=b, μ=x, T=log(N*corr(d_solve))/(b-d_solve))
                    ) - W_o, 0.01, 0.99)
        catch e
        end

        push!(dfits, d_solve)
        push!(mufits, mu_solve)

        sleep(0.01)
    end
    return (ds = tumorinfo.d, dfits = dfits, mufits = mufits)
end

function plot_turnover_violin(ds, Wa, Wo; N, b=1., mu, usecorrection = true, plotargs...)
    scalex=5
    d_uni = unique(ds)
    bins = [findall(isequal(d), ds) for d in d_uni]
    
    corr = d -> usecorrection ? (1 - d/b) : 1 
    
    p = plot(; layout=(1,2), size=(600,250), legend=:none, margin=3Plots.mm, xlab=L"d", xticks=(0:scalex, 0.:1/scalex:1.), plotargs...)

    scatter!(p[1], ds*scalex, Wa, ylab=L"W_{a}", alpha=0.2, ms = 3, c=:darkblue, ylim=(0,0.3))
    violin!(p[1], ds*scalex, Wa, marker = (5, 0.2, :darkblue), alpha=0.4, c=:lightblue)
    scatter!(p[1], d_uni*scalex, [median(Wa[bin]) for bin in bins], marker=(:hline, 3*scalex, :darkblue) )
    plot!(p[1], (0.:0.01:b)*scalex, d -> Turnover.W_orphaned(d/scalex/b; N=N * corr(d/scalex/b)))

    scatter!(p[2], ds*scalex, Wo, ylab=L"W_{o}", alpha=0.5, ms = 3, c=:black, ylim=(0,1))
    violin!(p[2], ds*scalex, Wo, marker = (5, 0.2, :darkblue), alpha=0.4, c=:lightblue)
    scatter!(p[2], d_uni*scalex, [median(Wo[bin]) for bin in bins], marker=(:hline, 3*scalex, :darkblue))
    
    plot!(p[2], (0.:0.01:b)*scalex, d -> min(1,
            Turnover.W_estranged(d/scalex/b; b=b, μ=mu, T=log(N * corr(d/scalex) )/(b-d/scalex)))
    )
end

function plot_infresult_violin(ds, dfits, mufits; mu, size=(600,300), plotargs...)
    scalex=5
    d_uni = unique(ds)
    bins = [findall(isequal(d), ds) for d in d_uni]
    
    p = plot(; layout=(1,2), size=size, legend=:none, margin=3Plots.mm, yguidefontrotation=-90,
             aspect = scalex, xlab=L"d", yaxis=((0,1),0:0.2:1), xticks=(0:scalex, 0.:1/scalex:1.), plotargs...)

    scatter!(p[1], ds*scalex, dfits, ylab=L"d_\mathrm{fit}", marker = (5, 0.2, :darkblue), alpha=0.4)
    violin!(p[1], ds*scalex, dfits, marker = (5, 0.2, :darkblue), alpha=0.4, c=:lightblue)
    scatter!(p[1], d_uni*scalex, [median(dfits[bin]) for bin in bins], marker=(:hline, 3*scalex, :darkblue) )
    plot!(p[1], [0, scalex], [0,1], lw=2, c=:red)
    
    scatter!(p[2], ds*scalex, mufits, ylab=L"\mu_\mathrm{fit}", marker = (5, 0.2, :darkblue), alpha=0.4)
    violin!(p[2], ds*scalex, mufits, marker = (5, 0.2, :darkblue), alpha=0.4, c=:lightblue)
    scatter!(p[2], d_uni*scalex, [median(mufits[bin]) for bin in bins], marker=(:hline, 3*scalex, :darkblue))
    hline!(p[2], [mu], lw=2, c=:red)
end




function plot_series_violin(ds, dfits, mufits; mu, size=(800,300))
    bins = sort!(unique(ds))
    dfits_binned = [ [dfits[findall(isequal(d), ds)]...] for d in bins]
    mufits_binned = [ [mufits[findall(isequal(d), ds)]...] for d in bins]
    length.(dfits_binned) |> println
    
    scalex = 5
    p = plot(layout=(1,3), legend=:none, yguidefontrotation=-90, size=size,
             aspect_ratio=scalex, xaxis=(L"b", (0,scalex)), yaxis=((0,1),0:0.2:1), guidefontsize=15, tickfont=15, xticks=(0:scalex, range(0,1,length=scalex+1)))

    plot!(p[1], 0:scalex,range(0,1,length=scalex+1), c=:red, lw=2.)
    violin!(p[1], ds*scalex, dfits, ylab=L"b_{inf}", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[1], ds*scalex, dfits, ylab=L"b_{inf}", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[1], bins*scalex, mean.(dfits_binned), yerror = std.(dfits_binned), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)
    
    hline!(p[2], [mu], c=:red, lw=2.)
    violin!(p[2], ds*scalex, mufits, ylab=L"\mu_{inf}", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[2], ds*scalex, mufits, ylab=L"\mu_{inf}", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[2], bins*scalex, mean.(mufits_binned), yerror = std.(mufits_binned), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)
    
    mubeta = [mufits_binned[i]./(1 .- dfits_binned[i]) for i =1:length(bins)]
    plot!(p[3], 0:scalex,range(0,1,length=scalex+1), c=:red, lw=2., xlab=L"\mu/(1-b/a)", yguidefontrotation=-0)
    violin!(p[3], mu./(1 .-ds)*scalex, mufits./(1 .-dfits), ylab=L"\mu_{inf}/(1-b_{inf}/a)", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[3], mu./(1 .-ds)*scalex, mufits./(1 .-dfits), ylab=L"\mu_{inf}/(1-b_{inf}/a)", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[3], mu./(1 .-bins)*scalex, mean.(mubeta), yerror = std.(mubeta), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue) 
end

function plot_series_scatter(ds, dfits, mufits; mu, size=(800,300))
    bins = sort!(unique(ds))
    dfits_binned = [ [dfits[findall(isequal(d), ds)]...] for d in bins]
    mufits_binned = [ [mufits[findall(isequal(d), ds)]...] for d in bins]
    length.(dfits_binned) |> println
    
    p = plot(layout=(1,3), legend=:none, yguidefontrotation=-90, size=size,
             aspect_ratio=1, xaxis=(L"b", (0,1), 0.:0.2:1.), yaxis=((0,1),0:0.2:1), guidefontsize=15, tickfont=15)

    plot!(p[1], 0:1,0:1, c=:red, lw=2.)
    scatter!(p[1], ds, dfits, ylab=L"b_{inf}", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[1], bins, mean.(dfits_binned), yerror = std.(dfits_binned), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)
    
    hline!(p[2], [mu], c=:red, lw=2.)
    scatter!(p[2], ds, mufits, ylab=L"\mu_{inf}", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[2], bins, mean.(mufits_binned), yerror = std.(mufits_binned), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)
    
    mubeta = [mufits_binned[i]./(1 .- dfits_binned[i]) for i =1:length(bins)]
    plot!(p[3], 0:1, 0:1, c=:red, lw=2., xlab=L"\mu/(1-b/a)", yguidefontrotation=-0)
    scatter!(p[3], mu./(1 .-ds), mufits./(1 .-dfits), ylab=L"\mu_{inf}/(1-b_{inf}/a)", marker = (5, 0.2, :darkblue), alpha=0.4)
    scatter!(p[3], mu./(1 .-bins), mean.(mubeta), yerror = std.(mubeta), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)
end