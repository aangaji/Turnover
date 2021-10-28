function sampletumor_mfreqs(sampletumor)
    mutations = unique( vcat(sampletumor.mutations...) )
    m_ind = Dict(mutations .=> 1:length(mutations))
    freqs = zeros(length(mutations))
    for row in eachrow(sampletumor)
        for (m,f) in zip(row.mutations, row.frequencies)
            freqs[m_ind[m]] += f
        end
    end
    freqs./=nrow(sampletumor)
    return DataFrame(mutation = mutations, frequency=freqs)
end

function get_turnover(tumorinfo; useknown_N = false, useknown_T = false, Nthresh_orph, Nthresh_estr, tumor_sample_func = df -> df, mut_freqs_func = TumorGrowth.mutation_freqs)
    clade_turnover = Float64[]
    clone_turnover = Float64[]

    @showprogress for sim in eachrow(tumorinfo)

        tumor, mutations, b, d= tumor_sample_func( sim.tumor ), sim.mutations, sim.b, sim.d

        htypes = unique(tumor.mutations)
        freqs = if useknown_N
            Dict( 1:nrow(mutations) .=> 1 ./ mutations.N_birth )
        elseif useknown_T
            Dict( 1:nrow(mutations) .=> @. 1 / exp( (b-d)*mutations.t_birth ) )
        else
            mut_freqs_func(tumor) |> seq -> Dict(seq.mutation .=> seq.frequency)
        end

        orphaned_tumor = DataFrame( mutations = unique( filter.(m-> freqs[m] > 1/Nthresh_orph, htypes) ) )
        estranged_tumor = DataFrame(mutations = filter( muts -> all(  freqs[m] > 1/Nthresh_estr for m in muts), htypes) )

        W_a = Turnover.orphaned_red_treeless(orphaned_tumor) |> df -> sum(df.isorphaned)/sum(df.isgreen)

        W_o = Turnover.estranged_treeless(estranged_tumor) |> df -> sum(df.isestranged)/sum(df.isgreen)

        push!(clade_turnover, W_a)
        push!(clone_turnover, W_o)

        sleep(0.01)
    end
    return (ds = tumorinfo.d, Wa = clade_turnover, Wo = clone_turnover)
end

function infer_params( tumorinfo; Nthresh_orph, Nthresh_estr, Wa, Wo)
    dfits = []
    mufits = []

    @showprogress for i in 1:nrow(tumorinfo)
        b, mu = tumorinfo.b[i], tumorinfo.μ[i]
        W_a, W_o = Wa[i], Wo[i]
        
        d_solve = min(1., 2*log(Nthresh_orph)*W_a)*b
        
        mu_solve = missing
        try
#             mu_solve, n = bisection(x -> min(1., Turnover.W_estranged(d_solve; b=b, μ=x, T=log(Nthresh_estr)/(b-d_solve)))- W_o, [0.01,0.99], 100)
            mu_solve = fzero(x -> min(1., Turnover.W_estranged(d_solve; b=b, μ=x, T=log(Nthresh_estr)/(b-d_solve)))- W_o, 0.01, 0.99)
            catch e
        end

        push!(dfits, d_solve)
        push!(mufits, mu_solve)

        sleep(0.01)
    end
    return (ds = tumorinfo.d, dfits = dfits, mufits = mufits)
end

# function plot_infresult(ds, dfits, mufits; mu, size=(800,400))
    
#     p = plot(layout=(1,2), size=size, legend=:none, margin=3Plots.mm, yguidefontrotation=-90,
#              aspect_ratio=1, xaxis=(L"d", (0,1), 0.:0.2:1.), yaxis=((0,1),0:0.2:1), guidefontsize=15, tickfont=15, grid=false)

#     scatter!(p[1], ds, dfits, ylab=L"d_\mathrm{fit}", marker = (5, 0.2, :darkblue), alpha=0.4)
#     plot!(p[1], 0:1,0:1, lw=2, c=:red)
#     scatter!(p[1], unique(ds), [mean(filter(!ismissing, dfits[ds .== d])) for d in unique(ds)], yerror=[std(filter(!ismissing, dfits[ds .== d])) for d in unique(ds)], marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)

#     scatter!(p[2], ds, mufits, ylab=L"\mu_\mathrm{fit}", marker = (5, 0.2, :darkblue), alpha=0.4)
#     hline!(p[2], [mu], lw=2, c=:red)

#     mus = [filter(!ismissing, mufits[ds .== d]) for d in unique(ds)]
#     mask = .!isempty.(mus)

#     scatter!(p[2], unique(ds)[mask], mean.(mus[mask]), 
#         yerror= std.(mus[mask]), marker = (:hline, 14, :darkblue), markerstrokecolor=:darkblue)
# end

# function plot_turnover(ds, Wa, Wo; Nthresh_orph, Nthresh_estr, b=1., mu, plotargs...)
    
#     p = plot(layout=(1,2), size=(600,250), legend=:none, margin=3Plots.mm, xlab=L"d", xlim=(0,1))

#     scatter!(p[1], ds, Wa, ylab=L"W_{a}", alpha=0.2, ms = 3, c=:black, ylim=(0,0.3))
#     scatter!(p[1], unique(ds), [mean(filter(!ismissing, Wa[ds .== d])) for d in unique(ds)], 
#         yerror=[std(filter(!ismissing, Wa[ds .== d])) for d in unique(ds)], ms=10., marker=:hline, c=:black)
#     plot!(p[1], 0.:0.01:b, d-> Turnover.W_orphaned(d/b; N=Nthresh_orph))

#     scatter!(p[2], ds, Wo, ylab=L"W_{o}", alpha=0.5, ms = 3, c=:black, ylim=(0,1))

#     Wo_bins = [filter(!ismissing, Wo[ds .== d]) for d in unique(ds)]
#     mask = .!isempty.(Wo_bins)

#     scatter!(p[2], unique(ds)[mask], mean.(Wo_bins[mask]), 
#         yerror= std.(Wo_bins[mask]), ms=10., marker=:hline, c=:black)
#     plot!(p[2], 0.:0.01:b, d-> min( 1, Turnover.W_estranged(d; b=b, μ=mu, T=log(Nthresh_estr)/(b-d))) )
# end

function plot_turnover_violin(ds, Wa, Wo; Nthresh_orph, Nthresh_estr, b=1., mu, plotargs...)
    scalex=5
    d_uni = unique(ds)
    bins = [findall(isequal(d), ds) for d in d_uni]
    
    p = plot(layout=(1,2), size=(600,250), legend=:none, margin=3Plots.mm, xlab=L"d", xticks=(0:scalex, 0.:1/scalex:1.), plotargs...)

    scatter!(p[1], ds*scalex, Wa, ylab=L"W_{a}", alpha=0.2, ms = 3, c=:darkblue, ylim=(0,0.3))
    violin!(p[1], ds*scalex, Wa, marker = (5, 0.2, :darkblue), alpha=0.4, c=:lightblue)
    scatter!(p[1], d_uni*scalex, [median(Wa[bin]) for bin in bins], marker=(:hline, 3*scalex, :darkblue) )
    plot!(p[1], (0.:0.01:b)*scalex, Turnover.W_orphaned.((0.:0.01:b)/b; N=Nthresh_orph))

    scatter!(p[2], ds*scalex, Wo, ylab=L"W_{o}", alpha=0.5, ms = 3, c=:black, ylim=(0,1))
    violin!(p[2], ds*scalex, Wo, marker = (5, 0.2, :darkblue), alpha=0.4, c=:lightblue)
    scatter!(p[2], d_uni*scalex, [median(Wo[bin]) for bin in bins], marker=(:hline, 3*scalex, :darkblue))
    plot!(p[2], (0.:0.01:b)*scalex, d -> min(1, Turnover.W_estranged(d/scalex/b; b=b, μ=mu, T=log(Nthresh_estr)/(b-d/scalex))) )
end

function plot_infresult_violin(ds, dfits, mufits; mu, size=(600,300), plotargs...)
    scalex=5
    d_uni = unique(ds)
    bins = [findall(isequal(d), ds) for d in d_uni]
    
    p = plot(layout=(1,2), size=size, legend=:none, margin=3Plots.mm, yguidefontrotation=-90,
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