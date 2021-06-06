using Pkg; Pkg.activate("C:/Users/Arman Angaji/OneDrive - Universität zu Köln/Dokumente/Uni-Köln/Masterarbeit/Workspace/Julia_Master/MasterProject_Julia")

using DataFrames, StatsBase, Plots, Statistics, LaTeXStrings, Interact, LsqFit

include("Turnover.jl")
include("test.jl")
using .Turnover

using TumorGrowth: DataFrame, clones_by_mutations, nonspatial, nonspatial!

#######################################

sim = include("turnover_data/estranged_turnover_v2.jl")
expect = include("turnover_data/expected_estranged_turnover_v2.jl")
sim[1]


#######################################
begin data = expect
    ds = data[1][:d]
    mus = getindex.(data, :μ)
    t = [ [ res[:turnover][i] for res in data] for i=1:length(ds)]
    stdev = [ [ res[:std][i] for res in data] for i=1:length(ds)]
end

let p=plot(xlim=(0., 0.5), ylim=(0,1))
    for i=1:4:length(ds)
        plot!(mus, t[i], yerror=stdev[i], lab="")
    end
    display(p)
end

#########################################

let dslider = slider(ds)

    function model(mus, p)
        d = p[1]
        map( mus ) do mu
            min(1., W_estranged(d; b=1., μ=mu, T=log(200)/(1-d)))
        end
    end

    @manipulate for d in dslider
        i = dslider[:index].val
        fit = curve_fit(model, mus, t[i], [0.5])

        plot(ylim=(0,1), xlab=" μ")
        scatter!(mus, t[i], yerror=stdev[i], lab="")
        plot!(mus, model(mus, d), lab="d$d", style=:dash)
        d_fit = round(first(fit.param), digits=3)
        plot!(mus, model(mus, fit.param), lab="fit d$d_fit")
    end
end


let dslider = slider(ds, label="d")

    function model(Ls, p)
        d, mu = p
        (d<0. || d>=1. || mu<=0. || mu>=1.) && return fill(Inf, length(Ls))
        map( Ls ) do L
            min(1., W_estranged(d; b=1., μ=mu*L, T=log(200)/(1-d)))
        end
    end

    mup = first(mus)
    Ls = mus./mup

    fits = map( tp-> curve_fit(model, Ls, tp, [0.5, 0.5]), t)
    d_err = abs.(getindex.(getfield.(fits, :param),1) .- ds)
    mu_err = abs.(getindex.(getfield.(fits, :param),2) .- mup)

    @manipulate for d in dslider
        i = dslider[:index].val
        fit = fits[i]

        p = plot(layout=(2,1), size=(400,600), legend=:topleft, margin=5Plots.mm)
        plot!(p[1], xlab="L", ylab=:turnover, ylims=(0,1) )
        scatter!(p[1], Ls, t[i], yerror=stdev[i], lab="")
        plot!(p[1], Ls, model(Ls, [ds[i], first(mus)]), lab= "d$(ds[i]) μ$mup", style=:dash)
        d_fit, mu_fit = round.(fit.param, digits=3)
        plot!(p[1], Ls, model(Ls, fit.param), lab="fit d$d_fit μ$mu_fit")

        plot!(p[2], xlab="d", ylab=:error )
        plot!(p[2], ds, d_err, lab="", c=:blue)
        plot!(p[2], ds, mu_err, lab="", c=:red)
        scatter!(p[2], [ds[i]], [d_err[i]], lab="d$d_fit", c=:blue)
        scatter!(p[2], [ds[i]], [mu_err[i]], lab="μ$mu_fit", c=:red)
    end
end

#################################

b, d, mu, Nthresh = params = (b=1., d=0.5, μ=0.4, Nthresh=2000)
out = neutral_growth(Nthresh; params..., return_obs=false)
cutoff = length(out[:tumor])
mu*b/(b-d)*Nthresh

out = neutral_growth!(out[:tumor], 1000000; params...)

tumor = DataFrame(out[:tumor][1:cutoff])

function reduced_μ!(htumor, x)
    reduced = vcat(htumor.mutations...) |> unique! |> sort! |>
        muts -> filter!(muts) do m
            rand()<=x
        end
    for (i,muts) in enumerate(htumor.mutations)
        filter!(in(reduced), muts)
        n = htumor.n[i]
        htumor.n[i] = 0
        if isempty(muts)
            htumor.n[1] += n
        else
            htumor.n[last(muts)+1] += n
        end
    end
    htumor
end
reduced_μ(htumor, x) = reduced_μ!(deepcopy(htumor), x)

Ls = 0.3:0.02:0.5
reps = 50
tumors = map( Ls ) do L
    [filter(h->!iszero(h.n), reduced_μ(tumor, L) ) for _=1:reps]
end

t = map(tumors) do set
        filter!(!isnan, map(set) do tumor
            res = estranged_treeless(tumor)
            sum(res.isestranged)/sum(res.isgreen)
        end )
    end


scatter(Ls, mean.(t), yerror=std.(t), lab="", xlab=:L, ylab=:turnover)
plot!(Ls, L -> W_estranged(d; b=b, μ=L*mu, T=log(Nthresh)/(b-d) ), lab="" )
plot!(Ls, L -> W_estranged(d; b=b, μ=L*mu, T=Inf ), lab="" )

### theory + low noise
t = [ W_estranged(d; b=b, μ=L*mu, T=log(Nthresh)/(b-d) ) .+ 0.1*randn(50) for L in Ls]

let
    function model(Ls, p)
        d, mu = p
        (d<0. || d>=1. || mu<=0. || mu>=1.) && return fill(Inf, length(Ls))
        map( Ls ) do L
            min(1., W_estranged(d; b=b, μ=mu*L, T=log(Nthresh)/(b-d)))
        end
    end
    fit = curve_fit(model, Ls, mean.(t), [0.5, 0.5])
    d_err, mu_err = abs.(fit.param .- (d, mu) )

    p = scatter(Ls, mean.(t), yerror=std.(t), lab="", xlab=:L, ylab=:turnover)
    plot!(Ls, model(Ls, [d, mu]), lab="" )
    plot!(Ls, model(Ls, fit.param), lab="" )
    println(fit.param)
    display(p)
end

let
    function model(Ls, p)
        d, mu = p
        (d<0. || d>=1. || mu<=0. || mu>=1.) && return fill(Inf, length(Ls))
        map( Ls ) do L
            min(1., W_estranged(d; b=b, μ=mu*L, T=log(Nthresh)/(b-d)))
        end
    end

    tvecs = [ rand.(t) for _=1:100 ]

    fits = map( tp-> curve_fit(model, Ls, tp, [0.5, 0.5]), tvecs)
    d_fits, mu_fits = getindex.(getfield.(fits, :param),1), getindex.(getfield.(fits, :param),2)
    d_err, mu_err = abs.(d_fits .- d), abs.(mu_fits .- mu)

    println(" ",minimum(d_fits), " ", maximum(mu_fits))

    @manipulate for i in slider(1:length(fits))
        fit = fits[i]

        p = plot(layout=(2,1), size=(400,600), legend=:topleft, margin=5Plots.mm)
        plot!(p[1], xlab="L", ylab=:turnover, ylims=(0,1) )
        scatter!(p[1], Ls, tvecs[i], lab="")
        plot!(p[1], Ls, model(Ls, [d, mu]), lab= "d$d μ$mu", style=:dash)
        d_fit, mu_fit = round.(fit.param, digits=3)
        plot!(p[1], Ls, model(Ls, fit.param), lab="fit d$d_fit μ$mu_fit")

        plot!(p[2], xlab="d", ylab=:error )
        plot!(p[2], 1:length(fits), d_err, lab="", c=:blue)
        plot!(p[2], 1:length(fits), mu_err, lab="", c=:red)
        scatter!(p[2], [i], [d_err[i]], lab="d$d_fit", c=:blue)
        scatter!(p[2], [i], [mu_err[i]], lab="μ$mu_fit", c=:red)
    end
end
