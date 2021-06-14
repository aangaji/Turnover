projectdir = "C:/Users/Arman Angaji/OneDrive - Universität zu Köln/Dokumente/Uni-Köln/Masterarbeit/Workspace/Julia_Master/MasterProject_Julia"
using Pkg; Pkg.activate(projectdir)

using DataFrames, StatsBase, Plots, Statistics, LaTeXStrings, Interact, LsqFit, CSV

include("Turnover.jl")
include("test.jl")
using .Turnover

using TumorGrowth: DataFrame, clones_by_mutations, nonspatial, nonspatial!, data_import

######################################

sim = include("turnover_data/estranged_turnover_v2.jl")
expect = include("turnover_data/expected_estranged_turnover_v2.jl")
sim[1]


#######################################

let data = sim
    ds = data[1][:d]
    mus = getindex.(data, :μ)
    t = [ [ res[:turnover][i] for res in data] for i=1:length(ds)]
    stdev = [ [ res[:std][i] for res in data] for i=1:length(ds)]

    p=plot(xlim=(0., 0.5), ylim=(0,1), xlab=:μ)
    for i=1:4:length(ds)
        plot!(mus, t[i], yerror=stdev[i], lab="d$(ds[i])")
    end
    display(p)
end

#########################################

let data = sim
    ds = data[1][:d]
    mus = getindex.(data, :μ)
    t = [ [ res[:turnover][i] for res in data] for i=1:length(ds)]
    stdev = [ [ res[:std][i] for res in data] for i=1:length(ds)]

    dslider = slider(ds)
    pyplot()
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


let data = expect
    ds = data[1][:d]
    mus = getindex.(data, :μ)
    t = [ [ res[:turnover][i] for res in data] for i=1:length(ds)]
    stdev = [ [ res[:std][i] for res in data] for i=1:length(ds)]

    dslider = slider(ds, label="d")

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

#########################
### Single Simulation ###
#########################

### run
Nthresh = 2000
b, d, mu = params = (b=1., d=0.2, μ=0.3)
out = neutral_growth(Nthresh; params..., return_obs=false)
cutoff = length(out[:tumor])

orphtumor = neutral_growth!(deepcopy(out[:tumor]), 200000; Nthresh=Nthresh, params...)[:tumor] |> DataFrame
out = neutral_growth!(deepcopy(out[:tumor]), 200000; params...)
tumor = out[:tumor][1:cutoff] |> DataFrame


### save
# CSV.write("turnover_data/simulated_tumors/tumor_N1000000_b$(b)_d$(d)_μ$(mu)_Nthresh$(Nthresh).csv", tumor)
#
#
# ### load
# b, d, mu, Nthresh = params = (b=1., d=0.5, μ=0.4, Nthresh=2000)
# tumor = data_import("turnover_data/simulated_tumors/tumor_N1000000_b$(b)_d$(d)_μ$(mu)_Nthresh$(Nthresh).csv", delim=",")


Ls = 0.1:0.1:1.
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
# t = [ W_estranged(d; b=b, μ=L*mu, T=log(Nthresh)/(b-d) ) .+ 0.1*randn(50) for L in Ls]

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
        global b, Nthresh
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


##########################
### known death rate d ###
##########################

let d = 0.25
    function model(Ls, p)
        global b, Nthresh
        mu = p[1]
        (mu<=0. || mu>=1.) && return fill(Inf, length(Ls))
        map( Ls ) do L
            min(1., W_estranged(d; b=b, μ=mu*L, T=log(Nthresh)/(b-d)))
        end
    end
    fit = curve_fit(model, Ls, mean.(t), [0.5])
    mu_fit = fit.param[1]

    p1 = scatter(Ls, mean.(t), ribbon=std.(t), yerror = std.(t)./sqrt.(length.(t)), lab="", xlab=:L, ylab=:turnover, legend=:bottomright)
    plot!(Ls, model(Ls, [mu]), lab="" )
    plot!(Ls, model(Ls, fit.param), lab="fit : μ$(round(mu_fit,digits=3))" )

    tvecs = [ rand.(t) for _=1:1000 ]
    fits = map( tp-> curve_fit(model, Ls, tp, [0.5]), tvecs)
    mu_fits = getindex.(getfield.(fits, :param),1)
    p2 = histogram(mu_fits, xlab="μ", lab="")
    vline!([mu_fit], lab="")

    println(" ")
    println("μ : ", mu_fit)
    println("μ : ", mean(mu_fits), " ± ", 1.96*std(mu_fits)/sqrt(length(mu_fits)))
    plot(p1, p2, layout=(2,1), size=(400,600))
end

let tumors = map( Ls ) do L
        [filter(h->!iszero(h.n), reduced_μ(orphtumor, L) ) for _=1:reps]
    end

    t = map(tumors) do set
        filter!(!isnan, map(set) do tumor
            res = orphaned_red_treeless(tumor)
            sum(res.isorphaned)/sum(res.isgreen)
        end )
    end

    qvals = map(t) do Ws
            2*log(Nthresh).*Ws
        end

    global d, b, Nthresh
    p = scatter(Ls, mean.(t), ribbon=std.(t), lab="", xlab=:L, ylab=:turnover)
    plot!(Ls, L-> W_orphaned(d/b; N=Nthresh),
        lab="W$(round(W_orphaned(d/b; N=Nthresh),digits=3)) d$d" )
    plot!(Ls, mean.(t), yerror=std.(t)./sqrt.(length.(t)), lab="" )

    println("  ")
    println("W : ", mean(mean.(t)), " ± ", 1.96*std(mean.(t))/sqrt(length(t)) )
    println("q : ", mean(mean.(qvals)), " ± ", 1.96*std(mean.(qvals))/sqrt(length(t)) )
    display(p)
end


using Distributions, NLsolve, Roots
z = nlsolve(x->cdf.(Normal(), x) .- 0.975, [1.]).zero[1]

let W_l = orphaned_red_treeless(filter!(h->!iszero(h.n), orphtumor)) |>
                df -> sum(df.isorphaned)/sum(df.isgreen)

    qval = 2*log(Nthresh)*W_l
    println(qval)
end

function bisection(f, interval, iter; precision=1e-9)
    xlow, xup = interval
    sign(f(xlow)) == sign(f(xup)) && error("no zero in interval $interval")
    any((isnan.(f.(interval)))) && error("f not defined on interval $interval")
    for i=1:iter
        xm = (xlow+xup)/2
        abs(f(xm)) < precision && return (xm, i)
        if sign(f(xm)) != sign(f(xlow))
            xup = xm
        else
            xlow = xm
        end
    end
    println(xup)
    error("no convergence within $iter steps")
end

let d = 0.25

    W_c = estranged_treeless(filter!(h->!iszero(h.n), tumor)) |> df -> sum(df.isestranged)/sum(df.isgreen)

    function model(mu)
        global b, Nthresh
        min(1., W_estranged(d; b=b, μ=mu, T=log(Nthresh)/(b-d)))- W_c
    end

    mu_solve, n = bisection(model, [0.01,0.99], 100)
    # mu_solve = fzero(model, 0.01, 1.)

    println(" ")
    println("μ : ", mu_solve)
    p1 = plot(0.:0.001:1., model,
        size=(500,400), legend=:none, c=:black, xlab=:μ, ylab=L"W_c")
    hline!([0.], c=:black); vline!([mu_solve], c=:red)

    p2 = scatter(Ls, mean.(t), ribbon=std.(t), yerror = std.(t)./sqrt.(length.(t)), lab="", xlab=:L, ylab=:turnover, legend=:bottomright)
    plot!(Ls, L-> min(1., W_estranged(d; b=b, μ=L*mu_solve, T=log(Nthresh)/(b-d))), lab="" )

    plot(p1, p2)
end
