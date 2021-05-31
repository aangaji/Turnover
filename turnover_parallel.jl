projectdir = "C:/Users/Arman Angaji/OneDrive - Universität zu Köln/Dokumente/Uni-Köln/Masterarbeit/Workspace/Julia_Master/MasterProject_Julia"
using Pkg; Pkg.activate(projectdir)
using Distributed, SharedArrays, SentinelArrays

addprocs(6)
@everywhere using Pkg
@everywhere projectdir = "C:/Users/Arman Angaji/OneDrive - Universität zu Köln/Dokumente/Uni-Köln/Masterarbeit/Workspace/Julia_Master/MasterProject_Julia"
@everywhere Pkg.activate(projectdir)
workers()
# rmprocs(workers())

@everywhere using StatsBase, LinearAlgebra, DataFrames
@everywhere using TumorGrowth
@everywhere include("treeless_turnover.jl")

include(projectdir * "/pipeline_wrapper_parallel.jl")

# @everywhere function orphaned(tumor, params; Nthresh)
#     N, b, d, μ, ρ, dim, id = params
#
#     orphcount = 0
#     mtypes, mutations = clones_by_mutations(tumor; res = 1/Nthresh)
#     freqs = nrow.(mtypes)
#
#     k = 0
#     for (mind,redclone) in enumerate(mtypes)
#         n_clone = freqs[mind]
#         # n_clone < nrow(tumor)/200 && continue
#         greenmuts = filter!(mg -> mg>mutations[mind], unique(vcat(redclone.mutations...)))
#         indices = [ findfirst(isequal(m_green), mutations) for m_green in greenmuts ]
#         n_green = getindex(freqs, filter!(!isnothing, indices) )
#         # filter!(n -> n> nrow(tumor)/200, n_green)
#         isempty(n_green) && continue
#         k += length(n_green)
#         orphcount += count(iszero, n_clone .- n_green ) #/ length(n_green)
#     end
#     println("ρ$ρ d$d turnover: $(orphcount/k)")
#     return orphcount/k
# end

# @everywhere function estranged(tumor, params; Nthresh)
#     N, b, d, μ, ρ, dim, id = params
#
#     countestr = 0
#     mtypes, typemuts = haplotypes(tumor)
#     mutations, freqs = mutation_freqs(tumor) |> df -> (df.mutation, df.frequency)
#
#     countgreens = 0
#     for (mind,greentype) in enumerate(mtypes)
#         length(typemuts[mind]) < 2 && continue
#         green_m = last(typemuts[mind])
#         freqs[findfirst(isequal(green_m), mutations)] < 1/Nthresh && continue
#
#         countgreens += 1
#
#         redtype = typemuts[mind][1:end-1]
#         if isnothing( findfirst(isequal(redtype), typemuts) )
#             countestr += 1
#         end
#     end
#     println("ρ$ρ d$d turnover: $(countestr/countgreens)")
#     return countestr/countgreens
# end

@everywhere function orphaned_wrapper(tumor, params; Nthresh, based_on)

    N, b, d, μ, ρ, dim, id = params
    htumor = DataFrame(mutations = unique(tumor.mutations))
    freqs = mutation_freqs(tumor) |> seq -> Dict(seq.mutation .=> seq.frequency)
    filter!.(m-> freqs[m] > 1/Nthresh, htumor.mutations)
    res = if based_on == :green
        orphaned_green_treeless(htumor)
    elseif based_on == :red
        orphaned_red_treeless(htumor)
    else
        error("based_on is :green or :red")
    end

    turnover = sum(res.isorphaned) / sum(res.isgreen)
    println("ρ$ρ d$d turnover: $turnover")
    return turnover
end

@everywhere function estranged_wrapper(tumor, params; Nthresh)
    N, b, d, μ, ρ, dim, id = params
    htumor = DataFrame(mutations = unique(tumor.mutations))
    freqs = mutation_freqs(tumor) |> seq -> Dict(seq.mutation .=> seq.frequency)
    res = estranged_treeless(htumor)
    filter!(m -> freqs[m.mutation] > 1/Nthresh, res)

    turnover = sum(res.isestranged) / sum(res.isgreen)
    println("ρ$ρ d$d turnover: $turnover")
    return turnover
end

### RUN
[true,false,true] .& [true,false,false]

ρ, d = nothing, nothing
Nthresh, based_on = 200, :green
@time pipeline_wrapper( (t,p) -> orphaned_wrapper(t,p; Nthresh=Nthresh, based_on=based_on);
                        input_dir = projectdir*"/simulated_tumors/2d/bulk", ρ = ρ, d = d,
                        output_path = projectdir*"/turnover/turnover_data/orphaned_N$(Nthresh)_$(based_on).csv")

ρ, d = [0.8, 1.0, 1.3, 1.6, 1.8, 2.0, 2.2, 2.3], nothing
Nthresh = 200
@time pipeline_wrapper( (t,p) -> estranged_wrapper(t,p; Nthresh=Nthresh);
                        input_dir = projectdir*"/simulated_tumors/2d", ρ = ρ, d = d,
                        output_path = projectdir*"/turnover/turnover_data/estranged_rho_less2.4_N$(Nthresh).csv")

# @time pipeline_wrapper( (t,p) -> estranged(t,p; Nthresh=Nthresh);
#                         input_dir = projectdir*"/simulated_tumors/2d/bulk", ρ = ρ, d = d,
#                         output_path = projectdir*"/turnover/turnover_data/estranged_N$(Nthresh)_check.csv")

### LOAD RESULTS

# red
orphdata = DataFrame(CSV.File(projectdir*"/turnover/turnover_data/orphaned_rho_less2.4_N100_red.csv", delim=";"))
orphdata = DataFrame(CSV.File(projectdir*"/turnover/turnover_data/orphaned_N100_red.csv", delim=";"))
#green
orphdata = DataFrame(CSV.File(projectdir*"/turnover/turnover_data/orphaned_N100_green.csv", delim=";"))

b = unique(orphdata.b)[1]
scatter( unique(orphdata.d) ./ b, [ mean(filter(t -> t.d == d, orphdata).result) for d=unique(orphdata.d) ],
        xlims=(0,1), ylims=(0,1), xlab=:d, ylab=:orphaned, lab="N?")

scatter!(orphdata.d ./ b, orphdata.result, xlims=(0.,1.), ylims=(0.,1.), xlab=:d, ylab=:orphaned, lab="", alpha=0.2)



estrdata = DataFrame(CSV.File(projectdir*"/turnover/turnover_data/estranged_N500_new.csv", delim=";"))

b = unique(orphdata.b)[1]
scatter!( unique(estrdata.d) ./ b,
        [ mean(  filter(t -> !isnan(t.result) && t.d == d, estrdata).result  ) for d=unique(estrdata.d) ],
        xlims=(0,1), ylims=(0.,1.), xlab=:d, ylab=:estranged,
        lab="", legend=:topleft)

scatter!(estrdata.d ./ b, estrdata.result, xlims=(0.,0.69), lab="", alpha=0.1)

# savefig(projectdir*"/turnover/turnover_plots/estranged_rho_less2.4_N200.png")
