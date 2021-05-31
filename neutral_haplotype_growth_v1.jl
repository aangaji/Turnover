using StaticArrays, Distributions, ProgressMeter

"""
    Haplotype(mutations, n, t_birth, parent)

Type used by the `neutral_growth` simulation. Fields are `mutations` an ordered list of mutations defining the haplotype,
`n` the number of cells of that type, `t_birth` the time of birth of the first cell of this type and `parent` the rowindex of the parental haplotype
within the tumor Array.
"""
mutable struct Haplotype
    mutations :: SVector{L, Int} where L
    n :: Int
    t_birth :: Float64
    parent :: Int
end

"""  
    neutral_growth!(tumor, obs, endsize, T=Inf; b, d, μ, t=0., Nthresh=endsize, Tthresh = Inf )

Neutral (exponential) growth simulation of `Haplotypes`. Receives a Vector of `Haplotypes`, an according `obs` Matrix, an `endsize` and final time `T` value.
`obs` records the following observables (**currently disabled for performance**): <br>
* obs[m,1] records the parents haplotype n at birth of mutation m <br>
* obs[m,2] records the sum of n of all haplotypes containing the parents last mutation while not
         containing m <br>
         this includes estranged parents while avoiding offspring of m <br>
* obs[m,3] records tumor size at birth of m <br>
Kwargs are `b` birth-rate, `d` death-rate, `μ` mutation-rate, `t` current time, `Nthresh` and `Tthresh` size and time after which `μ` is set to zero.

Output is a named Tuple `(:tumor, :obs, :N, :t)` containing the final `Haplotype` Vector, observables, endsize and final time.
"""
function neutral_growth!(tumor, obs, endsize, T=Inf; b, d, μ, t=0., Nthresh=endsize, Tthresh = Inf )

    mutID = length(tumor)-1
    obsarr = Vector{Vector{Int}}()
    obstemp = Vector{Int}(undef, 3)
    F = cumsum(getfield.(tumor,:n))
    N = last(F)

    prog = ProgressUnknown("Progress: N ")
    while 0 < N < endsize && t < T
        μp = N < Nthresh && t < Tthresh ? μ : 0.

        p = rand()*(b+d)
        m = rand(1:N)
        row = findfirst(i-> i >= m, F)
        parent = tumor[row]

        t += rand(Exponential(1.))/((b+d)*N)
        if p > d

            if rand() < μp

                mutID += 1
                push!(tumor, Haplotype(push(parent.mutations, mutID), 1, t, row) )
                push!(F, N+1)

#                 obstemp[1] = parent.n
#                 obstemp[2] = isone(row) ? N : sum( getfield.(filter(h-> row-1 in h.mutations, tumor[row : mutID]),:n) )
#                 obstemp[3] = N
#                 push!(obsarr, copy(obstemp) )
            else
                parent.n += 1
                F[row:end] .+= 1
            end
            N += 1

        else
            parent.n -= 1
            F[row:end] .-= 1
            N -= 1
        end

        ProgressMeter.update!(prog, N)
    end
    obsnew = Matrix{Int}(undef, length(obsarr), 3)
#     for i=1:3
#         obsnew[:, i] = getindex.(obsarr, i)
#     end
    return (tumor=tumor, obs=vcat(obs, obsnew), N=N, t=t)
end

"""
    neutral_growth(endsize, T; args... )

Neutral (exponential) growth simulation of `Haplotypes`. Receives `endsize` and final time `T`. See documentation of `neutral_growth!`
for kwargs.

Output is a named Tuple `(:tumor, :obs, :N, :t)` containing the final `Haplotype` Vector, observables, endsize and final time.
"""
function neutral_growth(endsize, T=Inf; args... )
    mutID = 0
    tumor = [Haplotype(SVector{0, Int}(), 1, 0., 0 )]
    obs = Matrix{Int}(undef, 0, 3)

    out =  neutral_growth!(tumor, obs, endsize, T;  args...)
    return iszero(out[:N]) ? neutral_growth(endsize, T; args... ) : out
end

"""
    q(d; b, μ )
Function returning the extinction probability for given `d` death-rate, `b` birth-rate and `μ` mutation-rate. <br>
**Branchingprocess `v1`**: only mutation of single randomly picked offspring possible at division. `μ` is the probability of mutation per cell division.

"""
q(d; b, μ ) = d / (b*(1-μ))

;
