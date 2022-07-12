using StaticArrays, Distributions, ProgressMeter, DataFrames

export neutral_growth_clones, neutral_growth_clones!

function neutral_growth_clones!(tumor::Vector{Haplotype}, clonesizes::Vector{Int}, obs::Vector{Vector{Int}}, endsize::Int, T=Inf; b, d, μ, t=0., Nthresh=endsize, Tthresh = T, showprogress=true )

    mutID = length(tumor)-1
    F = cumsum(getfield.(tumor,:n))
    N = last(F)

    if showprogress prog = ProgressUnknown("Progress: N ") end
    while 0 < N < endsize && t < T
        μp = N < Nthresh && t < Tthresh ? μ : 0.

        p = rand()*(b+d)
        m = rand(1:N)
        row = findfirst(i-> i >= m, F)
        parent = tumor[row]

        t += rand(Exponential(1.))/((b+d)*N)
        if p > d

            mnew =  (rand()<μp/2) + (rand()<μp/2)
            for n=1:mnew
                mutID += 1
                push!(tumor, Haplotype(push(parent.mutations, mutID), 1, t, row) )
                push!(F, N+n)

                push!(clonesizes, 1)

                push!(obs, clonesizes[parent.mutations])
            end
            if mnew == 0
                parent.n += 1
                F[row:end] .+= 1
            elseif mnew == 2
                parent.n -= 1
                F[row:end] .-= 1
            end
            clonesizes[parent.mutations] .+= 1
            N += 1

        else
            parent.n -= 1
            F[row:end] .-= 1
            N -= 1

            clonesizes[parent.mutations] .-= 1
        end

        if showprogress ProgressMeter.update!(prog, N) end
    end

    return (tumor=tumor, clonesizes=clonesizes, obs=obs, N=N, t=t)
end

function neutral_growth_clones(endsize, T=Inf; args... )
    mutID = 0
    tumor = [Haplotype(SVector{0, Int}(), 1, 0., 0 )]
    obs = Vector{Vector{Int}}()
    clonesizes = Int[]

    out = neutral_growth_clones!(tumor, clonesizes, obs, endsize, T; args... )
    return iszero(out[:N]) ? neutral_growth_clones(endsize, T; args... ) : out
end

# out = neutral_growth(100; b=1., d=0.2, μ=0.3)
#
# out.clonesizes
#
# tumor = out.tumor |> DataFrame
