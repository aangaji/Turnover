module Turnover

    using DataFrames
    using TumorGrowth: clones_by_mutations

    export mfreqs, reduced_μ, reduced_μ!, remove_mutations, remove_mutations!
    export W_orphaned, W_estranged
    export estranged, estranged_treeless, estranged_expected
    export orphaned_red, orphaned_green, orphaned_red_treeless, orphaned_green_treeless, orphaned_red_expected, orphaned_green_expected

    ##########################
    ### Neutral Simulation ###
    ##########################

    include("neutral_haplotype_growth_v2.jl")
#     include("neutral_haplotype_growth_v2_removedead.jl")
    include("neutral_haplotype_growth_clones_v2.jl")

    """
        remove_mutations!(htumor, keep)
    Removes all mutations that are not in the collection `keep` and deletes their haplotypes from the haplotype tumor. Remaining identical haplotypes are recombined. Mutation indices will generally not coincide with row indices anymore!
    """
    function remove_mutations!(htumor, keep)
        removed = setdiff(1:nrow(htumor)-1, keep)
        for (i,muts) in enumerate(htumor.mutations)
            filter!(in(keep), muts)
            n = htumor.n[i]
            htumor.n[i] = 0
            if isempty(muts)
                htumor.n[1] += n
            else
                htumor.n[last(muts)+1] += n
            end
        end
        delete!(htumor, removed.+1)
        return htumor
    end
    remove_mutations(htumor, keep) = remove_mutations!(deepcopy(htumor), keep)

    """
        reduced_μ!(htumor, x)
    Artificially reduces mutation rate of a **complete (!) haplotype-tumor** by sampling from its set of mutations, where each mutation is selected with probability x. <br>
    Remaining identical haplotypes are recombined and mutation indices mapped to `1:nrow(htumor)-1`.
    """
    function reduced_μ!(htumor, x)
        keep = findall(rand(nrow(htumor)-1).<=x)
        remove_mutations!(htumor, keep)
        rebrand = (keep .=> 1:nrow(htumor)-1)
        for muts in htumor.mutations
            replace!(muts, rebrand...)
        end
        return htumor
    end
    reduced_μ(htumor, x) = reduced_μ!(deepcopy(htumor), x)

    ##########################
    ### Theoretical result ###
    ##########################

    """
        W_orphaned(q; N)
    Analytical solution to (red) **orphaned turnover**. Under this definition the probability of being orphaned
    depends on mutation extinction probability `q = d/b` and (threshold-) size `N` only.
    """
    W_orphaned(q; N) = q / (2*log(N)) * (1-N^(-2))

    """
        W_estranged(q; N)
    Analytical solution to **estranged turnover**. A haplotype clone's probability of becoming estranged, i.e. losing it's haplotype parent,
    depends on birth rate `b`, death rate `d`, mutation rate `μ` and (threshold-) time `T`. <br>
    Assumes **Branchingprocess `v2`**: independent mutation of offspring pair at division. `μ` is the probability of mutation per cell division.
    The mutation probability per cell division per cell is `μ/2`.
    """
    W_estranged(d; b, μ, T) = (μ/2)*(d + b*(μ/2)^2) / (1-(μ/2))^2 / (b*(1-μ/2)-d) * (1-exp(-2*(b-d-b*μ/2)*T)) / (1-exp(-2*b*μ/2*T))

    #########################
    ### METHODS WITH TREE ###
    #########################

    """
        mfreqs(htumor)

    Get mutation frequencies for a **complete (!) haplotype-tumor**,
    `complete` meaning that haplotypes with size `n` = 0 are listed.
    """
    function mfreqs(htumor)
        count = zeros(Int, nrow(htumor)-1)
        for (i, muts) in enumerate(htumor.mutations[2:end])
            count[muts] .+= htumor.n[i+1]
        end
        count ./ sum(htumor.n)
    end

    """
        estranged(htumor)

    Estranged turnover on **complete (!) haplotype-tumors**,
    *complete* meaning that haplotypes with size *n* = 0 are listed.
    The function reports for all offspring mutations wether their parental
    haplotype went extinct in form of

        DataFrame( mutation :: Vector{Int}, isestranged :: Vector{Bool})
    """
    function estranged(htumor)
        # freqs = mfreqs(htumor)
        # daughtertypes = htumor.mutations[ [false, (freqs .> 0)...] .& (length.(htumor.mutations) .> 1)  ]

        daughtertypes = htumor.mutations[ (htumor.n .> 0) .& (length.(htumor.mutations) .> 1)  ]

        parents_n = [htumor.n[muts[end-1]+1] for muts=daughtertypes]
        return DataFrame(  mutation = last.(daughtertypes), isestranged = iszero.(parents_n))
    end

    ### red-clones

    """
        orphaned_red(htumor)

    Orphaned turnover based on parent clones on **complete (!) haplotype-tumors**,
    `complete` meaning that haplotypes with size `n` = 0 are listed.
    The function returns the orphan-status of all offspring mutations of every clone as

        DataFrame( mutation :: Vector{Int}, isorphaned :: Vector{Bool})

    This result reduces to the method based on offspring clones by counting unique
    mutations in `mutation` to get total offspring count and by counting a mutation as
    orphaned if it is orphaned with respect to any of the parental clones

        total count: orphaned_red(htumor).mutation |> unique |> length
        orphan count: unique(orphaned_red(htumor)).orphaned |> count

    """
    function orphaned_red(htumor)
        mtypes = clones_by_mutations(htumor)[1]
        freqs = sum.(getproperty.(mtypes,:n))

        orphans = DataFrame(mutation = Int[], isorphaned = Bool[])
        for (m,redclone) in enumerate(mtypes)
            n_clone = freqs[m]
            m_green = last.(redclone.mutations[Not(1)])
            filter!( m-> !iszero(freqs[m]), m_green )
            n_green = freqs[m_green]
            isempty(n_green) && continue

            append!(orphans.mutation, m_green)
            append!(orphans.isorphaned, iszero.(n_clone .- n_green))
        end
        return orphans
    end

    # To recover "green-clones" method:
    #     - total count: orphaned_red(htumor).mutation |> unique |> length
    #     - orphan count: unique(orphaned_red(htumor)).orphaned |> count


    ### green-clones
    """
        orphaned_green(htumor)

    Orphaned turnover based on offspring clones on **complete (!) haplotype-tumors**,
    `complete` meaning that haplotypes with size `n` = 0 are listed.
    The function returns the orphan-status of all offspring mutations as

        DataFrame( mutation :: Vector{Int}, isorphaned :: Vector{Bool})
    """
    function orphaned_green(tumor)
        mtypes, mutations = clones_by_mutations(tumor)
        n = sum.(getproperty.(mtypes,:n))

        typeone = tumor.mutations[Not(1)] .|> first |> unique
        m_green = filter(m -> !iszero(n[m]), setdiff(mutations, typeone))
        m_red = [ first(mtypes[m].mutations)[end-1] for m=m_green]

        countorph = count(iszero, n[m_red] .- n[m_green])

        orphans = DataFrame(mutation = m_green, isorphaned = iszero.(n[m_red] .- n[m_green]) )
        return orphans
    end


    #########################
    ### TREE-LESS METHODS ###
    #########################

    """
        estranged_treeless(htumor)

    Tree-less version of the `estranged` algorithm on **incomplete (!) haplotype-tumors**. In an incomplete
    haplotype-tumor no haplotype has size `n` = 0 (those are to be removed beforehand).
    The function returns

        DataFrame( mutation :: Vector{Int}, isestranged :: Vector{Float64}, isgreen :: Vector{Float64} )

    Summing over `isestranged` gives the number of estranged mutations and summing over `isgreen` the total number
    of offspring haplotypes.
    """
    function estranged_treeless(htumor)
        # htumor = filter(h -> !iszero(h.n), tumor)

        mtypes, mutations = clones_by_mutations(htumor)

        truncs = [ intersect(Set.(mtype.mutations)...) for mtype in mtypes]
        htypes = Set.(htumor.mutations)

        mutdict = Dict(mutations .=> 1:length(mutations))
        estrangedmuts = DataFrame(  mutation = mutations,
                                    isestranged = zeros(length(mutations)),
                                    isgreen = zeros(length(mutations))
                                    )
        for (i,m) in enumerate(mutations)
            trunc = truncs[i]
            truncp = filter( mp -> truncs[mutdict[mp]] == trunc, trunc)
            t = length(truncp)

            lastparentestranged = !(setdiff(trunc, truncp) in htypes)
            lastmutalive = trunc in htypes

            if t > 1
                estrangedmuts.isestranged[i] = lastmutalive/t
                estrangedmuts.isgreen[i] = lastmutalive/t
            elseif isone(length(trunc))
                estrangedmuts.isestranged[i] = 0.
                estrangedmuts.isgreen[i] = 0.
            else
                estrangedmuts.isestranged[i] = lastparentestranged && lastmutalive
                estrangedmuts.isgreen[i] = lastmutalive
            end

        end
        return estrangedmuts
    end



    ### orphaned tree-less

    """
        orphaned_green_treeless(htumor)

    Tree-less version of the `orphaned_green` algorithm on **incomplete (!) haplotype-tumors**. In an incomplete
    haplotype-tumor no haplotype has size `n` = 0 (those are to be removed beforehand).
    The function returns

        DataFrame( mutation :: Vector{Int}, isorphaned :: Vector{Float64}, isgreen :: Vector{Float64} )

    Summing over `isorphaned` gives the number of orphans and summing over `isgreen` the total number
    of offspring mutations.
    """
    function orphaned_green_treeless(htumor)
        # htumor = filter(h -> !iszero(h.n), tumor)

        mtypes, mutations = clones_by_mutations(htumor)

        truncs = [ intersect(Set.(mtype.mutations)...) for mtype in mtypes]

        mutdict = Dict(mutations .=> 1:length(mutations))
        orphanedmuts = DataFrame(  mutation = mutations,
                                    isorphaned = zeros(length(mutations)),
                                    isgreen = zeros(length(mutations))
                                    )
        for (i,m) in enumerate(mutations)
            trunc = truncs[i]
            truncp = filter( mp -> truncs[mutdict[mp]] == trunc, trunc)
            t = length(truncp)
            orphanedmuts.isorphaned[i] = 1 - 1/t
            orphanedmuts.isgreen[i] = truncp == trunc ? 1 - 1/t : 1
        end
        return orphanedmuts
    end

    """
        orphaned_red_treeless(htumor)

    Tree-less version of the `orphaned_red` algorithm on **incomplete (!) haplotype-tumors**. In an incomplete
    haplotype-tumor no haplotype has size `n` = 0 (those are to be removed beforehand).
    The function returns

        DataFrame( mutation :: Vector{Int}, isorphaned :: Vector{Float64}, isgreen :: Vector{Float64} )

    Summing over `isorphaned` gives the number of orphans and summing over `isgreen` the total number
    of offspring mutations.
    """
    function orphaned_red_treeless(htumor)
        # htumor = filter(h -> !iszero(h.n), tumor)

        mtypes, mutations = clones_by_mutations(htumor)

        truncs = [ intersect(Set.(mtype.mutations)...) for mtype in mtypes]

        mutdict = Dict(mutations .=> 1:length(mutations))
        orphanedmuts = DataFrame(  mutation = mutations,
                                    isorphaned = zeros(length(mutations)),
                                    isgreen = zeros(length(mutations))
                                    )
        for (i,m) in enumerate(mutations)
            trunc = truncs[i]
            twigs = setdiff( union(mtypes[i].mutations...), trunc )
            truncp = filter( mp -> truncs[mutdict[mp]] == trunc, trunc)
            t = length(truncp)
            orphanedmuts.isorphaned[i] = (t-1)/2
            orphanedmuts.isgreen[i] = length(twigs) + (t-1)/2
        end
        return orphanedmuts
    end

    """
        orphaned_red_expected(q, clonesizes, nbirth)

    Expected orphaned turnover considering full ancestry of a mutant clone. Given extinction probability `q` and observed ancestral clone sizes at birth `nbirth` an offspring mutation `m` is expected to be orphaned with respect to an ancestor of size `n` at birth of `m` with probability `q^n`.
    The function returns the probability of being orphaned for all offspring mutations as

        DataFrame( mutation :: Vector{Int}, isorphaned :: Vector{Float64})
    """
    function orphaned_red_expected(q, nbirth)
        extant = filter(m -> !isempty(nbirth[m]), 1:length(nbirth))
        orphans = DataFrame(mutation = Int[], isorphaned = Float64[])
        for m_green in extant
            parentsizes = nbirth[m_green]
            append!( orphans.mutation, fill(m_green, length(parentsizes)) )
            append!( orphans.isorphaned, q .^ parentsizes)
        end
        return orphans
    end

    """
        orphaned_green_expected(htumor, q, nbirth)

    Expected orphaned turnover based on offspring mutations on **complete (!) haplotype-tumors**,
    `complete` meaning that haplotypes with size `n` = 0 are listed. Given extinction probability `q`
    and observed parent clone sizes at birth `nbirth` a (green) offspring mutation is expected to be orphaned
    with probability `q^n`.
    The function returns the probability of being orphaned for all offspring mutations as

        DataFrame( mutation :: Vector{Int}, isorphaned :: Vector{Float64})
    """
    function orphaned_green_expected(htumor, q, nbirth)
        mtypes = clones_by_mutations(htumor)[1]
        freqs = sum.(getproperty.(mtypes,:n))

        orphans = DataFrame(mutation = Int[], isorphaned = Float64[])
        for (m,redclone) in enumerate(mtypes)
            m_green = last.(redclone.mutations[Not(1)])
            filter!( m-> !iszero(freqs[m]), m_green )
            isempty(m_green) && continue

            append!(orphans.mutation, m_green)
            append!(orphans.isorphaned, q .^ nbirth[m_green] )
        end
        return orphans
    end

    """
        estranged(htumor, q, nbirth)

    Estranged turnover on **complete (!) haplotype-tumors**,
    *complete* meaning that haplotypes with size *n* = 0 are listed. Given extinction probability `q`
    and observed parent clone sizes at birth `nbirth`, an offspring mutation is expected to be estranged
    with probability `q^n`.
    The function returns the probability of being orphaned for all offspring mutations as

        DataFrame( mutation :: Vector{Int}, isestranged :: Vector{Float64})
    """
    function estranged_expected(htumor, q, nbirth)
        freqs = mfreqs(htumor)
        daughters = filter(h -> h.n >= 1 && length(h.mutations)>1, htumor)

        mutations = last.(daughters.mutations)
        return DataFrame(  mutation = mutations, isestranged = q .^ nbirth[mutations])
    end
end  # module Turnover

# function estranged_treeless_v2(htumor)
#     clones, mutations = clones_by_mutations(htumor)
#     m_indices = Dict(mutations .=> 1:length(mutations))
#     htypes = filter(h-> length(h)>1, htumor.mutations)
#     out = DataFrame(  clone = htypes, isestranged = Vector{Bool}(undef, length(htypes)) )
#
#     for (i,h) in enumerate(htypes)
#         n = length(h)
#         isestranged = true
#         for m in h
#             isestranged || break
#             for potential_parent in clones[m_indices[m]].mutations
#                 if length(potential_parent) == n-1 && isone(length(setdiff(h, potential_parent)))
#                     isestranged = false
#                     break
#                 end
#             end
#         end
#         out.isestranged[i] = isestranged
#     end
#     return out
# end
