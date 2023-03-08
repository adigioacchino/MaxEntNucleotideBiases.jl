# for now, MEVB code deals only with all k-mers with k up to 3. 
# a more general scenario, where an arbitrary number of parameters associated to k-mers is fixed to 0, can be in principle done.

"""
    GaugeMaskVariables(motifs::Vector{String})
This function returns a Array{Bool} corresponding to motifs that have to be inferred after
all motifs containing Ts are fixed to 0 thanks to gauge transformations.
"""
function GaugeMaskVariables(motifs::Vector{String})
    return [!occursin("T", m) for m in motifs]
end


"""
    zerosum_gauge(model_pars::Dict{String, Float64})
This function takes as input a dictionary describing the parmeters of a model
and changes the gauge so that one and two point parameters are "zero-sum". It
does not modify three point parameters.
"""
function zerosum_gauge(model_pars::Dict{String, Float64})
    new_pars = copy(model_pars)
    kmax = maximum(length.(keys(new_pars)))
    ##############
    #### k=2
    ##############
    if kmax >= 2
        mot2s = [x for x in keys(model_pars) if length(x)==2]    
        # compute Cs & modify
        C1s = [- mean([new_pars[x] for x in mot2s if (x[1] == m2[1])]) for m2 in mot2s]
        [new_pars[mot2s[i]] += c1 for (i, c1) in enumerate(C1s)]
        [new_pars[string(mot2s[i][1])] -= c1/4 for (i, c1) in enumerate(C1s)]
        C2s = [- mean([new_pars[x] for x in mot2s if (x[2] == m2[2])]) for m2 in mot2s]
        [new_pars[mot2s[i]] += c2 for (i, c2) in enumerate(C2s)]
        [new_pars[string(mot2s[i][2])] -= c2/4 for (i, c2) in enumerate(C2s)]
    end
    ##############
    #### k=1
    ##############
    C = - mean([new_pars[nt] for nt in dna_alphabet])
    [new_pars[nt] += C for nt in dna_alphabet]
    return new_pars
end