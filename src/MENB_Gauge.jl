# for now, MEVB code deals only with all k-mers with k up to 3. 
# a more general scenario, where an arbitrary number of parameters associated to k-mers is fixed to 0, can be in principle done.

"""
    gauge_mask_variables(motifs::Vector{String})

This function returns a Array{Bool} corresponding to motifs that have to be inferred after
all motifs containing Ts are fixed to 0 thanks to gauge transformations.
"""
function gauge_mask_variables(motifs::Vector{String})
    return [(!occursin("T", m) | (length(m)==3 && m[1] != 'T' && m[2] == 'T' && m[3] != 'T')) for m in motifs]
end


"""
    gauge_zerosum(model_pars::NucleotideModel)
    
This function takes as input a dictionary describing the parmeters of a model
and changes the gauge into the "zero-sum" gauge.
"""
function gauge_zerosum(model::NucleotideModel)
    model_pars = get_forces_dict(model)
    new_pars = get_forces_dict(model)
    kmax = model.Lmotifs
    ##############
    #### k=3
    ##############
    if kmax == 3
        mot3s = [x for x in keys(model_pars) if length(x)==3]
        C1s = [- mean([new_pars[x] for x in mot3s if (x[2:3] == m3[2:3])]) for m3 in mot3s]
        [new_pars[mot3s[i]] += c1 for (i, c1) in enumerate(C1s)]
        [new_pars[mot3s[i][2:3]] -= c1/4 for (i, c1) in enumerate(C1s)]
        C3s = [- mean([new_pars[x] for x in mot3s if (x[1:2] == m3[1:2])]) for m3 in mot3s]
        [new_pars[mot3s[i]] += c3 for (i, c3) in enumerate(C3s)]
        [new_pars[mot3s[i][1:2]] -= c3/4 for (i, c3) in enumerate(C3s)]
    end
    ##############
    #### k=2
    ##############
    if kmax >= 2
        mot2s = [x for x in keys(model_pars) if length(x)==2]    
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
    return NucleotideModel(model.motifs, [new_pars[m] for m in model.motifs])
end