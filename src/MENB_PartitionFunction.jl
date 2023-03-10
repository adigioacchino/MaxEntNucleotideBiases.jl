"""
    GenerateTM(model_pars::Dict{String, Float64})
Return the trasfer matrix used for the computation of the
partition function. The first index correspond to the left-most nucleotide,
the last index to the right-most one.
"""
function GenerateTM(model_pars::Dict{String, Float64})
    ks = keys(model_pars) 
    if maximum(length.(ks)) == 2
        k1 = [k for k in ks if length(k)==1]
        k2 = [k for k in ks if length(k)==2]
        kmers = [a*b for a in dna_alphabet, b in dna_alphabet]
        p1t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k1], dims=1)[1] ./ 2
        p2t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k2], dims=1)[1]        
        TM = exp.(p1t .+ p2t)
    elseif maximum(length.(ks)) == 3
        k1 = [k for k in ks if length(k)==1]
        k2 = [k for k in ks if length(k)==2]
        k3 = [k for k in ks if length(k)==3]
        pre_kmers = [a*b for a in k2, b in k2]
        zeromask = [x[2] == x[3] ? 1 : 0 for x in pre_kmers]
        kmers = [x[1]*x[2]*x[4] for x in pre_kmers]
        p1t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k1]) ./ 3
        p2t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k2]) ./ 2
        p3t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k3])        
        pre_TM = exp.(p1t .+ p2t .+ p3t)
        TM = pre_TM .* zeromask
    else
        println("Error, maximum length of motifs used to generate the TM should be 2 or 3!")
    end
    return TM
end


"""
    GenerateTMLast(model_pars::Dict{String, Float64})
Return the last transfer matrix used for the computation of the
partition function. The first index correspond to the left-most nucleotide,
the last index to the right-most one.
"""
function GenerateTMLast(model_pars::Dict{String, Float64})
    ks = keys(model_pars) 
    if maximum(length.(ks)) == 2
        kmers = [a*b for a in dna_alphabet, b in dna_alphabet]
        k1mers = [b  for a in dna_alphabet, b in dna_alphabet]
        k1 = [k for k in ks if length(k)==1]
        k2 = [k for k in ks if length(k)==2]
        p1t = sum([(count.(k, kmers, overlap=true)
                    + count.(k, k1mers, overlap=true)) .* model_pars[k] for k in k1], dims=1)[1] ./ 2
        p2t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k2], dims=1)[1]        
        TM = exp.(p1t .+ p2t)
    elseif maximum(length.(ks)) == 3 # I could also always create this kind of matrices (with 0 parameters for all 3-point motifs), the code will be simpler but slower
        k1 = [k for k in ks if length(k)==1]
        k2 = [k for k in ks if length(k)==2]
        k3 = [k for k in ks if length(k)==3]        
        pre_kmers = [a*b for a in k2, b in k2]
        zeromask = [x[2] == x[3] ? 1 : 0 for x in pre_kmers]
        kmers = [x[1]*x[2]*x[4] for x in pre_kmers]
        k2mers = [x[2]*x[3]     for x in pre_kmers]
        k12mers = [string(x[2]) for x in pre_kmers]
        k13mers = [string(x[3]) for x in pre_kmers]        
        p1t = sum([(count.(k, kmers, overlap=true) 
                    + count.(k, k12mers, overlap=true)
                    + 2 * count.(k, k13mers, overlap=true)) .* model_pars[k] for k in k1], dims=1)[1] ./ 3
        p2t = sum([(count.(k, kmers, overlap=true)
                    + count.(k, k2mers, overlap=true)) .* model_pars[k] for k in k2], dims=1)[1] ./ 3
        p3t = sum([count.(k, kmers, overlap=true) .* model_pars[k] for k in k3], dims=1)[1]
        pre_TM = exp.(p1t .+ p2t .+ p3t)
        TM = pre_TM .* zeromask
    else
        println("Error, maximum length of motifs used to generate the TM should be 2 or 3!")
    end
    return TM
end


"""
    EvalLogZ(model_pars::Dict{String, Float64}, L::Int)
Compute the partition function of a model of length L through the transfer matrix
method. model_pars is a dict of motif => parameter in the Hamiltonian.
"""
function EvalLogZ(model_pars::Dict{String, Float64}, L::Int)
    TM = GenerateTM(model_pars)
    TM_last = GenerateTMLast(model_pars)
    log_factors = 0
    tP = copy(TM_last)
    for i in 1:(L-3)
        if i%10 == 0 
            f = norm(tP)
            log_factors += log(f)
            tP = TM * tP * (1/f)
        else
            tP = TM * tP
        end
    end
    tP = TM * tP
    return log(tr(tP)) + log_factors
end


"""
    EvalLogZFast(model_pars::Dict{String, Float64}, L::Int)
Compute the partition function of a model of length L by taking the real part of
the largest eigenvalue of the transfer matrix. model_pars is a dict of 
motif => parameter in the Hamiltonian.
"""
function EvalLogZFast(model_pars::Dict{String, Float64}, L::Int)
    TM = GenerateTM(model_pars)
    return L * log(maximum(real.(eigvals(TM))))
end


"""
    EvalLogZFast(model_pars::Dict{String, Float64}, L::Int, fast::Bool=false)
    Compute the partition function of a model of length L.
    If fast==false, use the transfer matrix method. 
    If fast==true
        Compute the partition function of a model of length L by taking the real part of
the largest eigenvalue of the transfer matrix. model_pars is a dict of 
motif => parameter in the Hamiltonian.
"""
function EvalLogZ(model::NucleotideModel, L::Int, fast::Bool=false)
    model_pars = ForcesDict(model)
    if fast
        res = EvalLogZFast(model_pars, L)
    else
        res = EvalLogZ(model_pars, L)
    end
    return res
end