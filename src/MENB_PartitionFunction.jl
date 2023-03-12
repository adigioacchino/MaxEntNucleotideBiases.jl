"""
    generate_transfer_matrix(model::NucleotideModel)

Return the trasfer matrix used for the computation of the
partition function. The first index correspond to the left-most nucleotide,
the last index to the right-most one.
"""
function generate_transfer_matrix(model::NucleotideModel)
    if model.Lmotifs == 2
        kmers = [a*b for a in dna_alphabet, b in dna_alphabet]
        p1t = sum(
                [count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==1],
                dims=1)[1] ./ 2
        p2t = sum(
                [count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==2],
                dims=1)[1]
        TM = exp.(p1t .+ p2t)
    elseif model.Lmotifs == 3
        k2 = [m for m in model.motifs if length(m)==2]
        pre_kmers = [a*b for a in k2, b in k2]
        zeromask = [x[2] == x[3] ? 1 : 0 for x in pre_kmers]
        kmers = [x[1]*x[2]*x[4] for x in pre_kmers]
        p1t = sum(
                [count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==1]
                 ) ./ 3
        p2t = sum(
                [count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==2]
                 ) ./ 2
        p3t = sum(
                [count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==3]
                 )        
        pre_TM = exp.(p1t .+ p2t .+ p3t)
        TM = pre_TM .* zeromask
    else
        println("Error, maximum length of motifs used to generate the TM should be 2 or 3!")
    end
    return TM
end


"""
    generate_transfer_matrix_last(model::NucleotideModel)

Return the last transfer matrix used for the computation of the
partition function. The first index correspond to the left-most nucleotide,
the last index to the right-most one.
"""
function generate_transfer_matrix_last(model::NucleotideModel)
    if model.Lmotifs == 2
        kmers = [a*b for a in dna_alphabet, b in dna_alphabet]
        k1mers = [b  for _ in dna_alphabet, b in dna_alphabet]
        p1t = sum(
                [(count.(m, kmers, overlap=true)
                    + count.(m, k1mers, overlap=true)).* f for (m, f) in zip(model.motifs, model.forces) if length(m)==1], 
                 dims=1)[1] ./ 2
        p2t = sum(
                [count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==2], 
                 dims=1)[1]        
        TM = exp.(p1t .+ p2t)
    elseif model.Lmotifs == 3
        k2 = [m for m in model.motifs if length(m)==2]
        pre_kmers = [a*b for a in k2, b in k2]
        zeromask = [x[2] == x[3] ? 1 : 0 for x in pre_kmers]
        kmers = [x[1]*x[2]*x[4] for x in pre_kmers]
        k2mers = [x[2]*x[3]     for x in pre_kmers]
        k12mers = [string(x[2]) for x in pre_kmers]
        k13mers = [string(x[3]) for x in pre_kmers]        
        p1t = sum([(count.(m, kmers, overlap=true) 
                    + count.(m, k12mers, overlap=true)
                    + 2 * count.(m, k13mers, overlap=true)) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==1], 
                  dims=1)[1] ./ 3
        p2t = sum([(count.(m, kmers, overlap=true)
                    + count.(m, k2mers, overlap=true)) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==2], 
                  dims=1)[1] ./ 3
        p3t = sum([count.(m, kmers, overlap=true) .* f for (m, f) in zip(model.motifs, model.forces) if length(m)==3], 
                  dims=1)[1]
        pre_TM = exp.(p1t .+ p2t .+ p3t)
        TM = pre_TM .* zeromask
    else
        println("Error, maximum length of motifs used to generate the TM should be 2 or 3!")
    end
    return TM
end


"""
    compute_logz(model::NucleotideModel, L::Int)

Compute the partition function of a model of length L through the transfer matrix
method.
"""
function compute_logz(model::NucleotideModel, L::Int)
    TM = generate_transfer_matrix(model)
    TM_last = generate_transfer_matrix_last(model)
    log_factors = 0.
    K = 10
    nK = norm(TM^K)
    lnK = log(nK)
    TMK = TM^K
    tP = copy(TM_last)
    for _ in 1:(L-2)÷K
        log_factors += lnK
        tP = TMK * tP * (1/nK)
    end
    for _ in 1:(L-2)%K
        tP = TM * tP
    end
    return log(tr(tP)) + log_factors
end


"""
    compute_logz_fast(model::NucleotideModel, L::Int)

Compute the partition function of a model of length L by taking the real part of
the largest eigenvalue of the transfer matrix.
"""
function compute_logz_fast(model::NucleotideModel, L::Int)
    TM = generate_transfer_matrix(model)
    return L * log(maximum(real.(eigvals(TM))))
end