"""
    compute_minus_energy(seq::String, model::NucleotideModel)

Given a sequence `seq` and the model parameters `model_pars`, 
compute (minus) the energy of this sequence.
"""
function compute_minus_energy(seq::AbstractString, model::NucleotideModel)
    Lmotifs = model.Lmotifs
    model_pars = get_forces_dict(model)
    minus_energy = zeros(Lmotifs * length(seq))
    for i in 1:length(seq)-Lmotifs+1
        @inbounds minus_energy[Lmotifs*(i-1)+1] = model_pars[seq[i:i+Lmotifs-1]]
        for j in 1:Lmotifs-1
            @inbounds minus_energy[Lmotifs*(i-1)+1+j] = model_pars[seq[i:i+Lmotifs-1][1:j]]
        end
    end
    if Lmotifs == 3
        @inbounds minus_energy[end-2] = model_pars[seq[end-1:end]]
        @inbounds minus_energy[end-1] = model_pars[seq[end-1:end-1]]
        @inbounds minus_energy[end] = model_pars[seq[end:end]]
    elseif Lmotifs == 2
        @inbounds minus_energy[end] = model_pars[seq[end:end]]
    end
    return sum(minus_energy)
end


"""
    compute_loglikelihood(seq::String, model::NucleotideModel; logZ=missing)

Given a sequence `seq` and the model parameters `model_pars`,
compute the log-likelihood (energy minus log of Z) of this sequence.
logZ can be passed directly if pre-computed, otherwise is it computed each time this function is called.
"""
function compute_loglikelihood(seq::String, model::NucleotideModel; logZ=missing)
    e = compute_minus_energy(seq, model)
    if ismissing(logZ)
        if model.Lmotifs == 1
            logz = 0
        else
            L = length(seq)
            logz = compute_logz(model, L)
        end
    else
            logz = logZ
    end
    return e - logz
end


"""
    _mcstep!(L::Int, mlk::Int, curr_seq::Vector{String}, 
                  model::NucleotideModel, 
                  beta::Float64)
              
Performs a MonteCarlo step; L is curr_seq length, mlk is the maximum length of the
motifs considered, curr_seq is the starting sequence, model specifies the model, 
and beta is the inverse temperature.
"""
function _mcstep!(L::Int, mlk::Int, curr_seq::Vector{String}, 
                  model::NucleotideModel, 
                  beta::Float64)
    pos = rand(1:L)
    curr_nt = curr_seq[pos]
    curr_nt2num = dna2num_alphabet[curr_nt]
    proposed_nt = dna_alphabet[((curr_nt2num-1) + rand(1:3)) % 4 + 1]    
    b0_pos = max(1, pos-mlk+1)
    b1_pos = pos - 1
    a0_pos = pos + 1
    a1_pos = min(L, pos+mlk-1)
    old_kmer = join([curr_seq[b0_pos:b1_pos]; [curr_nt];      curr_seq[a0_pos:a1_pos]])
    new_kmer = join([curr_seq[b0_pos:b1_pos]; [proposed_nt];  curr_seq[a0_pos:a1_pos]])    
    mEnew = compute_minus_energy(new_kmer, model)
    mEold = compute_minus_energy(old_kmer, model)
    minus_deltaE = mEnew - mEold
    if rand() < exp(beta * minus_deltaE)
        curr_seq[pos] = proposed_nt
    end
end


"""
    sample_metropolis(L::Int, model::NucleotideModel; 
                             beta::Float64=1.0, Nsamples::Int=1, Nsteps::Int=1, Ntherm::Int=L*10, 
                             startseq::Union{String,Missing}=missing)

Sample Nsample sequences of length L at inverse temperature beta using a model given in model. 
Ntherm is the number of MonteCarlo steps done before starting collecting sequences, 
startseq is the starting sequence.
"""
function sample_metropolis(L::Int, model::NucleotideModel; 
                             beta::Float64=1.0, Nsamples::Int=1, Nsteps::Int=1, Ntherm::Int=L*10, 
                             startseq::Union{String,Missing}=missing)
    if ismissing(startseq)
        startseq = rand(dna_alphabet, L)
    else
        startseq = string.(collect(startseq))
    end
    mlk = model.Lmotifs   
    # thermalization
    curr_seq = copy(startseq)
    for _ in 1:Ntherm
        _mcstep!(L, mlk, curr_seq, model, beta)
    end
    # collect samples
    samples = Vector{String}()
    Lchain = Nsamples * Nsteps
    for i in 1:Lchain
        if i % Nsteps == 0
            push!(samples, join(curr_seq))
            _mcstep!(L, mlk, curr_seq, model, beta)
        else
            _mcstep!(L, mlk, curr_seq, model, beta)
        end
    end
    return samples
end


"""
    compute_entropy(L::Int, model::NucleotideModel; fast::Bool=false)

Compute the entropy of the model with parameters `model_pars` and having length `L`.
The computation exploit the derivatives of the partition function, that can be
approximated to make it faster (using `fast=true`), although it should be very quick
in any case.
"""
function compute_entropy(L::Int, model::NucleotideModel; fast::Bool=false)
    function lambda_EvalLogZ(lambda::Float64)
        model_mod = NucleotideModel(model.motifs, model.forces .* lambda)
        if fast
            return compute_logz_fast(model_mod, L)
        else
            return compute_logz(model_mod, L)    
        end
    end
    ave_E = -FiniteDiff.finite_difference_derivative(lambda_EvalLogZ, 1.)
    if fast
        return ave_E + compute_logz_fast(model, L)
    else
        return ave_E + compute_logz(model, L)    
    end
end


"""
    compute_pressure(L::Int, model::NucleotideModel; fast::Bool=false)
    
The pressure on a genome can be quantified as how the genome is different from a random
uniform one. Therefore, this function compute the difference in the entropy of a random
uniform model and the model described with `model_pars`, with lenght `L`, rescaled with
the entropy of the uniform model. Setting `fast=true` allows for a quicker estimation of
the entropy of the model (it can be useful for extremely large sequences).
"""
function compute_pressure(L::Int, model::NucleotideModel; fast::Bool=false)
    return 1 - compute_entropy(L, model, fast=fast) / (L * log(length(dna_alphabet)))
end


"""
    compute_symmetrized_kl(L::Int, model1::NucleotideModel, 
                        model2::NucleotideModel; fast::Bool=false)

Compute the symmetrized version of the Kullback-Leibler divergence between the distributions
defined by `model1_pars` and `model2_pars`. Setting `fast=true` allows for a quicker estimation 
of the KL divergence (it can be useful for extremely large sequences).
"""
function compute_symmetrized_kl(L::Int, model1::NucleotideModel, 
                        model2::NucleotideModel; fast::Bool=false)
    function LambdaEvalLogZ12(lambda::Float64)
        model1_mod = NucleotideModel(model1.motifs, model1.forces .+ model2.forces .* lambda)
        if fast
            return compute_logz_fast(model1_mod, L)
        else
            return compute_logz(model1_mod, L)    
        end
    end
    function LambdaEvalLogZ21(lambda::Float64)
        model2_mod = NucleotideModel(model2.motifs, model2.forces .+ model1.forces .* lambda)
        if fast
            return compute_logz_fast(model2_mod, L)
        else
            return compute_logz(model2_mod, L)    
        end
    end        
    ave_E2_1 = -FiniteDiff.finite_difference_derivative(LambdaEvalLogZ12, 0.)
    ave_E1_2 = -FiniteDiff.finite_difference_derivative(LambdaEvalLogZ21, 0.)    
    if fast
        logZ1 = compute_logz_fast(model1, L)
        logZ2 = compute_logz_fast(model2, L)
        S1 = compute_entropy(L, model1, fast=true)
        S2 = compute_entropy(L, model2, fast=true)        
    else
        logZ1 = compute_logz(model1, L)
        logZ2 = compute_logz(model2, L)
        S1 = compute_entropy(L, model1, fast=false)
        S2 = compute_entropy(L, model2, fast=false)   
    end
    return (ave_E2_1 + ave_E1_2 + logZ1 + logZ2 - S1 - S2) / 2
end