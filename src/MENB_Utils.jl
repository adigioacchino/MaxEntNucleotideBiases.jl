"""
    compute_minus_energy(seq::String, model_pars::Dict{String, Float64})
Given a sequence `seq` and the model parameters `model_pars`, 
compute (minus) the energy of this sequence. It uses all available threads
to increase speedup.
"""
function compute_minus_energy(seq::String, model_pars::Dict{String, Float64})
    Lmotifs = maximum(length.(keys(model_pars)))
    minus_energy = zeros(Lmotifs * length(seq))
    Threads.@threads for i in 1:length(seq)-Lmotifs+1
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
    compute_minus_energy_Lmodel3(seq::String, model_pars::Dict{String, Float64})
Given a sequence `seq` and the model parameters `model_pars`, 
compute (minus) the energy of this sequence. It only works for 
motifs of length 3, and it is slightly faster than `compute_minus_energy`.
It uses all available threads to increase speedup.
"""
function compute_minus_energy_Lmodel3(seq::String, model_pars::Dict{String, Float64})
    minus_energy = zeros(3 * (length(seq) - 2) + 3)
    Threads.@threads for i in 1:length(seq)-2
        @inbounds minus_energy[3*(i-1)+1] = model_pars[seq[i:i+2]]
        @inbounds minus_energy[3*(i-1)+2] = model_pars[seq[i:i+1]]
        @inbounds minus_energy[3*(i-1)+3] = model_pars[seq[i:i]]
    end
    @inbounds minus_energy[end-2] = model_pars[seq[end-1:end]]
    @inbounds minus_energy[end-1] = model_pars[seq[end-1:end-1]]
    @inbounds minus_energy[end] = model_pars[seq[end:end]]
    return sum(minus_energy)
end


"""
    compute_loglikelihood(seq::String, model_pars::Dict{String, Float64}; logZ=missing)
Given a sequence `seq` and the model parameters `model_pars`,
compute the log-likelihood (energy minus log of Z) of this sequence.
logZ can be passed directly if pre-computed, otherwise is it computed each time this function is called.
"""
function compute_loglikelihood(seq::String, model_pars::Dict{String, Float64}; logZ=missing)
    e = compute_minus_energy(seq, model_pars)
    ks = keys(model_pars)
    mlk = maximum(length.(ks))
    if ismissing(logZ)
        if mlk == 1
            logz = 0
        else
            L = length(seq)
            logz = eval_log_Z(model_pars, L)
        end
    else
            logz = logZ
    end
    return e - logz
end


"""
    _MC_step!(L::Int, mlk::Int, curr_seq::Vector{String}, 
              model_pars::Dict{String, Float64}, 
              beta::Float64)
Performs a MonteCarlo step; L is curr_seq length, mlk is the maximum length of the
motifs considered, curr_seq is the starting sequence, model_pars are the model 
parameters (dict of motifs => floats), and beta is the inverse temperature.
"""
function _MC_step!(L::Int, mlk::Int, curr_seq::Vector{String}, 
                  model_pars::Dict{String, Float64}, 
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
    mEnew = compute_minus_energy(new_kmer, model_pars)
    mEold = compute_minus_energy(old_kmer, model_pars)
    minus_deltaE = mEnew - mEold
    if rand() < exp(beta * minus_deltaE)
        curr_seq[pos] = proposed_nt
    end
end


"""
    Metropolis_sampling(L::Int, model_pars::Dict{String, Float64}; 
                        beta::Float64=1.0, Nsamples::Int=1, Nsteps::Int=1, Ntherm::Int=L*10, 
                        startseq=missing)
Sample Nsample sequences of length L at inverse temperature beta using a model defined by the
model parameters given in model_pars (dict of motifs => floats). Ntherm is the number of 
MonteCarlo steps done before starting collecting sequences, startseq is the starting sequence.
"""
function Metropolis_sampling(L::Int, model_pars::Dict{String, Float64}; 
                             beta::Float64=1.0, Nsamples::Int=1, Nsteps::Int=1, Ntherm::Int=L*10, 
                             startseq::Union{String,Missing}=missing)
    if ismissing(startseq)
        startseq = rand(dna_alphabet, L)
    else
        startseq = string.(collect(startseq))
    end
    ks = keys(model_pars)
    mlk = maximum(length.(ks))    
    # thermalization
    curr_seq = copy(startseq)
    for _ in 1:Ntherm
        _MC_step!(L, mlk, curr_seq, model_pars, beta)
    end
    # collect samples
    samples = Vector{String}()
    Lchain = Nsamples * Nsteps
    for i in 1:Lchain
        if i % Nsteps == 0
            push!(samples, join(curr_seq))
            _MC_step!(L, mlk, curr_seq, model_pars, beta)
        else
            _MC_step!(L, mlk, curr_seq, model_pars, beta)
        end
    end
    return samples
end


"""
    compute_entropy(L::Int, model_pars::Dict{String, Float64}; fast::Bool=false)
Compute the entropy of the model with parameters `model_pars` and having length `L`.
The computation exploit the derivatives of the partition function, that can be
approximated to make it faster (using `fast=true`), although it should be very quick
in any case.
"""
function compute_entropy(L::Int, model_pars::Dict{String, Float64}; fast::Bool=false)
    function lambda_eval_log_Z(lambda::Float64)
        mp = copy(model_pars)
        [mp[k] *= lambda for k in keys(mp)]
        if fast
            return eval_log_Zfast(mp, L)
        else
            return eval_log_Z(mp, L)    
        end
    end
    ave_E = -FiniteDiff.finite_difference_derivative(lambda_eval_log_Z, 1.)
    if fast
        return ave_E + eval_log_Zfast(model_pars, L)
    else
        return ave_E + eval_log_Z(model_pars, L)    
    end
end


"""
    compute_pressure(L::Int, model_pars::Dict{String, Float64}; fast::Bool=false)
The pressure on a genome can be quantified as how the genome is different from a random
uniform one. Therefore, this function compute the difference in the entropy of a random
uniform model and the model described with `model_pars`, with lenght `L`, rescaled with
the entropy of the uniform model. Setting `fast=true` allows for a quicker estimation of
the entropy of the model (it can be useful for extremely large sequences).
"""
function compute_pressure(L::Int, model_pars::Dict{String, Float64}; fast::Bool=false)
    return 1 - compute_entropy(L, model_pars) / (L * log(length(dna_alphabet)))
end


"""
    compute_symmetrized_KL(L::Int, model1_pars::Dict{String, Float64}, 
                           model2_pars::Dict{String, Float64}; fast::Bool=false)
Compute the symmetrized version of the Kullback-Leibler divergence between the distributions
defined by `model1_pars` and `model2_pars`. Setting `fast=true` allows for a quicker estimation 
of the KL divergence (it can be useful for extremely large sequences).
"""
function compute_symmetrized_KL(L::Int, model1_pars::Dict{String, Float64}, 
                        model2_pars::Dict{String, Float64}; fast::Bool=false)
    function lambda_eval_log_Z_12(lambda::Float64)
        mp = copy(model1_pars)
        [mp[k] += lambda*model2_pars[k] for k in keys(mp)]
        if fast
            return MaxEntViralBiases.eval_log_Zfast(mp, L)
        else
            return MaxEntViralBiases.eval_log_Z(mp, L)    
        end
    end
    function lambda_eval_log_Z_21(lambda::Float64)
        mp = copy(model2_pars)
        [mp[k] += lambda*model1_pars[k] for k in keys(mp)]
        if fast
            return MaxEntViralBiases.eval_log_Zfast(mp, L)
        else
            return MaxEntViralBiases.eval_log_Z(mp, L)    
        end
    end        
    ave_E2_1 = -FiniteDiff.finite_difference_derivative(lambda_eval_log_Z_12, 0.)
    ave_E1_2 = -FiniteDiff.finite_difference_derivative(lambda_eval_log_Z_21, 0.)    
    if fast
        logZ1 = MaxEntViralBiases.eval_log_Zfast(model1_pars, L)
        logZ2 = MaxEntViralBiases.eval_log_Zfast(model2_pars, L)
        S1 = compute_entropy(L, model1_pars, fast=true)
        S2 = compute_entropy(L, model2_pars, fast=true)        
    else
        logZ1 = MaxEntViralBiases.eval_log_Z(model1_pars, L)
        logZ2 = MaxEntViralBiases.eval_log_Z(model2_pars, L)
        S1 = compute_entropy(L, model1_pars, fast=false)
        S2 = compute_entropy(L, model2_pars, fast=false)   
    end
    return (ave_E2_1 + ave_E1_2 + logZ1 + logZ2 - S1 - S2) / 2
end