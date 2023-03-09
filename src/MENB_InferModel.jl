# For now only allow to infer: fields-only, pairs, 3points (by giving an integer to the main function!)

"""
    _PreprocessSeqs(seqs_in::Vector{String})
Pre-process the vector of sequences `seqs_in`, so that sequences only contains letters A, C, G, T.
The pre-processed sequence vector is the returned.
"""
function _PreprocessSeqs(seqs_in::Vector{String})
    seqs_out = [replace(s, "U" => "T") for s in seqs_in]
    for s in seqs_out
        unique_nt_counts = sum([count(m, s) for m in dna_alphabet])
        if length(s) != unique_nt_counts
            println("Warning: one or more sequences have a non-ACGT/U nucleotide symbol.")
            break
        end
    end
    return seqs_out
end
    
    
"""
    _ComputeNObs(seqs::Vector{String}, independent_motifs::Vector{String}, L::Int)
For each motif in `independent_motifs`, compute the number of observed motifs in each sequence
in `seqs`, then divide by the sequence length, take the average of these intensive fractions
over the sequences, and multiply by the model length L.
"""
function _ComputeNObs(seqs::Vector{String}, independent_motifs::Vector{String}, L::Int)
    Lseqs = length.(seqs)
    return L .* [mean(count.(m, seqs, overlap=true) ./ Lseqs) for m in independent_motifs]    
end
    
    
"""
    ModelFit(seqs::Vector{String}, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; 
             add_pseudocount::Bool=false, tolerance::Float64=0.01, max_iter::Int=100, 
             verbose::Bool=true)
    ModelFit(seq::String, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; 
             add_pseudocount::Bool=false, tolerance::Float64=0.01, max_iter::Int=100, 
             verbose::Bool=true)

Fit the model parameters, which are:
- only 1-point functions (fields) if Lmotifs==1;
- 1-point and 2-point functions (fields, 2-mer forces) if Lmotifs==2;
- 1-point, 2-point and 3-point functions (fields, 2-mer and 3-mer forces) if Lmotifs==3.

`Lmodel` is the number of nucleotides used for the inference. If `seqs` is made of 
sequences of constant length (or if it is String and not a vector), the default option
will use the sequences length as `Lmodel`; otherwise, the default behaviour is to take
`Lmodel` = 5000.

If `add_pseudocount`, a single pseudocount is added for each observed number of
nucleotides and dinucleotides.

If `fast`, the partition function is estimated through the top eigenvalue of the
transfer matrix alone (much faster, but slightly less precise, expecially
for short sequences).

Finally, `tolerance` and `max_iter` are parameters for the Newton-Raphson algorithm 
used to solve the system of equations.
"""
function ModelFit(seqs::Vector{String}, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; 
                  add_pseudocount::Bool=false, tolerance::Float64=0.01, max_iter::Int=100, 
                  verbose::Bool=false, fast::Bool=false)
    seqs_dna = _PreprocessSeqs(seqs)
    # compute sequence length
    if ismissing(Lmodel)
        if all(length.(seqs_dna) .== length(seqs_dna[1]))
            L = length(seqs_dna[1])
        else
            L = 5000
            println("Warning: Lmodel is not specified and the sequence set provided contains sequences",
                    " of different lengths. Using an arbitrary value of L=5000.")
        end
    else
        L = Lmodel
    end
    
    # compute all motifs and gaugemaks array
    if Lmotifs == 1 # only in this case, the gauge will be set so that the sum of the exp of the fields is 1
        Lseqs = length.(seqs)
        n_freqs = [mean(count.(m, seqs_dna, overlap=true) ./ Lseqs) for m in dna_alphabet]
        res = Dict(zip(dna_alphabet, log.(n_freqs)))
        return res
    elseif Lmotifs == 2
        all_motifs = [[a for a in dna_alphabet]; 
                      [a*b for a in dna_alphabet for b in dna_alphabet]
                     ]
    elseif Lmotifs == 3
        all_motifs = [[a for a in dna_alphabet]; 
                      [a*b for a in dna_alphabet for b in dna_alphabet]; 
                      [a*b*c for a in dna_alphabet for b in dna_alphabet for c in dna_alphabet]
                     ]
    else
        println("max motif length is 3, exiting...")
        return
    end
    gaugemask = GaugeMaskVariables(all_motifs)
    independent_motifs = all_motifs[gaugemask]
    dependent_motifs = all_motifs[.!gaugemask]
    mp_dep = Dict(zip(dependent_motifs, zeros(length(dependent_motifs))))
    #println("Indep motifs:", independent_motifs)

    # compute n_obs for each non-masked motif, add pseudocounts if add_pseudocount
    n_obs = _ComputeNObs(seqs_dna, independent_motifs, L)
    if add_pseudocount
        n_obs = [x+1 for x in n_obs]
    end   
    
    # prepare for inference: "closure" of eval_logZ -> here the sorting of x is important, 
    #    it is the one defined at the definition of independent_motifs 
    function ClosedEvalLogZ(x::Vector{Float64})
        mp_ind = Dict(zip(independent_motifs, x))
        mp = merge(mp_ind, mp_dep)
        if fast
            return EvalLogZFast(mp, L)
        else
            return EvalLogZ(mp, L)    
        end
    end
    
    # inference: starting point
    if Lmotifs == 2 # use as starting point the solution with Lmotifs == 1
        Lseqs = length.(seqs)
        n_obs_T = L * mean(count.("T", seqs_dna, overlap=true) ./ Lseqs)
        vars = zeros(sum(gaugemask))
        [vars[i] = log(n_obs[i] / L) - log(n_obs_T / L) for i in 1:3]
    elseif Lmotifs == 3 # use as starting point the solution with Lmotifs == 2
        L2ress = ModelFit(seqs_dna, 2, L; add_pseudocount=true,
                          tolerance=0.01, max_iter=100, verbose=false, fast=true)
        independent_motifs2 = [m for m in independent_motifs if length(m) < 3]
        vars = [[L2ress[m] for m in independent_motifs2];
                zeros(length(independent_motifs)-length(independent_motifs2))]
    end

    # inferece: Newton-Raphson algorithm
    for l in 1:max_iter
        if verbose
            println("Starting iteration $(l)...")
            flush(stdout)
        end
        ns = FiniteDiff.finite_difference_gradient(ClosedEvalLogZ, vars)
        dn = FiniteDiff.finite_difference_hessian(ClosedEvalLogZ, vars)
        delta = inv(dn) * (n_obs .- ns)
        if maximum(abs.(n_obs .- ns)) <= tolerance
            break
        end
        vars .+= delta
    end
    
    # format result as a Dict{String, Float64}, including the motifs set to 0 via gauge
    res_1 = Dict(zip(independent_motifs, vars))
    return merge(res_1, mp_dep)
end

function ModelFit(seq::String, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; 
                  add_pseudocount::Bool=false, tolerance::Float64=0.01, max_iter::Int=100, 
                  verbose::Bool=true)
    return ModelFit([seq], Lmotifs, Lmodel; add_pseudocount, tolerance, max_iter, verbose)
end