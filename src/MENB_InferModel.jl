"""
    _preprocess_seqs(seqs_in::Vector{String})

Pre-process the vector of sequences `seqs_in`, so that sequences only contains letters A, C, G, T.
The pre-processed sequence vector is the returned.
"""
function _preprocess_seqs(seqs_in::Vector{String})
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
    _compute_freqs(seqs::Vector{String}, motifs::Vector{String})
    
For each motif in `motifs`, compute the frequency of observed motifs in each sequence
in `seqs`, then divide by the sequence length and take the average of these intensive 
fractions over the sequences.
"""
function _compute_freqs(seqs::Vector{String}, motifs::Vector{String})
    Lseqs = length.(seqs)
    return [mean(count.(m, seqs, overlap=true) ./ Lseqs) for m in motifs]    
end
    
"""
    fitmodel(seqs::Union{Vector{String}, String}, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; 
                  pseudocount_param::Float64=0.0, tolerance::Float64=0.01, max_iter::Int=100, 
                  verbose::Bool=false, fast::Union{Bool,String}="auto", ZS_gauge::Bool=true)

Fit the model parameters, which are:
- only 1-point functions (fields) if Lmotifs==1;
- 1-point and 2-point functions (fields, 2-mer forces) if Lmotifs==2;
- 1-point, 2-point and 3-point functions (fields, 2-mer and 3-mer forces) if Lmotifs==3.

`Lmodel` is the number of nucleotides used for the inference. If `seqs` is made of 
sequences of constant length (or if it is String and not a vector), the default option
will use the sequences length as `Lmodel`; otherwise, the default behaviour is to take
`Lmodel` = 5000.

`pseudocount_param` is the fraction of weigth coming from random uniform sequences
in the computation of the motif's frequences from the data.

If `fast`, the partition function is estimated through the top eigenvalue of the
transfer matrix alone (much faster, but slightly less precise, expecially
for short sequences). The default value, "auto", automatically uses the fast 
evaluation for long sequences.

`tolerance` and `max_iter` are parameters for the Newton-Raphson algorithm 
used to solve the system of equations.

`ZS_gauge` specifies whether the result has to be put in the zero sum gauge
before being returned.
"""
function fitmodel(seqs::Union{Vector{String}, String}, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; 
                  pseudocount_param::Float64=0.0, tolerance::Float64=0.01, max_iter::Int=100, 
                  verbose::Bool=false, fast::Union{Bool,String}="auto", ZS_gauge::Bool=true)
    (typeof(seqs) == String) ? (pre_seqs_dna = [seqs]) : (pre_seqs_dna = seqs)
    seqs_dna = _preprocess_seqs(pre_seqs_dna)
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
    
    mots_list = vcat([join.([collect(Iterators.product([dna_alphabet for _ in 1:k]...))...]) for k in 1:Lmotifs]...)
    freqs_obs = _compute_freqs(seqs_dna, mots_list)

    return fitmodel(freqs_obs, mots_list, L; pseudocount_param, tolerance, max_iter, 
                    verbose, fast, ZS_gauge)
end


"""
    fitmodel(freqs::Vector{<:Number}, motifs::Vector{String}, L::Int=5000; 
    pseudocount_param::Float64=0.0, tolerance::Float64=0.01, max_iter::Int=100, 
    verbose::Bool=false, fast::Union{Bool,String}="auto", ZS_gauge::Bool=true)

Fit the model parameters directly from the frequencies `freqs` of the motifs `motifs`.
The other parameters are the same as for the function that takes as input the sequences
directly, and are documented above.
"""
function fitmodel(freqs::Vector{<:Number}, motifs::Vector{String}, L::Int=5000; 
    pseudocount_param::Float64=0.0, tolerance::Float64=0.01, max_iter::Int=100, 
    verbose::Bool=false, fast::Union{Bool,String}="auto", ZS_gauge::Bool=true)
    
    Lmotifs = maximum(length.(motifs))

    # deal with pre_fast=auto
    if typeof(fast) == String && fast != "auto"
        println("`fast` must be a Bool, or the string 'auto'.")
        return nothing
    elseif fast == "auto"
        (L > 5000) ? (fast_bool = true) : (fast_bool = false)
    else
        fast_bool = fast
    end

    # compute all motifs and gaugemaks array
    if Lmotifs == 1 # only in this case, the gauge will be set so that the sum of the exp of the fields is 1
        return NucleotideModel(motifs, log.(freqs))
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
    gaugemask = gauge_mask_variables(all_motifs)
    independent_motifs = all_motifs[gaugemask]
    dependent_motifs = all_motifs[.!gaugemask]

    # compute n_obs for each non-masked motif, add pseudocounts if add_pseudocount
    freqs_dict = Dict(zip(motifs, freqs))
    n_obs = [freqs_dict[m] .* L for m in independent_motifs] 

    # add pseudocounts
    for i in 1:length(n_obs)
        len_motif = length(independent_motifs[i])
        n_obs[i] = (n_obs[i] + L * pseudocount_param / 4^len_motif) / (1 + pseudocount_param)
    end   
    
    # prepare for inference: "closure" of eval_logZ -> here the sorting of x is important, 
    #    it is the one defined at the definition of independent_motifs 
    function closed_compute_logz(x::Vector{Float64})
        mots = [independent_motifs; dependent_motifs]
        fors = [x; zeros(length(dependent_motifs))]
        if fast_bool
            return compute_logz_fast(NucleotideModel(mots, fors), L)
        else
            return compute_logz(NucleotideModel(mots, fors), L)    
        end
    end
    
    # inference: starting point
    if Lmotifs == 2 # use as starting point the solution with Lmotifs == 1
        n_obs_T = freqs_dict["T"] * L
        vars = zeros(sum(gaugemask))
        [vars[i] = log(n_obs[i] / L) - log(n_obs_T / L) for i in 1:3]
    elseif Lmotifs == 3 # use as starting point the solution with Lmotifs == 2
        motifs2pts = vcat([join.([collect(Iterators.product([dna_alphabet for _ in 1:k]...))...]) for k in 1:2]...)
        freqs2pts = [freqs_dict[m] for m in motifs2pts]
        L2ress = get_forces_dict(fitmodel(freqs2pts, motifs2pts, L; pseudocount_param=pseudocount_param,
                            tolerance=0.01, max_iter=100, verbose=false, fast=true))
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
        ns = FiniteDiff.finite_difference_gradient(closed_compute_logz, vars)
        if maximum(abs.(n_obs .- ns)) <= tolerance
            break
        end
        dn = FiniteDiff.finite_difference_hessian(closed_compute_logz, vars)
        delta = inv(dn) * (n_obs .- ns)
        vars .+= delta
    end

    # format result as a NucleotideModel, including the motifs set to 0 via gauge
    mots = [independent_motifs; dependent_motifs]
    fors = [vars; zeros(length(dependent_motifs))]
    if ZS_gauge
        res = gauge_zerosum(NucleotideModel(mots, fors))
    else
        res = NucleotideModel(mots, fors)
    end

    return res
end