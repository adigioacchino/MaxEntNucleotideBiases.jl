var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MaxEntNucleotideBiases","category":"page"},{"location":"#MaxEntNucleotideBiases","page":"Home","title":"MaxEntNucleotideBiases","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MaxEntNucleotideBiases.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is not registered. Install with:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(url=\"https://github.com/adigioacchino/MaxEntNucleotideBiases.jl\")","category":"page"},{"location":"#Introduction-and-references","page":"Home","title":"Introduction and references","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TBD","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"#Basic-usage","page":"Home","title":"Basic usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We will take sub-sequence from the Influenza H5N1 PB2 segment (strain used: A/Anhui/1/2005) and we will use it as working example to test the function of the package. We start by loading the package and the sequence:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using MaxEntNucleotideBiases\nexample_seq = \"ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using MaxEntNucleotideBiases\nexample_seq = \"ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using MaxEntNucleotideBiases\nexample_seq = \"ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT\"\nmod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2)\nmod3 = MaxEntNucleotideBiases.fitmodel(example_seq, 3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package allows to fit up to 1-, 2- or 3-mers (i.e. single nucleotide, dinucleotide or trinucleotide) forces.  To do so, we will use the fitmodel function, whose second argument is the 'k' of the maximum k-mers we want to consider (1, 2 or 3).","category":"page"},{"location":"","page":"Home","title":"Home","text":"We will start by fitting the single nucleotide biases:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mod1 = MaxEntNucleotideBiases.fitmodel(example_seq, 1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The result is a NucleotideModel object, which contains the motifs, the fitted parameters and the value of 'k'. Similarly, we can fit dinucleotide and trinucleotide biases:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2)\nmod3 = MaxEntNucleotideBiases.fitmodel(example_seq, 3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"From a model we can obtain a dictionary (motif => force) with the get_forces_dict function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"MaxEntNucleotideBiases.get_forces_dict(mod2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we can compute several interesting quantities:","category":"page"},{"location":"","page":"Home","title":"Home","text":"the energy of any sequence with a given model, with the compute_minus_energy function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L = 1000;\nseq = join(rand([\"A\",\"C\",\"G\",\"T\"], L));\nMaxEntNucleotideBiases.compute_minus_energy(mod3, seq)","category":"page"},{"location":"","page":"Home","title":"Home","text":"the log-likelihood of any sequence with a given model, with the compute_loglikelihood function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L = 1000;\nseq = join(rand([\"A\",\"C\",\"G\",\"T\"], L));\nMaxEntNucleotideBiases.compute_loglikelihood(mod2, seq)","category":"page"},{"location":"#Gauge-choices","page":"Home","title":"Gauge choices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are several gauge choices that can be used to fit the model. The default choice is the 'zero-sum' gauge, that is described in the paper. The other notable choice is the 'lattice gas gauge', which is again described in the paper and is implemented in this package. In particular, the results of the inference can be directly returned in the lattice gas gauge with the ZS_gauge argument passed as false to the fitmodel function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2, ZS_gauge=false)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Moreover, the zero-sum gauge can be restored after the inference with the gauge_zerosum function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2, ZS_gauge=false) # hide\nMaxEntNucleotideBiases.gauge_zerosum(mod2)","category":"page"},{"location":"#Entropy-and-related-quantities","page":"Home","title":"Entropy and related quantities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The simple structure of the models fitted with this package allows to compute exactly several quantities of interest that, in most other cases, require (uncontrolled) approximations. The first example is the entropy of the probability distribution defined with the model.  It can be computed with the compute_entropy function, and since this quantity depends on the length of the sequences considered, a reference length must be passed as argument:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L = 1000;\nMaxEntNucleotideBiases.compute_entropy(L, mod3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The entropy of the model can be used to quantify how different the probability distribution defined by the model is from the uniform distribution, which in turn is a measure of the 'pressure' that is exerted on the nucleotides to be used in a certain way. In this package, the pressure ranges from 0 (unconstrained model) to 1 (fully constrained model), and it is computed with the compute_pressure function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L = 1000;\nMaxEntNucleotideBiases.compute_pressure(L, mod3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, it can be interesting to quantify how much two models are different from each other.  This package address this problem with the compute_symmetrized_kl function, which computes the symmetrized Kullback-Leibler divergence between two models:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L = 1000;\nmod3_unif = MaxEntNucleotideBiases.fitmodel(join(rand([\"A\",\"C\",\"G\",\"T\"], L)), 3)\nMaxEntNucleotideBiases.compute_symmetrized_kl(L, mod3, mod3_unif)","category":"page"},{"location":"#Sampling","page":"Home","title":"Sampling","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Another interseting feature of this package is the possibility to sample sequences from an inferred model. This package implements a simple Metropolis routine to sample sequences from a model:","category":"page"},{"location":"","page":"Home","title":"Home","text":"L = 1000;\nMaxEntNucleotideBiases.sample_metropolis(L, mod3, \n                                         Nsamples=10, Nsteps=3000, Ntherm=10_000, beta=1.)","category":"page"},{"location":"#Model-I/O","page":"Home","title":"Model I/O","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The model inference can take some time, so it can be useful to save the model on disk and load it later. This can be done with the save_model and load_model functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"MaxEntNucleotideBiases.writemodel(mod3, \"mod3.menb\")\nmod3 = MaxEntNucleotideBiases.readmodel(\"mod3.menb\")","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MaxEntNucleotideBiases]","category":"page"},{"location":"#MaxEntNucleotideBiases._compute_nobs-Tuple{Vector{String}, Vector{String}, Int64}","page":"Home","title":"MaxEntNucleotideBiases._compute_nobs","text":"_compute_nobs(seqs::Vector{String}, independent_motifs::Vector{String}, L::Int)\n\nFor each motif in independent_motifs, compute the number of observed motifs in each sequence in seqs, then divide by the sequence length, take the average of these intensive fractions over the sequences, and multiply by the model length L.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases._mcstep!-Tuple{Int64, Int64, Vector{String}, MaxEntNucleotideBiases.NucleotideModel, Float64}","page":"Home","title":"MaxEntNucleotideBiases._mcstep!","text":"_mcstep!(L::Int, mlk::Int, curr_seq::Vector{String}, \n              model::NucleotideModel, \n              beta::Float64)\n\nPerforms a MonteCarlo step; L is currseq length, mlk is the maximum length of the motifs considered, currseq is the starting sequence, model specifies the model,  and beta is the inverse temperature.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases._preprocess_seqs-Tuple{Vector{String}}","page":"Home","title":"MaxEntNucleotideBiases._preprocess_seqs","text":"_preprocess_seqs(seqs_in::Vector{String})\n\nPre-process the vector of sequences seqs_in, so that sequences only contains letters A, C, G, T. The pre-processed sequence vector is the returned.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_entropy-Tuple{Int64, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.compute_entropy","text":"compute_entropy(L::Int, model::NucleotideModel; fast::Bool=false)\n\nCompute the entropy of the model with parameters model_pars and having length L. The computation exploit the derivatives of the partition function, that can be approximated to make it faster (using fast=true), although it should be very quick in any case.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_loglikelihood-Tuple{String, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.compute_loglikelihood","text":"compute_loglikelihood(seq::String, model::NucleotideModel; logZ=missing)\n\nGiven a sequence seq and the model parameters model_pars, compute the log-likelihood (energy minus log of Z) of this sequence. logZ can be passed directly if pre-computed, otherwise is it computed each time this function is called.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_logz-Tuple{MaxEntNucleotideBiases.NucleotideModel, Int64}","page":"Home","title":"MaxEntNucleotideBiases.compute_logz","text":"compute_logz(model::NucleotideModel, L::Int)\n\nCompute the partition function of a model of length L through the transfer matrix method.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_logz_fast-Tuple{MaxEntNucleotideBiases.NucleotideModel, Int64}","page":"Home","title":"MaxEntNucleotideBiases.compute_logz_fast","text":"compute_logz_fast(model::NucleotideModel, L::Int)\n\nCompute the partition function of a model of length L by taking the real part of the largest eigenvalue of the transfer matrix.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_minus_energy-Tuple{AbstractString, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.compute_minus_energy","text":"compute_minus_energy(seq::String, model::NucleotideModel)\n\nGiven a sequence seq and the model parameters model_pars,  compute (minus) the energy of this sequence.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_pressure-Tuple{Int64, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.compute_pressure","text":"compute_pressure(L::Int, model::NucleotideModel; fast::Bool=false)\n\nThe pressure on a genome can be quantified as how the genome is different from a random uniform one. Therefore, this function compute the difference in the entropy of a random uniform model and the model described with model_pars, with lenght L, rescaled with the entropy of the uniform model. Setting fast=true allows for a quicker estimation of the entropy of the model (it can be useful for extremely large sequences).\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.compute_symmetrized_kl-Tuple{Int64, MaxEntNucleotideBiases.NucleotideModel, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.compute_symmetrized_kl","text":"compute_symmetrized_kl(L::Int, model1::NucleotideModel, \n                    model2::NucleotideModel; fast::Bool=false)\n\nCompute the symmetrized version of the Kullback-Leibler divergence between the distributions defined by model1_pars and model2_pars. Setting fast=true allows for a quicker estimation  of the KL divergence (it can be useful for extremely large sequences).\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.fitmodel","page":"Home","title":"MaxEntNucleotideBiases.fitmodel","text":"fitmodel(seqs::Union{Vector{String}, String}, Lmotifs::Int, Lmodel::Union{Int,Missing}=missing; \n              pseudocount_param::Float64=0.0, tolerance::Float64=0.01, max_iter::Int=100, \n              verbose::Bool=false, fast::Union{Bool,String}=\"auto\", ZS_gauge::Bool=true)\n\nFit the model parameters, which are:\n\nonly 1-point functions (fields) if Lmotifs==1;\n1-point and 2-point functions (fields, 2-mer forces) if Lmotifs==2;\n1-point, 2-point and 3-point functions (fields, 2-mer and 3-mer forces) if Lmotifs==3.\n\nLmodel is the number of nucleotides used for the inference. If seqs is made of  sequences of constant length (or if it is String and not a vector), the default option will use the sequences length as Lmodel; otherwise, the default behaviour is to take Lmodel = 5000.\n\npseudocount_param is the fraction of weigth coming from random uniform sequences in the computation of the motif's frequences from the data.\n\nIf fast, the partition function is estimated through the top eigenvalue of the transfer matrix alone (much faster, but slightly less precise, expecially for short sequences). The default value, \"auto\", automatically uses the fast  evaluation for long sequences.\n\ntolerance and max_iter are parameters for the Newton-Raphson algorithm  used to solve the system of equations.\n\nZS_gauge specifies whether the result has to be put in the zero sum gauge before being returned.\n\n\n\n\n\n","category":"function"},{"location":"#MaxEntNucleotideBiases.gauge_mask_variables-Tuple{Vector{String}}","page":"Home","title":"MaxEntNucleotideBiases.gauge_mask_variables","text":"gauge_mask_variables(motifs::Vector{String})\n\nThis function returns a Array{Bool} corresponding to motifs that have to be inferred after all motifs containing Ts are fixed to 0 thanks to gauge transformations.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.gauge_zerosum-Tuple{MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.gauge_zerosum","text":"gauge_zerosum(model_pars::NucleotideModel)\n\nThis function takes as input a dictionary describing the parmeters of a model and changes the gauge into the \"zero-sum\" gauge.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.generate_transfer_matrix-Tuple{MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.generate_transfer_matrix","text":"generate_transfer_matrix(model::NucleotideModel)\n\nReturn the trasfer matrix used for the computation of the partition function. The first index correspond to the left-most nucleotide, the last index to the right-most one.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.generate_transfer_matrix_last-Tuple{MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.generate_transfer_matrix_last","text":"generate_transfer_matrix_last(model::NucleotideModel)\n\nReturn the last transfer matrix used for the computation of the partition function. The first index correspond to the left-most nucleotide, the last index to the right-most one.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.get_forces_dict-Tuple{MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.get_forces_dict","text":"get_forces_dict(model::NucleotideModel)\n\nReturn a dictionary of the form motif => force.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.readmodel-Tuple{AbstractString}","page":"Home","title":"MaxEntNucleotideBiases.readmodel","text":"readmodel(fname::AbstractString)\n\nRead file (that should have been written by writemodel)  to get a NucleotideModel, which is returned.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.sample_metropolis-Tuple{Int64, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.sample_metropolis","text":"sample_metropolis(L::Int, model::NucleotideModel; \n                         beta::Float64=1.0, Nsamples::Int=1, Nsteps::Int=1, Ntherm::Int=L*10, \n                         startseq::Union{String,Missing}=missing)\n\nSample Nsample sequences of length L at inverse temperature beta using a model given in model.  Ntherm is the number of MonteCarlo steps done before starting collecting sequences,  startseq is the starting sequence.\n\n\n\n\n\n","category":"method"},{"location":"#MaxEntNucleotideBiases.writemodel-Tuple{AbstractString, MaxEntNucleotideBiases.NucleotideModel}","page":"Home","title":"MaxEntNucleotideBiases.writemodel","text":"writemodel(fname::AbstractString, model::NucleotideModel)\n\nWrite model to fname.\n\n\n\n\n\n","category":"method"}]
}
