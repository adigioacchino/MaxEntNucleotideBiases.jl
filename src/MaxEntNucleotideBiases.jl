module MaxEntNucleotideBiases

    using LinearAlgebra
    using Statistics
    using FiniteDiff

    const dna_alphabet = ["A", "C", "G", "T"]
    const dna2num_alphabet = Dict(["A" => 1, "C" => 2, "G" => 3, "T" => 4])
    const len_alphabet = 4
    
    include("MENB_PartitionFunction.jl") # functions to compute the partition function
    include("MENB_Gauge.jl") # functions to deal with the gauge degrees of freedom
    include("MENB_InferModel.jl") # main functions to infer the model
    include("MENB_Utils.jl") # additional functions to compute energy, marginals, sample sequences, etc...
    
end
