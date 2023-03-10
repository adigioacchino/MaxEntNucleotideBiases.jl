module MaxEntNucleotideBiases

    using LinearAlgebra
    using Statistics
    using FiniteDiff

    const dna_alphabet = ["A", "C", "G", "T"]
    const dna2num_alphabet = Dict(["A" => 1, "C" => 2, "G" => 3, "T" => 4])
    const len_alphabet = 4
 
    struct NucleotideModel
        motifs::Vector{String}
        forces::Vector{Float64}
        Lmotifs::Int
        function NucleotideModel(mots::Vector{<:AbstractString}, fors::Vector{<:Number}, Lmots::Int)
            (maximum(length.(mots)) != Lmots) && error("inconsistent Lmotifs")
            (length(mots) != length(fors)) && error("inconsistent number of motifs and forces")
            return new(mots, fors, Lmots)
        end
    end
    NucleotideModel(mots::Vector{<:AbstractString}, fors::Vector{<:Number}) = NucleotideModel(mots, fors, maximum(length.(mots)))

    """
        ForcesDict(model::NucleotideModel)

    Return a dictionary of the form motif => force.
    """
    function ForcesDict(model::NucleotideModel)
        Dict(zip(model.motifs, model.forces))
    end


    include("MENB_PartitionFunction.jl") # functions to compute the partition function
    include("MENB_Gauge.jl") # functions to deal with the gauge degrees of freedom
    include("MENB_InferModel.jl") # main functions to infer the model
    include("MENB_Utils.jl") # additional functions to compute energy, marginals, sample sequences, etc...
    
end
