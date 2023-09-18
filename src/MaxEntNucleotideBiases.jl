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
            # set to zero the forces of the motifs of length up to Lmots that are not in mots
            all_mots = copy(mots)
            all_fors = copy(fors)
            for l in 1:Lmots
                mots_mat = collect(Iterators.product([dna_alphabet for _ in 1:l]...))
                mots_list = join.([mots_mat...])
                for m in mots_list
                    if !(m in mots)
                        push!(all_mots, m)
                        push!(all_fors, 0.0)
                    end
                end
            end
            return new(all_mots, all_fors, Lmots)
        end
    end
    NucleotideModel(mots::Vector{<:AbstractString}, fors::Vector{<:Number}) = NucleotideModel(mots, fors, maximum(length.(mots)))

    """
        get_forces_dict(model::NucleotideModel)

    Return a dictionary of the form motif => force.
    """
    function get_forces_dict(model::NucleotideModel)
        Dict(zip(model.motifs, model.forces))
    end

    """
        writemodel(fname::AbstractString, model::NucleotideModel)

    Write `model` to `fname`.
    """
    function writemodel(fname::AbstractString, model::NucleotideModel)
        open(fname,"w") do f
            for m in model.motifs
                write(f, "$(m) ")
            end
            write(f, "\n")
            for p in model.forces
                write(f, "$(p) ")
            end
        end
    end

    """
        readmodel(fname::AbstractString)

    Read `file` (that should have been written by `writemodel`) 
    to get a NucleotideModel, which is returned.
    """
    function readmodel(fname::AbstractString)
        lines = readlines(fname)
        motifs = split(lines[1])
        forces = parse.(Float64, split(lines[2]))
        return NucleotideModel(motifs, forces)
    end


    include("MENB_PartitionFunction.jl") # functions to compute the partition function
    include("MENB_Gauge.jl") # functions to deal with the gauge degrees of freedom
    include("MENB_InferModel.jl") # main functions to infer the model
    include("MENB_Utils.jl") # additional functions to compute energy, marginals, sample sequences, etc...
    
end
