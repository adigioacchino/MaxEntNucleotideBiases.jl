using MaxEntNucleotideBiases
using Test
using Random
using Statistics
using Scratch

########################################################
# useful variables
########################################################
# this is an Influenza H5N1 PB2 segment, strain used: A/Anhui/1/2005
testseq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCTGAAAGGCTAAAACATGGAACCTTCGGTCCCGTCCATTTTCGAAACCAAGTTAAAATACGCCGCCGAGTTGATATAAATCCTGGCCATGCAGATCTCAGTGCTAAAGAAGCACAAGATGTCATCATGGAGGTCGTTTTCCCAAATGAAGTGGGAGCTAGAATATTGACATCAGAGTCACAATTGACAATAACGAAAGAGAAGAAAGAAGAGCTCCAAGATTGTAAGATTGCTCCCTTAATGGTTGCATACATGTTGGAAAGGGAACTGGTCCGCAAAACCAGATTCCTACCGGTAGCAAGCGGAACAAGCAGTGTGTACATTGAGGTATTGCATTTGACTCAAGGGACCTGCTGGGAACAGATGTACACTCCAGGCGGAGAAGTGAGAAACGACGATGTTGACCAGAGTTTGATCATCGCTGCCAGAAACATTGTTAGGAGAGCAACGGTATCAGCGGATCCACTGGCATCACTGCTGGAGATGTGTCACAGCACACAAATTGGTGGGATAAGGATGGTGGACATCCTTAGGCAAAACCCAACTGAGGAACAAGCTGTGGGTATATGCAAAGCAGCAATGGGTCTGAGGATCAGTTCATCCTTTAGCTTTGGAGGCTTCACTTTCAAAAGAACAAGTGGATCATCCGTCACGAAGGAAGAGGAAGTGCTTACAGGCAACCTCCAAACATTGAAAATAAGAGTACATGAGGGGTATGAAGAGTTCACAATGGTTGGACGGAGGGCAACAGCTATCCTGAGGAAAGCAACTAGAAGGCTGATTCAGTTGATAGTAAGTGGAAGAGACGAACAATCAATCGCTGAGGCAATCATTGTAGCAATGGTGTTCTCACAGGAGGATTGCATGATAAAGGCAGTCCGGGGCGATTTGAATTTCGTAAACAGAGCAAACCAAAGATTAAACCCCATGCATCAACTCCTGAGACATTTTCAAAAGGACGCAAAAGTGCTATTTCAGAATTGGGGAATTGAACCCATTGATAATGTCATGGGGATGATCGGAATATTACCTGACCTGACTCCCAGCACAGAAATGTCACTGAGAAGAGTAAGAGTTAGTAAAGTGGGAGTGGATGAATATTCCAGCACTGAGAGAGTAATTGTAAGTATTGACCGTTTCTTAAGGGTTCGAGATCAGCGGGGGAACGTACTCTTATCTCCCGAAGAGGTCAGCGAAACCCAGGGAACAGAGAAATTGACAATAACATATTCATCATCAATGATGTGGGAAATCAACGGTCCTGAGTCAGTGCTTGTTAACACCTATCAATGGATCATCAGAAACTGGGAAACTGTGAAGATTCAATGGTCTCAAGACCCCACGATGCTGTACAATAAGATGGAGTTTGAACCGTTCCAATCCTTGGTACCTAAGGCTGCCAGAGGTCAATACAGTGGATTTGTGAGAACACTATTCCAACAAATGCGTGACGTACTGGGGACATTTGATACTGTCCAGATAATAAAGCTGCTACCATTTGCAGCAGCCCCACCAGAGCAGAGCAGAATGCAGTTTTCTTCTCTAACTGTGAATGTGAGAGGCTCAGGAATGAGAATACTCGTAAGGGGCAATTCCCCTGTGTTCAACTACAATAAGGCAACCAAAAGGCTTACCGTTCTTGGAAAGGACGCAGGTGCATTAACAGAGGATCCAGATGAGGGGACAACCGGAGTGGAGTCTGCAGTACTGAGGGAATTCCTAATTCTAGGCAAGGAGGACAAAAGATATGGACCAGCATTGAGTATCAATGAACTGAGCAACCTTGCGAAAGGGGAGAAAGCTAATGTGCTGATAGGACAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTCTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAG"
nts = ["A","C","G","T"]
########################################################
########################################################
@testset verbose=true "MaxEntNucleotideBiases tests" begin
    ########################################################
    println("Running preliminary computations...")
    println("(using $(Threads.nthreads()) threads for sampling)") 
    testseq5000 = (testseq^3)[1:5000] # if a sequence too short is used the fact that the gauge freedoom is exact only 
                                      # for L->\infty creates some issues for the i=3 case during testing
                                      # (things that should be >0 are about -0.005 or so).
    max_loglik_mods = [MaxEntNucleotideBiases.fitmodel(testseq5000, i, ZS_gauge=false) for i in 1:3]
    
    samples = Vector{String}[]
    for i in 1:3
        t_sample = [String[], String[]]
            Threads.@threads for j in 1:2
                t_sample[j] = MaxEntNucleotideBiases.sample_metropolis(5000, max_loglik_mods[i], 
                                                Nsamples=2500, Nsteps=5000, Ntherm=50000, beta=1.)
            end
        t_sample = vcat(t_sample...)
        push!(samples, t_sample)
    end

    println("Done! Starting tests.")
    ########################################################
    
    # test that the constructors work as expected
    @testset "checks on NucleotideModel constructor" begin
        for i in 1:3
            mod = max_loglik_mods[i]
            all_motifs = mod.motifs
            all_forces = mod.forces
            some_motifs = all_motifs[1:end-2] 
            some_forces = all_forces[1:end-2]
            new_mod = MaxEntNucleotideBiases.NucleotideModel(some_motifs, some_forces)
            @test sort(new_mod.motifs) == sort(all_motifs)
            new_fors_dict = MaxEntNucleotideBiases.get_forces_dict(new_mod)
            orig_fors_dict = MaxEntNucleotideBiases.get_forces_dict(mod)
            for m in some_motifs
                @test new_fors_dict[m] == orig_fors_dict[m]
            end
            for m in all_motifs[end-1:end]
                @test new_fors_dict[m] == 0.
            end
        end
    end

    # simple test to see that main functions run correctly
    @testset "checks on random uniform sequences inference" begin
        randseq = join(rand(Xoshiro(0), nts, 5_000_000))
        mod1 = MaxEntNucleotideBiases.get_forces_dict(MaxEntNucleotideBiases.fitmodel(randseq, 1, 5000))
        for k in keys(mod1)
            @test exp(mod1[k]) ≈ 0.25 atol=0.02
        end
        mod2 = MaxEntNucleotideBiases.get_forces_dict(MaxEntNucleotideBiases.fitmodel(randseq, 2, 5000))
        for k in keys(mod2)
            @test mod2[k] ≈ 0. atol=0.02
        end
        mod3 = MaxEntNucleotideBiases.get_forces_dict(MaxEntNucleotideBiases.fitmodel(randseq, 3, 5000))
        for k in keys(mod3)
            @test mod3[k] ≈ 0. atol=0.05
        end
    end

    # test for model fitting from frequences directly
    @testset "checks fit from frequences" begin
        for k in 1:3
            forcedict = MaxEntNucleotideBiases.get_forces_dict(max_loglik_mods[k])
            all_mots = vcat([join.([collect(Iterators.product([nts for _ in 1:kk]...))...]) for kk in 1:k]...)
            all_freqs = MaxEntNucleotideBiases._compute_freqs([testseq5000], all_mots)
            t_mod = MaxEntNucleotideBiases.fitmodel(all_freqs, all_mots, 5000, ZS_gauge=false)
            t_forcedict = MaxEntNucleotideBiases.get_forces_dict(t_mod)
            for m in all_mots
                @test forcedict[m] ≈ t_forcedict[m] atol=0.001
            end
        end
    end

    # take a sequence, infer a model, do variation of all parameters, check that loglik is lower
    @testset "checks on loglikelihoods" begin
        for i in 1:3
            max_loglik_model = max_loglik_mods[i]
            loglik = MaxEntNucleotideBiases.compute_loglikelihood(testseq5000, max_loglik_model)
            epsilon = 0.05
            for n in max_loglik_model.motifs
                alt_model = MaxEntNucleotideBiases.get_forces_dict(max_loglik_model)
                alt_model[n] = alt_model[n] + epsilon
                alt_loglik = MaxEntNucleotideBiases.compute_loglikelihood(testseq5000, 
                    MaxEntNucleotideBiases.NucleotideModel(max_loglik_model.motifs, [alt_model[m] for m in max_loglik_model.motifs]))
                @test loglik > alt_loglik skip = (i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
                alt_model = MaxEntNucleotideBiases.get_forces_dict(max_loglik_model)
                alt_model[n] = alt_model[n] - epsilon
                alt_loglik = MaxEntNucleotideBiases.compute_loglikelihood(testseq5000, 
                    MaxEntNucleotideBiases.NucleotideModel(max_loglik_model.motifs, [alt_model[m] for m in max_loglik_model.motifs]))
                @test loglik > alt_loglik skip = (i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
            end
        end
    end


    # sample sequences and check number of motifs
    @testset "checks on sampling and inference" begin
        for i in 1:3
            max_loglik_model = max_loglik_mods[i]
            sample = samples[i]
            for m in max_loglik_model.motifs
                c1 = count(m, testseq5000, overlap=true)
                c2 = mean(count.(m, sample, overlap=true))
                @test (c1-c2)/c1 ≈ 0. atol=0.025
            end
        end
    end

    # test multiseq: take a sequence, split and check that the result is consistent
    @testset "multiseq inference" begin
        for i in 1:3
            # test if using many times the same sequences give the same result
            max_loglik_model = MaxEntNucleotideBiases.get_forces_dict(max_loglik_mods[i])
            multiseq_model = MaxEntNucleotideBiases.get_forces_dict(MaxEntNucleotideBiases.fitmodel([testseq5000 for _ in 1:10], i, 5000, ZS_gauge = false))
            @testset "same sequence multiple times" for m in keys(max_loglik_model)
                @test multiseq_model[m] ≈ max_loglik_model[m] atol=0.001
            end
            # test if using sample sequences from the same model gives a similar result
            sample = samples[i]
            multiseq_model = MaxEntNucleotideBiases.get_forces_dict(MaxEntNucleotideBiases.fitmodel([testseq5000 for _ in 1:10], i, 5000, ZS_gauge = false))
            @testset "multiple sampled sequences" for m in keys(max_loglik_model)
                @test multiseq_model[m] ≈ max_loglik_model[m] atol=0.05
            end
        end
    end

    # test zero-sum gauge 
    @testset "change gauge into zerosum" begin
        for i in 1:3
            mod = max_loglik_mods[i]
            mod_zg = MaxEntNucleotideBiases.gauge_zerosum(mod)
            ll1 = MaxEntNucleotideBiases.compute_loglikelihood(testseq5000, mod) 
            ll2 = MaxEntNucleotideBiases.compute_loglikelihood(testseq5000, mod_zg)
            @test (ll1 - ll2) / ll1 ≈ 0. atol=1e-3 skip=(i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
            mot1s = [x for x in mod.motifs if length(x)==1] 
            mot2s = [x for x in mod.motifs if length(x)==2]
            mot3s = [x for x in mod.motifs if length(x)==3]
            mod_zg = MaxEntNucleotideBiases.get_forces_dict(mod_zg)
            @test sum([mod_zg[x] for x in mot1s]) ≈ 0. atol=1e-9
            if i >= 2
                @test sum([mod_zg[x] for x in mot2s]) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot2s if (x[1] == m2[1])]) for m2 in mot2s])) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot2s if (x[2] == m2[2])]) for m2 in mot2s])) ≈ 0. atol=1e-9
            end
            if i >= 3
                @test sum([mod_zg[x] for x in mot3s]) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot3s if (x[1:2] == m3[1:2])]) for m3 in mot3s])) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot3s if (x[2:3] == m3[2:3])]) for m3 in mot3s])) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot3s if (x[1] == m3[1])]) for m3 in mot3s])) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot3s if (x[2] == m3[2])]) for m3 in mot3s])) ≈ 0. atol=1e-9
                @test maximum(abs.([mean([mod_zg[x] for x in mot3s if (x[3] == m3[3])]) for m3 in mot3s])) ≈ 0. atol=1e-9
            end
        end
    end

    # test KL divergence 
    @testset "KL divergence with gauge changes" begin
        mod = max_loglik_mods[3]
        mod_zg = MaxEntNucleotideBiases.gauge_zerosum(mod)
        @test MaxEntNucleotideBiases.compute_symmetrized_kl(5000, mod, mod_zg) ≈ 0. atol=0.1
    end

    # test fast evaluation
    @testset "Fast computation of compute_logz" begin
        for i in 2:3
            mod = max_loglik_mods[i]
            mod_f = MaxEntNucleotideBiases.fitmodel(testseq5000, i, fast=true)
            @test MaxEntNucleotideBiases.compute_symmetrized_kl(5000, mod, mod_f) ≈ 0. atol=0.1
        end
    end

    # test model write/read
    @testset "model IO" begin
        dir = get_scratch!("test_IO")
        for i in 1:3
            mod = max_loglik_mods[i]
            fpath = joinpath(dir, "tempfile.txt")
            MaxEntNucleotideBiases.writemodel(fpath, mod)
            @test readdir(dir) == ["tempfile.txt"]
            mod2 = MaxEntNucleotideBiases.readmodel(fpath)
            @test all(mod.motifs .== mod2.motifs)
            @test all(mod.forces .== mod2.forces)
            @test all(mod.Lmotifs .== mod2.Lmotifs)
        end
        delete_scratch!("test_IO")
    end
end