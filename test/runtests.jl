using MaxEntNucleotideBiases
using Test
using Random
using Statistics

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
    max_loglik_mods = [MaxEntNucleotideBiases.ModelFit(testseq5000, i, ZS_gauge=false) for i in 1:3]
    
    samples = Vector{String}[]
    for i in 1:3
        t_sample = [String[], String[]]
            Threads.@threads for j in 1:2
                t_sample[j] = MaxEntNucleotideBiases.MetropolisSampling(5000, max_loglik_mods[i], 
                                                Nsamples=2500, Nsteps=5000, Ntherm=50000, beta=1.)
            end
        t_sample = vcat(t_sample...)
        push!(samples, t_sample)
    end

    println("Done! Starting tests.")
    ########################################################

    # simple test to see that main functions run correctly
    @testset "random uniform sequences" begin
        randseq = join(rand(Xoshiro(0), nts, 5_000_000))
        mod1 = MaxEntNucleotideBiases.ForcesDict(MaxEntNucleotideBiases.ModelFit(randseq, 1, 5000))
        for k in keys(mod1)
            @test exp(mod1[k]) ≈ 0.25 atol=0.02
        end
        mod2 = MaxEntNucleotideBiases.ForcesDict(MaxEntNucleotideBiases.ModelFit(randseq, 2, 5000))
        for k in keys(mod2)
            @test mod2[k] ≈ 0. atol=0.02
        end
        mod3 = MaxEntNucleotideBiases.ForcesDict(MaxEntNucleotideBiases.ModelFit(randseq, 3, 5000))
        for k in keys(mod3)
            @test mod3[k] ≈ 0. atol=0.05
        end
    end

    # take a sequence, infer a model, do variation of all parameters, check that loglik is lower
    @testset "checks on loglikelihoods" begin
        for i in 1:3
            max_loglik_model = max_loglik_mods[i]
            loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, max_loglik_model)
            epsilon = 0.05
            for n in max_loglik_model.motifs
                alt_model = MaxEntNucleotideBiases.ForcesDict(max_loglik_model)
                alt_model[n] = alt_model[n] + epsilon
                alt_loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, 
                    MaxEntNucleotideBiases.NucleotideModel(max_loglik_model.motifs, [alt_model[m] for m in max_loglik_model.motifs]))
                @test loglik > alt_loglik skip = (i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
                alt_model = MaxEntNucleotideBiases.ForcesDict(max_loglik_model)
                alt_model[n] = alt_model[n] - epsilon
                alt_loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, 
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
            max_loglik_model = MaxEntNucleotideBiases.ForcesDict(max_loglik_mods[i])
            multiseq_model = MaxEntNucleotideBiases.ForcesDict(MaxEntNucleotideBiases.ModelFit([testseq5000 for _ in 1:10], i, 5000, ZS_gauge = false))
            @testset "same sequence multiple times" for m in keys(max_loglik_model)
                @test multiseq_model[m] ≈ max_loglik_model[m] atol=0.001
            end
            # test if using sample sequences from the same model gives a similar result
            sample = samples[i]
            multiseq_model = MaxEntNucleotideBiases.ForcesDict(MaxEntNucleotideBiases.ModelFit([testseq5000 for _ in 1:10], i, 5000, ZS_gauge = false))
            @testset "multiple sampled sequences" for m in keys(max_loglik_model)
                @test multiseq_model[m] ≈ max_loglik_model[m] atol=0.05
            end
        end
    end

    # test zero-sum gauge 
    @testset "change gauge into zerosum" begin
        for i in 1:3
            mod = max_loglik_mods[i]
            mod_zg = MaxEntNucleotideBiases.ZerosumGauge(mod)
            ll1 = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, mod) 
            ll2 = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, mod_zg)
            @test (ll1 - ll2) / ll1 ≈ 0. atol=1e-3 skip=(i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
            mot1s = [x for x in mod.motifs if length(x)==1] 
            mot2s = [x for x in mod.motifs if length(x)==2]
            mot3s = [x for x in mod.motifs if length(x)==3]
            mod_zg = MaxEntNucleotideBiases.ForcesDict(mod_zg)
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
        mod_zg = MaxEntNucleotideBiases.ZerosumGauge(mod)
        @test MaxEntNucleotideBiases.ComputeSymmetrizedKL(5000, mod, mod_zg) ≈ 0. atol=0.1
    end

    # test fast evaluation
    @testset "Fast computation of EvalLogZ" begin
        for i in 2:3
            mod = max_loglik_mods[i]
            mod_f = MaxEntNucleotideBiases.ModelFit(testseq5000, i, fast=true)
            @test MaxEntNucleotideBiases.ComputeSymmetrizedKL(5000, mod, mod_f) ≈ 0. atol=0.1
        end
    end

end