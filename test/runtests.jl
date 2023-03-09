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
    println("Running preliminary computations...")
    testseq5000 = (testseq^3)[1:5000] # if a sequence too short is used the fact that the gauge freedoom is exact only 
                                      # for L->\infty creates some issues for the i=3 case during testing
                                      # (things that should be >0 are about -0.005 or so).
    max_loglik_mods = [MaxEntNucleotideBiases.ModelFit(testseq5000, i) for i in 1:3]
    println("Done! Starting tests.")

    # simple test to see that main functions run correctly
    @testset "random uniform sequences" begin
        randseq = join(rand(Xoshiro(0), nts, 5_000_000))
        mod1 = MaxEntNucleotideBiases.ModelFit(randseq, 1, 5000)
        for k in keys(mod1)
            @test exp(mod1[k]) ≈ 0.25 atol=0.02
        end
        mod2 = MaxEntNucleotideBiases.ModelFit(randseq, 2, 5000)
        for k in keys(mod2)
            @test mod2[k] ≈ 0. atol=0.02
        end
        mod3 = MaxEntNucleotideBiases.ModelFit(randseq, 3, 5000)
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
            for n in keys(max_loglik_model)
    #            println(n)
                alt_model = copy(max_loglik_model)
                alt_model[n] = alt_model[n] + epsilon
                alt_loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, alt_model)
                @test loglik > alt_loglik skip = (i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
                alt_model = copy(max_loglik_model)
                alt_model[n] = alt_model[n] - epsilon
                alt_loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, alt_model)
                @test loglik > alt_loglik skip = (i==1) # i=1 does not work because Z=1 is hardcoded for 1-pt models, and only holds for a specific gauge
            end
        end
    end


    # sample sequences and check number of motifs
    @testset "checks on sampling and inference" begin
        for i in 1:3
            max_loglik_model = max_loglik_mods[i]
            sample = MaxEntNucleotideBiases.MetropolisSampling(5000, max_loglik_model, Nsamples=5000, Nsteps=5000, Ntherm=50000, beta=1.)
            for m in keys(max_loglik_model)
                c1 = count(m, testseq5000, overlap=true)
                c2 = mean(count.(m, sample, overlap=true))
                @test (c1-c2)/c1 ≈ 0. atol=0.02
            end
        end
    end

    # test multiseq: take a sequence, split and check that the result is consistent
    @testset verbose=true "multiseq inference" begin
        for i in 1:3
            # test if using many times the same sequences give the same result
            max_loglik_model = max_loglik_mods[i]
            multiseq_model = MaxEntNucleotideBiases.ModelFit([testseq5000 for _ in 1:10], i, 5000)
            @testset "same sequence multiple times" for m in keys(max_loglik_model)
                @test multiseq_model[m] ≈ max_loglik_model[m] atol=0.001
            end

            sample = MaxEntNucleotideBiases.MetropolisSampling(5000, max_loglik_model, Nsamples=100, Nsteps=5000, Ntherm=50000, beta=1.)
            multiseq_model = MaxEntNucleotideBiases.ModelFit([testseq5000 for _ in 1:10], i, 5000)
            @testset "multiple sampled sequences" for m in keys(max_loglik_model)
                @test multiseq_model[m] ≈ max_loglik_model[m] atol=0.05
            end
        end
    end

end