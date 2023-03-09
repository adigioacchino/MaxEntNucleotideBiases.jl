using MaxEntNucleotideBiases
using Test
using Random

########################################################
# useful variables
########################################################
# this is an Influenza H5N1 PB2 segment, strain used: A/Anhui/1/2005
testseq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCTGAAAGGCTAAAACATGGAACCTTCGGTCCCGTCCATTTTCGAAACCAAGTTAAAATACGCCGCCGAGTTGATATAAATCCTGGCCATGCAGATCTCAGTGCTAAAGAAGCACAAGATGTCATCATGGAGGTCGTTTTCCCAAATGAAGTGGGAGCTAGAATATTGACATCAGAGTCACAATTGACAATAACGAAAGAGAAGAAAGAAGAGCTCCAAGATTGTAAGATTGCTCCCTTAATGGTTGCATACATGTTGGAAAGGGAACTGGTCCGCAAAACCAGATTCCTACCGGTAGCAAGCGGAACAAGCAGTGTGTACATTGAGGTATTGCATTTGACTCAAGGGACCTGCTGGGAACAGATGTACACTCCAGGCGGAGAAGTGAGAAACGACGATGTTGACCAGAGTTTGATCATCGCTGCCAGAAACATTGTTAGGAGAGCAACGGTATCAGCGGATCCACTGGCATCACTGCTGGAGATGTGTCACAGCACACAAATTGGTGGGATAAGGATGGTGGACATCCTTAGGCAAAACCCAACTGAGGAACAAGCTGTGGGTATATGCAAAGCAGCAATGGGTCTGAGGATCAGTTCATCCTTTAGCTTTGGAGGCTTCACTTTCAAAAGAACAAGTGGATCATCCGTCACGAAGGAAGAGGAAGTGCTTACAGGCAACCTCCAAACATTGAAAATAAGAGTACATGAGGGGTATGAAGAGTTCACAATGGTTGGACGGAGGGCAACAGCTATCCTGAGGAAAGCAACTAGAAGGCTGATTCAGTTGATAGTAAGTGGAAGAGACGAACAATCAATCGCTGAGGCAATCATTGTAGCAATGGTGTTCTCACAGGAGGATTGCATGATAAAGGCAGTCCGGGGCGATTTGAATTTCGTAAACAGAGCAAACCAAAGATTAAACCCCATGCATCAACTCCTGAGACATTTTCAAAAGGACGCAAAAGTGCTATTTCAGAATTGGGGAATTGAACCCATTGATAATGTCATGGGGATGATCGGAATATTACCTGACCTGACTCCCAGCACAGAAATGTCACTGAGAAGAGTAAGAGTTAGTAAAGTGGGAGTGGATGAATATTCCAGCACTGAGAGAGTAATTGTAAGTATTGACCGTTTCTTAAGGGTTCGAGATCAGCGGGGGAACGTACTCTTATCTCCCGAAGAGGTCAGCGAAACCCAGGGAACAGAGAAATTGACAATAACATATTCATCATCAATGATGTGGGAAATCAACGGTCCTGAGTCAGTGCTTGTTAACACCTATCAATGGATCATCAGAAACTGGGAAACTGTGAAGATTCAATGGTCTCAAGACCCCACGATGCTGTACAATAAGATGGAGTTTGAACCGTTCCAATCCTTGGTACCTAAGGCTGCCAGAGGTCAATACAGTGGATTTGTGAGAACACTATTCCAACAAATGCGTGACGTACTGGGGACATTTGATACTGTCCAGATAATAAAGCTGCTACCATTTGCAGCAGCCCCACCAGAGCAGAGCAGAATGCAGTTTTCTTCTCTAACTGTGAATGTGAGAGGCTCAGGAATGAGAATACTCGTAAGGGGCAATTCCCCTGTGTTCAACTACAATAAGGCAACCAAAAGGCTTACCGTTCTTGGAAAGGACGCAGGTGCATTAACAGAGGATCCAGATGAGGGGACAACCGGAGTGGAGTCTGCAGTACTGAGGGAATTCCTAATTCTAGGCAAGGAGGACAAAAGATATGGACCAGCATTGAGTATCAATGAACTGAGCAACCTTGCGAAAGGGGAGAAAGCTAATGTGCTGATAGGACAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTCTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAG"
nts = ["A","C","G","T"]
########################################################
########################################################

# simple test to see that main functions run correctly
@test_skip @testset "random uniform sequences" begin
    randseq = join(rand(Xoshiro(0), nts, 5_000_000))
    mod1 = MaxEntNucleotideBiases.ModelFit(randseq,1)
    for k in keys(mod1)
        @test exp(mod1[k]) ≈ 0.25 atol=0.02
    end
    mod2 = MaxEntNucleotideBiases.ModelFit(randseq, 2, 1000)
    for k in keys(mod2)
        @test mod2[k] ≈ 0. atol=0.02
    end
    mod3 = MaxEntNucleotideBiases.ModelFit(randseq, 3, 1000)
    for k in keys(mod3)
        @test mod3[k] ≈ 0. atol=0.05
    end
end

# take a sequence, infer a model, do variation of all parameters, check that loglik is lower
@testset "checks on loglikelihoods" begin
    for i in 1:3
        testseq5000 = (testseq^3)[1:5000] # if a sequence too short is used the fact that the gauge freedoom is exact for L->\infty
                                          # creates some issues for the i=3 case.
        max_loglik_model = MaxEntNucleotideBiases.ModelFit(testseq5000, i)
        loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, max_loglik_model)
        epsilon = 0.05
        for n in keys(max_loglik_model)
#            println(n)
            alt_model = copy(max_loglik_model)
            alt_model[n] = alt_model[n] + epsilon
            alt_loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, alt_model)
            if i == 1
                @test loglik > alt_loglik skip = true # does not work because Z=1 is hardcoded, and only holds for a specific gauge
            else
#                println(loglik - alt_loglik)
                @test loglik > alt_loglik 
            end
            alt_model = copy(max_loglik_model)
            alt_model[n] = alt_model[n] - epsilon
            alt_loglik = MaxEntNucleotideBiases.ComputeLoglikelihood(testseq5000, alt_model)
            if i == 1
                @test loglik > alt_loglik skip = true # does not work because Z=1 is hardcoded, and only holds for a specific gauge
            else
#                println(loglik - alt_loglik)
                @test loglik > alt_loglik 
            end
        end
    end
end


# sample sequences and check number of motifs