```@meta
CurrentModule = MaxEntNucleotideBiases
```

# MaxEntNucleotideBiases

Documentation for [MaxEntNucleotideBiases](https://github.com/adigioacchino/MaxEntNucleotideBiases.jl).

## Installation
This package is not registered. Install with:

```julia
import Pkg
Pkg.add(url="https://github.com/adigioacchino/MaxEntNucleotideBiases.jl")
```
## Introduction and references
TBD
## Examples

### Basic usage
We will take sub-sequence from the Influenza H5N1 PB2 segment (strain used: A/Anhui/1/2005) and we will use it as working example to test the function of the package.
We start by loading the package and the sequence:
```julia
using MaxEntNucleotideBiases
example_seq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT"
```

```@setup load_module_sequence
using MaxEntNucleotideBiases
example_seq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT"
```

```@setup fit_models
using MaxEntNucleotideBiases
example_seq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT"
mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2)
mod3 = MaxEntNucleotideBiases.fitmodel(example_seq, 3)
```

The package allows to fit up to 1-, 2- or 3-mers (i.e. single nucleotide, dinucleotide or trinucleotide) forces. 
To do so, we will use the `fitmodel` function, whose second argument is the 'k' of the maximum k-mers we want to consider (1, 2 or 3).

We will start by fitting the single nucleotide biases:
```@repl load_module_sequence
mod1 = MaxEntNucleotideBiases.fitmodel(example_seq, 1)
```
The result is a `NucleotideModel` object, which contains the motifs, the fitted parameters and the value of 'k'.
Similarly, we can fit dinucleotide and trinucleotide biases:
```@repl load_module_sequence
mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2)
mod3 = MaxEntNucleotideBiases.fitmodel(example_seq, 3)
```

From a model we can obtain a dictionary (motif => force) with the `get_forces_dict` function:
```@repl fit_models
MaxEntNucleotideBiases.get_forces_dict(mod2)
```

Finally, we can compute several interesting quantities:
 - the energy of any sequence with a given model, with the `compute_minus_energy` function:
```@repl fit_models
L = 1000;
seq = join(rand(["A","C","G","T"], L));
MaxEntNucleotideBiases.compute_minus_energy(mod3, seq)
```
- the log-likelihood of any sequence with a given model, with the `compute_loglikelihood` function:
```@repl fit_models
L = 1000;
seq = join(rand(["A","C","G","T"], L));
MaxEntNucleotideBiases.compute_loglikelihood(mod2, seq)
```

### Gauge choices
There are several gauge choices that can be used to fit the model.
The default choice is the 'zero-sum' gauge, that is described in the paper.
The other notable choice is the 'lattice gas gauge', which is again described in the paper and is implemented in this package.
In particular, the results of the inference can be directly returned in the lattice gas gauge with the `ZS_gauge` argument passed as `false` to the `fitmodel` function:
```@repl load_module_sequence
mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2, ZS_gauge=false)
```
Moreover, the zero-sum gauge can be restored after the inference with the `gauge_zerosum` function:
```@repl load_module_sequence
mod2 = MaxEntNucleotideBiases.fitmodel(example_seq, 2, ZS_gauge=false) # hide
MaxEntNucleotideBiases.gauge_zerosum(mod2)
```

### Entropy and related quantities
The simple structure of the models fitted with this package allows to compute exactly several quantities of interest that, in most other cases, require (uncontrolled) approximations.
The first example is the entropy of the probability distribution defined with the model. 
It can be computed with the `compute_entropy` function, and since this quantity depends on the length of the sequences considered, a reference length must be passed as argument:
```@repl fit_models
L = 1000;
MaxEntNucleotideBiases.compute_entropy(L, mod3)
```

The entropy of the model can be used to quantify how different the probability distribution defined by the model is from the uniform distribution, which in turn is a measure of the 'pressure' that is exerted on the nucleotides to be used in a certain way.
In this package, the pressure ranges from 0 (unconstrained model) to 1 (fully constrained model), and it is computed with the `compute_pressure` function:
```@repl fit_models
L = 1000;
MaxEntNucleotideBiases.compute_pressure(L, mod3)
```

Finally, it can be interesting to quantify how much two models are different from each other. 
This package address this problem with the `compute_symmetrized_kl` function, which computes the symmetrized Kullback-Leibler divergence between two models:
```@repl fit_models
L = 1000;
mod3_unif = MaxEntNucleotideBiases.fitmodel(join(rand(["A","C","G","T"], L)), 3)
MaxEntNucleotideBiases.compute_symmetrized_kl(L, mod3, mod3_unif)
```

### Sampling
Another interseting feature of this package is the possibility to sample sequences from an inferred model.
This package implements a simple Metropolis routine to sample sequences from a model:

```@repl fit_models
L = 1000;
MaxEntNucleotideBiases.sample_metropolis(L, mod3, 
                                         Nsamples=10, Nsteps=3000, Ntherm=10_000, beta=1.)
```

### Model I/O
The model inference can take some time, so it can be useful to save the model on disk and load it later.
This can be done with the `save_model` and `load_model` functions:
```julia
MaxEntNucleotideBiases.writemodel(mod3, "mod3.menb")
mod3 = MaxEntNucleotideBiases.readmodel("mod3.menb")
```

## API

```@index
```

```@autodocs
Modules = [MaxEntNucleotideBiases]
```
