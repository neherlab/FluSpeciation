## Phylodynamics of rapidly adapting pathogens: extinction and speciation of a Red Queen.
### Le Yan, Richard A. Neher, and Boris I Shraiman

#### [prerint on biorxiv](https://www.biorxiv.org/content/early/2018/10/29/455444)
#### Abstract:
Rapidly evolving pathogens like influenza viruses can persist by accumulating antigenic novelty fast enough to evade the adaptive immunity of the host population, yet without continuous accumulation of genetic diversity. This dynamical state is often compared to the Red Queen evolving as fast as it can just to maintain its foothold in the host population: Accumulation of antigenic novelty is balanced by the build-up of host immunity. Such Red Queen States (RQS) of continuous adaptation in large rapidly mutating populations are well understood in terms of Traveling Wave (TW) theories of population genetics. Here we shall make explicit the mapping of the established Multi-strain Susceptible-Infected-Recovered (SIR) model onto the TW theory and demonstrate that a pathogen can persist in RQS if cross-immunity is long-ranged and its population size is large populations allowing for rapid adaptation.
We then investigate the stability of this state focusing on the rate of extinction and the rate of *speciation* defined as antigenic divergence of viral strains beyond the range of cross-inhibition. RQS states are transient, but in a certain range of evolutionary parameters can exist for the time long compared to the typical time to the most recent common ancestor.
In this range the steady TW is unstable and the antigenic advance of the lead strains relative to the typical co-circulating viruses tends to oscillate. This results in large fluctuations in prevalence that facilitate extinction. We shall demonstrate that the rate of TW fission into antigenically uncoupled viral populations is related to fluctuations of diversity and construct a *phase diagram* identifying different regimes of viral phylodynamics as a function of evolutionary parameters.

### Contents
 * the directory `flu_data` contains alignments and `Snakefile` to run `augur`. The output was visualized in `auspice` to generate Figure 1. The script `plot_all_Tmrca.py` generates the saw-tooth graphs in Fig.~1.
 * `flu_data/split_trees.py`: script that splits a global influenza B tree into the sublineages Victoria and Yamagata
 * `flu_data/plot_tmrca.py`: script that uses the output of the `augur` pipeline to produce the tabular files `{lineage}_tmrca_trajectory.dat`. The latter are tab separated files with year in the first and Tmrca in the second columns.
 * the MATLAB file `mfdist.m` implements the stochastic fitness class simulation of a large population in a one dimensional landscape. It the was used to produce simulation data on extinction times.
 * the file `pop2tw.m` is an analogous fitness class simulation of two coupled populations to investigate the speciation behavior.
 * the file `FluEpiTreeNM.m` implements a tree structured viral population use to simulate the SIR model.

