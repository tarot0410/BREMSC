# BREMSC
An R package for joint clustering with single cell transcriptomics (scRNA-seq) and proteomics (CITE-seq) data.

Version: 0.1.0 (Date: 2019-11-26)

See [Homepage @ Wei Chen's Lab](http://www.pitt.edu/~wec47/singlecell.html)

## Introduction
**BREMSC** is an R package (with core function *BREMSC*) for joint clustering droplet-based scRNA-seq and CITE-seq data. As an extension of DIMMSC, it uses separate Dirichlet mixture priors to characterize variations across cell types for each data source, and uses random effects to incorporate the two data sources. A Bayesian framework with Gibbs-sampling is used for parameter estimation. This package can directly work on raw count data from droplet-based scRNA-seq and CITE-seq experiments without any data transformation, and it can provide clustering uncertainty for each cell.

## Installation

Install DIMMSC from Github
```
library(devtools)
install_github("tarot0410/BREMSC")
```
Or terminal command (first download *BREMSC* source file from Wei Chen's Lab website)
```
R CMD INSTALL BREMSC_0.1.0.tar
```

## Usage
```
BREMSC(dataProtein, dataRNA, K, nCores = 1, nMCMC = 500, useGene = 100, diagnostic = T, sd_alpha = 1, sd_b = 0.5, sigmaB = 0.8)

```

## Arguments
* *dataProtein* : a D*C matrix with D proteins and C cells
* *dataRNA* : a G*C matrix with G genes and C cells
* *K* : number of desired clusters
* *nCores* : number of CPUs to use
* *nMCMC* : number of MCMCs
* *useGene* : number of genes to keep
* *diagnostic* : whether store a dignostic plot
* *sd_alpha* : random walk depth for alpha
* *sd_b* : random walk depth for b
* *sigmaB* : initial estimate of SD of random effects

## Values
BREMSC returns a list object containing:
* *clusterID* : estimated label for each cell
* *posteriorProb* : a C*K matrix with probability that each cell belongs to each cluster
* *sdRF* : estimated SD of random effects
* *logLik* : estimated log likelihood
* *vecB* : estimated random effects for each cell
* *alphaMtxProtein* : a K*D matrix of alpha estimates for protein source
* *alphaMtxRNA* : a K*G matrix of alpha estimates for RNA source

## Example:
```
# First load pacakge
library(BREMSC)

# Next load the example simulated data (dataADT: protein data; dataRNA: RNA data)
data("dataADT")
data("dataRNA")

# Test run of BREMSC (using small number of MCMC here to save time)
# nMCMC>200 is suggested in real application
# Parallel computing is highly recommended. nCores>=40 is recommended for ~10,000 cells.
testRun <- BREMSC(dataADT, dataRNA, K=4, nCores=5, nMCMC=20)

```

## Publications
* **Xinjun Wang**, Zhe Sun, Yanfu Zhang, Heng Huang, Kong Chen, Ying Ding, Wei Chen. BREM-SC: A Bayesian Random Effects Mixture Model for Joint Clustering Single Cell Multi-omics Data. *Submitted* 2019.

## Contact
Xinjun Wang (xiw119@pitt.edu), [Wei Chen](http://www.pitt.edu/~wec47/index.html).
