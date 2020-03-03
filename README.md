# BREMSC
**BREMSC** is an R package (with core functions *jointDIMMSC* and *BREMSC*) for joint clustering droplet-based scRNA-seq and CITE-seq data. *jointDIMMSC* is developed as a direct extension of DIMMSC, which assumes full indenpendency between single cell RNA and surface protein data. To take the correlation between two data sources into consideration, we further develop *BREMSC*, which uses random effects to incorporate the two data sources. This package can directly work on raw count data from droplet-based scRNA-seq and CITE-seq experiments without any data transformation, and it can provide clustering uncertainty for each cell.

Version: 0.1.0 (Date: 2019-11-26)

See [Homepage @ Wei Chen's Lab](http://www.pitt.edu/~wec47/singlecell.html)

## Installation

Install BREMSC from Github
```
install.packages("devtools")
library(devtools)
install_github("tarot0410/BREMSC")
```
Or terminal command (first download *BREMSC* source file from Wei Chen's Lab website)
```
R CMD INSTALL BREMSC_0.1.0.tar
```

# Function *jointDIMMSC*
## Introduction
*jointDIMMSC* is developed as an extension of DIMMSC, which assumes full indenpendency between single cell RNA and surface protein data. We construct the joint likelihood of the two data sources as their product, and use EM algorithm for parameter inference. In practice, the computational speed for *jointDIMMSC* is much faster than *BREMSC*, but the model assumption is more stringent.

## Usage
```
jointDIMMSC(dataProtein, dataRNA, K, useGene = 100, maxiter = 100, tol = 1e-04, lik.tol = 0.01)
```

## Arguments
* *dataProtein* : a D*C matrix with D proteins and C cells
* *dataRNA* : a G*C matrix with G genes and C cells
* *K* : number of desired clusters
* *useGene* : number of genes to keep
* *maxiter* : maximum number of iterations
* *tol* : a convergence tolerance for the difference of vector pie between iterations
* *lik.tol* : a convergence tolerance for the difference of log-likelihoods between iterations

## Values
jointDIMMSC returns a list object containing:
* *clusterID* : estimated label for each cell
* *posteriorProb* : a C*K matrix with probability that each cell belongs to each cluster
* *logLik* : estimated log likelihood
* *alphaMtxProtein* : a K*D matrix of alpha estimates for protein source
* *alphaMtxRNA* : a K*G matrix of alpha estimates for RNA source

## Example:
```
# First load BREMSC R package
library(BREMSC)

# Next load the example simulated data (dataADT: protein data; dataRNA: RNA data)
data("dataADT")
data("dataRNA")

# Test run of jointDIMMSC
testRun <- jointDIMMSC(dataADT, dataRNA, K=4)
```

# Function *BREMSC*
## Introduction
Similar to *jointDIMMSC*, *BREMSC* uses separate Dirichlet mixture priors to characterize variations across cell types for each data source, but it further uses random effects to incorporate the two data sources. A Bayesian framework with Gibbs-sampling is used for parameter estimation. The computational speed for *BREMSC* is much slower than *jointDIMMSC*. In practice, nMCMC>500 is suggested in real application, and running with more than 3 chains (set as a parameter) are strongly recommended for better stability.

## Usage
```
BREMSC(dataProtein, dataRNA, K, nCores = 1, nMCMC = 500, useGene = 100, diagnostic = T, sd_alpha = 1, sd_b = 0.5, sigmaB = 0.8)

```

## Arguments
* *dataProtein* : a D*C matrix with D proteins and C cells
* *dataRNA* : a G*C matrix with G genes and C cells
* *K* : number of desired clusters
* *nChains* : number of chains to run in parallel (need to have >= nChains number of threads available), 3-5 is recommended, default 3
* *nMCMC* : number of MCMCs, default 1000
* *sd_alpha_range* : range of random walk depth for alpha, each chain uses a different depth, default c(0.5, 1.5)
* *sd_b_range* : range of random walk depth for b, each chain uses a different depth, default c(0.2, 1)
* *sigmaB* : initial estimate of SD of random effects, default 0.8

## Values
BREMSC returns a list object containing:
* *clusterID* : estimated label for each cell
* *posteriorProb* : a C*K matrix with probability that each cell belongs to each cluster
* *sdRF_final* : estimated SD of random effects from last MCMC
* *logLik_final* : estimated log likelihood from last MCMC
* *vecB* : estimated random effects for each cell
* *alphaMtxProtein* : a K*D matrix of alpha estimates for protein source
* *alphaMtxRNA* : a K*G matrix of alpha estimates for RNA source
* *vecLogLik* : a vector of log likelihood for each MCMC (used for convergence check!)
* *vecSDRF* : a vector of estimated SD of random effects from each MCMC

## Example:
```
# First load BREMSC R package
library(BREMSC)

# Next load the example simulated data (dataADT: protein data; dataRNA: RNA data)
data("dataADT")
data("dataRNA")

# Test run of BREMSC (using small number of MCMC here to save time)
testRun <- BREMSC(dataADT, dataRNA, K=4, nCores=2, nMCMC=100)
```

## Publications
* **Xinjun Wang**, Zhe Sun, Yanfu Zhang, Zhongli Xu, Hongyi Xin, Heng Huang, Richard H Duerr, Kong Chen, Ying Ding, Wei Chen. BREM-SC: A Bayesian Random Effects Mixture Model for Joint Clustering Single Cell Multi-omics Data. *Submitted* 2020.

## Contact
Xinjun Wang (xiw119@pitt.edu), [Wei Chen](http://www.pitt.edu/~wec47/index.html).
