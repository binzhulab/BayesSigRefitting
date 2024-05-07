# BayesSigRefitting
A Bayesian Approach to Select Mutational Signatures


### Introduction
A Bayesian approach to select and refit an optimal set of signatures from a
collection of reference signatures. BayesSigRefitting specifies a Bayesian hierarchical
model for signature selection, imposes model sparsity using Laplace prior,
and utilizes the Shotgun Stochastic Search algorithm to select reference signatures.
For more information please refer to the [user guide](https://github.com/binzhulab/BayesSigRefitting/blob/main/BayesSigRefitting-manual.pdf).
<br/>

### Installation
To install from Github, use the devtools R package:
```r
if (!requireNamespace("devtools", quietly = TRUE))  
	install.packages("devtools")
devtools::install_github("binzhulab/BayesSigRefitting/source")
```
Alternatively, download the package and follow the steps below. Download BayesSigRefitting_0.0.3.tar.gz (for Unix) or BayesSigRefitting_0.0.3.zip (for Windows, R version >= 4.0). To install BayesSigRefitting on Unix, enter the command (without the quotes) from a Unix prompt:
```
R CMD INSTALL BayesSigRefitting_0.0.3.tar.gz -l path_to_install_package
```
Alternatively, BayesSigRefitting_0.0.3.tar.gz (for Unix) or BayesSigRefitting_0.0.3.zip (for Windows, R version >= 4.0) from the [Github page](https://github.com/binzhulab/BayesSigRefitting) are available and one may use the following commands:
```
install.packages("BayesSigRefitting_0.0.3.tar.gz", repose = NULL, type = "source")
install.packages("BayesSigRefitting_0.0.3.zip", repose = NULL, type = "win.binary")
```
Once the installation is successful, it can be loaded on **R** by calling 
```
library(BayesSigRefitting)
```
<br/>

