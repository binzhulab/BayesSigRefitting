
\name{BayesSigRefitting-package}
\alias{BayesSigRefitting}
\docType{package}
\title{
    A Bayesian Approach to Select Mutational Signatures
}
\description{
    A Bayesian approach to select and refit an optimal set of signatures from a 
    collection of reference signatures. BayesSigRefitting specifies a Bayesian hierarchical 
    model for signature selection, imposes model sparsity using Laplace prior, 
    and utilizes the Shotgun Stochastic Search algorithm to select reference signatures.
}
\details{
Mutational signatures are the patterns of somatic mutations within the cancer genome. 
Analysis of mutational signatures helps understand causes of
cancer and inform treatment strategies. 
When the sample size is limited, signature refitting is commonly used. 
It requires selecting a set of reference signatures to be refitted, 
which has not been sufficiently explored. 
This package utilizes a Bayesian approach to select and refit an optimal
set of signatures from a collection of reference signatures. 
\bold{BayesSigRefitting} specifies a Bayesian hierarchical model for signature selection, 
imposes model sparsity using Laplace prior, and utilizes the 
Shotgun Stochastic Search algorithm to select reference signatures.

The two main functions in this package are \code{\link{BayesSigRef}} and 
\code{\link{SigRefitting}}.
}

\references{

}

\author{Min Hua <min.hua@nih.gov> and Bin Zhu <bin.zhu@nih.gov>}
\keyword{package}

