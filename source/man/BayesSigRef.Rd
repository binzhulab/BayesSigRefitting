
\name{BayesSigRef}
\alias{BayesSigRef}
\title{BayesSigRef}
\description{
    A Bayesian approach to select and refit an optimal set of signatures from a 
    collection of reference signatures.
}
\usage{ BayesSigRef(V, Sig.W="SBS", n.iter=3000, tol=0.05, n.core=1, sig.select=0.8)
}
\arguments{
    \item{V}{A numeric matrix with 96 rows and n columns with the rows representing
    the 96 mutation types and the columns representing the tumors. This is
    the Single Base Substitution (SBS) mutation counts matrix.
    The order of the mutation types should match the order of mutation types in 
    COSMIC SBS signatures, users need to rearrange the data if the order does not match.}

    \item{Sig.W}{The string "SBS" or a numeric matrix with 96 rows where
    the rows represent the 96 mutation types, and the columns represent the signatures.
    This is the reference signatures from which signatures for refitting will be selected. 
    The default setting is "SBS". If set to "SBS", then the COSMIC SBS 96 V3.1 
    signatures will be used. Users have the freedom to set it to the 
    reference signatures of their own choice, however, the order of the mutation types 
    should match that of the COSMIC SBS 96 category.}

    \item{n.iter}{The number of iterations the algorithm
    will repeat to conduct model searching. 
    A larger value gives better convergence, however, requires longer run time. 
    The default is 3000.}

    \item{tol}{The convergence tolerance for the algorithm. 
    A smaller value gives better convergence, however requires a 
    longer run time.
    The default is 0.05.}

    \item{n.core}{The number of cpu cores to use. The default is 1. }
    \item{sig.select}{The cut-off of the probability the signature will be selected. 
        The default is 0.8, meaning a signature will be selected only if its 
        posterior probability is greater than 0.8. 
        Users can change the value according to their own domain knowledge 
        and/or tolerance. }
}

\details{
The algorithm consistents of two steps, signature selection and 
signature refitting.
In the first step, \code{BayesSigRef} specifies a Bayesian hierarchical 
model for signature selection, 
imposes model sparsity using Laplace prior, 
and utilizes the Shotgun Stochastic Search algorithm
to select reference signatures.
The second step estimates the signature activities for each tumor 
using the selected signatures.
}

\value{ A list containing the objects:
    \itemize{
    \item \code{selection}: A named vector of the selected signatures, where the names
          are the signatures and the values are the locations of the selected signatures
          in the reference signatures. For example, \code{(SBS13, 19)} 
          means the selected signature is SBS13 and it is the 19th signature in the 
          reference signature set.

    \item \code{refitting}: A list containing the results of the signature
          selection and refitting, code{W} and \code{H} matrices. 
          \code{W} is the corresponding signature matrix, 
          and \code{H} is the corresponding activity matrix.

    \item \code{objects}: A list containing objects needed for calling 
        \code{\link{SigRefitting}} and the results of the signature selection step.
    }
}
%\references{}
\author{
Min Hua <min.hua@nih.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
     # Importing data
    data(SimData_V, package="BayesSigRefitting")

    # Signature selection and refitting (the default options could take considerable time)
    \donttest{ref <- BayesSigRef(SimData_V, Sig.W="SBS", n.iter=3000, tol=0.05, 
                                 n.core=5, sig.select=0.8)

    # Show the selected signatures
    ref$selection

    # Show the W matrix
    ref$refitting$W

    # Show the H matrix 
    ref$refitting$H
   }
}
\keyword{ mutational signatures } 
\seealso{ \code{\link{SigRefitting}} }


