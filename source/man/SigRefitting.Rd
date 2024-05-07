
\name{SigRefitting}
\alias{SigRefitting}
\title{Signature Refitting}
\description{
    Estimate the signature activities for each tumor using the selected signatures.
}
\usage{ SigRefitting(BSR.obj, sig.select=0.8)
}
\arguments{
    \item{BSR.obj}{The returned object from \code{\link{BayesSigRef}}. }
    \item{sig.select}{The lower bound to select signatures. The default is 0.8. }
}

\details{This function allows users to change the \code{sig.select} option 
    without re-running the signature selection. 
}

\value{ A list containing the objects:
    \itemize{
    \item \code{selection}: A named vector of the selected signatures, where the names
          are the signatures and the values are the locations of the selected signatures
          in the reference signatures. For example, \code{(SBS13, 19)} 
          means the selected signature is SBS13 and it is the 19th signature in the 
          reference signature set.

    \item \code{W}: A numeric matrix with 96 rows and k columns containing
         the selected signatures, where k is the number of the 
         selected signatures.

    \item \code{H}: A numeric matrix with k rows containing the
         signature activities of the selected signatures in W.
    }
}
%\references{}
\author{
Min Hua <min.hua@nih.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
    data(BSR.obj, package="BayesSigRefitting")

    \donttest{ret <- SigRefitting(BSR.obj, sig.select=0.95)}
}
\keyword{ mutational signatures } 
\seealso{ \code{\link{BayesSigRef}} }


