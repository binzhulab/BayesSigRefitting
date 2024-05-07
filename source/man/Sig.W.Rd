% --- Source file: data.Rd ---
\name{Sig.W}
\alias{Sig.W}

\title{
    COSMIC SBS Signatures
}
\description{
    The set of 96 COSMIC Single Base Substitution (SBS) Signatures 
}
\details{ 
    A 96 x 60 matrix, which contains 60 known reference SBS signatures
    from COSMIC SBS Mutational Signatures V3.3.1.
}
\examples{

    data(Sig.W, package="BayesSigRefitting")

    # Display a subset of data 
    Sig.W[1:5, 1:5]

}
\keyword{data}
\seealso{ \code{\link{BayesSigRef}} }

