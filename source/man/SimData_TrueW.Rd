% --- Source file: data.Rd ---
\name{SimData_TrueW}
\alias{SimData_TrueW}

\title{
    Profile Matrix of Signatures
}
\description{
    A 96 x 12 profile matrix of predefined ground truth signatures.
}
\details{ 
The 96 x 12 profile matrix of the predefined ground
truth signatures (SBS1, SBS2, SBS3, SBS5, SBS6, SBS8, SBS13, SBS17b, SBS18, SBS20,
SBS26, SBS30) used in the construction of \code{\link{SimData_V}}.
}
\examples{

    data(SimData_TrueW, package="BayesSigRefitting")

    # Display a subset of data 
    SimData_TrueW[1:5, 1:5]

}
\keyword{data}
\seealso{ \code{\link{SimData_V}}, \code{\link{SimData_TrueH}} }

