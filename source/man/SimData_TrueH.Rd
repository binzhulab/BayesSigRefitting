% --- Source file: data.Rd ---
\name{SimData_TrueH}
\alias{SimData_TrueH}

\title{
    Signature Activity Matrix
}
\description{
    A signature activity matrix of dimension 12 x 200
    used in the contruction of  \code{\link{SimData_V}}.
}
%\details{ }
\examples{

    data(SimData_TrueH, package="BayesSigRefitting")

    # Display a subset of data 
    SimData_TrueH[1:5, 1:5]

}
\keyword{data}
\seealso{ \code{\link{SimData_TrueW}}, \code{\link{SimData_V}} }

