% --- Source file: data.Rd ---
\name{SimData_V}
\alias{SimData_V}

\title{
    Mutation Catalog Matrix V
}
\description{
    A 96 x 200 matrix of simulated mutation counts.
}
\details{ 
The matrix contains 200 tumors for 96 single base substitution (SBS) types. The
simulated mutation counts mimic the real data in the BRCA560 breast cancer dataset. Each
element in the matrix is generated from the Poisson distribution with the mean corresponding
to each element of \eqn{W_{t} H_{t}}{W_tH_t}, 
where \eqn{W_{t}}{W_t} is the 96 x 12 profile matrix of the predefined ground
truth signatures (SBS1, SBS2, SBS3, SBS5, SBS6, SBS8, SBS13, SBS17b, SBS18, SBS20,
SBS26, SBS30), and \eqn{H_{t}}{H_t} is the 12 x 200 signature activity matrix. 
The matrices \eqn{W_{t}}{W_t} and \eqn{H_{t}}{H_t} are also included 
in the \code{data} folder as
\code{\link{SimData_TrueW}} and \code{\link{SimData_TrueH}}.
}
\examples{

    data(SimData_V, package="BayesSigRefitting")

    # Display a subset of data 
    SimData_V[1:5, 1:5]

}
\keyword{data}
\seealso{ \code{\link{SimData_TrueW}}, \code{\link{SimData_TrueH}} }

