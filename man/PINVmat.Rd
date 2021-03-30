\name{PINVmat}
\alias{PINVmat}

\title{
 Generalized Inverse of A Big Matrix Whit MATLAB
}
\description{
Calucated invesre of the  generilzed big matrix with MATLAB
}
\usage{
PINVmat(x, maxTriess = 3000, intervall = 30)
}

\arguments{
  \item{x}{	a numeric matrix }
  \item{maxTriess}{The maximum number of times the connection is check for an answer from the MATLAB server before giving up. Default values is 3000 times.}
  \item{intervall}{The interval in seconds between each poll for an answer.  Default interval is 30 (second).}
}
\details{
see pinv function in MATLAB.
}
\value{
a inverse generalized matrix
}
\references{
[1]Alavian, S. R. (2019). Creating LRQMM package for predicting the breeding value of animals
by corrected mixed quantile regression (Unpublished master's thesis). Ferdowsi University Of Mashhad. Mashhad. Iran.[Persian].
}
\author{Sayyed Reza Alavian }

\examples{
M <- rbind(
  c(20, 10, 15,  0,  2),
  c(10,  5,  8,  1,  0),
  c( 0,  1,  2,  6,  3))
#before running this code, be sure for Matlab installation in your system.
# >PINVmat(M)
}

\keyword{big Matrix }
\keyword{ Generalized Inverse }
\keyword{MATLAB}
