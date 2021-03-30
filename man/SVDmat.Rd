\name{SVDmat}
\alias{SVDmat}

\title{
Calculates SVD of Matrix in MATLAB
}
\description{
This function Calculates SVD of Matrix in MATLAB and produces the "economy size" decomposition.
}
\usage{
SVDmat(E,maxTriess = 3000,intervall = 30)
}
\arguments{
  \item{E}{a numeric matrix}
  \item{maxTriess}{The maximum number of times the connection is check for an answer from the MATLAB server before giving up. Default values is 3000 times.}
  \item{intervall}{The interval in seconds between each poll for an answer.  Default interval is 30 (second).}
}
\details{
This function use R working directory to biulte and read files. So there should be enough space.  All function's files remove after finishing calucation. This function is commonly used in big data.
}
\value{
  \item{d}{a vector containing the positive singular values}
  \item{u}{a matrix with the corresponding left singular vectors}
  \item{v}{a matrix with the corresponding right singular vectors}
}
\references{
[1]Alavian, S. R. (2019). Creating LRQMM package for predicting the breeding value of animals
by corrected mixed quantile regression (Unpublished master's thesis). Ferdowsi University Of Mashhad. Mashhad. Iran.[Persian].
}
\author{Sayyed Reza Alavian}
\note{
When this function stops abnormally (due an error or warning in MATLAB), you should close the MATLAB software window and disconnect the software. By performing this function again, the connection will be established. When more times need to the connection check for an answer from the MATLAB server before giving up, "maxTries" can be increase. When more times need to increase seconds between each poll for an answer, "interval" can be increase.
}


\examples{
M <- rbind(
  c(20, 10, 15,  0,  2),
  c(10,  5,  8,  1,  0),
  c( 0,  1,  2,  6,  3))
#before running this code, be sure for Matlab installation in your system.
# >SVDmat(M)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SVD }
\keyword{ MATLAB }
