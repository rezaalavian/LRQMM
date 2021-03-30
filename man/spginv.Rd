\name{spginv}
\alias{spginv}

\title{
 Generalized Inverse of a Sparse Matrix
}
\description{
Calucated invesre of the  generilzed sparse matrix with sparsesvd function in sparcesvd package and ginv function in MASS package.
}
\usage{
spginv(x)
}

\arguments{
  \item{x}{	a sparse real matrix in Matrix package format}
}
\details{
see sparsesvd function in sparcesvd package and ginv function in MASS package.
}
\value{
a inverse generalized sparse matrix
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
M <- Matrix::Matrix(M, sparse=TRUE)
spginv (M)
}

\keyword{Sparse Matrix }
\keyword{ Generalized Inverse }
