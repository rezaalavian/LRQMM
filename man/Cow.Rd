\name{Cow}
\alias{Cow}
\docType{data}
\title{Herd life Data of Iranian Holstein Cows}
\description{ Cow data include 100 cows with id records, father's record, mother' record, number of HYS, age of first calving and Herd life.}
\usage{data("Cow")}
\format{
  A data frame with 100 observations on the following 6 variables.
  \describe{
    \item{\code{REGNO}}{ The number form animal record as vector or column matrix}
    \item{\code{FREG}}{The number form father's animal record as vector or column matrix}
    \item{\code{MREG}}{The number form mother's animal record as vector or column matrix}
    \item{\code{HYS}}{a numeric vector levels of Herd, Year, Season}
    \item{\code{AGECAL}}{a numeric vector of age of first calving}
    \item{\code{HL}}{a numeric vector of Herd Life }
  }
}

\examples{
data(Cow)

}
\keyword{datasets}
