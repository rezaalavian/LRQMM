\name{lrqmm}
\alias{lrqmm}

\title{
Fitting Linear Quantile Regression Mixed Models With Relationship Matrix
}
\description{
Fit a quantile regression mixed model involved Relationship Matrix using a sparse implementation of the
  Frisch-Newton interior-point algorithm.
}
\usage{
lrqmm(id, sire, dam, X, Y,cova=NULL , alpha = 0 , tau = 0.5 , Factor = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{ The number form animal record as column matrix}
  \item{sire}{The number form father's animal record as column matrix}
  \item{dam}{The number form mother's animal record as column matrix}
  \item{X}{fixed effect(s) as column matrix}
  \item{Y}{a response column matrix}
  \item{cova}{{covariate effect(s) column matrix}}
  \item{alpha}{a parameter for raite error's varince to variance of random effects dependent on statistical model (Animal model, Sire model, etc.)  }
  \item{tau}{desired quantile}
  \item{Factor}{type of fixed effect that "TRUE" as factor variable and "FALSE" as quantitative variable}
}
\details{
The function computes an estimate on the tau-th quantile effects of the linear mixed model. This is a sparse implementation of the Frisch-Newton algorithm for quantile regression described in Portnoy and Koenker (1997).

We used "GeneticsPed", "Matrix", "kinship2", "MCMCglmm", "rsvd", "SparseM" and "quantreg" packages in this function. befor using "lrqmm" function be sure from installation this packages.


"GeneticsPed" available in

<https://bioconductor.org/packages/release/bioc/src/contrib/GeneticsPed_1.46.0.tar.gz> or
orders in <http://bioconductor.org/packages/release/bioc/html/GeneticsPed.html>.

other packages are available in CRAN.
}
\value{
\item{Fixed effects}{estimate for fixed effect(s) from linear quantile regression mixed model with its standard error}
\item{cova effects}{estimate for covariate effect(s) from linear quantile regression mixed model with its standard error}
\item{Random effects}{estimate for random effect(s) from linear quantile regression mixed model with its standard error}
\item{residuals}{estimate for model residuals from linear quantile regression mixed model}
\item{Time_between_start_to_end}{execution time of linear quantile regression mixed model}
\item{summary}{reporting quantile for  effects estimation,  mean absolute error for fitted model,  variance of response variable, variance of pedigree's random.effect, variance of record's random.effect, number of observations, pedigree's length, fix effect lavels and random effect lavels}
}
\references{
[1]Alavian, S. R. (2019). Creating LRQMM package for predicting the breeding value of animals
by corrected mixed quantile regression (Unpublished master's thesis). Ferdowsi University Of Mashhad. Mashhad. Iran.[Persian].

[2]Koenker, R. and S. Portnoy (1997). The Gaussian Hare and the Laplacean Tortoise: Computability of Squared-error vs Absolute Error Estimators, (with discussion). Statistical Science, 12, 279-300.
<https://www.jstor.org/stable/2246216>

[3]Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press. ISBN: 0521608279.


[4]Mrode, R. A. (2005). Linear Models for the Prediction of Animal Breeding Values. 3rd edition.
CABI International. ISBN: 9781780643915.
}
\author{Sayyed Reza Alavian}



\examples{
#Start(not run)


data(Cow)
with(lrqmm(id=REGNO,sire=FREG,dam=MREG,X=HYS,Y=HL,cova=AGECAL,alpha=1,tau=0.5,Factor=TRUE)
,data=Cow)


#End(not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{quantile regression}% use one of  RShowDoc("KEYWORDS")
\keyword{mixed models }% __ONLY ONE__ keyword per line
\keyword{Frisch-Newton interior-point }
