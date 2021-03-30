STDE<-function (coef,Y=Y,E=E,SVD=SVD,tau = tau,n = n)
{# This function writed in "summary.rq" in "quantreg" package but in below used and changed for lrqmm function.
  y<-Y
  x<-SVD$u
  eps <- .Machine$double.eps^(1/2)
  p <- length(coef)
  h <- quantreg::bandwidth.rq(tau, n, hs = TRUE)
  if (tau + h > 1) stop("tau + h > 1:  error in STDE")
  if (tau - h < 0) stop("tau - h < 0:  error in STDE")
  bhi <- quantreg::rq.fit(SparseM::as.matrix.csr(x),Y,method="sfn",tau=tau+h)$coef
  bhi<-SVD$v%*%(solve(diag(SVD$d))%*%as.matrix(bhi))
  blo <- quantreg::rq.fit(SparseM::as.matrix.csr(x),Y,method="sfn",tau=tau-h)$coef
  blo<-SVD$v%*%(solve(diag(SVD$d))%*%as.matrix(blo))
  dyhat <- E %*% (bhi - blo)
  f <- pmax(1e-1, (2 * h)/(dyhat - eps))
  fII <- diag(p)
  fII<-backsolve(qr(sqrt(f) * x)$qr[drop = FALSE],fII)
  fII <- fII %*% t(fII)
  cov <- tau * (1 - tau) * fII %*% crossprod(x) %*% fII
  SE <- sqrt(diag(cov))
  coef <- array(coef, c(p, 2))
  dimnames(coef) <- list(dimnames(x)[[2]], c("Value", "Std. Error"))
  coef[, 2] <- abs(SVD$v%*%(solve(diag(SVD$d))%*%as.matrix(SE)))
  class(coef) <- "STDE"
  return(coef)
}
