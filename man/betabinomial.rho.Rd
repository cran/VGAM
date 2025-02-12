\name{betabinomial.rho}
\alias{betabinomial.rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Beta-binomial Distribution Family Function
(with known rho)}
\description{
  Fits a beta-binomial distribution by maximum
  likelihood estimation, where the
  correlation coefficient rho is inputted.
  The parameter estimated is the mean.


}
\usage{
betabinomial.rho(lmu = "logitlink", imethod = 1, ishrinkage = 0.95)
}
%- maybe also 'usage' for other objects documented here.
%            ishrinkage = 0.95, nsimEIM = NULL, zero = 2
\arguments{
  \item{lmu, imethod, ishrinkage}{
    Same as \code{\link{betabinomial}}.
  
  }
}
\details{
  This family function conducts a logistic-like regression
  where the correlation parameter \eqn{\rho} of
  a betabinomial distribution is
  inputted by the user.
  The family function is somewhat like a simplified
  \code{\link{betabinomial}}.
  The argument \code{form2}
  (see \code{\link{vglm}})
  is used to input the
  \eqn{\rho} values, which must lie in \eqn{[0, 1]}.

  
  The default model has \eqn{\eta_1 = logit(\mu)}{eta1 =logit(mu)}.


}
\value{
  Same as \code{\link{betabinomial}}.

}
\author{ T. W. Yee }
%\note{
%}
%\section{Warning }{


%}
\seealso{
  \code{\link{betabinomial}},
  \code{\link{extbetabinomial}},
  \code{\link{betabinomialff}},
  \code{\link{Betabinom}},
  \code{\link{vglm}},
  \code{\link{binomialff}},
  \code{\link{betaff}}.


}
\examples{
\dontrun{
# Example 1
nn <- 10000; NN <- 100
bdata <- data.frame(N = NN, x2 = rnorm(nn),
                    x3 = rnorm(nn))
bdata <-
    transform(bdata,
              mu1  = logitlink(-0.5, inverse = TRUE),
              rho1 = logitlink( 0.5, inverse = TRUE),
              mu2  = logitlink(-0.5 + x2, inverse = TRUE),
              rho2 = logitlink(-0.5 + x3, inverse = TRUE))
bdata <- transform(bdata,
         y1 = rbetabinom(nn, size = N, prob = mu1, rho = rho1),
         y2 = rbetabinom(nn, size = N, prob = mu2, rho = rho2))
fit1 <- vglm(cbind(y1, N - y1) ~ 1, betabinomial.rho,
         form2 = ~ rho1, crit = "c", bdata, trace = TRUE)
coef(fit1, matrix = TRUE)
head(fit1@extra$rho)
max(abs(fitted(fit1) - with(bdata, mu1)))  # Should be 0

# Example 2
fit2 <- vglm(cbind(y2, N - y2) ~ x2, form2 = ~ rho2,
         betabinomial.rho, crit = "c",
         bdata, trace = TRUE)
coef(fit2, matrix = TRUE)
max(abs(fit2@extra$rho - with(bdata, rho2)))  # Should be 0
max(abs(fitted(fit2) - with(bdata, mu2)))  # Should be 0
}}
\keyword{models}
\keyword{regression}
