\name{biclaytoncop}
\alias{biclaytoncop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Clayton Copula (Bivariate) Family Function }
\description{
  Estimate the correlation parameter of
  the (bivariate) Clayton copula
  distribution by maximum likelihood estimation.

}
\usage{
biclaytoncop(lapar = "loglink", iapar = NULL, imethod = 1,
             parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar, iapar, imethod}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more link function choices.


  }
  \item{parallel, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.
  If \code{parallel = TRUE} then the constraint is also applied
  to the intercept.


  }
}
\details{

  The cumulative distribution function is
\deqn{P(u_1, u_2;\alpha) = (u_1^{-\alpha} +
                            u_2^{-\alpha}-1)^{-1/\alpha}}{%
      P(u1,u2,alpha) = (u1^(-alpha) + u2^(-alpha)-1)^(-1/alpha)}
  for \eqn{0 \leq \alpha }{0 <= alpha}.
  Here, \eqn{\alpha}{alpha} is the association parameter.
  The support of the function is the interior of the unit square;
  however, values of 0 and/or 1 are not allowed (currently).
  The marginal distributions are the standard uniform distributions.
  When \eqn{\alpha = 0}{alpha=0} the random variables are independent.


  This \pkg{VGAM} family function can handle multiple responses,
  for example, a six-column matrix where the first 2 columns
  is the first out of three responses,
  the next 2 columns being the next response, etc.



}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{

%A Model for Association in Bivariate Survival Data.


Clayton, D. (1982).
A model for association in bivariate survival data.
\emph{Journal of the Royal Statistical Society,
      Series B, Methodological},
\bold{44}, 414--422.


Schepsmeier, U. and Stober, J. (2014).
Derivatives and Fisher information of bivariate copulas.
\emph{Statistical Papers}
\bold{55}, 525--542.



}
\author{ R. Feyter and T. W. Yee }
\note{
  The response matrix must have a multiple of two-columns.
  Currently, the fitted
  value is a matrix with the same number of columns
  and values equal to 0.5.
  This is because each marginal distribution
  corresponds to a standard uniform distribution.


  This \pkg{VGAM} family function is fragile;
  each response must be in the interior of the unit square.
%  Setting \code{crit = "coef"} is sometimes a good idea because
%  inaccuracies in \code{\link{pbinorm}} might mean
%  unnecessary half-stepping will occur near the solution.


}

\seealso{
  \code{\link{rbiclaytoncop}},
  \code{\link{dbiclaytoncop}},
  \code{\link{kendall.tau}}.


}
\examples{
ymat <- rbiclaytoncop(n = (nn <- 1000), apar = exp(2))
bdata <- data.frame(y1 = ymat[, 1], y2 = ymat[, 2],
                    y3 = ymat[, 1], y4 = ymat[, 2], x2 = runif(nn))
summary(bdata)
\dontrun{ plot(ymat, col = "blue") }
fit1 <-
  vglm(cbind(y1, y2, y3, y4) ~ 1,  # 2 responses, e.g., (y1,y2) is the 1st
       biclaytoncop, data = bdata,
       trace = TRUE, crit = "coef")  # Sometimes a good idea
coef(fit1, matrix = TRUE)
Coef(fit1)
head(fitted(fit1))
summary(fit1)

# Another example; apar is a function of x2
bdata <- transform(bdata, apar = exp(-0.5 + x2))
ymat <- rbiclaytoncop(n = nn, apar = with(bdata, apar))
bdata <- transform(bdata, y5 = ymat[, 1], y6 = ymat[, 2])
fit2 <- vgam(cbind(y5, y6) ~ s(x2), data = bdata,
             biclaytoncop(lapar = "loglink"), trace = TRUE)
\dontrun{plot(fit2, lcol = "blue", scol = "orange", se = TRUE) }}
\keyword{models}
\keyword{regression}

% for real \eqn{\alpha}{alpha} in (-1,1).
