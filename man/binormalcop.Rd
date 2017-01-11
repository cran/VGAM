\name{binormalcop}
\alias{binormalcop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Copula (Bivariate) Family Function }
\description{
  Estimate the correlation parameter of
  the (bivariate) Gaussian copula
  distribution by maximum likelihood estimation.

}
\usage{
binormalcop(lrho = "rhobit", irho = NULL, imethod = 1,
            parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
%  apply.parint = TRUE,
\arguments{
  \item{lrho, irho, imethod}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more link function choices.


  }
  \item{parallel, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.
  If \code{parallel = TRUE} then the constraint is applied to the
  intercept too.


  }
}
\details{
  The cumulative distribution function is
  \deqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = \Phi_2
             ( \Phi^{-1}(y_1), \Phi^{-1}(y_2); \rho ) }{%
        P(Y1 <= y1, Y2 <= y2) =
        Phi_2(\Phi^(-1)(y_1), \Phi^(-1)(y_2); \rho)}
  for \eqn{-1 < \rho < 1}{-1 < rho < 1},
  \eqn{\Phi_2}{Phi_2} is the cumulative distribution function
  of a standard bivariate normal
  (see \code{\link{pbinorm}}),
  and \eqn{\Phi}{Phi} is the cumulative distribution function
  of a standard univariate normal
  (see \code{\link[stats]{pnorm}}).


  The support of the function is the interior of the unit square;
  however, values of 0 and/or 1 are not allowed.
  The marginal distributions are the standard uniform distributions.
  When \eqn{\rho = 0}{rho=0} the random variables are
  independent.


  This \pkg{VGAM} family function can handle multiple responses,
  for example, a six-column matrix where the first 2 columns
  is the first out of three responses,
  the next 2 columns being the next response, etc.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{

Schepsmeier, U. and Stober, J. (2013)
Derivatives and Fisher information of bivariate copulas.
\emph{Statistical Papers}.


}
\author{ T. W. Yee }
\note{
  The response matrix must have a multiple of two-columns.
  Currently, the fitted
  value is a matrix with the same number of columns and values equal to 0.5.
  This is because each marginal distribution corresponds to a standard
  uniform distribution.


  This \pkg{VGAM} family function is fragile;
  each response must be in the interior of the unit square.
  Setting \code{crit = "coef"} is sometimes a good idea because
  inaccuracies in \code{\link{pbinorm}} might mean
  unnecessary half-stepping will occur near the solution.

}

\seealso{
  \code{\link{rbinormcop}},
  \code{\link[stats]{pnorm}},
  \code{\link{kendall.tau}}.


}
\examples{
nn <- 1000
ymat <- rbinormcop(n = nn, rho = rhobit(-0.9, inverse = TRUE))
bdata <- data.frame(y1 = ymat[, 1],
                    y2 = ymat[, 2],
                    y3 = ymat[, 1],
                    y4 = ymat[, 2],
                    x2 = runif(nn))

summary(bdata)
\dontrun{ plot(ymat, col = "blue") }
fit1 <- vglm(cbind(y1, y2, y3, y4) ~ 1,  # 2 responses, e.g., (y1,y2) is the first
             fam = binormalcop,
             crit = "coef",  # Sometimes a good idea
             data = bdata, trace = TRUE)

coef(fit1, matrix = TRUE)
Coef(fit1)
head(fitted(fit1))
summary(fit1)

# Another example; rho is a linear function of x2
bdata <- transform(bdata, rho = -0.5 + x2)
ymat <- rbinormcop(n = nn, rho = with(bdata, rho))
bdata <- transform(bdata, y5 = ymat[, 1],
                          y6 = ymat[, 2])
fit2 <- vgam(cbind(y5, y6) ~ s(x2), data = bdata,
             binormalcop(lrho = "identitylink"), trace = TRUE)
\dontrun{ plot(fit2, lcol = "blue", scol = "orange", se = TRUE, las = 1) }
}
\keyword{models}
\keyword{regression}

% for real \eqn{\rho}{rho} in (-1,1).
