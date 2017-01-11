\name{bistudentt}
\alias{bistudentt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Student-t Family Function }
\description{
  Estimate the degrees of freedom and correlation parameters of
  the (bivariate) Student-t
  distribution by maximum likelihood estimation.

}
\usage{
bistudentt(ldf = "loglog", lrho = "rhobit",
           idf = NULL, irho = NULL, imethod = 1,
           parallel = FALSE, zero = "rho")
}
%- maybe also 'usage' for other objects documented here.
%apply.parint = TRUE,
\arguments{
  \item{ldf, lrho, idf, irho, imethod}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more link function choices.


  }
  \item{parallel, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The density function is
  \deqn{f(y_1, y_2; \nu, \rho) =
        \frac{1}{2\pi\sqrt{1-\rho^2}}
        (1 + y_1^2 + y_2^2 - 2\rho y_1 y_2) / (\nu (1-\rho^2))^{(\nu+2)/2} }{%
        f(y1, y2; nu, rho) =
        (1/(2*pi*sqrt(1-\rho^2))) *
        (1 + y1^2 + y_2^2 - 2*rho*y1*y2) / (nu*(1-rho^2))^((\nu+2)/2) }
  for \eqn{-1 < \rho < 1}{-1 < rho < 1},
  and real \eqn{y_1}{y1} and \eqn{y_2}{y2}.


% The support of the function is the interior of the unit square;
% however, values of 0 and/or 1 are not allowed.
% The marginal distributions are the standard uniform distributions.
% When \eqn{\rho = 0}{rho=0} the random variables are
% independent.


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
\author{ T. W. Yee,
with help from Thibault Vatter.


}
\note{
  The response matrix must have a multiple of two-columns.
  Currently, the fitted
  value is a matrix with the same number of columns and values equal to 0.0.


}
\section{Warning }{


  The working weight matrices have not been fully checked.


}

\seealso{
  \code{\link{dbistudentt}},
  \code{\link{binormal}},
  \code{\link[stats]{pt}}.


}
\examples{
nn <- 1000
mydof <- loglog(1, inverse = TRUE)
ymat <- cbind(rt(nn, df = mydof), rt(nn, df = mydof))
bdata <- data.frame(y1 = ymat[, 1], y2 = ymat[, 2],
                    y3 = ymat[, 1], y4 = ymat[, 2], x2 = runif(nn))
summary(bdata)
\dontrun{ plot(ymat, col = "blue") }
fit1 <- vglm(cbind(y1, y2, y3, y4) ~ 1,  # 2 responses, e.g., (y1,y2) is the 1st
             fam = bistudentt,  # crit = "coef",  # Sometimes a good idea
             data = bdata, trace = TRUE)

coef(fit1, matrix = TRUE)
Coef(fit1)
head(fitted(fit1))
summary(fit1)
}
\keyword{models}
\keyword{regression}

%


