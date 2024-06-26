\name{binormal}
\alias{binormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Normal Distribution Family Function }
\description{
  Maximum likelihood estimation of the five parameters of a
  bivariate normal distribution.

}
\usage{
binormal(lmean1 = "identitylink", lmean2 = "identitylink",
         lsd1   = "loglink",     lsd2   = "loglink",
         lrho   = "rhobitlink",
         imean1 = NULL,       imean2 = NULL,
         isd1   = NULL,       isd2   = NULL,
         irho   = NULL,       imethod = 1,
         eq.mean = FALSE,     eq.sd = FALSE,
         zero = c("sd", "rho"), rho.arg = NA)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean1, lmean2, lsd1, lsd2, lrho}{
  Link functions applied to the means, standard deviations and
  \code{rho} parameters.
  See \code{\link{Links}} for more choices.
  Being positive quantities, a log link is the default for the
  standard deviations.


  }
  \item{imean1, imean2, isd1, isd2, irho, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{eq.mean, eq.sd}{
  Logical or formula.
  Constrains the means or the standard deviations to be equal.




% 20150530; FALSE now; they work separately:
% Only one of these arguments may be assigned a value.



  }

% 20210923; rho.arg implemented.
  \item{rho.arg}{
    If \eqn{\rho}{rho} is known then this argument may
    be assigned the (scalar) value lying in \eqn{(-1, 1)}.
    The default is to estimate that parameter so that \eqn{M=5}.
    If known, then other arguments such as \code{lrho}
    and \code{irho} are ignored,
    and \code{"rho"} is removed from \code{zero}.


%   and arguments such as \code{eq.sd} might not work,
  }


}
\details{
  For the bivariate normal distribution,
  this fits a linear model (LM) to the means, and
  by default,
  the other parameters are intercept-only.
  The response should be a two-column matrix.
  The correlation parameter is \code{rho},
  which lies between \eqn{-1} and \eqn{1}
  (thus the \code{\link{rhobitlink}} link is
  a reasonable choice).
  The fitted means are returned as the fitted
  values, which is in
  the form of a two-column matrix.
  Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\section{Warning}{
  This function may be renamed to
  \code{normal2()} or something like that at
  a later date.


}


%\references{
%
%}

\author{ T. W. Yee }
\note{
  If both equal means and equal standard
  deviations are desired then use something
  like
\code{constraints = list("(Intercept)" = }
\code{matrix(c(1,1,0,0,0, 0,0,1,1,0 ,0,0,0,0,1), 5, 3))}
  and maybe
  \code{zero = NULL}
  etc.


}

\seealso{
  \code{\link{uninormal}},
  \code{\link{trinormal}},
  \code{\link{pbinorm}},
  \code{\link{bistudentt}},
  \code{\link{rhobitlink}}.


% \code{\link{gaussianff}},
% \code{\link{pnorm2}},


}
\examples{
set.seed(123); nn <- 1000
bdata <- data.frame(x2 = runif(nn), x3 = runif(nn))
bdata <- transform(bdata, y1 = rnorm(nn, 1 + 2 * x2),
                          y2 = rnorm(nn, 3 + 4 * x2))
fit1 <- vglm(cbind(y1, y2) ~ x2,
             binormal(eq.sd = TRUE), bdata, trace = TRUE)
coef(fit1, matrix = TRUE)
constraints(fit1)
summary(fit1)

# Estimated P(Y1 <= y1, Y2 <= y2) under the fitted model
var1  <- loglink(2 * predict(fit1)[, "loglink(sd1)"], inv = TRUE)
var2  <- loglink(2 * predict(fit1)[, "loglink(sd2)"], inv = TRUE)
cov12 <- rhobitlink(predict(fit1)[, "rhobitlink(rho)"], inv = TRUE)
head(with(bdata, pbinorm(y1, y2,
                         mean1 = predict(fit1)[, "mean1"],
                         mean2 = predict(fit1)[, "mean2"],
                         var1 = var1, var2 = var2, cov12 = cov12)))
}
\keyword{models}
\keyword{regression}
