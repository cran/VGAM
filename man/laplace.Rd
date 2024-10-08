\name{laplace}
\alias{laplace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Laplace Regression Family Function }
\description{
   Maximum likelihood estimation of
   the 2-parameter classical Laplace distribution.

}
\usage{
laplace(llocation = "identitylink", lscale = "loglink",
  ilocation = NULL, iscale = NULL, imethod = 1, zero = "scale")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{ Character.
  Parameter link functions for location parameter \eqn{a} and
  scale parameter \eqn{b}.
  See \code{\link[VGAM]{Links}} for more choices.


  }
  \item{ilocation, iscale}{
  Optional initial values.
  If given, it must be numeric and values are recycled to the
  appropriate length.
  The default is to choose the value internally.


  }
  \item{imethod}{
  Initialization method.
  Either the value 1 or 2.


  }
  \item{zero}{
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.


  }
}
\details{
  The Laplace distribution is often known as the
  \emph{double-exponential} distribution and,
  for modelling, has heavier tail than the normal distribution.
  The Laplace density function is
  \deqn{f(y) = \frac{1}{2b} \exp \left( - \frac{|y-a|}{b}
                    \right) }{%
        f(y) =  (1/(2b)) exp( -|y-a|/b ) }
  where \eqn{-\infty<y<\infty}{-Inf<y<Inf},
  \eqn{-\infty<a<\infty}{-Inf<a<Inf} and
  \eqn{b>0}.
  Its mean is \eqn{a} and its variance is \eqn{2b^2}.
  This parameterization is called the \emph{classical Laplace
  distribution} by Kotz et al. (2001), and the density is symmetric
  about \eqn{a}.


  For \code{y ~ 1} (where \code{y} is the response)
  the maximum likelihood estimate (MLE) for the location
  parameter is the sample median, and the MLE for \eqn{b} is
  \code{mean(abs(y-location))} (replace location by its MLE
  if unknown).



}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link[VGAM]{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link[VGAM]{vglm}}
  and \code{\link[VGAM]{vgam}}.



}
\references{

Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001).
\emph{The Laplace distribution and generalizations:
a revisit with applications to communications,
economics, engineering, and finance},
Boston: Birkhauser.


}
\author{ T. W. Yee }
\section{Warning}{
  This family function has not been fully tested.
  The MLE regularity conditions do \emph{not} hold for this
  distribution, therefore misleading inferences may result, e.g.,
  in the \code{summary} and \code{vcov} of the object.  Hence this
  family function might be withdrawn from \pkg{VGAM} in the future.



}
\note{
  This family function uses Fisher scoring.
  Convergence may be slow for non-intercept-only models;
  half-stepping is frequently required.


}

\seealso{
\code{\link{rlaplace}},
\code{\link[VGAMdata]{alaplace2}}
(which differs slightly from this parameterization),
\code{\link[VGAM]{exponential}},
\code{\link[stats]{median}}.


}

\examples{
ldata <- data.frame(y = rlaplace(nn <- 100, 2, scale = exp(1)))
fit <- vglm(y  ~ 1, laplace, ldata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
with(ldata, median(y))

ldata <- data.frame(x = runif(nn <- 1001))
ldata <- transform(ldata, y = rlaplace(nn, 2, scale = exp(-1 + 1*x)))
coef(vglm(y ~ x, laplace(iloc = 0.2, imethod = 2, zero = 1), ldata,
          trace = TRUE), matrix = TRUE)
}
\keyword{models}
\keyword{regression}

