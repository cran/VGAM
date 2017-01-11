\name{triangle}
\alias{triangle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Triangle Distribution Family Function }
\description{
  Estimating the parameter of the triangle distribution by maximum
  likelihood estimation.

}
\usage{
triangle(lower = 0, upper = 1,
         link = extlogit(min = 0, max = 1), itheta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower, upper}{lower and upper limits of the distribution.
     Must be finite.
     Called \eqn{A} and \eqn{B} respectively below.


   }

  \item{link}{
  Parameter link function applied to the parameter \eqn{\theta}{theta},
  which lies in \eqn{(A,B)}.
  See \code{\link{Links}} for more choices.
  The default constrains the estimate to lie in the interval.


  }
  \item{itheta}{
  Optional initial value for the parameter.
  The default is to compute the value internally.


  }
}
\details{
  The triangle distribution
  has a probability density function that consists of two lines
  joined at \eqn{\theta}{theta}, which is the location of the mode.
  The lines intersect the
  \eqn{y = 0} axis at \eqn{A} and \eqn{B}.
  Here, Fisher scoring is used.



  On fitting, the \code{extra} slot has components called \code{lower}
  and \code{upper} which contains the values of the above arguments
  (recycled to the right length).
  The fitted values are the mean of the distribution, which is
  \eqn{(A + B + \theta)/3}{(A + B + theta)/3}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{


Kotz, S. and van Dorp, J. R. (2004)
Beyond Beta: Other Continuous Families of Distributions
with Bounded Support and Applications.
Chapter 1.
World Scientific: Singapore.


Nguyen, H. D. and McLachlan, G. J. (2016)
Maximum likelihood estimation of triangular and polygon distributions.
\emph{Computational Statistics and Data Analysis},
\bold{102}, 23--36.




}
\author{ T. W. Yee }
\section{Warning}{
  The MLE regularity conditions do not hold for this
  distribution
  (e.g., the first derivative evaluated at the mode
  does not exist because it is not continuous)
  so that misleading inferences may result, e.g.,
  in the \code{summary} and \code{vcov} of the object.
  Additionally, convergence to the MLE often appears to fail.


}

\note{
  The response must contain values in \eqn{(A, B)}.
  For most data sets (especially small ones) it is very common for
  half-stepping to occur.


% 20130603
  Arguments \code{lower} and \code{upper} and \code{link} must match.
  For example, setting
  \code{lower = 0.2} and \code{upper = 4} and
  \code{link = extlogit(min = 0.2, max = 4.1)} will result in an error.
  Ideally \code{link = extlogit(min = lower, max = upper)}
  ought to work but it does not (yet)!
  Minimal error checking is done for this deficiency.



}
\seealso{
  \code{\link{Triangle}},
  \code{\link{Topple}},
  \code{\link{simulate.vlm}}.


}
\examples{
# Example 1
tdata <- data.frame(y = rtriangle(n <- 3000, theta = 3/4))
fit <- vglm(y ~ 1, triangle(link = "identitylink"), data = tdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
head(fit@extra$lower)
head(fitted(fit))
with(tdata, mean(y))

# Example 2; Kotz and van Dorp (2004), p.14
rdata <- data.frame(y = c(0.1, 0.25, 0.3, 0.4, 0.45, 0.6, 0.75, 0.8))
fit <- vglm(y ~ 1, triangle(link = "identitylink"), data = rdata, trace = TRUE,
            crit = "coef", maxit = 1000)
Coef(fit)  # The MLE is the 3rd order statistic, which is 0.3.
fit <- vglm(y ~ 1, triangle(link = "identitylink"), data = rdata, trace = TRUE,
            crit = "coef", maxit = 1001)
Coef(fit)  # The MLE is the 3rd order statistic, which is 0.3.
}
\keyword{models}
\keyword{regression}



% 20130603: yettodo: fix up so ideally
%        link = extlogit(min = lower, max = upper), itheta = NULL)
% works.
