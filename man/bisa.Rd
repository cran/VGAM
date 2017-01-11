\name{bisa}
\alias{bisa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Birnbaum-Saunders Distribution Family Function }
\description{
  Estimates the shape and scale parameters of the
  Birnbaum-Saunders distribution by maximum likelihood estimation.

}
\usage{
bisa(lscale = "loge", lshape = "loge", iscale = 1,
     ishape = NULL, imethod = 1, zero = "shape", nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nowarning}{ Logical. Suppress a warning?
  Ignored for \pkg{VGAM} 0.9-7 and higher.


  }


  \item{lscale, lshape}{
  Parameter link functions applied to the shape and scale parameters
  (\eqn{a} and \eqn{b} below).
  See \code{\link{Links}} for more choices.
  A log link is the default for both because they are positive.


  }
  \item{iscale, ishape}{
  Initial values for \eqn{a} and \eqn{b}.
  A \code{NULL} means an initial value is chosen internally using
  \code{imethod}.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} or \code{3} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for
  \code{ishape} and/or \code{iscale}.


  }
  \item{zero}{
  Specifies which linear/additive predictor is modelled as intercept-only.
  If used, choose one value from the set \{1,2\}.
  See \code{\link{CommonVGAMffArguments}} for more details.



% The default is none of them.



  }
}
\details{
  The (two-parameter) Birnbaum-Saunders distribution
  has a cumulative distribution function that can be written as
  \deqn{F(y;a,b) = \Phi[ \xi(y/b)/a] }{%
        F(y;a,k) = pnorm[xi(y/b)/a] }
  where \eqn{\Phi(\cdot)}{pnorm()} is the
  cumulative distribution function of a standard normal
  (see \code{\link[stats:Normal]{pnorm}}),
  \eqn{\xi(t) = \sqrt{t} - 1 / \sqrt{t}}{xi(t) = t^(0.5) - t^(-0.5)},
  \eqn{y > 0},
  \eqn{a>0} is the shape parameter,
  \eqn{b>0} is the scale parameter.
  The mean of \eqn{Y} (which is the fitted value) is
  \eqn{b(1 + a^2/2)}{b*(1 + a*a/2)}.
  and the variance is
  \eqn{a^2 b^2 (1 + \frac{5}{4}a^2)}{a^2 b^2 (1 + (5/4)*a^2)}.
  By default, \eqn{\eta_1 = \log(a)}{eta1 = log(a)} and
  \eqn{\eta_2 = \log(b)}{eta2 = log(b)} for this family function.


  Note that \eqn{a} and \eqn{b} are orthogonal,
  i.e., the Fisher information matrix is diagonal.
  This family function implements Fisher scoring, and
  it is unnecessary to compute any integrals numerically.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Lemonte, A. J. and Cribari-Neto, F. and Vasconcellos, K. L. P. (2007)
Improved statistical inference for the two-parameter
Birnbaum-Saunders distribution.
\emph{Computational Statistics \& Data Analysis}, \bold{51}, 4656--4681.


Birnbaum, Z. W. and Saunders, S. C. (1969)
A new family of life distributions.
\emph{Journal of Applied Probability}, \bold{6}, 319--327.


Birnbaum, Z. W. and Saunders, S. C. (1969)
Estimation for a family of life distributions with applications to fatigue.
\emph{Journal of Applied Probability}, \bold{6}, 328--347.


Engelhardt, M. and Bain, L. J. and Wright, F. T. (1981)
Inferences on the parameters of the Birnbaum-Saunders fatigue
life distribution based on maximum likelihood estimation.
\emph{Technometrics}, \bold{23}, 251--256.


Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 2,
New York: Wiley.


}
\author{ T. W. Yee }
%\note{
%
%}
%\section{Warning }{
%}

\seealso{
  \code{\link{pbisa}},
  \code{\link{inv.gaussianff}},
  \code{\link{CommonVGAMffArguments}}.


}
\examples{
bdata1 <- data.frame(x2 = runif(nn <- 1000))
bdata1 <- transform(bdata1, shape = exp(-0.5 + x2), scale = exp(1.5))
bdata1 <- transform(bdata1, y = rbisa(nn, scale, shape))
fit1 <- vglm(y ~ x2, bisa(zero = 1), data = bdata1, trace = TRUE)
coef(fit1, matrix = TRUE)

\dontrun{
bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
bdata2 <- transform(bdata2, y = rbisa(nn, scale, shape))
fit <- vglm(y ~ 1, bisa, data = bdata2, trace = TRUE)
with(bdata2, hist(y, prob = TRUE, ylim = c(0, 0.5), col = "lightblue"))
coef(fit, matrix = TRUE)
with(bdata2, mean(y))
head(fitted(fit))
x <- with(bdata2, seq(0, max(y), len = 200))
lines(dbisa(x, Coef(fit)[1], Coef(fit)[2]) ~ x, data = bdata2,
      col = "orange", lwd = 2) }
}
\keyword{models}
\keyword{regression}

