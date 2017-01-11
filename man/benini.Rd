\name{benini1}
\alias{benini1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Benini Distribution Family Function }
\description{
  Estimating the 1-parameter Benini distribution by maximum
  likelihood estimation.

}
\usage{
benini1(y0 = stop("argument 'y0' must be specified"), lshape = "loge",
        ishape = NULL, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y0}{
  Positive scale parameter.

  }
  \item{lshape}{
  Parameter link function and extra argument of the parameter \eqn{b},
  which is the shape parameter.
  See \code{\link{Links}} for more choices.
  A log link is the default because \eqn{b} is positive.

  }
  \item{ishape}{
  Optional initial value for the shape parameter.
  The default is to compute the value internally.

  }
  \item{imethod, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
  The Benini distribution
  has a probability density function that can be written
  \deqn{f(y) = 2 s \exp(-s[(\log(y/y_0))^2]) \log(y/y_0) / y }{%
        f(y) = 2*s*exp(-s * [(log(y/y0))^2]) * log(y/y0) / y}
  for \eqn{0 < y_0 < y}{0 < y0 < y}, and shape \eqn{s > 0}.
  The cumulative distribution function for \eqn{Y} is
  \deqn{F(y) = 1 - \exp(-s[(\log(y/y_0))^2]).}{%
        F(y) = 1 - exp(-s * [(log(y / y0))^2]). }
  Here, Newton-Raphson and Fisher scoring coincide.
  The median of \eqn{Y} is now returned as the fitted values, by default.
  This \pkg{VGAM} family function can handle a multiple
  responses, which is inputted as a matrix.


  On fitting, the \code{extra} slot has a component called
  \code{y0} which contains the value of the \code{y0}
  argument.


}
%\section{Warning}{
%
%
% The median of \eqn{Y}, which are returned as the fitted values,
% may be incorrect.
%
%
%}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


% Section 7.1, pp.235--8


}
\author{ T. W. Yee }
\note{
  Yet to do: the 2-parameter Benini distribution estimates another
  shape parameter \eqn{a}{a} too.
  Hence, the code may change in the future.


}
\seealso{
  \code{\link{Benini}}.


}
\examples{
y0 <- 1; nn <- 3000
bdata <- data.frame(y  = rbenini(nn, y0 = y0, shape = exp(2)))
fit <- vglm(y ~ 1, benini1(y0 = y0), data = bdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
fit@extra$y0
c(head(fitted(fit), 1), with(bdata, median(y)))  # Should be equal
}
\keyword{models}
\keyword{regression}

