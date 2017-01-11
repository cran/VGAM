\name{logF}
\alias{logF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Natural Exponential Family Generalized Hyperbolic Secant Distribution
        Family Function

}
\description{
   Maximum likelihood estimation of
   the 1-parameter log F distribution.


}
\usage{
 logF(lshape1 = "loge", lshape2 = "loge",
      ishape1 = NULL, ishape2 = 1, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1, lshape2}{
% Character.
  Parameter link functions for
  the shape parameters.
  Called \eqn{\alpha}{alpha} and \eqn{\beta}{beta} respectively.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape1, ishape2}{
  Optional initial values for the shape parameters.
  If given, it must be numeric and values are recycled to the
  appropriate length.
  The default is to choose the value internally.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{imethod}{
  Initialization method.
  Either the value 1, 2, or \ldots.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The density for this distribution is
 \deqn{f(y; \alpha, \beta) = \exp(\alpha y) / [B(\alpha,\beta)
                             (1 + e^y)^{\alpha + \beta}] }{%
       f(y;  alpha,  beta) =  exp(\alpha y) / [B(\alpha,\beta) *
                             (1 + exp(y))^(\alpha + \beta)] }
  where \eqn{y} is real,
  \eqn{\alpha > 0},
  \eqn{\beta > 0},
  \eqn{B(., .)} is the beta function \code{\link[base:Special]{beta}}.




}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

  Jones, M. C. (2008).
  On a class of distributions with simple exponential tails.
  \emph{Statistica Sinica},
  \bold{18}(3), 1101--1110.
% Section 3.2.


}
\author{ Thomas W. Yee }
%\section{Warning}{
%
%}

%\note{
%}

\seealso{
  \code{\link{dlogF}},
  \code{\link{logff}}.


% \code{\link{simulate.vlm}}.


}

\examples{
nn <- 1000
ldata <- data.frame(y1 = rnorm(nn, m = +1, sd = exp(2)),  # Not proper data
                    x2 = rnorm(nn, m = -1, sd = exp(2)),
                    y2 = rnorm(nn, m = -1, sd = exp(2)))  # Not proper data
fit1 <- vglm(y1 ~ 1 , logF, data = ldata, trace = TRUE)
fit2 <- vglm(y2 ~ x2, logF, data = ldata, trace = TRUE)
coef(fit2, matrix = TRUE)
summary(fit2)
vcov(fit2)

head(fitted(fit1))
with(ldata, mean(y1))
max(abs(head(fitted(fit1)) - with(ldata, mean(y1))))
}
\keyword{models}
\keyword{regression}

