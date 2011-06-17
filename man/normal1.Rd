\name{normal1}
\alias{normal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate Normal Distribution }
\description{
  Maximum likelihood estimation of the two parameters of a univariate
  normal distribution.

}
\usage{
normal1(lmean = "identity", lsd = "loge",
        emean = list(), esd = list(), imethod = 1, zero = -2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean, lsd}{
  Link functions applied to the mean and standard deviation.
  See \code{\link{Links}} for more choices.
  Being a positive quantity, a log link is the default for the
  standard deviation.

  }
  \item{emean, esd}{
  List. Extra argument for the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{imethod, zero}{ 
  See \code{\link{CommonVGAMffArguments}} for more information.
  If \code{lmean = loge} then try \code{imethod = 2}.

  }

}
\details{
  This fits a linear model (LM) as the first linear/additive predictor.
  So, by default, this is just the mean.
  By default,
  the log of the standard deviation is the second linear/additive predictor.
  The Fisher information matrix is diagonal.
  This \pkg{VGAM} family function can handle multiple responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

}

\author{ T. W. Yee }
\note{
    Yet to do: allow an argument such as \code{sameSD} that enables the
    standard devations to be the same.
    And a \code{parallel} argument.

}
\seealso{
    \code{\link{gaussianff}},
    \code{\link{posnormal1}},
    \code{\link{mix2normal1}},
    \code{\link{tobit}},
    \code{\link{cennormal1}},
    \code{\link{fnormal1}},
    \code{\link{skewnormal1}},
    \code{\link{dcennormal1}},
    \code{\link{huber}},
    \code{\link{studentt}},
    \code{\link{binormal}},
    \code{\link[stats:Normal]{dnorm}}.

}
\examples{
ndata <- data.frame(x2 = rnorm(nn <- 200))
ndata <- transform(ndata, y = rnorm(nn, mean = 1-3*x2, sd = exp(1+0.2*x2)))
fit <- vglm(y ~ x2, normal1(zero = NULL), ndata, trace = TRUE)
coef(fit, matrix = TRUE)

# Generate data from N(mu = theta = 10, sigma = theta) and estimate theta.
theta <- 10
ndata <- data.frame(y = rnorm(100, m = theta, sd = theta))
fit <- vglm(y ~ 1, normal1(lsd = "identity"), ndata,
            constraints = list("(Intercept)" = rbind(1, 1)))
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
