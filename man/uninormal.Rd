\name{uninormal}
\alias{uninormal}
\alias{normal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate Normal Distribution }
\description{
  Maximum likelihood estimation of the two parameters of a univariate
  normal distribution.


}
\usage{
uninormal(lmean = "identitylink", lsd = "loge", lvar = "loge",
          var.arg = FALSE, imethod = 1, isd = NULL, parallel = FALSE,
          smallno = 1e-05, zero = "sd")
}
%- maybe also 'usage' for other objects documented here.
%         apply.parint = FALSE,
\arguments{
  \item{lmean, lsd, lvar}{
  Link functions applied to the mean and standard deviation/variance.
  See \code{\link{Links}} for more choices.
  Being positive quantities, a log link is the default for the
  standard deviation and variance (see \code{var.arg}).


  }


% \item{emean, esd, evar}{
% List. Extra argument for the links.
% See \code{earg} in \code{\link{Links}} for general information.
%       emean = list(), esd = list(), evar = list(),
% }


  \item{var.arg}{
  Logical.
  If \code{TRUE} then the second parameter is the variance and
  \code{lsd} and \code{esd} are ignored,
  else the standard deviation is used
  and \code{lvar} and \code{evar} are ignored.


  }
  \item{smallno}{
  Numeric, positive but close to 0.
  Used specifically for quasi-variances; if the link for the
  mean is \code{\link{explink}} then any non-positive value
  of \code{eta} is replaced by this quantity (hopefully,
  temporarily and only during early iterations).


  }
  \item{imethod, parallel, isd, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  If \code{lmean = loge} then try \code{imethod = 2}.
  If \code{parallel = TRUE} then the parallelism constraint
  is not applied to the intercept.


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

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}

\author{ T. W. Yee }
\section{Warning}{
  \code{uninormal()} is the new name;
  \code{normal1()} is old and will be decommissioned soon.


}

\note{
  Yet to do: allow an argument such as \code{eq.sd} that enables
  the standard devations to be the same.


}
\seealso{
    \code{\link{gaussianff}},
    \code{\link{posnormal}},
    \code{\link{mix2normal}},
    \code{\link{normal.vcm}},
    \code{\link{Qvar}},
    \code{\link{tobit}},
    \code{\link{cens.normal}},
    \code{\link{foldnormal}},
    \code{\link{skewnormal}},
    \code{\link{double.cens.normal}},
    \code{\link{SURff}},
    \code{\link{AR1}},
    \code{\link{huber2}},
    \code{\link{studentt}},
    \code{\link{binormal}},
    \code{\link[stats:Normal]{dnorm}},
    \code{\link{simulate.vlm}}.


}
\examples{
udata <- data.frame(x2 = rnorm(nn <- 200))
udata <- transform(udata,
                   y1  = rnorm(nn, m = 1 - 3*x2, sd = exp(1 + 0.2*x2)),
                   y2a = rnorm(nn, m = 1 + 2*x2, sd = exp(1 + 2.0*x2)^0.5),
                   y2b = rnorm(nn, m = 1 + 2*x2, sd = exp(1 + 2.0*x2)^0.5))
fit1 <- vglm(y1 ~ x2, uninormal(zero = NULL), data = udata, trace = TRUE)
coef(fit1, matrix = TRUE)
fit2 <- vglm(cbind(y2a, y2b) ~ x2, data = udata, trace = TRUE,
             uninormal(var = TRUE, parallel = TRUE ~ x2,
                       zero = NULL))
coef(fit2, matrix = TRUE)

# Generate data from N(mu = theta = 10, sigma = theta) and estimate theta.
theta <- 10
udata <- data.frame(y3 = rnorm(100, m = theta, sd = theta))
fit3a <- vglm(y3 ~ 1, uninormal(lsd = "identitylink"), data = udata,
             constraints = list("(Intercept)" = rbind(1, 1)))
fit3b <- vglm(y3 ~ 1, uninormal(lsd = "identitylink", parallel = TRUE ~ 1,
                                zero = NULL), data = udata)
coef(fit3a, matrix = TRUE)
coef(fit3b, matrix = TRUE)  # Same as fit3a
}
\keyword{models}
\keyword{regression}
