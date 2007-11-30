\name{fnormal1}
\alias{fnormal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Folded Normal Distribution Family Function }
\description{
  Fits a (generalized) folded (univariate) normal distribution.
}
\usage{
fnormal1(lmean="identity", lsd="loge", emean=list(), esd=list(),
         imean=NULL, isd=NULL, a1=1, a2=1, nsimEIM=500,
         method.init=1, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  See \code{\link{CommonVGAMffArguments}} for general information about
  many of these arguments.

  \item{lmean, lsd}{
  Link functions for the mean and standard
  deviation parameters of the usual univariate normal distribution.
  They are \eqn{\mu}{mu} and \eqn{\sigma}{sigma} respectively.
  See \code{\link{Links}} for more choices.

  }
  \item{emean, esd}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{imean, isd}{
  Optional initial values for \eqn{\mu}{mu} and \eqn{\sigma}{sigma}.
  A \code{NULL} means a value is computed internally.

  }
  \item{a1, a2}{
  Positive weights, called \eqn{a_1}{a1} and \eqn{a_2}{a2} below.
  Each must be of length 1.

  }
  \item{nsimEIM, method.init, zero}{
  See \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
  If a random variable has an ordinary univariate normal distribution then
  the absolute value of that random variable has an ordinary \emph{folded
  normal distribution}. That is, the sign has not been recorded; only
  the magnitude has been measured.

  More generally, suppose \eqn{X} is normal with mean \code{mean} and
  standard deviation \code{sd}.
  Let \eqn{Y=\max(a_1 X, -a_2 X)}{Y=max(a1*X, -a2*X)}
  where \eqn{a_1}{a1} and \eqn{a_2}{a2} are positive weights.
  This means that \eqn{Y = a_1 X}{Y = a1*X} for \eqn{X > 0}, and
  \eqn{Y = a_2 X}{Y = a2*X} for \eqn{X < 0}.
  Then \eqn{Y} is said to have a \emph{generalized folded normal distribution}.
  The ordinary folded normal distribution corresponds to the
  special case \eqn{a_1 = a_2 = 1}{a1 = a2 = 1}.

  The probability density function of the ordinary folded normal distribution
  can be written
  \code{dnorm(y, mean, sd) + dnorm(y, -mean, sd)} for
  \eqn{y \ge 0}.
  By default, \code{mean} and \code{log(sd)} are the linear/additive
  predictors.
  Having \code{mean=0} and \code{sd=1} results in the
  \emph{half-normal} distribution.
  The mean of an ordinary folded normal distribution is
  \deqn{E(Y) = \sigma \sqrt{2/\pi} \exp(-\mu^2/(2\sigma^2)) +
               \mu [1-2\Phi(-\mu/\sigma)] }{%
        E(Y) = sigma*sqrt(2/pi)*exp(-mu^2/(2*sigma^2)) +
               mu*[1-2*Phi(-mu/sigma)] }
  and these are returned as the fitted values.
  Here, \eqn{\Phi()}{Phi} is the cumulative distribution function of a
  standard normal (\code{\link[stats]{pnorm}}).

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
  Lin, P. C. (2005)
  Application of the generalized folded-normal distribution to the process
  capability measures.
  \emph{International Journal of Advanced Manufacturing Technology},
  \bold{26}, 825--830.

}
\author{ Thomas W. Yee }
\note{
  The response variable for this family function is the same as
  \code{\link{normal1}} except positive values are required.
  Reasonably good initial values are needed.
  Fisher scoring using simulation is implemented.

}

\section{Warning }{
  Under- or over-flow may occur if the data is ill-conditioned.
  It is recommended that several different initial values be used
  to help avoid local solutions.

}
\seealso{ 
    \code{\link{rfnorm}},
    \code{\link{normal1}},
    \code{\link[stats]{dnorm}},
    \code{\link{skewnormal1}}.
}

\examples{
m =  2; SD = exp(1)
y = rfnorm(n <- 1000, m=m, sd=SD)
\dontrun{hist(y, prob=TRUE, main=paste("fnormal1(m=",m,", sd=",round(SD,2),")"))}
fit = vglm(y ~ 1, fam=fnormal1, trace=TRUE)
coef(fit, mat=TRUE)
(Cfit = Coef(fit))
mygrid = seq(min(y), max(y), len=200) # Add the fit to the histogram
\dontrun{lines(mygrid, dfnorm(mygrid, Cfit[1], Cfit[2]), col="red")}
}
\keyword{models}
\keyword{regression}
