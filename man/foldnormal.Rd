\name{foldnormal}
\alias{foldnormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Folded Normal Distribution Family Function }
\description{
  Fits a (generalized) folded (univariate) normal distribution.

}
\usage{
foldnormal(lmean = "identitylink", lsd = "loge", imean = NULL, isd = NULL,
           a1 = 1, a2 = 1, nsimEIM = 500, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{lmean, lsd}{
  Link functions for the mean and standard
  deviation parameters of the usual univariate normal distribution.
  They are \eqn{\mu}{mu} and \eqn{\sigma}{sigma} respectively.
  See \code{\link{Links}} for more choices.

  }


% \item{emean, esd}{
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
%        emean=list(), esd=list(),
% }


  \item{imean, isd}{
  Optional initial values for \eqn{\mu}{mu} and \eqn{\sigma}{sigma}.
  A \code{NULL} means a value is computed internally.
  See \code{\link{CommonVGAMffArguments}}.


  }
  \item{a1, a2}{
  Positive weights, called \eqn{a_1}{a1} and \eqn{a_2}{a2} below.
  Each must be of length 1.


  }
  \item{nsimEIM, imethod, zero}{
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
  standard normal (\code{\link[stats:Normal]{pnorm}}).


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


  Johnson, N. L. (1962)
  The folded normal distribution:
  accuracy of estimation by maximum likelihood.
  \emph{Technometrics},
  \bold{4}, 249--256.




}
\author{ Thomas W. Yee }
\note{
  The response variable for this family function is the same as
  \code{\link{uninormal}} except positive values are required.
  Reasonably good initial values are needed.
  Fisher scoring using simulation is implemented.


  See \code{\link{CommonVGAMffArguments}} for general information about
  many of these arguments.


  Yet to do: implement the results of Johnson (1962) which gives
  expressions for the EIM, albeit, under a different parameterization.
  Also, one element of the EIM appears to require numerical integration.



}

\section{Warning }{
  Under- or over-flow may occur if the data is ill-conditioned.
  It is recommended that several different initial values be used
  to help avoid local solutions.


}
\seealso{
    \code{\link{rfoldnorm}},
    \code{\link{uninormal}},
    \code{\link[stats:Normal]{dnorm}},
    \code{\link{skewnormal}}.


}

\examples{
\dontrun{ m <-  2; SD <- exp(1)
fdata <- data.frame(y = rfoldnorm(n <- 1000, m = m, sd = SD))
hist(with(fdata, y), prob = TRUE, main = paste("foldnormal(m = ", m,
     ", sd = ", round(SD, 2), ")"))
fit <- vglm(y ~ 1, foldnormal, data = fdata, trace = TRUE)
coef(fit, matrix = TRUE)
(Cfit <- Coef(fit))
# Add the fit to the histogram:
mygrid <- with(fdata, seq(min(y), max(y), len = 200))
lines(mygrid, dfoldnorm(mygrid, Cfit[1], Cfit[2]), col = "orange")
}
}
\keyword{models}
\keyword{regression}
