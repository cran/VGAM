\name{borel.tanner}
\alias{borel.tanner}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Borel-Tanner Distribution Family Function}
\description{
  Estimates the parameter of a Borel-Tanner distribution
  by maximum likelihood estimation.


}
\usage{
borel.tanner(Qsize = 1, link = "logit", imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Qsize}{
  A positive integer.
  It is called \eqn{Q} below and is the initial queue size.
  The default value \eqn{Q = 1} corresponds to the Borel distribution.


  }
  \item{link}{
  Link function for the parameter;
  see \code{\link{Links}} for more choices and for general information.


  }
  \item{imethod}{
  See \code{\link{CommonVGAMffArguments}}.
  Valid values are 1, 2, 3 or 4.


  }
}
\details{
  The Borel-Tanner distribution (Tanner, 1953) describes the
  distribution of the total number of customers served before a queue
  vanishes given a single queue with random arrival times of
  customers (at a constant rate \eqn{r} per unit time,
  and each customer taking a constant time \eqn{b} to be served).
  Initially the queue has \eqn{Q} people and the first one starts to
  be served.
  The two parameters appear in the density only in the form of the product
  \eqn{rb}, therefore we use \eqn{a=rb}, say, to denote the single
  parameter to be estimated.
  The density function is
  \deqn{f(y;a) =
  \frac{ Q }{(y-Q)!} y^{y-Q-1} a^{y-Q}  \exp(-ay)
  }{%
  f(y;a) =
  (Q / (y-Q)!) * y^(y-Q-1) * a^(y-Q) * exp(-ay)}
  where \eqn{y=Q,Q+1,Q+2,\ldots}{y=Q,Q+1,Q+2,...}.
  The case \eqn{Q=1} corresponds to the \emph{Borel} distribution
  (Borel, 1942).
  For the \eqn{Q=1} case it is necessary for \eqn{0 < a < 1} for the
  distribution to be proper.
  The Borel distribution is a basic Lagrangian distribution of the
  first kind.
  The Borel-Tanner distribution is an \eqn{Q}-fold convolution of the
  Borel distribution.


  The mean is \eqn{Q/(1-a)} (returned as the fitted values) and the
  variance is \eqn{Q a / (1-a)^3}{Q*a/(1-a)^3}.
  The distribution has a very long tail unless \eqn{a} is small.
  Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
Tanner, J. C. (1953)
A problem of interference between two queues.
\emph{Biometrika}, \bold{40}, 58--69.


Borel, E. (1942)
Sur l'emploi du theoreme de Bernoulli pour faciliter le
calcul d'une infinite de coefficients.
Application au probleme de l'attente a un guichet.
\emph{Comptes Rendus, Academie des Sciences, Paris, Series A},
\bold{214}, 452--456.


    Page 328 of
    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
    \emph{Univariate Discrete Distributions},
    3rd edition,
    Hoboken, New Jersey: Wiley.


Consul, P. C. and Famoye, F. (2006)
\emph{Lagrangian Probability Distributions},
Boston, MA, USA: Birkhauser.


}
\author{ T. W. Yee }
%\note{
%
%}

\seealso{
  \code{\link{rbort}},
  \code{\link{poissonff}},
  \code{\link{felix}}.


}
\examples{
bdata <- data.frame(y = rbort(n <- 200))
fit <- vglm(y ~ 1, borel.tanner, data = bdata, trace = TRUE, crit = "c")
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

