\name{hzeta}
\alias{hzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Haight's Zeta Family Function }
\description{
  Estimating the parameter of Haight's zeta distribution


}
\usage{
hzeta(lshape = "loglog", ishape = NULL, nsimEIM = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape}{
  Parameter link function for the parameter,
  called \eqn{\alpha}{alpha} below.
  See \code{\link{Links}} for more choices.
  Here, a log-log link keeps the parameter greater than one, meaning
  the mean is finite.


  }
  \item{ishape,nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
    The probability function is
    \deqn{f(y) = (2y-1)^{(-\alpha)} - (2y+1)^{(-\alpha)},}{%
          f(y) = (2y-1)^(-alpha) - (2y+1)^(-alpha),}
    where the parameter \eqn{\alpha>0}{alpha>0}
    and \eqn{y=1,2,\ldots}{y=1,2,...}.
    The function \code{\link{dhzeta}} computes this probability function.
    The mean of \eqn{Y}, which is returned as fitted values, is
    \eqn{(1-2^{-\alpha}) \zeta(\alpha)}{(1-2^(-alpha))*zeta(alpha)}
    provided \eqn{\alpha > 1}{alpha > 1}, where \eqn{\zeta}{zeta} is
    Riemann's zeta function.
    The mean is a decreasing function of \eqn{\alpha}{alpha}.
    The mean is infinite if \eqn{\alpha \leq 1}{alpha <= 1}, and
    the variance is infinite if \eqn{\alpha \leq 2}{alpha <= 2}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

    Pages 533--4 of
    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
    \emph{Univariate Discrete Distributions},
    3rd edition,
    Hoboken, New Jersey: Wiley.

}
\author{ T. W. Yee }
%\note{
%}


\seealso{
  \code{\link{Hzeta}},
  \code{\link{zeta}},
  \code{\link{zetaff}},
  \code{\link{loglog}},
  \code{\link{simulate.vlm}}.


}
\examples{
shape <- exp(exp(-0.1))  # The parameter
hdata <- data.frame(y = rhzeta(n = 1000, shape))
fit <- vglm(y ~ 1, hzeta, data = hdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)  # Useful for intercept-only models; should be same as shape
c(with(hdata, mean(y)), head(fitted(fit), 1))
summary(fit)
}
\keyword{models}
\keyword{regression}




%# Generate some hzeta random variates
%set.seed(123)
%nn <- 400
%x <- 1:20
%shape <- 1.1  # The parameter
%probs <- dhzeta(x, shape)
%\dontrun{
%plot(x, probs, type="h", log="y")}
%cs <- cumsum(probs)
%tab <- table(cut(runif(nn), brea = c(0,cs,1)))
%index <- (1:length(tab))[tab>0]
%y <- rep(index, times=tab[index])



