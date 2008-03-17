\name{yulesimon}
\alias{yulesimon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Yule-Simon Family Function }
\description{
  Estimating the parameter of the Yule-Simon distribution.

}
\usage{
yulesimon(link="loge", earg=list(), irho=NULL, nsimEIM=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, earg}{
  Link function and extra argument for the \eqn{\rho}{rho} parameter.
  See \code{\link{Links}} for more choices and for general information.

  }
  \item{irho}{
  Optional initial value for the (positive) parameter. 
  See \code{\link{CommonVGAMffArguments}} for more information.
  The default is to obtain an initial value internally. Use this argument
  if the default fails.

  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
    The probability function is
    \deqn{f(y;\rho) = rho*beta(y,rho+1),}{%
          f(y;rho) = rho*beta(y,rho+1),}
    where the parameter \eqn{\rho>0}{rho>0}
    and \eqn{y=1,2,\ldots}{y=1,2,...}.
    The function \code{\link{dyules}} computes this probability function.
    The mean of \eqn{Y}, which is returned as fitted values, is
    \eqn{\rho/(\rho-1)}{rho/(rho-1)}
    provided \eqn{\rho > 1}{rho > 1}.
    The variance of \eqn{Y} is
    \eqn{\rho^2/((\rho-1)^2 (\rho-2))}{rho^2/((rho-1)^2 (rho-2))}
    provided \eqn{\rho > 2}{rho > 2}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{ 

    Simon, H. A. (1955)
    On a class of skew distribution functions.
    \emph{Biometrika},
    \bold{42},
    425--440.

}
\author{ T. W. Yee }
%\note{ 
%}

\seealso{
    \code{\link{ryules}}.

}
\examples{
x = runif(n <- 1000)
y = ryules(n, rho=exp(1.5-x))
table(y)
fit = vglm(y ~ x, yulesimon, trace=TRUE)
coef(fit, matrix=TRUE)
summary(fit)
}
\keyword{models}
\keyword{regression}


%# Generate some yulesimon random variates
%set.seed(123)
%nn = 400
%x = 1:20
%alpha = 1.1  # The parameter
%probs = dyulesimon(x, alpha)
%\dontrun{
%plot(x, probs, type="h", log="y")}
%cs = cumsum(probs)
%tab = table(cut(runif(nn), brea = c(0,cs,1)))
%index = (1:length(tab))[tab>0]
%y = rep(index, times=tab[index]) 



