\name{fff}
\alias{fff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ F Distribution Family Function }
\description{
   Maximum likelihood estimation of the (2-parameter) F distribution.
}
\usage{
fff(link="loge", earg=list(), idf1=NULL, idf2=NULL,
    method.init=1, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function for both parameters.
  See \code{\link{Links}} for more choices.
  The default keeps the parameters positive.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{idf1, idf2}{
  Numeric and positive. 
  Initial value for the parameters.
  The default is to choose each value internally.

  }
  \item{method.init}{
  Initialization method. Either the value 1 or 2.
  If both fail try setting values for \code{idf1} and \code{idf2}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The value must be from the set \{1,2\}, corresponding
  respectively to \eqn{df1} and \eqn{df2}.
  By default all linear/additive predictors are modelled as
  a linear combination of the explanatory variables.

  }
}
\details{
  The F distribution is named after Fisher and has a density function 
  that has two parameters, called \code{df1} and \code{df2} here.
  This function treats these degrees of freedom as \emph{positive reals}
  rather than integers. 
  The mean of the distribution is
  \eqn{df2/(df2-2)} provided \eqn{df2>2}, and its variance is
  \eqn{2 df2^2 (df1+df2-2)/(df1 (df2-2)^2 (df2-4))}{2*df2^2*(df1+df2-2)/
       (df1*(df2-2)^2*(df2-4))} provided \eqn{df2>4}.
  The estimated mean is returned as the fitted values.
  Although the F distribution can be defined to accommodate a
  non-centrality parameter \code{ncp}, it is assumed zero here.
}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

}
\author{ T. W. Yee }
\section{Warning}{
  Numerical problems will occur when the estimates of the parameters
  are too low. 
}

\note{ 
  This family function uses the BFGS quasi-Newton update formula for the
  working weight matrices.  Consequently the estimated variance-covariance
  matrix may be inaccurate or simply wrong! The standard errors must be
  therefore treated with caution; these are computed in functions such
  as \code{vcov()} and \code{summary()}.

}
\seealso{
  \code{\link[stats:Fdist]{FDist}}.
}
\examples{
x = runif(n <- 4000)
df1 = exp(2+0.5*x)
df2 = exp(2-0.5*x)
y = rf(n, df1, df2)
fit = vglm(y  ~ x, fff, trace=TRUE)
fit = vglm(y  ~ x, fff(link="logoff", earg=list(offset=0.5)), trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
vcov(fit)   # caution needed!
}
\keyword{models}
\keyword{regression}



%# Another simpler example
%set.seed(123)
%df1 = exp(2)
%df2 = exp(2)
%y = rf(n <- 1000, df1, df2)
%fit1 = vglm(y  ~ 1, fff, trace=TRUE)
%fit2 = vglm(y  ~ 1, fff(link="logoff", earg=list(offset=0.5)), trace=TRUE)
%vcov(fit1)   # caution needed!
%vcov(fit2)   # caution needed!

