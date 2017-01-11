\name{fff}
\alias{fff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ F Distribution Family Function }
\description{
   Maximum likelihood estimation of the (2-parameter) F distribution.
}
\usage{
fff(link = "loge", idf1 = NULL, idf2 = NULL, nsimEIM = 100,
    imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function for both parameters.
  See \code{\link{Links}} for more choices.
  The default keeps the parameters positive.

  }
  \item{idf1, idf2}{
  Numeric and positive.
  Initial value for the parameters.
  The default is to choose each value internally.

  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{imethod}{
  Initialization method. Either the value 1 or 2.
  If both fail try setting values for \code{idf1} and \code{idf2}.


  }
% \item{zero}{
% An integer-valued vector specifying which
% linear/additive predictors are modelled as intercepts only.
% The value must be from the set \{1,2\}, corresponding
% respectively to \eqn{df1} and \eqn{df2}.
% By default all linear/additive predictors are modelled as
% a linear combination of the explanatory variables.
%
%
% }

}
\details{
  The F distribution is named after Fisher and has a density function
  that has two parameters, called \code{df1} and \code{df2} here.
  This function treats these degrees of freedom as \emph{positive reals}
  rather than integers.
  The mean of the distribution is
  \eqn{df2/(df2-2)} provided \eqn{df2>2},
  and its variance is
  \eqn{2 df2^2 (df1+df2-2)/(df1 (df2-2)^2 (df2-4))}{2*df2^2*(df1+df2-2)/
       (df1*(df2-2)^2*(df2-4))} provided \eqn{df2>4}.
  The estimated mean is returned as the fitted values.
  Although the F distribution can be defined to accommodate a
  non-centrality parameter \code{ncp}, it is assumed zero here.
  Actually it shouldn't be too difficult to handle any known
  \code{ncp}; something to do in the short future.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.



}
\author{ T. W. Yee }
\section{Warning}{
  Numerical problems will occur when the estimates of the parameters
  are too low or too high.


}

%\note{
%  This family function uses the BFGS quasi-Newton update formula for the
%  working weight matrices.  Consequently the estimated variance-covariance
%  matrix may be inaccurate or simply wrong! The standard errors must be
%  therefore treated with caution; these are computed in functions such
%  as \code{vcov()} and \code{summary()}.
%
%}
\seealso{
  \code{\link[stats:Fdist]{FDist}}.


}
\examples{
\dontrun{
fdata <- data.frame(x2 = runif(nn <- 2000))
fdata <- transform(fdata, df1 = exp(2+0.5*x2),
                          df2 = exp(2-0.5*x2))
fdata <- transform(fdata, y   = rf(nn, df1, df2))
fit <- vglm(y  ~ x2, fff, data = fdata, trace = TRUE)
coef(fit, matrix = TRUE)
}
}
\keyword{models}
\keyword{regression}



