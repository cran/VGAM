\name{paretoIV}
\alias{paretoIV}
\alias{paretoIII}
\alias{paretoII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pareto(IV/III/II) Distribution Family Functions }
\description{
  Estimates three of the parameters of the Pareto(IV) distribution
  by maximum likelihood estimation. Some special cases of this
  distribution are also handled.


}
\usage{
paretoIV(location = 0, lscale = "loge", linequality = "loge", lshape = "loge",
         iscale = 1, iinequality = 1, ishape = NULL, imethod = 1)
paretoIII(location = 0, lscale = "loge", linequality = "loge",
          iscale = NULL, iinequality = NULL)
paretoII(location = 0, lscale = "loge", lshape = "loge",
         iscale = NULL, ishape = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{location}{
  Location parameter, called \eqn{a} below.
  It is assumed known.


  }
  \item{lscale, linequality, lshape}{
  Parameter link functions for the
  scale parameter (called \eqn{b} below),
  inequality parameter (called \eqn{g} below), and
  shape parameter (called \eqn{s} below).
  See \code{\link{Links}} for more choices.
  A log link is the default for all because all these parameters are
  positive.


  }
  \item{iscale, iinequality, ishape}{
  Initial values for the parameters.
  A \code{NULL} value means that it is obtained internally.
  If convergence failure occurs, use these arguments to input
  some alternative initial values.


  }
  \item{imethod}{
  Method of initialization for the shape parameter.
  Currently only values 1 and 2 are available.
  Try the other value if convergence failure occurs.


  }
}
\details{
  The Pareto(IV) distribution, which is used in actuarial science,
  economics, finance and telecommunications,
  has a cumulative distribution function that can be written
  \deqn{F(y) = 1 - [1 + ((y-a)/b)^{1/g}]^{-s}}{%
        F(y) = 1 - [1 + ((y-a)/b)^(1/g)]^(-s)}
  for \eqn{y > a}, \eqn{b>0}, \eqn{g>0} and \eqn{s>0}.
  The \eqn{a}{a} is called the \emph{location} parameter,
  \eqn{b} the \emph{scale} parameter,
  \eqn{g} the \emph{inequality} parameter, and
  \eqn{s} the \emph{shape} parameter.


  The location parameter is assumed known otherwise the
  Pareto(IV) distribution will not be a regular family.
  This assumption is not too restrictive in modelling
  because in typical applications this parameter is known,
  e.g., in insurance and reinsurance it is pre-defined by
  a contract and can be represented as a deductible or a
  retention level.


  The inequality parameter is so-called because of its
  interpretation in the economics context. If we choose a
  unit shape parameter value and a zero location parameter
  value then the inequality parameter is the Gini index of
  inequality, provided \eqn{g \leq 1}{g<=1}.


  The fitted values are currently the median, e.g.,
  \code{\link{qparetoIV}} is used for \code{paretoIV()}.



% The fitted values are currently \code{NA} because I
% haven't worked out what the mean of \eqn{Y} is yet.


% The mean of \eqn{Y} is
% \eqn{\alpha k/(k-1)}{alpha*k/(k-1)} provided \eqn{k>1}.
% Its variance is
% \eqn{\alpha^2 k /((k-1)^2 (k-2))}{alpha^2 k /((k-1)^2 (k-2))}
% provided \eqn{k>2}.


% The maximum likelihood estimator for the location parameter is
% \code{min(y)}, i.e., the smallest response value.


  There are a number of special cases of the Pareto(IV) distribution.
  These include the Pareto(I), Pareto(II), Pareto(III), and Burr family
  of distributions.
  Denoting \eqn{PIV(a,b,g,s)} as the Pareto(IV) distribution,
  the Burr distribution \eqn{Burr(b,g,s)} is \eqn{PIV(a=0,b,1/g,s)},
  the Pareto(III) distribution \eqn{PIII(a,b,g)} is \eqn{PIV(a,b,g,s=1)},
  the Pareto(II) distribution \eqn{PII(a,b,s)} is \eqn{PIV(a,b,g=1,s)},
  and
  the Pareto(I) distribution \eqn{PI(b,s)} is \eqn{PIV(b,b,g=1,s)}.
  Thus the Burr distribution can be fitted using the
  \code{\link{negloge}} link
  function and using the default \code{location=0} argument.
  The Pareto(I) distribution can be fitted using \code{\link{paretoff}}
  but there is a slight change in notation: \eqn{s=k} and
  \eqn{b=\alpha}{b=alpha}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Johnson N. L., Kotz S., and Balakrishnan N. (1994)
\emph{Continuous Univariate Distributions, Volume 1},
2nd ed.
New York: Wiley.



Brazauskas, V. (2003)
Information matrix for Pareto(IV), Burr, and related distributions.
\emph{Comm. Statist. Theory and Methods}
\bold{32}, 315--325.


Arnold, B. C. (1983)
\emph{Pareto Distributions}.
Fairland, Maryland: International Cooperative Publishing House.


}
\author{ T. W. Yee }
\note{
  The \code{extra} slot of the fitted object has a component called
  \code{"location"} which stores the location parameter value(s).


}
\section{Warning }{
   The Pareto(IV) distribution is very general,
   for example, special cases include the Pareto(I), Pareto(II),
   Pareto(III), and Burr family of distributions.
   [Johnson et al. (1994) says on p.19 that fitting Type IV by ML is
   very difficult and rarely attempted].
   Consequently, reasonably good initial values are recommended,
   and convergence to a local solution may occur. For this
   reason setting \code{trace=TRUE} is a good idea for monitoring
   the convergence.  Large samples are ideally required to get
   reasonable results.




 }
\seealso{
  \code{\link{ParetoIV}},
  \code{\link{paretoff}},
  \code{\link{gpd}}.


}
\examples{
pdata <- data.frame(y = rparetoIV(2000, scale = exp(1),
                                  ineq = exp(-0.3), shape = exp(1)))
\dontrun{par(mfrow = c(2, 1))
with(pdata, hist(y)); with(pdata, hist(log(y))) }
fit <- vglm(y ~ 1, paretoIV, data = pdata, trace = TRUE)
head(fitted(fit))
summary(pdata)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

