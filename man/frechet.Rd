\name{frechet}
\alias{frechet}
%\alias{frechet2}
%\alias{frechet3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Frechet Distribution Family Function }
\description{
   Maximum likelihood estimation of the 2-parameter
   Frechet distribution.


%  and 3-parameter


}
\usage{
frechet(location = 0, lscale = "loge", lshape = logoff(offset = -2),
        iscale = NULL, ishape = NULL, nsimEIM = 250, zero = NULL)
}
%frechet3(anchor = NULL, ldifference = "loge", lscale = "loge",
%         lshape = "loglog",
%         ilocation = NULL, iscale = NULL, ishape = NULL,
%         zero = NULL, effpos = .Machine$double.eps^0.75)
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{location}{
  Numeric. Location parameter.
  It is called \eqn{a} below.


  }
  \item{lscale, lshape}{
  Link functions for the parameters;
  see \code{\link{Links}} for more choices.


  }
  \item{iscale, ishape, zero, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }


% \item{edifference}{   %
% Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }


% \item{anchor}{
% An ``anchor'' point for estimating the location parameter. This must
% be a value no greater than \code{min(y)} where \code{y} is the response.
% The location parameter is \eqn{A - D} where
% \eqn{A} is the anchor,
% \eqn{D} is the ``difference'' (default is to make this positive).
% The default value of \code{anchor} means \code{min(y)} is chosen.
%
% }
% \item{ldifference}{
% Parameter link function for the difference \eqn{D} between the anchor
% point and the location parameter estimate.
% The default keeps this difference positive so that numerical
% problems are less likely to occur.

% }
% \item{ilocation}{
% Optional initial value for the location parameter.
% A good choice can speed up the convergence rate markedly.
% A \code{NULL} means it is chosen internally.

% }
}
\details{
  The (3-parameter) Frechet distribution has a density function
  that can be written
  \deqn{f(y) = \frac{sb}{(y-a)^2} [b/(y-a)]^{s-1} \, \exp[-(b/(y-a))^s] }{%
        f(y) = ((s*b) / (y-a)^2) * exp[-(b/(y-a))^s] * [b/(y-a)]^(s-1)}
  for \eqn{y > a} and scale parameter \eqn{b > 0}.
  The positive shape parameter is \eqn{s}.
  The cumulative distribution function is
  \deqn{F(y) = \exp[-(b/(y-a))^s]. }{%
        F(y) =  exp[-(b/(y-a))^s].}
  The mean of \eqn{Y} is \eqn{a + b \Gamma(1-1/s)}{a + b*gamma(1-1/s)}
  for \eqn{s > 1} (these are returned as the fitted values).
  The variance of \eqn{Y} is
  \eqn{b^2 [ \Gamma(1-2/s) - \Gamma^2(1-1/s)]}{b^2 * [gamma(1 - 2/s) - gamma(1 - 1/s)^2]}
  for \eqn{s > 2}.


  Family \code{frechet} has \eqn{a} known, and
  \eqn{\log(b)}{log(b)} and
  \eqn{\log(s - 2)}{log(s - 2)} are the default linear/additive predictors.
  The working weights are estimated by simulated Fisher scoring.


% Note that the \code{\link{loglog}} link ensures \eqn{s > 1}.


% whereas \code{frechet3} estimates it.
% Estimating \eqn{a} well requires a lot of data and
% a good choice of \code{ilocation} will help speed up convergence.
% For \code{frechet3} the default linear/additive predictors are
% \eqn{\log(D)}{log(D)},
% It would be great if the first linear/additive predictor was a direct
% function of the location parameter, but this can run the risk that
% the estimate is out of range (i.e., greater than \code{min(y)}).



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications
      in Engineering and Science},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee }
\section{Warning}{
%  Convergence for \code{frechet3} can be very slow, especially if the
%  initial value for the location parameter is poor. Setting something
%  like \code{maxit = 200, trace = TRUE} is a good idea.


  Family function \code{frechet} may fail for low values of
  the shape parameter, e.g., near 2 or lower.


}

%\note{

% Family function \code{frechet3} uses
% the BFGS quasi-Newton update formula for the
% working weight matrices.  Consequently the estimated variance-covariance
% matrix may be inaccurate or simply wrong! The standard errors must be
% therefore treated with caution; these are computed in functions such
% as \code{vcov()} and \code{summary()}.


% If \code{fit} is a \code{frechet3} fit then \code{fit@extra$location}
% is the final estimate of the location parameter, and
% \code{fit@extra$LHSanchor} is the anchor point.


%}

\seealso{
  \code{\link{rfrechet}},
  \code{\link{gev}}.


}
\examples{
\dontrun{
set.seed(123)
fdata <- data.frame(y1 = rfrechet(nn <- 1000, shape = 2 + exp(1)))
with(fdata, hist(y1))
fit2 <- vglm(y1 ~ 1, frechet, data = fdata, trace = TRUE)
coef(fit2, matrix = TRUE)
Coef(fit2)
head(fitted(fit2))
with(fdata, mean(y1))
head(weights(fit2, type = "working"))
vcov(fit2)
}
}
\keyword{models}
\keyword{regression}

