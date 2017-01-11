\name{lognormal}
\alias{lognormal}
%\alias{lognormal3}
%%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lognormal Distribution }
\description{
  Maximum likelihood estimation of the (univariate)
  lognormal distribution.

}
\usage{
lognormal(lmeanlog = "identitylink", lsdlog = "loge", zero = "sdlog")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmeanlog, lsdlog}{
  Parameter link functions applied to the mean and (positive)
  \eqn{\sigma}{sigma} (standard deviation) parameter.
  Both of these are on the log scale.
  See \code{\link{Links}} for more choices.


  }


% \item{emeanlog, esdlog}{
%         emeanlog = list(), esdlog = list(),
%          emeanlog = list(), esdlog = list(),
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
% }


  \item{zero}{
  Specifies which
  linear/additive predictor is modelled as intercept-only.
  For \code{lognormal()},
  the values can be from the set \{1,2\} which correspond to
  \code{mu}, \code{sigma}, respectively.
  See \code{\link{CommonVGAMffArguments}} for more information.


% For \code{lognormal3()},
% the values must be from the set \{1,2,3\} where 3 is for
% \eqn{\lambda}{\lambda}.



  }


% \item{powers.try}{
% Numerical vector. The initial \eqn{lambda} is chosen
% as the best value from \code{min(y) - 10^powers.try} where
% \code{y} is the response.

% }


% \item{delta}{
% Numerical vector. An alternative method for
% obtaining an initial \eqn{lambda}. Here, \code{delta = min(y)-lambda}.
% If given, this supersedes the \code{powers.try} argument.
% The value must be positive.

% }
}
\details{
  A random variable \eqn{Y} has a 2-parameter lognormal distribution
  if \eqn{\log(Y)}{log(Y)}
  is distributed \eqn{N(\mu, \sigma^2)}{N(mu, sigma^2)}.
  The expected value of \eqn{Y}, which is
  \deqn{E(Y) = \exp(\mu + 0.5 \sigma^2)}{%
        E(Y) = exp(mu + 0.5 sigma^2)}
  and not \eqn{\mu}{mu}, make up the fitted values.
  The variance of \eqn{Y} is
  \deqn{Var(Y) = [\exp(\sigma^2) -1] \exp(2\mu + \sigma^2).}{%
        Var(Y) = [exp(sigma^2) -1] * exp(2 mu + sigma^2).}



% A random variable \eqn{Y} has a 3-parameter lognormal distribution
% if \eqn{\log(Y-\lambda)}{log(Y-lambda)}
% is distributed \eqn{N(\mu, \sigma^2)}{N(mu, sigma^2)}. Here,
% \eqn{\lambda < Y}{lambda < Y}.
% The expected value of \eqn{Y}, which is
% \deqn{E(Y) = \lambda + \exp(\mu + 0.5 \sigma^2)}{%
%       E(Y) =  lambda + exp(mu + 0.5 sigma^2)}
% and not \eqn{\mu}{mu}, make up the fitted values.



% \code{lognormal()} and \code{lognormal3()} fit the 2- and 3-parameter
% lognormal distribution respectively. Clearly, if the location
% parameter \eqn{\lambda=0}{lambda=0} then both distributions coincide.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee }
%\note{
%  The more commonly used 2-parameter lognormal distribution is the
%  3-parameter lognormal distribution with \eqn{\lambda}{lambda} equal
%  to zero---see \code{\link{lognormal3}}.
%
%
%}


%\section{Warning}{
%  Regularity conditions are not satisfied for the 3-parameter case:
%  results may be erroneous.
%  May withdraw it in later versions.
%
%
%}


\seealso{
  \code{\link[stats]{rlnorm}},
  \code{\link{uninormal}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.


% \code{\link{lognormal3}},


}

\examples{
ldata2 <- data.frame(x2 = runif(nn <- 1000))
ldata2 <- transform(ldata2, y1 = rlnorm(nn, mean = 1 + 2 * x2, sd = exp(-1)),
                            y2 = rlnorm(nn, mean = 1, sd = exp(-1 + x2)))
fit1 <- vglm(y1 ~ x2, lognormal(zero = 2), data = ldata2, trace = TRUE)
fit2 <- vglm(y2 ~ x2, lognormal(zero = 1), data = ldata2, trace = TRUE)
coef(fit1, matrix = TRUE)
coef(fit2, matrix = TRUE)
}
\keyword{models}
\keyword{regression}



%lognormal3(lmeanlog = "identitylink", lsdlog = "loge",
%           powers.try = (-3):3, delta = NULL, zero = 2)

%lambda <- 4
%ldata3 <- data.frame(y3 = lambda + rlnorm(1000, m = 1.5, sd = exp(-0.8)))
%fit3 <- vglm(y3 ~ 1, lognormal3, data = ldata3, trace = TRUE, crit = "c")
%coef(fit3, matrix = TRUE)
%summary(fit3)


%ldata <- data.frame(y1 = rlnorm(nn <- 1000, meanlog = 1.5, sdlog = exp(-0.8)))
%fit1 <- vglm(y1 ~ 1, lognormal, data = ldata, trace = TRUE, crit = "c")
%coef(fit1, matrix = TRUE)
%Coef(fit1)


