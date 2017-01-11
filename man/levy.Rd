\name{levy}
\alias{levy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Levy Distribution Family Function }
\description{
Estimates the scale parameter of the Levy distribution
by maximum likelihood estimation.


}
\usage{
levy(location = 0, lscale = "loge", iscale = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{location}{
  Location parameter. Must have a known value.
  Called \eqn{a} below.

% otherwise it is estimated (the default).


  }
  \item{lscale}{
  Parameter link function for the (positive) scale parameter \eqn{b}.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale}{
  Initial value for the \eqn{b} parameter.
  By default, an initial value is chosen internally.


  }
}
\details{
  The Levy distribution is one of three stable distributions
  whose density function has a tractable form.
  The formula for the density is
 \deqn{f(y;b) = \sqrt{\frac{b}{2\pi}}
       \exp \left( \frac{-b}{2(y - a)}
            \right) / (y - a)^{3/2} }{%
  f(y;b) = sqrt(b / (2 pi))
       exp( -b / (2(y - a))) / (y - a)^{3/2} }
  where \eqn{a<y<\infty}{a<y<Inf} and \eqn{b>0}.
  Note that if \eqn{a} is very close to \code{min(y)}
  (where \code{y} is the response), then numerical problem will occur.
  The mean does not exist.
  The median is returned as the fitted values.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Nolan, J. P. (2005)
  \emph{Stable Distributions: Models for Heavy Tailed Data}.


% p.5



}
\author{ T. W. Yee }


%\note{
%  If \eqn{\delta}{delta} is given, then only one parameter is estimated
%  and the default is \eqn{\eta_1=\log(\gamma)}{eta1=log(gamma)}.
%  If \eqn{\delta}{delta} is not given, then \eqn{\eta_2=\delta}{eta2=delta}.
%
%
%}


\seealso{
  The Nolan article was at
  \code{http://academic2.american.edu/~jpnolan/stable/chap1.pdf}.



% \code{\link{dlevy}}.


}
\examples{
nn <- 1000; loc1 <- 0; loc2 <- 10
myscale <- 1  # log link ==> 0 is the answer
ldata <- data.frame(y1 = loc1 + myscale/rnorm(nn)^2,  # Levy(myscale, a)
                    y2 = rlevy(nn, loc = loc2, scale = exp(+2)))

# Cf. Table 1.1 of Nolan for Levy(1,0)
with(ldata, sum(y1 > 1) / length(y1))  # Should be 0.6827
with(ldata, sum(y1 > 2) / length(y1))  # Should be 0.5205

fit1 <- vglm(y1 ~ 1, levy(location = loc1), data = ldata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
summary(fit1)
head(weights(fit1, type = "work"))

fit2 <- vglm(y2 ~ 1, levy(location = loc2), data = ldata, trace = TRUE)
coef(fit2, matrix = TRUE)
Coef(fit2)
c(median = with(ldata, median(y2)), fitted.median = head(fitted(fit2), 1))
}
\keyword{models}
\keyword{regression}


%%\eqn{\delta + \gamma \Gamma(-0.5) / (2\sqrt{\pi})}{delta +
%%               gamma * gamma(-0.5) / (2*sqrt(pi))}
%%where \code{gamma} is a parameter but \code{gamma()} is the gamma function.

%%mygamma = exp(1)    # log link ==> 1 is the answer
%% alternative:
%%w = rgamma(n, shape=0.5)   # W ~ Gamma(0.5) distribution
%%mean(w)    # 0.5
%%mean(1/w)
%%y = delta + mygamma / (2 * w)    # This is Levy(mygamma, delta)
%%mean(y)
%%set.seed(123)

%%sum(y > 3) / length(y)  # Should be 0.4363
%%sum(y > 4) / length(y)  # Should be 0.3829
%%sum(y > 5) / length(y)  # Should be 0.3453

%fit <- vglm(y ~ 1, levy(idelta = delta, igamma = mygamma),
%           data = ldata, trace = TRUE)  # 2 parameters




