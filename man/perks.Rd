\name{perks}
\alias{perks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perks Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  Perks distribution.

}
\usage{
perks(lscale = "loge", lshape = "loge",
      iscale = NULL,   ishape = NULL,
      gscale = exp(-5:5), gshape = exp(-5:5),
      nsimEIM = 500, oim.mean = FALSE, zero = NULL, nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nowarning}{ Logical. Suppress a warning?
  Ignored for \pkg{VGAM} 0.9-7 and higher.


   }


  \item{lscale, lshape}{
  Parameter link functions applied to the
  shape parameter \code{shape},
  scale parameter \code{scale}.
  All parameters are treated as positive here
  See \code{\link{Links}} for more choices.


  }

% \item{eshape, escale}{
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
% }


  \item{iscale, ishape}{
  Optional initial values.
  A \code{NULL} means a value is computed internally.


  }
  \item{gscale, gshape}{
  See \code{\link{CommonVGAMffArguments}}.


  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}}.

  }
  \item{oim.mean}{
  To be currently ignored.

  }
}
\details{
The Perks distribution
has cumulative distribution function
\deqn{F(y; \alpha, \beta) =
1 -
\left\{
\frac{1 + \alpha}{1 + \alpha e^{\beta y}}
\right\}^{1 / \beta}
}{%
F(y; alpha, beta) = 1 - ((1 + \alpha)/(1 + alpha * e^(beta * y)))^(1 / beta)
}
which leads to a probability density function
\deqn{f(y; \alpha, \beta) =
\left[ 1 + \alpha \right]^{1 / \beta}
\alpha  e^{\beta y} / (1 + \alpha e^{\beta y})^{1 + 1 / \beta}
}{%
f(y; alpha, beta) = [ 1 + alpha]^(1 / \beta) * alpha * exp(beta * y) / (1 + alpha * exp(beta * y))^(1 + 1 / beta)
}
for \eqn{\alpha > 0}{alpha > 0},
\eqn{\beta > 0}{beta > 0},
\eqn{y > 0}.
Here, \eqn{\beta}{beta} is called the scale parameter \code{scale},
and \eqn{\alpha}{alpha} is called a shape parameter.
The moments for this distribution do
not appear to be available in closed form.



Simulated Fisher scoring is used and multiple responses are handled.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Perks, W. (1932)
On some experiments in the graduation of mortality statistics.
\emph{Journal of the Institute of Actuaries},
\bold{63}, 12--40.



Richards, S. J. (2012)
A handbook of parametric survival models for actuarial use.
\emph{Scandinavian Actuarial Journal}.
1--25.


}

\author{ T. W. Yee }
\section{Warning }{
A lot of care is needed because
this is a rather difficult distribution for parameter estimation.
If the self-starting initial values fail then try experimenting
with the initial value arguments, especially \code{iscale}.
Successful convergence depends on having very good initial values.
Also, monitor convergence by setting \code{trace = TRUE}.


}

\seealso{
  \code{\link{dperks}},
  \code{\link{simulate.vlm}}.


}

\examples{
\dontrun{ set.seed(123)
pdata <- data.frame(x2 = runif(nn <- 1000))  # x2 unused
pdata <- transform(pdata, eta1  = -1,
                          ceta1 =  1)
pdata <- transform(pdata, shape1 = exp(eta1),
                          scale1 = exp(ceta1))
pdata <- transform(pdata, y1 = rperks(nn, shape = shape1, scale = scale1))
fit1 <- vglm(y1 ~ 1, perks, data = pdata, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1)
}
}
\keyword{models}
\keyword{regression}


%# fit1 <- vglm(y1 ~ 1, perks, data = pdata, trace = TRUE)
%# fit2 <- vglm(y1 ~ 1, perks(imeth = 2), data = pdata, trace = TRUE)
% Argument \code{probs.y} is used only when \code{imethod = 2}.

