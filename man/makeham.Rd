\name{makeham}
\alias{makeham}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Makeham Distribution Family Function }
\description{
  Maximum likelihood estimation of the 3-parameter
  Makeham distribution.

}
\usage{
makeham(lscale = "loge", lshape = "loge", lepsilon = "loge",
        iscale = NULL,   ishape = NULL,   iepsilon = NULL,
        gscale = exp(-5:5),gshape = exp(-5:5), gepsilon = exp(-4:1),
        nsimEIM = 500, oim.mean = TRUE, zero = NULL, nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{nowarning}{ Logical. Suppress a warning?
  Ignored for \pkg{VGAM} 0.9-7 and higher.


   }


  \item{lshape, lscale, lepsilon}{
  Parameter link functions applied to the
  shape parameter \code{shape},
  scale parameter \code{scale}, and
  other parameter \code{epsilon}.
  All parameters are treated as positive here
  (cf. \code{\link{dmakeham}} allows \code{epsilon = 0}, etc.).
  See \code{\link{Links}} for more choices.


  }

% \item{eshape, escale, eepsilon}{
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{ishape, iscale, iepsilon}{
  Optional initial values.
  A \code{NULL} means a value is computed internally.
  A value must be given for \code{iepsilon} currently, and this
  is a sensitive parameter!


  }
  \item{gshape, gscale, gepsilon}{
  See \code{\link{CommonVGAMffArguments}}.
  }



  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}}.
  Argument \code{probs.y} is used only when \code{imethod = 2}.


  }
  \item{oim.mean}{
  To be currently ignored.

  }
}
\details{
The Makeham distribution, which adds another parameter
to the Gompertz distribution,
has cumulative distribution function
\deqn{F(y; \alpha, \beta, \varepsilon) =
1 - \exp
\left\{
-y \varepsilon + \frac {\alpha}{\beta}
\left[ 1 - e^{\beta y} \right]
\right\}
}{%
F(y; alpha, beta, epsilon) = 1 - exp(-y * epsilon + (alpha / beta) * [1 - e^(beta * y)])
}
which leads to a probability density function
\deqn{f(y; \alpha, \beta, \varepsilon) =
\left[
\varepsilon + \alpha e^{\beta y} \right]
\;
\exp
\left\{
-y \varepsilon + \frac {\alpha}{\beta}
\left[ 1 - e^{\beta y} \right]
\right\},
}{%
f(y; alpha, beta, epsilon) = (epsilon + alpha * e^(beta y) ) * exp(-y * epsilon + (alpha / beta) * [1 - e^(beta * y)])
}
for \eqn{\alpha > 0}{alpha > 0},
\eqn{\beta > 0}{beta > 0},
\eqn{\varepsilon \geq 0}{epsilon >= 0},
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
%\references{
%
%}

\author{ T. W. Yee }
\section{Warning }{
A lot of care is needed because
this is a rather difficult distribution for parameter estimation,
especially when the shape parameter is large relative to the
scale parameter.
If the self-starting initial values fail then try experimenting
with the initial value arguments, especially \code{iepsilon}.
Successful convergence depends on having very good initial values.
More improvements could be made here.
Also, monitor convergence by setting \code{trace = TRUE}.


A trick is to fit a \code{\link{gompertz}} distribution and use
it for initial values; see below.
However, this family function is currently numerically fraught.


}

\seealso{
  \code{\link{dmakeham}},
  \code{\link{gompertz}},
  \code{\link{simulate.vlm}}.


}

\examples{
\dontrun{ set.seed(123)
mdata <- data.frame(x2 = runif(nn <- 1000))
mdata <- transform(mdata, eta1  = -1,
                          ceta1 =  1,
                          eeta1 = -2)
mdata <- transform(mdata, shape1 = exp(eta1),
                          scale1 = exp(ceta1),
                          epsil1 = exp(eeta1))
mdata <- transform(mdata,
         y1 = rmakeham(nn, shape = shape1, scale = scale1, eps = epsil1))

# A trick is to fit a Gompertz distribution first
fit0 <- vglm(y1 ~ 1, gompertz, data = mdata, trace = TRUE)
fit1 <- vglm(y1 ~ 1, makeham, data = mdata,
             etastart = cbind(predict(fit0), log(0.1)), trace = TRUE)

coef(fit1, matrix = TRUE)
summary(fit1)
}
}
\keyword{models}
\keyword{regression}


%# fit1 <- vglm(y1 ~ 1, makeham, data = mdata, trace = TRUE)
%# fit2 <- vglm(y1 ~ 1, makeham(imeth = 2), data = mdata, trace = TRUE)

