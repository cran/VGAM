% 20131230; adapted from simulate.Rd from R 3.0.2


% \newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}


\name{simulate.vlm}
\title{Simulate Responses for VGLMs and VGAMs}
\description{
  Simulate one or more responses from the distribution
  corresponding to a fitted model object.

}
\usage{
\method{simulate}{vlm}(object, nsim = 1, seed = NULL, \dots)
}
\alias{simulate.vlm}
\arguments{
  \item{object}{an object representing a fitted model.
  Usually an object of class
  \code{\link{vglm-class}}
  or
  \code{\link{vgam-class}}.



  }
  \item{nsim, seed}{
  Same as \code{\link[stats]{simulate}}.



  }
% \item{seed}{an object specifying if and how the random number
%   generator should be initialized (\sQuote{seeded}).\cr
%   For the "lm" method, either \code{NULL} or an integer that will be
%   used in a call to \code{set.seed} before simulating the response
%   vectors.  If set, the value is saved as the \code{"seed"} attribute
%   of the returned value.  The default, \code{NULL} will not change the
%   random generator state, and return \code{\link{.Random.seed}} as the
%   \code{"seed"} attribute, see \sQuote{Value}.
% }


  \item{\dots}{additional optional arguments.}



}
\value{
  Similar to \code{\link[stats]{simulate}}.
  Note that many \pkg{VGAM} family functions can handle multiple responses.
  This can result in a longer data frame with more rows
  (\code{nsim} multiplied by \code{n} rather than the
   ordinary \code{n}).
  In the future an argument may be available so that there
  is always \code{n} rows no matter how many responses were
  inputted.



% Typically, a list of length \code{nsim} of simulated responses.  Where
% appropriate the result can be a data frame (which is a special type of
% list).
% %% a *matrix* seems very natural and is more efficient
% %% for large-scale simulation, already for stats:::simulate.lm (in ../R/lm.R )



% For the \code{"lm"} method, the result is a data frame with an
% attribute \code{"seed"}.  If argument \code{seed} is \code{NULL}, the
% attribute is the value of \code{\link{.Random.seed}} before the
% simulation was started; otherwise it is the value of the argument with
% a \code{"kind"} attribute with value \code{as.list(\link{RNGkind}())}.



}
\details{
  This is a methods function for \code{\link[stats]{simulate}}
  and hopefully should behave in a very similar manner.
  Only \pkg{VGAM} family functions with a \code{simslot} slot
  have been implemented for \code{\link[stats]{simulate}}.




}
\seealso{
  Currently the \pkg{VGAM} family functions with a
  \code{simslot} slot are:
  \code{\link{alaplace1}},
  \code{\link{alaplace2}},
  \code{\link{betabinomial}},
  \code{\link{betabinomialff}},
  \code{\link{betaR}},
  \code{\link{betaff}},
  \code{\link{biamhcop}},
  \code{\link{bifrankcop}},
  \code{\link{bilogistic}},
  \code{\link{binomialff}},
  \code{\link{binormal}},
  \code{\link{binormalcop}},
  \code{\link{biclaytoncop}},
  \code{\link{cauchy}},
  \code{\link{cauchy1}},
  \code{\link{chisq}},
  \code{\link{dirichlet}},
  \code{\link{dagum}},
  \code{\link{erlang}},
  \code{\link{exponential}},
  \code{\link{bifgmcop}},
  \code{\link{fisk}},
  \code{\link{gamma1}},
  \code{\link{gamma2}},
  \code{\link{gammaR}},
  \code{\link{gengamma.stacy}},
  \code{\link{geometric}},
  \code{\link{gompertz}},
  \code{\link{gumbelII}},
  \code{\link{hzeta}},
  \code{\link{inv.lomax}},
  \code{\link{inv.paralogistic}},
  \code{\link{kumar}},
  \code{\link{lgamma1}},
  \code{\link{lgamma3}},
  \code{\link{lindley}},
  \code{\link{lino}},
  \code{\link{logff}},
  \code{\link{logistic1}},
  \code{\link{logistic}},
  \code{\link{lognormal}},
  \code{\link{lomax}},
  \code{\link{makeham}},
  \code{\link{negbinomial}},
  \code{\link{negbinomial.size}},
  \code{\link{paralogistic}},
  \code{\link{perks}},
  \code{\link{poissonff}},
  \code{\link{posnegbinomial}},
  \code{\link{posnormal}},
  \code{\link{pospoisson}},
  \code{\link{polya}},
  \code{\link{polyaR}},
  \code{\link{posbinomial}},
  \code{\link{rayleigh}},
  \code{\link{riceff}},
  \code{\link{simplex}},
  \code{\link{sinmad}},
  \code{\link{slash}},
  \code{\link{studentt}},
  \code{\link{studentt2}},
  \code{\link{studentt3}},
  \code{\link{triangle}},
  \code{\link{uninormal}},
  \code{\link{yulesimon}},
  \code{\link{zageometric}},
  \code{\link{zageometricff}},
  \code{\link{zanegbinomial}},
  \code{\link{zanegbinomialff}},
  \code{\link{zapoisson}},
  \code{\link{zapoissonff}},
  \code{\link{zigeometric}},
  \code{\link{zigeometricff}},
  \code{\link{zinegbinomial}},
  \code{\link{zipf}},
  \code{\link{zipoisson}},
  \code{\link{zipoissonff}}.



% \code{\link{logF}},
% \code{\link{tobit}},


  See also
  \code{\link{RNG}} about random number generation in \R,
  \code{\link{vglm}}, \code{\link{vgam}} for model fitting.


}

\section{Warning}{
  With multiple response and/or multivariate responses,
  the order of the elements may differ.
  For some \pkg{VGAM} families, the order is
  \eqn{n \times N \times F}{n x N x F},
  where \eqn{n} is the sample size,
  \eqn{N} is \code{nsim} and
  \eqn{F} is \code{ncol(fitted(vglmObject))}.
  For other \pkg{VGAM} families, the order is
  \eqn{n \times F \times N}{n x F x N}.
  An example of each is given below.


}


\examples{
nn <- 10; mysize <- 20; set.seed(123)
bdata <- data.frame(x2 = rnorm(nn))
bdata <- transform(bdata,
  y1   = rbinom(nn, size = mysize, p = logit(1+x2, inverse = TRUE)),
  y2   = rbinom(nn, size = mysize, p = logit(1+x2, inverse = TRUE)),
  f1   = factor(as.numeric(rbinom(nn, size = 1,
                                  p = logit(1+x2, inverse = TRUE)))))
(fit1 <- vglm(cbind(y1, aaa = mysize - y1) ~ x2,  # Matrix response (2-colns)
              binomialff, data = bdata))
(fit2 <- vglm(f1 ~ x2, binomialff, model = TRUE, data = bdata)) # Factor response

set.seed(123); simulate(fit1, nsim = 8)
set.seed(123); c(simulate(fit2, nsim = 3))  # Use c() when model = TRUE

# An n x N x F example
set.seed(123); n <- 100
bdata <- data.frame(x2 = runif(n), x3 = runif(n))
bdata <- transform(bdata, y1 = rnorm(n, 1 + 2 * x2),
                          y2 = rnorm(n, 3 + 4 * x2))
fit1 <- vglm(cbind(y1, y2) ~ x2, binormal(eq.sd = TRUE), data = bdata)
nsim <- 1000  # Number of simulations for each observation
my.sims <- simulate(fit1, nsim = nsim)
dim(my.sims)  # A data frame
aaa <- array(unlist(my.sims), c(n, nsim, ncol(fitted(fit1))))  # n by N by F
summary(rowMeans(aaa[, , 1]) - fitted(fit1)[, 1])  # Should be all 0s
summary(rowMeans(aaa[, , 2]) - fitted(fit1)[, 2])  # Should be all 0s

# An n x F x N example
n <- 100; set.seed(111); nsim <- 1000
zdata <- data.frame(x2 = runif(n))
zdata <- transform(zdata, lambda1 =  loge(-0.5 + 2 * x2, inverse = TRUE),
                          lambda2 =  loge( 0.5 + 2 * x2, inverse = TRUE),
                          pstr01  = logit( 0,            inverse = TRUE),
                          pstr02  = logit(-1.0,          inverse = TRUE))
zdata <- transform(zdata, y1 = rzipois(n, lambda = lambda1, pstr0 = pstr01),
                          y2 = rzipois(n, lambda = lambda2, pstr0 = pstr02))
zip.fit  <- vglm(cbind(y1, y2) ~ x2, zipoissonff, data = zdata, crit = "coef")
my.sims <- simulate(zip.fit, nsim = nsim)
dim(my.sims)  # A data frame
aaa <- array(unlist(my.sims), c(n, ncol(fitted(zip.fit)), nsim))  # n by F by N
summary(rowMeans(aaa[, 1, ]) - fitted(zip.fit)[, 1])  # Should be all 0s
summary(rowMeans(aaa[, 2, ]) - fitted(zip.fit)[, 2])  # Should be all 0s
}
\keyword{models}
\keyword{datagen}




