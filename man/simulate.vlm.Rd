% 20131230; adapted from simulate.Rd from R 3.0.2


\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}


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
  Note that many \pkg{VGAM} family functions can handle
  multiple responses. This can result in a longer data frame
  with more rows
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
  \code{\link{amh}},
  \code{\link{betabinomial}},
  \code{\link{betabinomial.ab}},
  \code{\link{beta.ab}},
  \code{\link{betaff}},
  \code{\link{bifrankcop}},
  \code{\link{bilogistic4}},
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
  \code{\link{fgm}},
  \code{\link{fisk}},
  \code{\link{gamma1}},
  \code{\link{gamma2}},
  \code{\link{gamma2.ab}},
  \code{\link{gengamma}},
  \code{\link{geometric}},
  \code{\link{gompertz}},
  \code{\link{gumbelII}},
  \code{\link{hzeta}},
  \code{\link{invlomax}},
  \code{\link{invparalogistic}},
  \code{\link{kumar}},
  \code{\link{lgammaff}},
  \code{\link{lgamma3ff}},
  \code{\link{lindley}},
  \code{\link{lino}},
  \code{\link{logff}},
  \code{\link{logistic1}},
  \code{\link{logistic2}},
  \code{\link{lognormal}},
  \code{\link{lognormal3}},
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
\examples{
nn <- 10; mysize <- 20
set.seed(123)
bdata <- data.frame(x2 = rnorm(nn))
bdata <- transform(bdata,
  y1   = rbinom(nn, size = mysize, p = logit(1+x2, inverse = TRUE)),
  y2   = rbinom(nn, size = mysize, p = logit(1+x2, inverse = TRUE)),
  f1   = factor(as.numeric(rbinom(nn, size = 1,
                                   p = logit(1+x2, inverse = TRUE)))))
(fit1 <- vglm(cbind(y1, aaa = mysize - y1) ~ x2,  # Matrix response (2-colns)
              binomialff, data = bdata))
(fit2 <- vglm(f1 ~ x2,  # Factor response
              binomialff, model = TRUE, data = bdata))

set.seed(123)
simulate(fit1, nsim = 8)

set.seed(123)
c(simulate(fit2, nsim = 3))  # Use c() when model = TRUE
}
\keyword{models}
\keyword{datagen}




