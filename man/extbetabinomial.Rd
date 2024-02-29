\name{extbetabinomial}
\alias{extbetabinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extended
  Beta-binomial Distribution Family Function }
\description{
  Fits an extended beta-binomial distribution
  by maximum
  likelihood estimation.  The two parameters
  here are the mean and correlation coefficient.


}
\usage{
extbetabinomial(lmu = "logitlink", lrho = "cloglink", 
      zero = "rho", irho = 0, grho = c(0, 0.05, 0.1, 0.2),
      vfl = FALSE, Form2 = NULL,
      imethod = 1, ishrinkage = 0.95)
}
%- maybe also 'usage' for other objects documented here.
%            ishrinkage = 0.95, nsimEIM = NULL, zero = 2
\arguments{
  \item{lmu, lrho}{
  Link functions applied to the two parameters.
  See \code{\link{Links}} for more choices.
  The first default ensure the mean remain
  in \eqn{(0, 1)},
  while the second allows for a slightly negative
  correlation parameter: you could say it
  lies in
  \eqn{(\max(-\mu/(N-\mu-1), -(1 - \mu)/(N-(1-\mu)-1)), 1)}
  where \eqn{\mu} is the mean (probability) and
  \eqn{N} is \code{size}.
  See below for details.
  For \code{lrho}, 
  \code{\link{cloglink}} is a good choice
  because it handles parameter values from 1
  downwards.
  Other good choices include
  \code{logofflink(offset = 1)} and
  \code{\link{rhobitlink}}.

  
  }
  \item{irho, grho}{
  The first is similar to \code{\link{betabinomial}}
  and it is a good idea to use this argument
  because to conduct a
  grid search based on \code{grho} is expensive.
  The default is effectively a binomial distribution.
  Set \code{irho = NULL} to perform a grid search
  which is more reliable but slow.



  }
  \item{imethod}{
  Similar to \code{\link{betabinomial}}.


  }
  \item{zero}{
  Similar to \code{\link{betabinomial}}.
  Also,
  see \code{\link{CommonVGAMffArguments}}
  for more information.
  Modelling \code{rho} with covariates requires
  large samples.



  }
  \item{ishrinkage}{
    See \code{\link{betabinomial}} and
    \code{\link{CommonVGAMffArguments}} for
    information.


  }
  \item{vfl, Form2}{
    See \code{\link{CommonVGAMffArguments}}.
    If \code{vfl = TRUE} then \code{Form2}
    should be a formula specifying the terms
    for \eqn{\eta_2} and all others are
    used for \eqn{\mu}.
    It is similar to \code{\link{uninormal}}.
    If these arguments are used then
    \code{cbind(0, log(size1 / (size1 - 1)))}
    should be used as an offset, and
    set \code{zero = NULL} too.



    
%   \code{\link{negbinomial}}.


  }
}
\details{
  The \emph{extended} beta-binomial
  distribution (EBBD) proposed
  by Prentice (1986)
  allows for a slightly negative correlation
  parameter whereas the ordinary BBD
  \code{\link{betabinomial}}
  only allows values in \eqn{(0, 1)} so it
  handles overdispersion only.
  When negative, the data is underdispersed
  relative to an ordinary binomial distribution.



  Argument \code{rho} is used here for the
  \eqn{\delta} used in Prentice (1986) because
  it is the correlation between the
  (almost) Bernoulli trials.
  (They are actually simple binary variates.)
  We use here  
  \eqn{N} for the number of trials
  (e.g., litter size),
  \eqn{T=NY} is the number of successes, and
  \eqn{p} (or \eqn{\mu})
  is the probability of a success
  (e.g., a malformation).
  That is, \eqn{Y} is the \emph{proportion}
  of successes. Like
  \code{\link{binomialff}}, the fitted values
  are the
  estimated probability
  of success
  (i.e., \eqn{E[Y]} and not \eqn{E[T]})
  and the prior weights \eqn{N} are attached
  separately on the object in a slot.





  The probability function is difficult
  to write but it involves three
  series of products.
  Recall \eqn{Y = T/N} is the real response
  being modelled, where \eqn{T} is the
  (total) sum of \eqn{N} correlated
  (almost) Bernoulli trials.



  The default model is
  \eqn{\eta_1 = logit(\mu)}{eta1 =logit(mu)}
  and \eqn{\eta_2 = clog(\rho)}{eta2 = clog(rho)}
  because the first
  parameter lies between 0 and 1.
  The second link is \code{\link{cloglink}}.
  The mean of \eqn{Y} is
  \eqn{p=\mu}{p = mu}
  and the variance of \eqn{Y} is
  \eqn{\mu(1-\mu)(1+(N-1)\rho)/N}{mu(1-mu)(1+(N-1)rho)/N}.
  Here, the correlation \eqn{\rho}{rho}
  may be slightly negative
  and is the correlation between the \eqn{N}
  individuals within a litter.
  A \emph{litter effect} is typically reflected
  by a positive value of \eqn{\rho}{rho} and
  corresponds to \emph{overdispersion} with
  respect to the binomial distribution.
  Thus an \emph{exchangeable} error structure
  is assumed between units within a litter
  for the EBBD.



  This family function uses Fisher scoring.
  Elements of the second-order expected
  derivatives 
  are computed numerically, which may
  fail for models very near the boundary of the
  parameter space.
  Usually, the computations
  are expensive for large \eqn{N} because of
  a \code{for} loop, so
  it may take a long time.


}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}}.


  Suppose \code{fit} is a fitted EBB
  model. Then \code{depvar(fit)}
  are the sample proportions \eqn{y},
  \code{fitted(fit)} returns
  estimates of \eqn{E(Y)},
  and \code{weights(fit, type = "prior")} returns
  the number of trials \eqn{N}.


}
\references{

Prentice, R. L. (1986).
Binary regression using an extended beta-binomial
distribution,
with discussion of correlation induced by
covariate measurement errors.
\emph{Journal of the American Statistical
  Association},
\bold{81}, 321--327.


}

\author{ T. W. Yee }
\note{

  
  This function is recommended over
  \code{\link{betabinomial}} and
  \code{\link{betabinomialff}}.
  It processes the input in the
  same way as \code{\link{binomialff}}.
  But it does not handle the case \eqn{N \leq 2}
  very well because there are two parameters to
  estimate, not one, for each row of the input.
  Cases where \eqn{N > 2}
  can be selected via the
  \code{subset} argument of \code{\link{vglm}}.


}
\section{Warning }{


  Modelling \code{rho} using covariates well
  requires much data
  so it is usually best to leave \code{zero}
  alone.
  It is good to set \code{trace = TRUE} and
  play around with \code{irho} if there are
  problems achieving convergence.
  Convergence problems will occur when the
  estimated \code{rho} is close to the
  lower bound,
  i.e., the underdispersion
  is almost too severe for the EBB to cope.


}
\seealso{
  \code{\link{Extbetabinom}},
  \code{\link{betabinomial}},
  \code{\link{betabinomialff}},
  \code{\link{binomialff}},
  \code{\link{dirmultinomial}},
  \code{\link{cloglink}},
  \code{\link{lirat}}.


}
\examples{
# Example 1
edata <- data.frame(N = 10, mu = 0.5, rho = 0.1)
edata <- transform(edata,
      y = rextbetabinom(100, N, mu, rho = rho))
fit1 <- vglm(cbind(y, N-y) ~ 1, extbetabinomial, edata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
head(cbind(depvar(fit1), weights(fit1, type = "prior")))

# Example 2: VFL model
\dontrun{N <- size1 <- 10; nn <- 2000; set.seed(1)
edata <-  # Generate the data set. Expensive.
    data.frame(x2 = runif(nn),
               ooo =  log(size1 / (size1 - 1)))
edata <- transform(edata, x1copy = 1, x2copy = x2,
  y2 = rextbetabinom(nn, size1,  # Expensive
         logitlink(1 + x2, inverse = TRUE),
         cloglink(ooo + 1 - 0.5 * x2, inv = TRUE)))
fit2 <- vglm(data = edata,
        cbind(y2, N - y2) ~ x2 + x1copy + x2copy,
        extbetabinomial(zero = NULL, vfl = TRUE,
                 Form2 = ~ x1copy + x2copy - 1),
        offset = cbind(0, ooo), trace = TRUE)
coef(fit2, matrix = TRUE)
wald.stat(fit2, values0 = c(1, 1, -0.5))
}}
\keyword{models}
\keyword{regression}
