\name{posnegbinomial}
\alias{posnegbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Negative Binomial Distribution Family Function }
\description{
  Maximum likelihood estimation of the two parameters of a positive
  negative binomial distribution.

}
\usage{
posnegbinomial(zero = "size", type.fitted = c("mean", "munb", "prob0"),
               mds.min = 0.001, nsimEIM = 500, cutoff.prob = 0.999,
               eps.trig = 1e-07, max.support = 4000, max.chunk.MB = 30,
               lmunb = "loge", lsize = "loge", imethod = 1,
               imunb = NULL, iprobs.y = NULL,
               gprobs.y = ppoints(8), isize = NULL,
               gsize.mux = exp(c(-30, -20, -15, -10, -6:3)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmunb}{
  Link function applied to the \code{munb} parameter, which is the mean
  \eqn{\mu_{nb}}{munb} of an ordinary negative binomial distribution.
  See \code{\link{Links}} for more choices.


  }
  \item{lsize}{
  Parameter link function applied to the dispersion parameter,
  called \code{k}.
  See \code{\link{Links}} for more choices.


  }
  \item{isize}{
  Optional initial value for \code{k}, an index parameter.
  The value \code{1/k} is known as a dispersion parameter.
  If failure to converge occurs try different values (and/or use
  \code{imethod}).
  If necessary this vector is recycled to length equal to the number
  of responses.
  A value \code{NULL} means an initial value for each response is
  computed internally using a range of values.


  }
  \item{nsimEIM, zero, eps.trig}{
  See \code{\link{CommonVGAMffArguments}}.


  }
  \item{mds.min, iprobs.y, cutoff.prob}{
    Similar to \code{\link{negbinomial}}.


  }
  \item{imunb, max.support}{
    Similar to \code{\link{negbinomial}}.


  }
  \item{max.chunk.MB, gsize.mux}{
    Similar to \code{\link{negbinomial}}.


  }
  \item{imethod, gprobs.y}{
  See \code{\link{negbinomial}}.


  }
  \item{type.fitted}{
    See \code{\link{CommonVGAMffArguments}} for details.


  }
}
\details{
  The positive negative binomial distribution is an ordinary negative
  binomial distribution but with the probability of a zero response
  being zero. The other probabilities are scaled to sum to unity.


  This family function is based on \code{\link{negbinomial}} and most
  details can be found there. To avoid confusion, the parameter
  \code{munb} here corresponds to the mean of an ordinary negative
  binomial distribution \code{\link{negbinomial}}. The mean of
  \code{posnegbinomial} is
  \deqn{\mu_{nb} / (1-p(0))}{%
            munb / (1-p(0))}
  where \eqn{p(0) = (k/(k + \mu_{nb}))^k}{p(0) = (k/(k + munb))^k} is the
  probability an ordinary negative binomial distribution has a zero value.


  The parameters \code{munb} and \code{k} are not independent in the
  positive negative binomial distribution, whereas they are in the
  ordinary negative binomial distribution.


  This function handles \emph{multiple} responses, so that a matrix
  can be used as the response. The number of columns is the number
  of species, say, and setting \code{zero = -2} means that \emph{all}
  species have a \code{k} equalling a (different) intercept only.


}
\section{Warning}{
  This family function is fragile;
  at least two cases will lead to numerical problems.
  Firstly,
  the positive-Poisson model corresponds to \code{k} equalling infinity.
  If the data is positive-Poisson or close to positive-Poisson,
  then the estimated \code{k} will diverge to \code{Inf} or some
  very large value.
  Secondly, if the data is clustered about the value 1 because
  the \code{munb} parameter is close to 0
  then numerical problems will also occur.
  Users should set \code{trace = TRUE} to monitor convergence.
  In the situation when both cases hold, the result returned
  (which will be untrustworthy) will depend on the initial values.



The negative binomial distribution (NBD) is a strictly unimodal distribution.
Any data set that does not exhibit a mode (in the middle) makes the
estimation problem difficult.
The positive NBD inherits this feature.
Set \code{trace = TRUE} to monitor convergence.



  See the example below of a data set where \code{posbinomial()}
  fails; the so-called solution is \emph{extremely} poor.
  This is partly due to a lack of a
  unimodal shape because the number of counts decreases only.
  This long tail makes it very difficult to estimate the mean
  parameter with any certainty. The result too is that the
  \code{size} parameter is numerically fraught.


%  Then trying a \code{\link{loglog}} link might help
%  handle this problem.


  This \pkg{VGAM} family function inherits the same warnings as
  \code{\link{negbinomial}}.
  And if \code{k} is much less than 1 then the estimation may
  be slow.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.


}
\references{
  Barry, S. C. and Welsh, A. H. (2002)
  Generalized additive modelling and zero inflated count data.
  \emph{Ecological Modelling},
  \bold{157},
  179--188.


  Williamson, E. and Bretherton, M. H. (1964)
  Tables of the logarithmic series distribution.
  \emph{Annals of Mathematical Statistics},
  \bold{35},
  284--297.


}
\author{ Thomas W. Yee }
\note{
  If the estimated \eqn{k} is very large then fitting a
  \code{\link{pospoisson}} model is a good idea.



  If both \code{munb} and \eqn{k} are large then it may be
  necessary to decrease \code{eps.trig} and increase
  \code{max.support} so that the EIMs are positive-definite,
  e.g.,
  \code{eps.trig = 1e-8} and \code{max.support = Inf}.



}

\seealso{
  \code{\link{rposnegbin}},
  \code{\link{pospoisson}},
  \code{\link{negbinomial}},
  \code{\link{zanegbinomial}},
  \code{\link[stats:NegBinomial]{rnbinom}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{corbet}},
  \code{\link{logff}},
  \code{\link{simulate.vlm}}.


% \code{\link[MASS]{rnegbin}}.


}

\examples{
pdata <- data.frame(x2 = runif(nn <- 1000))
pdata <- transform(pdata, y1 = rposnegbin(nn, munb = exp(0+2*x2), size = exp(1)),
                          y2 = rposnegbin(nn, munb = exp(1+2*x2), size = exp(3)))
fit <- vglm(cbind(y1, y2) ~ x2, posnegbinomial, data = pdata, trace = TRUE)
coef(fit, matrix = TRUE)
dim(depvar(fit))  # Using dim(fit@y) is not recommended


# Another artificial data example
pdata2 <- data.frame(munb = exp(2), size = exp(3)); nn <- 1000
pdata2 <- transform(pdata2, y3 = rposnegbin(nn, munb = munb, size = size))
with(pdata2, table(y3))
fit <- vglm(y3 ~ 1, posnegbinomial, data = pdata2, trace = TRUE)
coef(fit, matrix = TRUE)
with(pdata2, mean(y3))  # Sample mean
head(with(pdata2, munb/(1-(size/(size+munb))^size)), 1)  # Population mean
head(fitted(fit), 3)
head(predict(fit), 3)


# Example: Corbet (1943) butterfly Malaya data
fit <- vglm(ofreq ~ 1, posnegbinomial, weights = species, data = corbet)
coef(fit, matrix = TRUE)
Coef(fit)
(khat <- Coef(fit)["size"])
pdf2 <- dposnegbin(x = with(corbet, ofreq), mu = fitted(fit), size = khat)
print(with(corbet, cbind(ofreq, species, fitted = pdf2*sum(species))), dig = 1)
\dontrun{
with(corbet,
matplot(ofreq, cbind(species, fitted = pdf2*sum(species)), las = 1,
        xlab = "Observed frequency (of individual butterflies)",
        type = "b", ylab = "Number of species", col = c("blue", "orange"),
        main = "blue 1s = observe; orange 2s = fitted"))
}

\dontrun{
# This data (courtesy of Maxim Gerashchenko) causes posbinomial() to fail
pnbd.fail <- data.frame(
 y1 = c(1:16, 18:21, 23:28, 33:38, 42, 44, 49:51, 55, 56, 58,
 59, 61:63, 66, 73, 76, 94, 107, 112, 124, 190, 191, 244),
 ofreq = c(130, 80, 38, 23, 22, 11, 21, 14, 6, 7, 9, 9, 9, 4, 4, 5, 1,
           4, 6, 1, 3, 2, 4, 3, 4, 5, 3, 1, 2, 1, 1, 4, 1, 2, 2, 1, 3,
           1, 1, 2, 2, 2, 1, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1))
fit.fail <- vglm(y1 ~ 1, weights = ofreq, posnegbinomial,
               trace = TRUE, data = pnbd.fail)
}
}
\keyword{models}
\keyword{regression}


% bigN = with(corbet, sum(species))



%posnegbinomial(lmunb = "loge", lsize = "loge", imunb = NULL,
%               isize = NULL, zero = "size", nsimEIM = 250,
%               probs.y = 0.75, cutoff.prob = 0.999,
%               max.support = 2000, max.chunk.MB = 30,
%               gsize = exp((-4):4), ishrinkage = 0.95, imethod = 1)




