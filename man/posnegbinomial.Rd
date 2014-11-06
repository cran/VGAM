\name{posnegbinomial}
\alias{posnegbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Negative Binomial Distribution Family Function }
\description{
  Maximum likelihood estimation of the two parameters of a positive
  negative binomial distribution.

}
\usage{
posnegbinomial(lmunb = "loge", lsize = "loge",
               isize = NULL, zero = -2, nsimEIM = 250,
               ishrinkage = 0.95, imethod = 1)

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
  \item{nsimEIM, zero}{ 
  See \code{\link{CommonVGAMffArguments}}.


  }
  \item{ishrinkage, imethod}{
  See \code{\link{negbinomial}}.

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


  This function handles \emph{multivariate} responses, so that a matrix
  can be used as the response. The number of columns is the number
  of species, say, and setting \code{zero = -2} means that \emph{all}
  species have a \code{k} equalling a (different) intercept only.


}
\section{Warning}{
  The Poisson model corresponds to \code{k} equalling infinity.
  If the data is Poisson or close to Poisson, numerical problems may
  occur. Possibly a loglog link could be added in the future to try help
  handle this problem.


  This \pkg{VGAM} family function is computationally expensive
  and usually runs slowly;
  setting \code{trace = TRUE} is useful for monitoring convergence.


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
    This family function handles multiple responses.


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
\dontrun{
pdata <- data.frame(x2 = runif(nn <- 1000))
pdata <- transform(pdata, y1 = rposnegbin(nn, munb = exp(0+2*x2), size = exp(1)),
                          y2 = rposnegbin(nn, munb = exp(1+2*x2), size = exp(3)))
fit <- vglm(cbind(y1, y2) ~ x2, posnegbinomial, data = pdata, trace = TRUE)
coef(fit, matrix = TRUE)
dim(depvar(fit))  # dim(fit@y) is not as good


# Another artificial data example
pdata2 <- data.frame(munb = exp(2), size = exp(3)); nn <- 1000
pdata2 <- transform(pdata2, y3 = rposnegbin(nn, munb = munb, size = size))
with(pdata2, table(y3))
fit <- vglm(y3 ~ 1, posnegbinomial, pdata2, trace = TRUE)
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
print( with(corbet, cbind(ofreq, species, fitted = pdf2*sum(species))), digits = 1)
with(corbet,
matplot(ofreq, cbind(species, fitted = pdf2*sum(species)), las = 1,
        xlab = "Observed frequency (of individual butterflies)",
        type = "b", ylab = "Number of species", col = c("blue", "orange"),
        main = "blue 1s = observe; orange 2s = fitted"))
}
}
\keyword{models}
\keyword{regression}


% bigN = with(corbet, sum(species))



