\name{dirmul.old}
\alias{dirmul.old}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a Dirichlet-Multinomial Distribution }
\description{
  Fits a Dirichlet-multinomial distribution to a matrix of
  non-negative integers.
}
\usage{
dirmul.old(link = "loge", ialpha = 0.01, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to each of the \eqn{M} (positive)
  shape parameters \eqn{\alpha_j}{alpha_j} for \eqn{j=1,\ldots,M}.
  See \code{\link{Links}} for more choices.
  Here, \eqn{M} is the number of columns of the response matrix.

  }
  \item{ialpha}{
  Numeric vector. Initial values for the
  \code{alpha} vector. Must be positive.
  Recycled to length \eqn{M}.

  }
  \item{parallel}{
  A logical, or formula specifying which terms have equal/unequal
  coefficients.

  }
  \item{zero}{
  An integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\}.

  }
}
% formula is p.49 of Lange 2002.
\details{
  The Dirichlet-multinomial distribution, which is somewhat similar to a
  Dirichlet distribution, has probability function
  \deqn{P(Y_1=y_1,\ldots,Y_M=y_M) =
     {2y_{*} \choose {y_1,\ldots,y_M}}
        \frac{\Gamma(\alpha_{+})}{\Gamma(2y_{*}+\alpha_{+})}
        \prod_{j=1}^M \frac{\Gamma(y_j+\alpha_{j})}{\Gamma(\alpha_{j})}}{%
        P(Y_1=y_1,\ldots,Y_M=y_M) = C_{y_1,\ldots,y_M}^{2y_{*}}
         Gamma(alpha_+) / Gamma( 2y_* + alpha_+)
         prod_{j=1}^M [ Gamma( y_j+ alpha_j) / Gamma( alpha_j)]}
       for \eqn{\alpha_j > 0}{alpha_j  > 0},
       \eqn{\alpha_+ = \alpha_1 + \cdots + \alpha_M}{alpha_+ = alpha_1 + \cdots + alpha_M},
       and \eqn{2y_{*} = y_1 + \cdots + y_M}{2y_* = y_1 + \cdots + y_M}.
       Here, \eqn{a \choose b}{C_b^a} means ``\eqn{a} choose \eqn{b}'' and
       refers to combinations (see \code{\link[base]{choose}}).
The (posterior) mean is
  \deqn{E(Y_j) = (y_j + \alpha_j) / (2y_{*} + \alpha_{+})}{%
        E(Y_j) = (y_j +  alpha_j) / (2y_{*} +  alpha_+)}
  for \eqn{j=1,\ldots,M}{j=1,\ldots,M}, and these are returned as
  the fitted values as a \eqn{M}-column matrix.


%  One situation that arises for the Dirichlet-multinomial distribution
%   is a locus with M codominant alleles. If in a sample of y_* people,
%   allele i appears y_j times, then the maximum likelihood estimate of
%   the ith allele frequency is y_j / (2y_*).


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{
Lange, K. (2002)
\emph{Mathematical and Statistical Methods for Genetic Analysis},
2nd ed.  New York: Springer-Verlag.


Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


Paul, S. R., Balasooriya, U. and Banerjee, T. (2005)
Fisher information matrix of the Dirichlet-multinomial distribution.
\emph{Biometrical Journal}, \bold{47}, 230--236.


Tvedebrink, T. (2010)
Overdispersion in allelic counts and \eqn{\theta}-correction in forensic genetics.
\emph{Theoretical Population Biology}, \bold{78}, 200--210.



%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}
\author{ Thomas W. Yee }

\note{
  The response should be a matrix of non-negative values.
  Convergence seems to slow down if there are zero values.
  Currently, initial values can be improved upon.


  This function is almost defunct and may be withdrawn soon.
  Use \code{\link{dirmultinomial}} instead.


}

\seealso{
  \code{\link{dirmultinomial}},
  \code{\link{dirichlet}},
  \code{\link{betabinomialff}},
  \code{\link{multinomial}}.


}

\examples{
# Data from p.50 of Lange (2002)
alleleCounts <- c(2, 84, 59, 41, 53, 131, 2, 0,
       0, 50, 137, 78, 54, 51, 0, 0,
       0, 80, 128, 26, 55, 95, 0, 0,
       0, 16, 40, 8, 68, 14, 7, 1)
dim(alleleCounts) <- c(8, 4)
alleleCounts <- data.frame(t(alleleCounts))
dimnames(alleleCounts) <- list(c("White","Black","Chicano","Asian"),
                    paste("Allele", 5:12, sep = ""))

set.seed(123)  # @initialize uses random numbers
fit <- vglm(cbind(Allele5,Allele6,Allele7,Allele8,Allele9,
                  Allele10,Allele11,Allele12) ~ 1, dirmul.old,
             trace = TRUE, crit = "c", data = alleleCounts)

(sfit <- summary(fit))
vcov(sfit)
round(eta2theta(coef(fit), fit@misc$link, fit@misc$earg), digits = 2)  # not preferred
round(Coef(fit), digits = 2)  # preferred
round(t(fitted(fit)), digits = 4)  # 2nd row of Table 3.5 of Lange (2002)
coef(fit, matrix = TRUE)


pfit <- vglm(cbind(Allele5,Allele6,Allele7,Allele8,Allele9,
                   Allele10,Allele11,Allele12) ~ 1,
             dirmul.old(parallel = TRUE), trace = TRUE,
             data = alleleCounts)
round(eta2theta(coef(pfit, matrix = TRUE), pfit@misc$link,
                pfit@misc$earg), digits = 2)  # 'Right' answer
round(Coef(pfit), digits = 2)  # 'Wrong' answer due to parallelism constraint
}
\keyword{models}
\keyword{regression}


