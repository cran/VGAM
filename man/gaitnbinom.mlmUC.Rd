\name{Gaitnbinom.mlm}
\alias{Gaitnbinom.mlm}
\alias{dgaitnbinom.mlm}
\alias{pgaitnbinom.mlm}
\alias{qgaitnbinom.mlm}
\alias{rgaitnbinom.mlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated
  Negative Binomial Distribution
  (GAIT--NB--MLM--MLM)


% (multinomial logit model based; GAIT--NB--MLM--MLM)

}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated and --truncated
  negative binomial distribution,
  based on the multinomial logit model (MLM).
  This distribution is sometimes abbreviated as
  GAIT--NB--MLM--MLM.

  
}
\usage{
dgaitnbinom.mlm(x, size, prob = NULL, munb = NULL, alter = NULL,
                inflate = NULL, truncate = NULL, max.support = Inf,
                pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
                log.arg = FALSE)
pgaitnbinom.mlm(q, size, prob = NULL, munb = NULL, alter = NULL,
                inflate = NULL, truncate = NULL, max.support = Inf,
                pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
qgaitnbinom.mlm(p, size, prob = NULL, munb = NULL, alter = NULL,
                inflate = NULL, truncate = NULL, max.support = Inf,
                pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
rgaitnbinom.mlm(n, size, prob = NULL, munb = NULL, alter = NULL,
                inflate = NULL, truncate = NULL, max.support = Inf,
                pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log.arg}{
  Same meaning as in \code{\link[stats:NegBinomial]{NegBinomial}}.

  }

  \item{size, prob, munb}{
  Same meaning as in \code{\link[stats:NegBinomial]{NegBinomial}},
  i.e., for an ordinary negative binomial distribution.
  Argument \code{mu} has been renamed to \code{munb} to emphasize that
  it refers to the negative binomial mean.
  Only one of \code{prob} and \code{munb} should be specified.
  Short vectors are recycled.





  }
  \item{alter, inflate, truncate}{
  Se \code{\link{Gaitpois}}.


  
  
% Must be sorted and have unique values only.
  }
  \item{pstr.i, byrow.arg}{
  See \code{\link{Gaitpois}}.


  }
  \item{pobs.a, max.support}{
  See \code{\link{Gaitpois}}.

  
  
  }
}
\details{
  Largely similar to \code{\link{Gaitpois}} except a
  negative binomial distribution replaces the Poisson distribution.
  Special cases of these functions include
  \code{\link{dposnegbin}},
  \code{\link{dzanegbin}},
  \code{\link{dzinegbin}},
  \code{\link{pposnegbin}},
  \code{\link{pzanegbin}},
  \code{\link{pzinegbin}},
  \code{\link{qposnegbin}},
  \code{\link{qzanegbin}},
  \code{\link{qzinegbin}},
  \code{\link{rposnegbin}},
  \code{\link{rzanegbin}},
  \code{\link{rzinegbin}}
  corresponding to the special value 0.
  However, now arguments \code{alter}, \code{inflate} and \code{truncate}
  have been supplied in one function.




}
%\section{Warning }{
% See \code{\link{rgaitnbinom.mlm}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgaitnbinom.mlm(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgaitnbinom.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitnbinom.mlm} gives the density (PMF),
  \code{pgaitnbinom.mlm} gives the distribution function,
  \code{qgaitnbinom.mlm} gives the quantile function, and
  \code{rgaitnbinom.mlm} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:NegBinomial]{dnegbin}},
  \code{\link[stats:NegBinomial]{pnegbin}},
  \code{\link[stats:NegBinomial]{qnegbin}},
  \code{\link[stats:NegBinomial]{rnegbin}}
  respectively.
  


}
%\references{
%None.
%}

\author{ T. W. Yee. }
%\note{
% 20120405; no longer true to a superior method:
% For \code{rposnegbin}, the arguments of the function are fed
% into \code{\link[stats:NegBinomial]{rnegbin}} until \eqn{n} positive
% values are obtained.  This may take a long time if \code{lambda}
% has values close to 0.


%  The family function \code{\link{gaitnbinom.mlm}} estimates
%  the two parameters by maximum likelihood estimation.


%}




% \code{\link{gaitnegbinomial.mix}},
\seealso{
  \code{\link{gatnbinomial.mlm}},
  \code{\link{Gaitpois}},
  \code{\link[stats:NegBinomial]{NegBinomial}}.






% \code{\link{Posnegbin}},
% \code{\link{gaitnbinomial.mix}},
% \code{\link{Gtnbinom}},
% \code{\link{Ganbinom.mix}},
% \code{\link{Ganbinom.mlm}},

%  \code{\link{zanegbinson}},
%  \code{\link{zinegbinson}},


}
\examples{
avec <- c(5, 10)   # Alter these (special) values
ivec <- c(15, 20)  # Inflate these (special) values
tvec <- 0          # Truncate this (special) value
pobs.a <- c(0.1, 0.2)
pstr.i <- c(0.1, 0.1); size <- 3; munb <- 15; x <- 0:22
y <- rgaitnbinom.mlm(1000, size, munb = munb, alter = avec,
                     inflate = ivec, truncate = tvec, max.support = 20,
                     pstr.i = pstr.i, pobs.a = pobs.a, byrow = TRUE)
table(y)
(pmf <- dgaitnbinom.mlm(x, size, munb = munb, alter = avec,
                        inflate = ivec, truncate = tvec, max.support = 20,
                        pstr.i = pstr.i, pobs.a = pobs.a, byrow = TRUE))
\dontrun{ plot(x, pmf, type = "h", col = "blue", las = 1,
               main = "Heaped data that is also upper-truncated") }
}
\keyword{distribution}

