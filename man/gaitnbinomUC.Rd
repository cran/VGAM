\name{Gaitnbinom}
\alias{Gaitnbinom}
\alias{dgaitnbinom}
\alias{pgaitnbinom}
\alias{qgaitnbinom}
\alias{rgaitnbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated
  Negative Binomial Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated and --truncated
  negative binomial distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.
  Altogether it can be abbreviated as
  GAAIIT--NB(size.p, munb.p)--NB(size.a, munb.a)--MLM--NB(size.i,
  munb.i)--MLM,
  and it is also known as the GAIT-NB PNP combo where
  PNP stands for parametric and nonparametric.


}
\usage{
dgaitnbinom(x, size.p, prob.p = NULL, munb.p = NULL,
            alt.mix = NULL, alt.mlm = NULL,
            inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
            max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
            size.a = size.p, size.i = size.p,
            prob.a = prob.p, prob.i = prob.p,
            munb.a = munb.p, munb.i = munb.p,
            deflation = FALSE, log = FALSE)
pgaitnbinom(q, size.p, prob.p = NULL, munb.p = NULL,
            alt.mix = NULL, alt.mlm = NULL,
            inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
            max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
            size.a = size.p, size.i = size.p,
            prob.a = prob.p, prob.i = prob.p,
            munb.a = munb.p, munb.i = munb.p,
            lower.tail = TRUE)
qgaitnbinom(p, size.p, prob.p = NULL, munb.p = NULL,
            alt.mix = NULL, alt.mlm = NULL,
            inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
            max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
            size.a = size.p, size.i = size.p,
            prob.a = prob.p, prob.i = prob.p,
            munb.a = munb.p, munb.i = munb.p)
rgaitnbinom(n, size.p, prob.p = NULL, munb.p = NULL,
            alt.mix = NULL, alt.mlm = NULL,
            inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
            max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
            size.a = size.p, size.i = size.p,
            prob.a = prob.p, prob.i = prob.p,
            munb.a = munb.p, munb.i = munb.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link[stats]{rnbinom}}.

  }

  \item{size.p, prob.p, munb.p}{
  Same meaning as in \code{\link[stats]{rnbinom}}.
  See \code{\link[VGAM]{Gaitpois}} for generic information.
  
 
  }
  \item{size.a, prob.a, munb.a}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.
  
 
  }
  \item{size.i, prob.i, munb.i}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.
  
 
  }
  \item{truncate, max.support}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.

  }
  \item{alt.mix, inf.mix}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.
  

  }

  \item{alt.mlm, inf.mlm}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.


  }
  
  \item{pobs.mlm, pstr.mlm, byrow.ai}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.

  }
  
  \item{pobs.mix, pstr.mix}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.

  
  
  }
  \item{deflation}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.



  
  }
}
\details{
  These functions for the NBD are analogous to the Poisson,
  hence most details have been put in
  \code{\link[VGAM]{Gaitpois}}.
  The NBD has two possible parameterizations: one
  involving a probability (argument begins with \code{prob})
  and the other the mean (beginning with \code{mu}).
  Because \code{\link[stats]{NegBinomial}} only allows
  one of these arguments to be used, the functions here
  have the same behaviour.


}
%\section{Warning }{


%}
\value{
  \code{dgaitnbinom} gives the density,
  \code{pgaitnbinom} gives the distribution function,
  \code{qgaitnbinom} gives the quantile function, and
  \code{rgaitnbinom} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:NegBinomial]{dnbinom}},
  \code{\link[stats:NegBinomial]{pnbinom}},
  \code{\link[stats:NegBinomial]{qnbinom}},
  \code{\link[stats:NegBinomial]{rnbinom}}
  respectively.
  


}
%\references{

%Yee, T. W. and Ma, C. (2020).
%Generally--altered, --inflated and --truncated regression,
%with application to heaped and seeped counts.  % count data.
%\emph{In preparation}.


%}

\author{ T. W. Yee. }
\note{
  Four functions were moved from \pkg{VGAM} to \pkg{VGAMdata};
  they can be seen at \code{\link[VGAMdata]{Posnegbin}}.
  It is better to use
  \code{dgaitnbinom(x, size, munb.p = munb, truncate = 0)} instead of
  \code{dposnbinom(x, size, munb = munb)}, etc.





}



\seealso{
  \code{\link[VGAM]{Gaitpois}},
  \code{\link[VGAM]{multinomial}},
  \code{\link[VGAM]{Gaitbinom}},
  \code{\link[VGAM]{Gaitlog}},
  \code{\link[VGAM]{Gaitzeta}}.


% \code{\link[VGAMsecret]{Gaitgenpois1}},


% \code{\link{gaitpoisson}},
% \code{\link{Zapois}},
% \code{\link{Zipois}},
% \code{\link{Pospois}}
% \code{\link[stats:Poisson]{Poisson}};
  
%  \code{\link{Gaitlog.mix}}  and \code{\link{Gaitlog.mlm}},
%  \code{\link{Gaitpois.mix}} and \code{\link{Gaitpois.mlm}},
%  \code{\link{Gaitnbinom.mlm}},

% \code{\link{gaitpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},



}
\examples{size <- 10
ivec <- c(6, 14); avec <- c(8, 11); munb <- 10; xgrid <- 0:25
tvec <- 15; max.support <- 20; pobs.a <- 0.05; pstr.i <- 0.25
(ddd <- dgaitnbinom(xgrid, size, munb.p = munb, munb.a = munb + 5,
   truncate = tvec, max.support = max.support, pobs.mix = pobs.a,
   pobs.mlm = pobs.a, alt.mlm = avec,
   pstr.mix = pstr.i, inf.mix = ivec))
\dontrun{plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
              main = "GAIT PNP Combo PMF---NB Parent")
mylwd <- 1
abline(v = avec, col = 'blue', lwd = mylwd)
abline(v = ivec, col = 'purple', lwd = mylwd)
abline(v = tvec, col = 'tan', lwd = mylwd)
abline(v = max.support, col = 'magenta', lwd = mylwd)
abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
lines(xgrid, dnbinom(xgrid, size, mu = munb), col = 'gray', lty = 2)
lines(xgrid, ddd, type = "h", col = "pink", lwd = 7)  # GAIT PNP combo PMF
points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16)  }
}
\keyword{distribution}


% 20200815; checked identical results to [dpqr]gaitnbinom.mix() &
%   [dpqr]gaitnbinom.mlm().



