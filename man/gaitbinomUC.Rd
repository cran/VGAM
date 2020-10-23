\name{Gaitbinom}
\alias{Gaitbinom}
\alias{dgaitbinom}
\alias{pgaitbinom}
\alias{qgaitbinom}
\alias{rgaitbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated
  Binomial Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated and --truncated
  binomial distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.
  Altogether it can be abbreviated as
  GAAIIT--Binom(size.p, prob.p)--Binom(size.a, prob.a)--MLM--Binom(size.i,
  prob.i)--MLM,
  and it is also known as the GAIT-Binom PNP combo where
  PNP stands for parametric and nonparametric.


}
\usage{
dgaitbinom(x, size.p, prob.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
           max.support = NULL, pobs.mix = 0, pobs.mlm = 0,
           pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           deflation = FALSE, log = FALSE)
pgaitbinom(q, size.p, prob.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
           max.support = NULL, pobs.mix = 0, pobs.mlm = 0,
           pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           lower.tail = TRUE)
qgaitbinom(p, size.p, prob.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
           max.support = NULL, pobs.mix = 0, pobs.mlm = 0,
           pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p)
rgaitbinom(n, size.p, prob.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
           max.support = NULL, pobs.mix = 0, pobs.mlm = 0,
           pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link[stats]{Binomial}}.

  }

  \item{size.p, prob.p}{
  Same meaning as in \code{\link[stats]{Binomial}}.
  See \code{\link[VGAM]{Gaitpois}} for generic information.
  
 
  }
  \item{size.a, prob.a}{
  See \code{\link[VGAM]{Gaitpois}} for generic information.
  
 
  }
  \item{size.i, prob.i}{
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
  These functions for the GAIT binomial distribution
  are analogous to the GAIT Poisson,
  hence most details have been put in
  \code{\link[VGAM]{Gaitpois}}.



}
%\section{Warning }{


%}
\value{
  \code{dgaitbinom} gives the density,
  \code{pgaitbinom} gives the distribution function,
  \code{qgaitbinom} gives the quantile function, and
  \code{rgaitbinom} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:Binomial]{dbinom}},
  \code{\link[stats:Binomial]{pbinom}},
  \code{\link[stats:Binomial]{qbinom}},
  \code{\link[stats:Binomial]{rbinom}}
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
  Functions \code{\link[VGAMdata]{Posbinom}} have been moved
  to \pkg{VGAMdata}.
  It is better to use
  \code{dgaitbinom(x, size, prob, truncate = 0)} instead of
  \code{dposbinom(x, size, prob)}, etc.



  

}



\seealso{
  \code{\link[VGAM]{Gaitpois}},
  \code{\link[VGAM]{Gaitnbinom}},
  \code{\link[VGAM]{multinomial}},
  \code{\link[VGAM]{Gaitlog}},
  \code{\link[VGAM]{Gaitzeta}}.


% \code{\link[VGAMsecret]{Gaitgenpois1}},
% \code{\link{gaitpoisson.mlm}},
% \code{\link{gaitpoisson}},
% \code{\link{Zapois}},
% \code{\link{Zipois}},
% \code{\link{Pospois}}
% \code{\link[stats:Poisson]{Poisson}};
  
%  \code{\link{Gaitlog.mix}}  and \code{\link{Gaitlog.mlm}},
%  \code{\link{Gaitpois.mix}} and \code{\link{Gaitpois.mlm}},
%  \code{\link{Gaitbinom.mlm}},

% \code{\link{gaitpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},



}
\examples{ size <- 20
ivec <- c(6, 10); avec <- c(8, 11); prob <- 0.5; xgrid <- 0:25
tvec <- 14; pobs.a <- 0.05; pstr.i <- 0.15
(ddd <- dgaitbinom(xgrid, size, prob.p = prob, prob.a = prob + 5,
   truncate = tvec, pobs.mix = pobs.a,
   pobs.mlm = pobs.a, alt.mlm = avec,
   pstr.mix = pstr.i, inf.mix = ivec))
\dontrun{plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
              main = "GAIT PNP Combo PMF---Binomial Parent")
mylwd <- 1
abline(v = avec, col = 'blue', lwd = mylwd)
abline(v = ivec, col = 'purple', lwd = mylwd)
abline(v = tvec, col = 'tan', lwd = mylwd)
abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
lines(xgrid, dbinom(xgrid, size, prob), col = 'gray', lty = 2)
lines(xgrid, ddd, type = "h", col = "pink", lwd = 7)  # GAIT PNP combo PMF
points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16)  }
}
\keyword{distribution}


% 20200815; checked identical results to [dpqr]gaitbinom.mix() &
%   [dpqr]gaitbinom.mlm().



