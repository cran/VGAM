\name{Gaitdbinom}
\alias{Gaitdbinom}
\alias{dgaitdbinom}
\alias{pgaitdbinom}
\alias{qgaitdbinom}
\alias{rgaitdbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated, --Truncated and --Deflated
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
GAAIITDD--Binom(size.p, prob.p)--Binom(size.a, prob.a)--MLM--Binom(size.i, prob.i)--MLM--Binom(size.d, prob.d)--MLM.



% and it is also known as the GAIT-Binom PNP combo where
% PNP stands for parametric and nonparametric.


}
\usage{
dgaitdbinom(x, size.p, prob.p,
            a.mix = NULL, a.mlm = NULL,
            i.mix = NULL, i.mlm = NULL,
            d.mix = NULL, d.mlm = NULL, truncate = NULL,
            pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0,
            pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
            size.a = size.p, size.i = size.p, size.d = size.p,
            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
            log = FALSE, ...)
pgaitdbinom(q, size.p, prob.p,
            a.mix = NULL, a.mlm = NULL,
            i.mix = NULL, i.mlm = NULL,
            d.mix = NULL, d.mlm = NULL, truncate = NULL,
            pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0,
            pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
            size.a = size.p, size.i = size.p, size.d = size.p,
            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
            lower.tail = TRUE, ...)
qgaitdbinom(p, size.p, prob.p,
            a.mix = NULL, a.mlm = NULL,
            i.mix = NULL, i.mlm = NULL,
            d.mix = NULL, d.mlm = NULL, truncate = NULL,
            pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0,
            pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
            size.a = size.p, size.i = size.p, size.d = size.p,
            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p, ...)
rgaitdbinom(n, size.p, prob.p,
            a.mix = NULL, a.mlm = NULL,
            i.mix = NULL, i.mlm = NULL,
            d.mix = NULL, d.mlm = NULL, truncate = NULL,
            pobs.mix = 0, pobs.mlm = 0,
            pstr.mix = 0, pstr.mlm = 0,
            pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
            size.a = size.p, size.i = size.p, size.d = size.p,
            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p, ...)
}
%           max.support = NULL,  20220201: decided to make this unavailable.
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link[stats]{Binomial}}.

  }

  \item{size.p, prob.p}{
  Same meaning as in \code{\link[stats]{Binomial}}.
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{size.a, prob.a}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{size.i, prob.i}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{size.d, prob.d}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{truncate}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.


% max.support



  }
  \item{a.mix, i.mix, d.mix}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  

  }

  \item{a.mlm, i.mlm, d.mlm}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.


  }
  
  \item{pstr.mix, pstr.mlm, byrow.aid}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  }
  
  \item{pobs.mix, pobs.mlm}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  
  
  }
  \item{pdip.mix, pdip.mlm}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  
  
  }
% \item{deflation}{
% See \code{\link[VGAM]{Gaitdpois}} for generic information.
% }
  \item{\dots}{
  Arguments such as \code{max.support} that are ignored.
  This will occur internally within \code{\link[VGAM]{dgaitdplot}}.

  
  
  }
}
\details{
  These functions for the GAITD binomial distribution
  are analogous to the GAITD Poisson,
  hence most details have been put in
  \code{\link[VGAM]{Gaitdpois}}.



}
%\section{Warning }{


%}
\value{
  \code{dgaitdbinom} gives the density,
  \code{pgaitdbinom} gives the distribution function,
  \code{qgaitdbinom} gives the quantile function, and
  \code{rgaitdbinom} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:Binomial]{dbinom}},
  \code{\link[stats:Binomial]{pbinom}},
  \code{\link[stats:Binomial]{qbinom}},
  \code{\link[stats:Binomial]{rbinom}}
  respectively.
  


}
%\references{

%Yee, T. W. and Ma, C. (2022).
%Generally--altered, --inflated, --truncated and --deflated
%regression, with application to heaped and seeped data.
%\emph{In preparation}.


%}

\author{ T. W. Yee. }
\note{
  Functions \code{\link[VGAMdata]{Posbinom}} have been moved
  to \pkg{VGAMdata}.
  It is better to use
  \code{dgaitdbinom(x, size, prob, truncate = 0)} instead of
  \code{dposbinom(x, size, prob)}, etc.

  

}



\seealso{
  \code{\link[VGAM]{Gaitdpois}},
  \code{\link[VGAM]{Gaitdnbinom}},
  \code{\link[VGAM]{multinomial}},
  \code{\link[VGAM]{Gaitdlog}},
  \code{\link[VGAM]{Gaitdzeta}}.


% \code{\link[VGAMsecret]{Gaitgenpois1}},
% \code{\link{gaitpoisson.mlm}},
% \code{\link{gaitpoisson}},
% \code{\link{Zapois}},
% \code{\link{Zipois}},
% \code{\link{Pospois}}
% \code{\link[stats:Poisson]{Poisson}};
  
%  \code{\link{Gaitlog.mix}}  and \code{\link{Gaitlog.mlm}},
%  \code{\link{Gaitdpois.mix}} and \code{\link{Gaitdpois.mlm}},
%  \code{\link{Gaitdbinom.mlm}},

% \code{\link{gaitpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},



}
\examples{ size <- 20
ivec <- c(6, 10); avec <- c(8, 11); prob <- 0.25; xgrid <- 0:25
tvec <- 14; pobs.a <- 0.05; pstr.i <- 0.15
dvec <- 5; pdip.mlm <- 0.05
(ddd <- dgaitdbinom(xgrid, size, prob.p = prob, prob.a = prob + 0.05,
   truncate = tvec, pobs.mix = pobs.a,
   pdip.mlm = pdip.mlm, d.mlm = dvec,
   pobs.mlm = pobs.a, a.mlm = avec,
   pstr.mix = pstr.i, i.mix = ivec))
\dontrun{ dgaitdplot(c(size, prob), ylab = "Probability",
   xlab = "x", pobs.mix = pobs.mix,
   pobs.mlm = pobs.a, a.mlm = avec, all.lwd = 3,
   pdip.mlm = pdip.mlm, d.mlm = dvec, fam = "binom",
   pstr.mix = pstr.i, i.mix = ivec, deflation = TRUE,
   main = "GAITD Combo PMF---Binomial Parent")   }
}
\keyword{distribution}





%plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
%              main = "GAIT Combo PMF---Binomial Parent")
%mylwd <- 1
%abline(v = avec, col = 'blue', lwd = mylwd)
%abline(v = ivec, col = 'purple', lwd = mylwd)
%abline(v = tvec, col = 'tan', lwd = mylwd)
%abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
%lines(xgrid, dbinom(xgrid, size, prob), col = 'gray', lty = 2)
%lines(xgrid, ddd, type = "h", col = "pink", lwd = 7)  # GAIT combo PMF
%points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16)  }



% 20200815; checked identical results to [dpqr]gaitdbinom.mix() &
%   [dpqr]gaitdbinom.mlm().



