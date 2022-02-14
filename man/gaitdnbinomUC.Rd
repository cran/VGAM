\name{Gaitdnbinom}
\alias{Gaitdnbinom}
\alias{dgaitdnbinom}
\alias{pgaitdnbinom}
\alias{qgaitdnbinom}
\alias{rgaitdnbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated, --Truncated and --Deflated
  Negative Binomial Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated, --truncated and --deflated
  negative binomial (GAITD-NB) distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.
  Altogether it can be abbreviated as
  GAAIITDD--NB(size.p, munb.p)--NB(size.a, munb.a)--MLM--NB(size.i,
  munb.i)--MLM--NB(size.d, munb.d)--MLM.



% and it is also known as the GAITD-NB PNP combo
% GAITD-NB PNP combo PNP stands for parametric and nonparametric.
% whereas the GAITD-NB PNP combo



}
\usage{
dgaitdnbinom(x, size.p, munb.p,
             a.mix = NULL, a.mlm = NULL,
             i.mix = NULL, i.mlm = NULL,
             d.mix = NULL, d.mlm = NULL, truncate = NULL,
             max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
             pstr.mix = 0, pstr.mlm = 0,
             pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
             size.a = size.p, size.i = size.p, size.d = size.p,
             munb.a = munb.p, munb.i = munb.p, munb.d = munb.p,
             log = FALSE)
pgaitdnbinom(q, size.p, munb.p,
             a.mix = NULL, a.mlm = NULL,
             i.mix = NULL, i.mlm = NULL,
             d.mix = NULL, d.mlm = NULL, truncate = NULL,
             max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
             pstr.mix = 0, pstr.mlm = 0,
             pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
             size.a = size.p, size.i = size.p, size.d = size.p,
             munb.a = munb.p, munb.i = munb.p, munb.d = munb.p,
             lower.tail = TRUE)
qgaitdnbinom(p, size.p, munb.p,
             a.mix = NULL, a.mlm = NULL,
             i.mix = NULL, i.mlm = NULL,
             d.mix = NULL, d.mlm = NULL, truncate = NULL,
             max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
             pstr.mix = 0, pstr.mlm = 0,
             pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
             size.a = size.p, size.i = size.p, size.d = size.p,
             munb.a = munb.p, munb.i = munb.p, munb.d = munb.p)
rgaitdnbinom(n, size.p, munb.p,
             a.mix = NULL, a.mlm = NULL,
             i.mix = NULL, i.mlm = NULL,
             d.mix = NULL, d.mlm = NULL, truncate = NULL,
             max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
             pstr.mix = 0, pstr.mlm = 0,
             pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
             size.a = size.p, size.i = size.p, size.d = size.p,
             munb.a = munb.p, munb.i = munb.p, munb.d = munb.p)
}
%- maybe also 'usage' for other objects documented here.
%            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
%            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
%            prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link[stats]{rnbinom}}.

  }

  \item{size.p, munb.p}{
  Same meaning as in \code{\link[stats]{rnbinom}}.
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{size.a, munb.a}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{size.i, munb.i}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{size.d, munb.d}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  
 
  }
  \item{truncate, max.support}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  }
  \item{a.mix, i.mix, d.mix}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.
  

  }

  \item{a.mlm, i.mlm, d.mlm}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.


  }
  
  \item{pobs.mlm, pstr.mlm, byrow.aid}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  }
  
  \item{pobs.mix, pstr.mix}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  
  
  }
  \item{pdip.mix, pdip.mlm}{
  See \code{\link[VGAM]{Gaitdpois}} for generic information.

  
  
  }
% \item{deflation}{
% See \code{\link[VGAM]{Gaitdpois}} for generic information.
% }
}
\details{
  These functions for the NBD are analogous to the Poisson,
  hence most details have been put in
  \code{\link[VGAM]{Gaitdpois}}.
  The NBD has two possible parameterizations: one
  involving a probability (argument begins with \code{prob})
  and the other the mean (beginning with \code{mu}).
  Only the latter is supported here.

  

% Because \code{\link[stats]{NegBinomial}} only allows
% one of these arguments to be used, the functions here
% have the same behaviour.



  For now, arguments such as \code{prob.p} and \code{prob.a}
  are no longer supported.
  That's because \code{mu} is more likely to be
  used by most statisticians than \code{prob};
  see \code{\link[stats:NegBinomial]{dnbinom}}.



}
%\section{Warning }{


%}
\value{
  \code{dgaitdnbinom} gives the density,
  \code{pgaitdnbinom} gives the distribution function,
  \code{qgaitdnbinom} gives the quantile function, and
  \code{rgaitdnbinom} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:NegBinomial]{dnbinom}},
  \code{\link[stats:NegBinomial]{pnbinom}},
  \code{\link[stats:NegBinomial]{qnbinom}},
  \code{\link[stats:NegBinomial]{rnbinom}}
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
  Four functions were moved from \pkg{VGAM} to \pkg{VGAMdata};
  they can be seen at \code{\link[VGAMdata]{Posnegbin}}.
  It is preferable to use
  \code{dgaitdnbinom(x, size, munb.p = munb, truncate = 0)} instead of
  \code{dposnbinom(x, size, munb = munb)}, etc.





}



\seealso{
  \code{\link[VGAM]{gaitdnbinomial}},
  \code{\link[VGAM]{Gaitdpois}},
  \code{\link[VGAM]{multinomial}},
  \code{\link[VGAM]{Gaitdbinom}},
  \code{\link[VGAM]{Gaitdlog}},
  \code{\link[VGAM]{Gaitdzeta}}.


% \code{\link[VGAMsecret]{Gaitgenpois1}},


% \code{\link{gaitpoisson}},
% \code{\link{Zapois}},
% \code{\link{Zipois}},
% \code{\link{Pospois}}
% \code{\link[stats:Poisson]{Poisson}};
  
%  \code{\link{Gaitlog.mix}}  and \code{\link{Gaitlog.mlm}},
%  \code{\link{Gaitdpois.mix}} and \code{\link{Gaitdpois.mlm}},
%  \code{\link{Gaitdnbinom.mlm}},

% \code{\link{gaitpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},



}
\examples{size <- 10
ivec <- c(5, 6, 10, 14); avec <- c(8, 11); munb <- 10; xgrid <- 0:25
tvec <- 15; pobs.a <- 0.05; pstr.i <- 0.25
dvec <- 13; pdip.mlm <- 0.03; pobs.mlm <- 0.05
(ddd <- dgaitdnbinom(xgrid, size, munb.p = munb, munb.a = munb + 5,
   truncate = tvec, pobs.mix = pobs.a,
   pdip.mlm = pdip.mlm, d.mlm = dvec,
   pobs.mlm = pobs.a, a.mlm = avec,
   pstr.mix = pstr.i, i.mix = ivec))
\dontrun{dgaitdplot(c(size, munb), fam = "nbinom",
  ylab = "Probability", xlab = "x", xlim = c(0, 25),
  truncate = tvec, pobs.mix = pobs.mix,
  pobs.mlm = pobs.mlm, a.mlm = avec, all.lwd = 3,
  pdip.mlm = pdip.mlm, d.mlm = dvec,
  pstr.mix = pstr.i, i.mix = ivec, deflation = TRUE,
  main = "GAITD Combo PMF---NB Parent")   }
}
\keyword{distribution}


% 20200815; checked identical results to [dpqr]gaitdnbinom.mix() &
%   [dpqr]gaitdnbinom.mlm().




% 20211109; the following call doesnt work. I tried to get it
%   going but its quite involved. Abandoning it.
%   GenPois woulb be much easier than NBD because of NBDs arguments.
%dgaitdplot(c(size, NA, munb), fam = "nbinom",
%  baseparams.argnames = c("size", "munb"),
%  ylab = "Probability", xlab = "x",
%  truncate = tvec, max.support = max.support, pobs.mix = pobs.mix,
%  pobs.mlm = pobs.mlm, a.mlm = avec, all.lwd = 3,
%  pdip.mlm = pdip.mlm, d.mlm = dvec,
%  pstr.mix = pstr.i, i.mix = ivec, deflation = TRUE,
%  main = "GAITD PNP Combo PMF---Poisson Parent") 

