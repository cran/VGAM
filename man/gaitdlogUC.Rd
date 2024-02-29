\name{Gaitdlog}
\alias{Gaitdlog}
\alias{dgaitdlog}
\alias{pgaitdlog}
\alias{qgaitdlog}
\alias{rgaitdlog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally Altered, Inflated, Truncated
  and Deflated
  Logarithmic Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally altered, inflated, truncated
  and deflated
  logarithmic distribution.
  Both parametric and nonparametric variants are
  supported; these are based on finite mixtures
  of the parent with itself and the multinomial
  logit model (MLM) respectively.



% Altogether it can be abbreviated as
% GAAIITDD--Log(shape.p)--Log(shape.a)--MLM--
% Log(shape.i)--MLM--Log(shape.d)--MLM.



% and it is also known as the GAITD-Log PNP combo where
% PNP stands for parametric and nonparametric.


}
\usage{
dgaitdlog(x, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
          log = FALSE)
pgaitdlog(q, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
          lower.tail = TRUE)
qgaitdlog(p, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p)
rgaitdlog(n, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link{dlog}}.

  }

  \item{shape.p, shape.a, shape.i, shape.d}{
  Same meaning as \code{shape} for
  \code{\link{dlog}}, i.e., for an
  ordinary logarithmic distribution.
  See \code{\link{Gaitdpois}} for generic
  information.

  
 
  }
  \item{truncate, max.support}{
  See \code{\link{Gaitdpois}} for generic
  information.


  }
  \item{a.mix, i.mix, d.mix}{
  See \code{\link{Gaitdpois}} for generic
  information.


  

  }

  \item{a.mlm, i.mlm, d.mlm}{
  See \code{\link{Gaitdpois}} for generic
  information.



  }



  
  \item{pobs.mlm, pstr.mlm, pdip.mlm, byrow.aid}{
  See \code{\link{Gaitdpois}} for generic
  information.


  }
  
  \item{pobs.mix, pstr.mix, pdip.mix}{
  See \code{\link{Gaitdpois}} for generic information.
  
  
  }
% \item{deflation}{
% See \code{\link{Gaitdpois}} for generic information.
% }
}
\details{
  These functions for the logarithmic distribution
  are analogous to the Poisson,
  hence most details have been put in
  \code{\link[VGAM]{Gaitdpois}}.
  These functions do what
  \code{\link[VGAMdata]{Oalog}},
  \code{\link[VGAMdata]{Oilog}},
  \code{\link[VGAMdata]{Otlog}}
  collectively did plus much more.

  

%In the notation of Yee and Ma (2023)
%these functions allow for the special cases:
%(i) GAIT--Log(\code{shape.p})--Log(\code{shape.a},
%\code{a.mix}, \code{pobs.mix})--Log(\code{shape.i},
%\code{i.mix}, \code{pstr.mix});
%(ii) GAIT--Log(\code{shape.p})--MLM(\code{a.mlm},
%\code{pobs.mlm})--MLM(\code{i.mlm}, \code{pstr.mlm}).
%Model (i) is totally parametric while model (ii) is the most
%nonparametric possible.




}
\section{Warning }{
  See \code{\link[VGAM]{Gaitdpois}} about the dangers
  of too much inflation and/or deflation on
  GAITD PMFs, and the difficulties detecting such.



}
%\section{Warning }{
% See \code{\link{rgaitpois.mlm}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgaitpois.mlm(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgaitdlog.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitdlog} gives the density,
  \code{pgaitdlog} gives the distribution function,
  \code{qgaitdlog} gives the quantile function, and
  \code{rgaitdlog} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link{dlog}},
  \code{\link{plog}},
  \code{\link{qlog}},
  \code{\link{rlog}}
  respectively.
  


}
%\references{

%Yee, T. W. and Ma, C. (2024).
%Generally altered, inflated, truncated and deflated regression.
%\emph{Statistical Science}, \bold{39} (in press).


%}

\author{ T. W. Yee. }
\note{
  See \code{\link{Gaitdpois}} for general information also relevant
  to this parent distribution.
  

% Functions \code{\link[VGAMdata]{Posbinom}} have been moved
% to \pkg{VGAMdata}.
% It is better to use
% \code{dgaitdbinom(x, size, prob, truncate = 0)} instead of
% \code{dposbinom(x, size, prob)}, etc.






}




% \code{\link{gaitpoisson.mlm}},


\seealso{
  \code{\link{gaitdlog}},
  \code{\link{Gaitdpois}},
  \code{\link{dgaitdplot}},
  \code{\link{Gaitdzeta}},
  \code{\link{multinomial}},
  \code{\link[VGAMdata]{Oalog}},
  \code{\link[VGAMdata]{Oilog}},
  \code{\link[VGAMdata]{Otlog}}.





}
\examples{
ivec <- c(2, 10); avec <- ivec + 1; shape <- 0.995; xgrid <- 0:15
max.support <- 15; pobs.a <- 0.10; pstr.i <- 0.15
dvec <- 1; pdip.mlm <- 0.05
(ddd <- dgaitdlog(xgrid, shape,
   max.support = max.support, pobs.mix = pobs.a,
   pdip.mlm = pdip.mlm, d.mlm = dvec,
   a.mix = avec, pstr.mix = pstr.i, i.mix = ivec))
\dontrun{ dgaitdplot(shape, ylab = "Probability", xlab = "x",
   max.support = max.support, pobs.mix = 0,
   pobs.mlm = 0, a.mlm = avec, all.lwd = 3,
   pdip.mlm = pdip.mlm, d.mlm = dvec, fam = "log",
   pstr.mix = pstr.i, i.mix = ivec, deflation = TRUE,
   main = "GAITD Combo PMF---Logarithmic Parent")   }
}
\keyword{distribution}






