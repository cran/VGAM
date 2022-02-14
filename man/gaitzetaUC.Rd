\name{Gaitdzeta}
\alias{Gaitdzeta}
\alias{dgaitdzeta}
\alias{pgaitdzeta}
\alias{qgaitdzeta}
\alias{rgaitdzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated
  Zeta Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated and --truncated
  zeta distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.
  Altogether it can be abbreviated as
  GAAIIT--Zeta(shape.p)--Zeta(shape.a)--MLM--Zeta(shape.i)--MLM.



% and it is also known as the GAIT-Zeta PNP combo.


}
\usage{
dgaitdzeta(x, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0,
          byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
          log = FALSE)
pgaitdzeta(q, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0,
          byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
          lower.tail = TRUE)
qgaitdzeta(p, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0,
          byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p)
rgaitdzeta(n, shape.p, a.mix = NULL, a.mlm = NULL,
          i.mix = NULL, i.mlm = NULL,
          d.mix = NULL, d.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
          pstr.mix = 0, pstr.mlm = 0,
          pdip.mix = 0, pdip.mlm = 0,
          byrow.aid = FALSE,
          shape.a = shape.p, shape.i = shape.p, shape.d = shape.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link{dzeta}}.

  }

  \item{shape.p, shape.a, shape.i, shape.d}{
  Same meaning as \code{shape} for \code{\link{dzeta}},
  i.e., for an ordinary zeta distribution.
  See \code{\link{Gaitdpois}} for generic information.

  
 
  }
  \item{truncate, max.support}{
  See \code{\link{Gaitdpois}} for generic information.
  

  }
  \item{a.mix, i.mix, d.mix}{
  See \code{\link{Gaitdpois}} for generic information.


  

  }

  \item{a.mlm, i.mlm, d.mlm}{
  See \code{\link{Gaitdpois}} for generic information.



  }



  
  \item{pobs.mlm, pstr.mlm, pdip.mlm, byrow.aid}{
  See \code{\link{Gaitdpois}} for generic information.


  }
  
  \item{pobs.mix, pstr.mix, pdip.mix}{
  See \code{\link{Gaitdpois}} for generic information.
  
  
  }
% \item{deflation}{
% See \code{\link{Gaitdpois}} for generic information.
% }
}
\details{
  These functions for the zeta distribution are analogous to
  the Poisson, hence most details have been put in
  \code{\link[VGAM]{Gaitdpois}}.
  These functions do what
  \code{\link[VGAMdata]{Oazeta}},
  \code{\link[VGAMdata]{Oizeta}},
  \code{\link[VGAMdata]{Otzeta}}
  collectively did plus much more.

  

%In the notation of Yee and Ma (2022)
%these functions allow for the special cases:
%(i) GAIT--Zeta(\code{shape.p})--Zeta(\code{shape.a},
%\code{a.mix}, \code{pobs.mix})--Zeta(\code{shape.i},
%\code{i.mix}, \code{pstr.mix});
%(ii) GAIT--Zeta(\code{shape.p})--MLM(\code{a.mlm},
%\code{pobs.mlm})--MLM(\code{i.mlm}, \code{pstr.mlm}).
%Model (i) is totally parametric while model (ii) is the most
%nonparametric possible.




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
%  \code{rgaitdzeta.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitdzeta} gives the density,
  \code{pgaitdzeta} gives the distribution function,
  \code{qgaitdzeta} gives the quantile function, and
  \code{rgaitdzeta} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link{dzeta}},
  \code{\link{pzeta}},
  \code{\link{qzeta}},
  \code{\link{rzeta}}
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
  See \code{\link{Gaitdpois}} for general information also relevant
  to this parent distribution.
  

% Functions \code{\link[VGAMdata]{Poszeta}} have been moved
% to \pkg{VGAMdata}.
% It is better to use
% \code{dgaitdzeta(x, shape, truncate = 0)} instead of
% \code{dposzeta(x, shape)}, etc.





}




% \code{\link{gaitpoisson.mlm}},


\seealso{
  \code{\link{gaitdzeta}},
  \code{\link{Gaitdpois}},
  \code{\link{multinomial}},
  \code{\link[VGAMdata]{Oazeta}},
  \code{\link[VGAMdata]{Oizeta}},
  \code{\link[VGAMdata]{Otzeta}}.



%  \code{\link{Gaitdzeta}}.


}
\examples{
ivec <- c(2, 10); avec <- ivec + 4; shape <- 0.95; xgrid <- 0:29
tvec <- 15; max.support <- 25; pobs.a <- 0.10; pstr.i <- 0.15
(ddd <- dgaitdzeta(xgrid, shape, truncate = tvec,
   max.support = max.support, pobs.mix = pobs.a,
   a.mix = avec, pstr.mix = pstr.i, i.mix = ivec))
\dontrun{plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
              main = "GAIT PMF---Zeta Parent")
mylwd <- 0.5
abline(v = avec, col = 'blue', lwd = mylwd)
abline(v = ivec, col = 'purple', lwd = mylwd)
abline(v = tvec, col = 'tan', lwd = mylwd)
abline(v = max.support, col = 'magenta', lwd = mylwd)
abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
lines(xgrid, dzeta(xgrid, shape), col = 'gray', lty = "dashed")  # f_{\pi}
lines(xgrid, ddd, type = "h", col = "pink", lwd = 3)  # GAIT PMF
points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16, col = 'tan', cex = 2)  }
}
\keyword{distribution}






