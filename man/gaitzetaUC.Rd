\name{Gaitzeta}
\alias{Gaitzeta}
\alias{dgaitzeta}
\alias{pgaitzeta}
\alias{qgaitzeta}
\alias{rgaitzeta}
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
  GAAIIT--Zeta(shape.p)--Zeta(shape.a)--MLM--Zeta(shape.i)--MLM,
  and it is also known as the GAIT-Zeta PNP combo.


}
\usage{
dgaitzeta(x, shape.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, shape.a = shape.p,
          shape.i = shape.p, byrow.arg = FALSE, deflation = FALSE,
          log.arg = FALSE)
pgaitzeta(q, shape.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, shape.a = shape.p,
          shape.i = shape.p, byrow.arg = FALSE)
qgaitzeta(p, shape.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, shape.a = shape.p,
          shape.i = shape.p, byrow.arg = FALSE)
rgaitzeta(n, shape.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, shape.a = shape.p,
          shape.i = shape.p, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log.arg}{
  Same meaning as in \code{\link{dzeta}}.

  }

  \item{shape.p, shape.a, shape.i}{
  Same meaning as \code{shape} for \code{\link{dzeta}},
  i.e., for an ordinary zeta distribution.
  See \code{\link{Gaitpois}} for generic information.

  
 
  }
  \item{truncate, max.support}{
  See \code{\link{Gaitpois}} for generic information.
  

  }
  \item{alter.mix, inflate.mix}{
  See \code{\link{Gaitpois}} for generic information.


  

  }

  \item{alter.mlm, inflate.mlm}{
  See \code{\link{Gaitpois}} for generic information.



  }



  
  \item{pobs.mlm.a, pstr.mlm.i, byrow.arg}{
  See \code{\link{Gaitpois}} for generic information.


  }
  
  \item{pobs.mix.a, pstr.mix.i}{
  See \code{\link{Gaitpois}} for generic information.
  
  
  }
  \item{deflation}{
  See \code{\link{Gaitpois}} for generic information.


  
  }
}
\details{
  These functions allow any combination of 3 operator types:
  truncation, alteration and inflation.
  See \code{\link{Gaitpois}} for generic information.
  These functions do what
  \code{\link{Oazeta}},
  \code{\link{Oizeta}},
  \code{\link{Otzeta}}
  collectively did plus much more.

  

In the notation of Yee and Ma (2020)
these functions allow for the special cases:
(i) GAIT--Zeta(\code{shape.p})--Zeta(\code{shape.a},
\code{alter.mix}, \code{pobs.mix.a})--Zeta(\code{shape.i},
\code{inflate.mix}, \code{pstr.mix.i});
(ii) GAIT--Zeta(\code{shape.p})--MLM(\code{alter.mlm},
\code{pobs.mlm.a})--MLM(\code{inflate.mlm}, \code{pstr.mlm.i}).
Model (i) is totally parametric while model (ii) is the most
nonparametric possible.




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
%  \code{rgaitzeta.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitzeta} gives the density,
  \code{pgaitzeta} gives the distribution function,
  \code{qgaitzeta} gives the quantile function, and
  \code{rgaitzeta} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link{dzeta}},
  \code{\link{pzeta}},
  \code{\link{qzeta}},
  \code{\link{rzeta}}
  respectively.
  


}
%\references{

%Yee, T. W. and Ma, C. (2020).
%Generally--altered, --inflated and --truncated regression,
%with application to heaped and seeped count data.
%\emph{In preparation}.


%}

\author{ T. W. Yee. }
\note{
  See \code{\link{Gaitpois}} for general information also relevant
  to this parent distribution.
  



}




% \code{\link{gaitpoisson.mlm}},


\seealso{
  \code{\link{Gaitpois}},
  \code{\link{multinomial}},
  \code{\link{gaitzeta.mix}},
  \code{\link{Oazeta}},
  \code{\link{Oizeta}},
  \code{\link{Otzeta}}.



%  \code{\link{Gaitzeta}}.


}
\examples{
ivec <- c(2, 10); avec <- ivec + 4; shape <- 0.95; xgrid <- 0:29
tvec <- 15; max.support <- 25; pobs.a <- 0.10; pstr.i <- 0.15
(ddd <- dgaitzeta(xgrid, shape, truncate = tvec,
   max.support = max.support, pobs.mix.a = pobs.a,
   alter.mix = avec, pstr.mix.i = pstr.i, inflate.mix = ivec))
\dontrun{plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
              main = "GAIT PMF---Zeta Parent")
mylwd <- 0.5
abline(v = avec, col = 'green', lwd = mylwd)
abline(v = ivec, col = 'red', lwd = mylwd)
abline(v = tvec, col = 'tan', lwd = mylwd)
abline(v = max.support, col = 'magenta', lwd = mylwd)
abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
lines(xgrid, dzeta(xgrid, shape), col = 'gray', lty = "dashed")  # f_{\pi}
lines(xgrid, ddd, type = "h", col = "blue", lwd = 3)  # GAIT PMF
points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16)  }
}
\keyword{distribution}






