\name{Gaitlog}
\alias{Gaitlog}
\alias{dgaitlog}
\alias{pgaitlog}
\alias{qgaitlog}
\alias{rgaitlog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated
  Logarithmic Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated and --truncated
  logarithmic distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.
  Altogether it can be abbreviated as
  GAAIIT--Log(shape.p)--Log(shape.a)--MLM--Log(shape.i)--MLM,
  and it is also known as the GAIT-Log PNP combo.


}
\usage{
dgaitlog(x, shape.p, alt.mix = NULL, alt.mlm = NULL,
         inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
         max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
         pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
         shape.a = shape.p, shape.i = shape.p, deflation = FALSE,
         log = FALSE)
pgaitlog(q, shape.p, alt.mix = NULL, alt.mlm = NULL,
         inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
         max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
         pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
         shape.a = shape.p, shape.i = shape.p, lower.tail = TRUE)
qgaitlog(p, shape.p, alt.mix = NULL, alt.mlm = NULL,
         inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
         max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
         pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
         shape.a = shape.p, shape.i = shape.p)
rgaitlog(n, shape.p, alt.mix = NULL, alt.mlm = NULL,
         inf.mix = NULL, inf.mlm = NULL, truncate = NULL,
         max.support = Inf, pobs.mix = 0, pobs.mlm = 0,
         pstr.mix = 0, pstr.mlm = 0, byrow.ai = FALSE,
         shape.a = shape.p, shape.i = shape.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log, lower.tail}{
  Same meaning as in \code{\link{dlog}}.

  }

  \item{shape.p, shape.a, shape.i}{
  Same meaning as \code{shape} for \code{\link{dlog}},
  i.e., for an ordinary logarithmic distribution.
  See \code{\link{Gaitpois}} for generic information.

  
 
  }
  \item{truncate, max.support}{
  See \code{\link{Gaitpois}} for generic information.
  

  }
  \item{alt.mix, inf.mix}{
  See \code{\link{Gaitpois}} for generic information.


  

  }

  \item{alt.mlm, inf.mlm}{
  See \code{\link{Gaitpois}} for generic information.



  }



  
  \item{pobs.mlm, pstr.mlm, byrow.ai}{
  See \code{\link{Gaitpois}} for generic information.


  }
  
  \item{pobs.mix, pstr.mix}{
  See \code{\link{Gaitpois}} for generic information.
  
  
  }
  \item{deflation}{
  See \code{\link{Gaitpois}} for generic information.


  
  }
}
\details{
  These functions for the logarithmic distribution
  are analogous to the Poisson,
  hence most details have been put in
  \code{\link[VGAM]{Gaitpois}}.
  These functions do what
  \code{\link[VGAMdata]{Oalog}},
  \code{\link[VGAMdata]{Oilog}},
  \code{\link[VGAMdata]{Otlog}}
  collectively did plus much more.

  

%In the notation of Yee and Ma (2020)
%these functions allow for the special cases:
%(i) GAIT--Log(\code{shape.p})--Log(\code{shape.a},
%\code{alt.mix}, \code{pobs.mix})--Log(\code{shape.i},
%\code{inf.mix}, \code{pstr.mix});
%(ii) GAIT--Log(\code{shape.p})--MLM(\code{alt.mlm},
%\code{pobs.mlm})--MLM(\code{inf.mlm}, \code{pstr.mlm}).
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
%  \code{rgaitlog.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitlog} gives the density,
  \code{pgaitlog} gives the distribution function,
  \code{qgaitlog} gives the quantile function, and
  \code{rgaitlog} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link{dlog}},
  \code{\link{plog}},
  \code{\link{qlog}},
  \code{\link{rlog}}
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
  See \code{\link{Gaitpois}} for general information also relevant
  to this parent distribution.
  

% Functions \code{\link[VGAMdata]{Posbinom}} have been moved
% to \pkg{VGAMdata}.
% It is better to use
% \code{dgaitbinom(x, size, prob, truncate = 0)} instead of
% \code{dposbinom(x, size, prob)}, etc.






}




% \code{\link{gaitpoisson.mlm}},


\seealso{
  \code{\link{gaitlog}},
  \code{\link{Gaitpois}},
  \code{\link{multinomial}},
  \code{\link[VGAMdata]{Oalog}},
  \code{\link[VGAMdata]{Oilog}},
  \code{\link[VGAMdata]{Otlog}}.



%  \code{\link{Gaitzeta}}.


}
\examples{
ivec <- c(2, 10); avec <- ivec + 4; shape <- 0.95; xgrid <- 0:29
tvec <- 15; max.support <- 25; pobs.a <- 0.10; pstr.i <- 0.15
(ddd <- dgaitlog(xgrid, shape, truncate = tvec,
   max.support = max.support, pobs.mix = pobs.a,
   alt.mix = avec, pstr.mix = pstr.i, inf.mix = ivec))
\dontrun{plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
              main = "GAIT PMF---Logarithmic Parent")
mylwd <- 0.5
abline(v = avec, col = 'blue', lwd = mylwd)
abline(v = ivec, col = 'purple', lwd = mylwd)
abline(v = tvec, col = 'tan', lwd = mylwd)
abline(v = max.support, col = 'magenta', lwd = mylwd)
abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
lines(xgrid, dlog(xgrid, shape), col = 'gray', lty = "dashed")  # f_{\pi}
lines(xgrid, ddd, type = "h", col = "pink", lwd = 3)  # GAIT PMF
points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16, col = 'tan', cex = 2)  }
}
\keyword{distribution}






