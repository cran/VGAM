\name{Gaitnbinom.mix}
\alias{Gaitnbinom.mix}
\alias{dgaitnbinom.mix}
\alias{pgaitnbinom.mix}
\alias{qgaitnbinom.mix}
\alias{rgaitnbinom.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered, -Inflated and -Truncated
  Negative Binomial Distribution
  (GAIT--NB--NB--NB mixture)
}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally-altered, -inflated and -truncated
  negative binomial (NB) distribution,
  based on mixtures of NB distributions
  having different support.
  This distribution is sometimes abbreviated as
  GAIT--NB--NB--NB.

  
}
\usage{
dgaitnbinom.mix(x, size.p, prob.p = NULL, munb.p = NULL,
   alter = NULL, inflate = NULL, truncate = NULL, max.support = Inf,
   pobs.a = 0, pstr.i = 0, size.a = size.p, size.i = size.p,
   prob.a = prob.p, prob.i = prob.p,
   munb.a = munb.p, munb.i = munb.p, log.arg = FALSE)
pgaitnbinom.mix(q, size.p, prob.p = NULL, munb.p = NULL,
   alter = NULL, inflate = NULL, truncate = NULL, max.support = Inf,
   pobs.a = 0, pstr.i = 0, size.a = size.p, size.i = size.p,
   prob.a = prob.p, prob.i = prob.p,
   munb.a = munb.p, munb.i = munb.p)
qgaitnbinom.mix(p, size.p, prob.p = NULL, munb.p = NULL,
   alter = NULL, inflate = NULL, truncate = NULL, max.support = Inf,
   pobs.a = 0, pstr.i = 0, size.a = size.p, size.i = size.p,
   prob.a = prob.p, prob.i = prob.p,
   munb.a = munb.p, munb.i = munb.p)
rgaitnbinom.mix(n, size.p, prob.p = NULL, munb.p = NULL,
   alter = NULL, inflate = NULL, truncate = NULL, max.support = Inf,
   pobs.a = 0, pstr.i = 0, size.a = size.p, size.i = size.p,
   prob.a = prob.p, prob.i = prob.p,
   munb.a = munb.p, munb.i = munb.p)
}
\arguments{
  \item{x, p, q, n, log.arg}{
  Same meaning as in \code{\link[stats]{NegBinomial}}.

  }

  \item{size.p, size.a, size.i}{
  Same meaning as \code{size }in \code{\link[stats]{NegBinomial}},
  i.e., for an ordinary NB distribution.
  The first is for the main \emph{p}arent (inner) distribution, and
  the outer distribution(s) (usually spikes) concern
  the \emph{a}ltered and/or \emph{i}nflated values.
  Short vectors are recycled.

  }
  \item{prob.p, prob.a, prob.i}{
    Same meaning as \code{prob} in \code{\link[stats]{NegBinomial}}.

  }
  \item{munb.p, munb.a, munb.i}{
    Same meaning as \code{mu} in \code{\link[stats]{NegBinomial}}.
    Only one of the \code{prob}-type and \code{munb}-type arguments
    should be specified.
    Short vectors are recycled.


  
  }
  \item{alter, inflate, truncate}{
    See \code{\link{Gaitnbinom.mlm}}.
    The order of precedence is the same, viz. truncation first,
    then altering, and finally inflation.
    If \code{alter} and \code{pobs.a} are both of unit length
    then the default probability mass function (PMF)
    evaluated at \code{alter} is \code{pobs.a}.




  
% , but both \code{alter} and \code{inflate} cannot be
% specified together.
  
% Must be sorted and have unique values only.
  }
  \item{pobs.a, pstr.i}{
    See \code{\link{Gaitpois.mix}}.


  }
  \item{max.support}{
    Same as \code{\link{Gaitnbinom.mlm}}.

  
  }
}
\details{
  These functions are an alternative to \code{\link{Gaitnbinom.mlm}},
  and are to \code{\link{Gaitnbinom.mlm}} what
  \code{\link{Gaitpois.mix}} are to \code{\link{Gaitpois.mlm}}.
  See \code{\link{Gaitpois.mix}} for parallel information.
  This distribution may be suitable for heaped data whose parent
  distribution is overdispersed with respect to the Poisson.





In the notation of Yee and Ma (2019)
these functions concern the GAIT-NB-NB-NB distribution.
For information on the GAIT-NB-MLM-MLM distribution
see \code{\link{Gaitnbinom.mlm}}.





}
%\section{Warning }{
% See \code{\link{rgaitnegbin}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgaitnegbin(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgaitnegbin(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitnbinom.mix} gives the density (PMF),
  \code{pgaitnbinom.mix} gives the distribution function,
  \code{qgaitnbinom.mix} gives the quantile function, and
  \code{rgaitnbinom.mix} generates random deviates.


}
%\references{
%Yee, T. W. and Ma, C. C. (2019)
%Generally-altered, -inflated and -truncated count regression,
%with application to heaped and seeped data.
%\emph{In preparation}.


%, \bold{3}, 15--41.


%}

\author{ T. W. Yee. }
\note{
  It is intended that the minimal defaults of these functions
  are equivalent to \code{\link[stats]{NegBinomial}},
  e.g.,
  \code{dgaitnbinom.mix(x, size, munb.p = munb)} and
  \code{dnbinom(x, size, mu = munb)} should be identical, and
  ditto for the other 3 functions.



}



\seealso{
  \code{\link{Gaitnbinom.mlm}}
  for the GAIT-NB-MLM-MLM distribution,
  \code{\link[stats]{NegBinomial}}.


}
\examples{
ivec <- c(5, 15, 10); avec <- ivec; size = 10; munb <- 10
max.support <- 20; pobs.a <- 0.35; xvec <- 0:max.support
# GAT-NB-NB mixture
(pmf.a <- dgaitnbinom.mix(xvec, size.p = size, munb.p = munb,
     max.support = max.support, pobs.a = pobs.a, alter = avec))
sum(pmf.a)  # Should be 1
\dontrun{
ind4 <- match(xvec, avec, nomatch = 0) > 0  # xvec \%in\% avec
plot(xvec[ ind4], pmf.a[ ind4], type = "h", col = "orange", lwd = 1.1,
   las = 1, xlim = range(xvec), main = "GAT-NB-NB",
   ylim = c(0, max(pmf.a)), xlab = "y", ylab = "Probability")  # Spikes
lines(xvec[!ind4], pmf.a[!ind4], type = "h", col = "blue") }

# GIT-NB-NB mixture
pstr.i <- 0.15
(pmf.i <- dgaitnbinom.mix(xvec, size.p = size, munb.p = munb,
   max.support = max.support, pstr.i = pstr.i, inflate = ivec))
sum(pmf.i)  # Should be 1
\dontrun{  # Plot the components of pmf.i
plot(xvec, (1 - pstr.i) * dnbinom(xvec, size, mu  = munb), type = "h",
     col = "blue", las = 1, xlim = range(xvec),
     main = "GIT-NB-NB",  # The inner distribution
     ylim = c(0, max(pmf.i)), xlab = "y", ylab = "Probability")
spikes <- dnbinom(ivec, size, mu = munb) * pstr.i / sum(
          dnbinom(ivec, size, mu = munb))
start.pt <- dnbinom(ivec, size, mu = munb) *
            (1 - pstr.i) / pnbinom(max.support, size, mu = munb)
segments(ivec, start.pt,  # The outer distribution
         ivec, start.pt + spikes, col = "orange", lwd = 1.1) }
}
\keyword{distribution}

% 20191026



