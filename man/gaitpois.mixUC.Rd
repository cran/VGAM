\name{Gaitpois.mix}
\alias{Gaitpois.mix}
\alias{dgaitpois.mix}
\alias{pgaitpois.mix}
\alias{qgaitpois.mix}
\alias{rgaitpois.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered, -Inflated and -Truncated
  Poisson Distribution
  (GAIT--Pois--Pois--Pois mixture)
}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally-altered, -inflated and -truncated
  Poisson distribution,
  based on mixtures of Poisson distributions
  having different support.
  This distribution is sometimes abbreviated as
  GAIT--Pois--Pois--Pois.

}
\usage{
dgaitpois.mix(x, lambda.p, alter = NULL, inflate = NULL,
    truncate = NULL, max.support = Inf, pobs.a = 0, pstr.i = 0,
    lambda.a = lambda.p, lambda.i = lambda.p, log.arg = FALSE)
pgaitpois.mix(q, lambda.p, alter = NULL, inflate = NULL,
    truncate = NULL, max.support = Inf, pobs.a = 0, pstr.i = 0,
    lambda.a = lambda.p, lambda.i = lambda.p)
qgaitpois.mix(p, lambda.p, alter = NULL, inflate = NULL,
    truncate = NULL, max.support = Inf, pobs.a = 0, pstr.i = 0,
    lambda.a = lambda.p, lambda.i = lambda.p)
rgaitpois.mix(n, lambda.p, alter = NULL, inflate = NULL,
    truncate = NULL, max.support = Inf, pobs.a = 0, pstr.i = 0,
    lambda.a = lambda.p, lambda.i = lambda.p)
}
\arguments{
  \item{x, p, q, n, log.arg}{
  Same meaning as in \code{\link[stats]{Poisson}}.

  }

  \item{lambda.p, lambda.a, lambda.i}{
  Same meaning as in \code{\link[stats]{rpois}},
  i.e., for an ordinary Poisson distribution.
  The first is for the main \emph{p}arent (inner) distribution, and
  the outer distribution(s) (usually spikes) concern
  the \emph{a}ltered and/or \emph{i}nflated values.
  Short vectors are recycled.

  }
  \item{alter, inflate, truncate}{
    See \code{\link{Gaitpois.mix}}.
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
    Numeric, probabilities of an observed altered value
    or a structural inflated value.
    See \code{\link{Gaitpois.mlm}}, however these are
    ordinary vectors of the usual length
    rather than matrices.
    The probability of obtaining a value equal to
    an element of \code{alter} is \code{pobs.a}, while
    the probability of obtaining a value equal to
    an element of \code{inflate} \emph{structurally}
    is \code{pstr.i}. Actually, the probability of
    obtaining \code{inflate[i]} comes from a scaled Poisson
    distribution:
\code{pstr.i * dpois(inflate[i], lambda.p) / sum(dpois(inflate, lambda.p))}
    so that it isn't really a structural probability unless
    \code{length(inflate) == 1}.


    These two arguments are used only if \code{alter}
    and \code{inflate} are assigned values, respectively.


  }
  \item{max.support}{
    Same as \code{\link{Gaitpois.mlm}}.

  
  }
}
\details{
  These functions are an alternative to \code{\link{Gaitpois.mlm}}.
  The latter have an outer distribution that is
  based on the multinomial logit model (MLM)
  which allows the special values to be modelled nonparametrically
  in an unstructured manner.
  These functions here use the same parent distribution to
  model the special values, hence the special values are
  modelled in a more structured way; one could say they
  are more parametric than the MLM variant.  



  Jargonwise,
  the outer distribution concerns those special values which
  appear in \code{alter} or \code{inflate}, and
  the inner distribution concerns the remaining
  support points.
  General truncation is allowed throughout;
  the parent distribution has PMF that
  is merely scaled up to exclude those points
  (\code{truncate} and values beyond \code{max.support}).
  In particular,
  the GAT inner distribution has
  a parent that is truncated at \code{c(alter, truncated)}
  and any values beyond \code{max.support}.
  Similarly,
  the GIT inner distribution has
  a parent that is truncated at \code{truncated}
  and any values beyond \code{max.support}.


In the notation of Yee and Ma (2019)
these functions concern the GAIT-Pois-Pois-Pois distribution.
For the GAIT-Pois-MLM-MLM distribution
see \code{\link{Gaitpois.mlm}}.





}
%\section{Warning }{
% See \code{\link{rgaitpois}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgaitpois(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgaitpois(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitpois.mix} gives the density (PMF),
  \code{pgaitpois.mix} gives the distribution function,
  \code{qgaitpois.mix} gives the quantile function, and
  \code{rgaitpois.mix} generates random deviates.


}
\references{
Yee, T. W. and Ma, C. C. (2019)
Generally-altered, -inflated and -truncated count regression,
with application to heaped and seeped data.
\emph{In preparation}.


%, \bold{3}, 15--41.


}

\author{ T. W. Yee. }
%\note{
%}



\seealso{
  \code{\link{gatpoisson.mix}},
  \code{\link{Gaitpois.mlm}}
  for the GAIT-Pois-MLM-MLM distribution.


}
\examples{
ivec <- c(15, 10, 5); avec <- ivec; lambda <- 10
max.support <- 20; pobs.a <- 0.35; xvec <- 0:max.support
(pmf.a <- dgaitpois.mix(xvec, lambda,  # lambda.a = lambda,
     max.support = max.support, pobs.a = pobs.a, alter = avec))
sum(pmf.a)  # Should be 1
\dontrun{
ind4 <- match(xvec, avec, nomatch = 0) > 0  # xvec \%in\% avec
plot(xvec[ ind4], pmf.a[ ind4], type = "h", col = "orange", lwd = 1.1,
   las = 1, xlim = range(xvec), main = "GAT-Poisson-Poisson",
   ylim = c(0, max(pmf.a)), xlab = "y", ylab = "Probability")  # Spikes
lines(xvec[!ind4], pmf.a[!ind4], type = "h", col = "blue") }

# GIT-Poisson-Poisson mixture
pstr.i <- 0.20
(pmf.i <- dgaitpois.mix(xvec, lambda,  # lambda.a = lambda,
   max.support = max.support, pstr.i = pstr.i, inflate = ivec))
sum(pmf.i)  # Should be 1
\dontrun{  # Plot the components of pmf.i
spikes <- dpois(ivec, lambda) * pstr.i / sum(dpois(ivec, lambda))
start.pt <- dpois(ivec, lambda) * (1 - pstr.i) / ppois(max.support, lambda)
plot(xvec, (1 - pstr.i) * dpois(xvec, lambda), type = "h",
     col = "blue", las = 1, xlim = range(xvec),
     main = "GIT-Poisson-Poisson",  # The inner distribution
     ylim = c(0, max(pmf.i)), xlab = "y", ylab = "Probability")
segments(ivec, start.pt,  # The outer distribution
         ivec, start.pt + spikes, col = "orange", lwd = 1.1) }
}
\keyword{distribution}

% 20191019;
% \item{mixprob}{
%   Numeric, mixing probabilities.
% The parent distribution PMF is multiplied by
% \code{1-mixprob}, and then added to the product of
% \code{mixprob} and the normalized parent distribution PMF
% evaluated at the special values.
% Thus \code{mixprob} more directly modifies the size of the
% spikes and \code{1-mixprob} concerns more the inner parent
% distribution.
% See the examples figure below.
% }


