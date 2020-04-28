\name{Gaitpois}
\alias{Gaitpois}
\alias{dgaitpois}
\alias{pgaitpois}
\alias{qgaitpois}
\alias{rgaitpois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated
  Poisson Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated and --truncated
  Poisson distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.
  Altogether it can be abbreviated as
  GAAIIT--Pois(lambda.p)--Pois(lambda.a)--MLM--Pois(lambda.i)--MLM,
  and it is also known as the GAIT-Pois PNP combo.


}
\usage{
dgaitpois(x, lambda.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, lambda.a = lambda.p,
          lambda.i = lambda.p, byrow.arg = FALSE, deflation = FALSE,
          log.arg = FALSE)
pgaitpois(q, lambda.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, lambda.a = lambda.p,
          lambda.i = lambda.p, byrow.arg = FALSE)
qgaitpois(p, lambda.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, lambda.a = lambda.p,
          lambda.i = lambda.p, byrow.arg = FALSE)
rgaitpois(n, lambda.p, alter.mix = NULL, alter.mlm = NULL,
          inflate.mix = NULL, inflate.mlm = NULL, truncate = NULL,
          max.support = Inf, pobs.mix.a = 0, pobs.mlm.a = 0,
          pstr.mix.i = 0, pstr.mlm.i = 0, lambda.a = lambda.p,
          lambda.i = lambda.p, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log.arg}{
  Same meaning as in \code{\link[stats]{rpois}}.

  }

  \item{lambda.p, lambda.a, lambda.i}{
  Same meaning as in \code{\link[stats]{rpois}},
  i.e., for an ordinary Poisson distribution.
  The first is for the main \emph{p}arent (inner) distribution.
  The other two concern the parametric variant and
  these outer distributions (usually spikes) may be
  \emph{a}ltered and/or \emph{i}nflated.
  Short vectors are recycled.

  
 
  }
  \item{truncate, max.support}{
    numeric; these specify the set of truncated values.
  The default value of \code{NULL} means an empty set for the former.
  The latter is the
    maximum support value so that any value larger
  has been truncated (necessary because
  \code{truncate = (max.support + 1):Inf} is not allowed),
  hence is needed for truncating the upper tail of the distribution.
  Note that \code{max(truncate) < max.support} must be satisfied
  otherwise an error message will be issued.


  

  }
  \item{alter.mix, inflate.mix}{
  Vectors of nonnegative integers;
  the altered, inflated and truncated values for the
  parametric variant.
  Each argument must have unique values only.
  Assigning argument \code{alter.mix}
  means that \code{pobs.mix.a} will be used.
  Assigning argument \code{inflate.mix}
  means that \code{pstr.mix.i} will be used.
    If \code{alter.mix} is of unit length
    then the default probability mass function (PMF)
    evaluated at \code{alter.mix} will be \code{pobs.mix.a}.
    So having \code{alter.mix = 0} corresponds to the
    zero-inflated Poisson distribution (see \code{\link{Zipois}}).


  
% Must be sorted and have unique values only.
  }

  \item{alter.mlm, inflate.mlm}{
  Similar to the above, but for the nonparametric (MLM) variant.
  Assigning argument \code{alter.mlm}
  means that \code{pobs.mlm.a} will be used.
  Assigning argument \code{inflate.mlm}
  means that \code{pstr.mlm.i} will be used.
  Collectively, the above 6 arguments represent
  5 disjoint sets of
  special values and they are a proper subset of the support of the
  distribution.


  }



  
  \item{pobs.mlm.a, pstr.mlm.i, byrow.arg}{
    These arguments are coerced into a matrix of probabilities
    using \code{byrow.arg} to determine the order of the elements
 (similar to \code{\link[base]{matrix}}).
  The first argument is recycled if necessary to become
  \code{n x length(alter.mlm)}.
  The second argument becomes
  \code{n x length(inflate.mlm)}.
  Thus these arguments are not used unless
  \code{alter.mlm} and \code{inflate.mlm} are assigned.
  If \code{deflation = TRUE} then \code{pstr.mlm.i} may be
  negative.



%%% This paragraph only holds if inflation is the sole operator:
% One can think of this matrix as comprising of
% \emph{structural} probabilities.
% Then the matrix augmented with one more column on the RHS so it has
% dimension \code{n x (length(inflate) + 1)}
% and whose \code{\link[base]{rowSums}} is a vector of 1s.
% Finally,
% for \code{\link{rgaitpois.mlm}},
% a multinomial sample is taken and if it belongs to the final
% column then Poisson random variates are drawn.



  }
  
  \item{pobs.mix.a, pstr.mix.i}{
    Vectors of probabilities that are recycled if necessary to
    length \eqn{n}.
  The first  argument is used when \code{alter.mix}   is not \code{NULL}.
  The second argument is used when \code{inflate.mix} is not \code{NULL}.
  
  
  }
  \item{deflation}{
  Logical. If \code{TRUE} then \code{pstr.mlm.i} is allowed to have
  negative values,
  however, not too negative so that the final PMF becomes negative.
  Of course, if the values are negative then they cannot be
  interpreted as probabilities.
  In theory, one could also allow \code{pstr.mix.i} to be negative,
  but currently this is disallowed.


  
  }
}
\details{
  These functions allow any combination of 3 operator types:
  truncation, alteration and inflation.
The precedence is
truncation, then alteration and lastly inflation.
This order minimizes the potential interference among the operators.
Loosely, a set of probabilities is set to 0 by truncation
and the remaining probabilities are scaled up.
Then a different set of probabilities are set to some
values \code{pobs.mix.a} and/or \code{pobs.mlm.a}
and the remaining probabilities are rescaled up.
Then another different set of probabilities is inflated by
an amount \code{pstr.mlm.i} and/or proportional
to \code{pstr.mix.i}
so that individual elements in this set have two sources.
Then all the probabilities are
rescaled so that they sum to unity.




Both parametric and nonparametric variants are implemented.
They usually have arguments with suffix
\code{.mix} and \code{.mlm} respectively.
The MLM is a loose coupling that effectively separates
the parent distribution from the altered values.
Values inflated nonparametrically effectively have
their spikes shaved off.
The \code{.mix} variant has associated with it
\code{lambda.a} and \code{lambda.i}
because it is mixture of 3 Poisson distributions with
partitioned or nested support.




Any value of the support of the distribution that is
altered, inflated or truncated is called a \emph{special} value.
A special value that is altered may mean that its probability
increases or decreases relative to the parent distribution.
An inflated special value means that its probability has
increased, provided alteration elsewhere has not made it decrease
in the first case.



  Jargonwise,
  the outer distributions concern those special values which
  are altered or inflated, and
  the inner distribution concerns the remaining
  support points that correspond directly to
  the parent distribution.
  These functions do what
  \code{\link{Zapois}},
  \code{\link{Zipois}},
  \code{\link{Pospois}}
  collectively did plus much more.

  

In the notation of Yee and Ma (2020)
these functions allow for the special cases:
(i) GAIT--Pois(\code{lambda.p})--Pois(\code{lambda.a},
\code{alter.mix}, \code{pobs.mix.a})--Pois(\code{lambda.i},
\code{inflate.mix}, \code{pstr.mix.i});
(ii) GAIT--Pois(\code{lambda.p})--MLM(\code{alter.mlm},
\code{pobs.mlm.a})--MLM(\code{inflate.mlm}, \code{pstr.mlm.i}).
Model (i) is totally parametric while model (ii) is the most
nonparametric possible.






%  In particular,
%  the GAT inner distribution has
%  a parent that is truncated at \code{c(alter, truncated)}
%  and any values beyond \code{max.support}.
%  Similarly,
%  the GIT inner distribution has
%  a parent that is truncated at \code{truncated}
%  and any values beyond \code{max.support}.


  
%Altogether, the full distribution is 
%  based on mixtures of Poisson distributions
%  having different and/or nested support, and
%  mixtures of the multinomial logit models
%  having different and/or nested support.






% It can be considered a mixture of a multinomial distribution
% and an ordinary Poisson distribution.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-Poisson and Poisson
%  distributions become more similar.

% Unlike similar functions for the Poisson distribution, a zero value
% of \code{lambda} is not permitted here.


% These functions merge
% \code{dgaitpois.mix} and \code{dgaitpois.mlm},
% \code{pgaitpois.mix} and \code{pgaitpois.mlm},
% \code{qgaitpois.mix} and \code{qgaitpois.mlm},
% \code{rgaitpois.mix} and \code{rgaitpois.mlm}.
% These functions do what
% \code{dgapois},
% \code{dgipois},
% \code{dgtpois},
% \code{pgapois},
% \code{pgipois},
% \code{pgtpois},
% \code{qgapois},
% \code{qgipois},
% \code{qgtpois},
% \code{rgapois},
% \code{rgipois},
% \code{rgtpois}
% collectively did.
% The old functions may be relegated to \pkg{VGAMdata} in the
% future.



  

% because the arguments \code{alter}, \code{inflate} and \code{truncate}
% have been included in one function.

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
%  \code{rgaitpois.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitpois} gives the density,
  \code{pgaitpois} gives the distribution function,
  \code{qgaitpois} gives the quantile function, and
  \code{rgaitpois} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:Poisson]{dpois}},
  \code{\link[stats:Poisson]{ppois}},
  \code{\link[stats:Poisson]{qpois}},
  \code{\link[stats:Poisson]{rpois}}
  respectively.
  


}
\references{

Yee, T. W. and Ma, C. (2020).
Generally--altered, --inflated and --truncated regression,
with application to heaped and seeped count data.
\emph{In preparation}.


}

\author{ T. W. Yee. }
%\note{
% This is repeated info:
%  The MLM variant is called nonparametric while the
%  parent mixture variant is called parametric.
%  Values that are altered nonparametrically are effectively
%  removed from the main analysis because the MLM provides a
%  loose coupling.




%}




% \code{\link{gaitpoisson.mlm}},
% \code{\link{gaitpoisson.mix}},



\seealso{
  \code{\link{multinomial}},
  \code{\link{gaitpoisson.mlm}},
  \code{\link{gaitpoisson.mix}},
  \code{\link{Zapois}},
  \code{\link{Zipois}},
  \code{\link{Pospois}}
  \code{\link[stats:Poisson]{Poisson}};
  \code{\link{Gaitlog}},
  \code{\link{Gaitzeta}}.


  
%  \code{\link{Gaitlog.mix}}  and \code{\link{Gaitlog.mlm}},
%  \code{\link{Gaitpois.mix}} and \code{\link{Gaitpois.mlm}},
%  \code{\link{Gaitnbinom.mlm}},

% \code{\link{gaitpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},



}
\examples{
ivec <- c(6, 14); avec <- c(8, 11); lambda <- 10; xgrid <- 0:25
tvec <- 15; max.support <- 20; pobs.a <- 0.05; pstr.i <- 0.25
(ddd <- dgaitpois(xgrid, lambda, lambda.a = lambda + 5,
   truncate = tvec, max.support = max.support, pobs.mix.a = pobs.a,
   pobs.mlm.a = pobs.a, alter.mlm = avec,
   pstr.mix.i = pstr.i, inflate.mix = ivec))
\dontrun{plot(xgrid, ddd, type = "n", ylab = "Probability", xlab = "x",
              main = "GAIT PNP Combo PMF---Poisson Parent")
mylwd <- 0.5
abline(v = avec, col = 'green', lwd = mylwd)
abline(v = ivec, col = 'red', lwd = mylwd)
abline(v = tvec, col = 'tan', lwd = mylwd)
abline(v = max.support, col = 'magenta', lwd = mylwd)
abline(h = c(pobs.a, pstr.i, 0:1), col = 'gray', lty = "dashed")
lines(xgrid, dpois(xgrid, lambda), col = 'gray', lty = "dashed")  # f_{\pi}
lines(xgrid, ddd, type = "h", col = "blue", lwd = 3)  # GAIT PNP combo PMF
points(xgrid[ddd == 0], ddd[ddd == 0], pch = 16)  }
}
\keyword{distribution}


% 20200317; checked identical results to gaitpois.mix() & gaitpois.mlm().



