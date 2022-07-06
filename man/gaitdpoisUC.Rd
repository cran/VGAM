\name{Gaitdpois}
\alias{Gaitdpois}
\alias{dgaitdpois}
\alias{pgaitdpois}
\alias{qgaitdpois}
\alias{rgaitdpois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated, --Truncated
  and --Deflated Poisson Distribution



}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally--altered, --inflated, --truncated
  and --deflated
  Poisson distribution.
  Both parametric and nonparametric variants are supported;
  these are based on
  finite mixtures of the parent with itself
  and the multinomial logit model (MLM) respectively.



% Altogether it can be abbreviated as
% GAAIITDD--Pois(lambda.p)--Pois(lambda.a)--MLM--
% Pois(lambda.i)--MLM--Pois(lambda.d)--MLM.



% and it is also known as the GAITD-Pois PNP combo where
% PNP stands for parametric and nonparametric.


}
\usage{
dgaitdpois(x, lambda.p, a.mix = NULL, a.mlm = NULL, i.mix = NULL,
       i.mlm = NULL, d.mix = NULL, d.mlm = NULL, truncate = NULL,
       max.support = Inf, pobs.mix = 0, pobs.mlm = 0, pstr.mix = 0,
       pstr.mlm = 0, pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
       lambda.a = lambda.p, lambda.i = lambda.p,
       lambda.d = lambda.p, log = FALSE)
pgaitdpois(q, lambda.p, a.mix = NULL, a.mlm = NULL, i.mix = NULL,
       i.mlm = NULL, d.mix = NULL, d.mlm = NULL, truncate = NULL,
       max.support = Inf, pobs.mix = 0, pobs.mlm = 0, pstr.mix = 0,
       pstr.mlm = 0, pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
       lambda.a = lambda.p, lambda.i = lambda.p,
       lambda.d = lambda.p, lower.tail = TRUE)
qgaitdpois(p, lambda.p, a.mix = NULL, a.mlm = NULL, i.mix = NULL,
       i.mlm = NULL, d.mix = NULL, d.mlm = NULL, truncate = NULL,
       max.support = Inf, pobs.mix = 0, pobs.mlm = 0, pstr.mix = 0,
       pstr.mlm = 0, pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
       lambda.a = lambda.p, lambda.i = lambda.p, lambda.d = lambda.p)
rgaitdpois(n, lambda.p, a.mix = NULL, a.mlm = NULL, i.mix = NULL,
       i.mlm = NULL, d.mix = NULL, d.mlm = NULL, truncate = NULL,
       max.support = Inf, pobs.mix = 0, pobs.mlm = 0, pstr.mix = 0,
       pstr.mlm = 0, pdip.mix = 0, pdip.mlm = 0, byrow.aid = FALSE,
       lambda.a = lambda.p, lambda.i = lambda.p, lambda.d = lambda.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same meaning as in \code{\link[stats]{Poisson}}.

  }
  \item{log, lower.tail}{
  Same meaning as in \code{\link[stats]{Poisson}}.


  }

  \item{lambda.p, lambda.a, lambda.i, lambda.d}{
  Same meaning as in \code{\link[stats]{Poisson}},
  i.e., for an ordinary Poisson distribution.
  The first is for the main \emph{p}arent (or base) distribution.
  The next two concern the parametric variant and
  these distributions (usually spikes) may be
  \emph{a}ltered and/or \emph{i}nflated.
  The last one concerns the \emph{d}eflated variant.
  Short vectors are recycled.

  
 
  }
  \item{truncate, max.support}{
    numeric; these specify the set of truncated values.
    The default value of \code{NULL} means an empty set
    for the former.
  The latter is the
    maximum support value so that any value larger
  has been truncated (necessary because
  \code{truncate = (max.support + 1):Inf} is not allowed),
  hence is needed for truncating the upper tail of the
  distribution.  Note that \code{max(truncate) < max.support}
  must be satisfied otherwise an error message will be issued.


  

  }
  \item{a.mix, i.mix, d.mix}{
  Vectors of nonnegative integers;
  the altered, inflated and deflated values for the
  parametric variant.
  Each argument must have unique values only.
  Assigning argument \code{a.mix}
  means that \code{pobs.mix} will be used.
  Assigning \code{i.mix}
  means that \code{pstr.mix} will be used.
  Assigning \code{d.mix}
  means that \code{pdip.mix} will be used.
    If \code{a.mix} is of unit length
    then the default probability mass function (PMF)
    evaluated at \code{a.mix} will be \code{pobs.mix}.
    So having \code{a.mix = 0} corresponds to the
    zero-inflated Poisson distribution (see \code{\link{Zipois}}).


  
% Must be sorted and have unique values only.
  }

  \item{a.mlm, i.mlm, d.mlm}{
  Similar to the above, but for the nonparametric (MLM) variant.
  For example, assigning \code{a.mlm}
  means that \code{pobs.mlm} will be used.
  Collectively, the above 7 arguments represent
  7 disjoint sets of
  special values and they are a proper subset of the support
  of the distribution.


  }



  
  \item{pobs.mlm, pstr.mlm, pdip.mlm, byrow.aid}{
    The first three arguments are coerced into a matrix of
    probabilities
    using \code{byrow.aid} to determine the order of the elements
 (similar to \code{byrow} in \code{\link[base]{matrix}}, and
  the \code{.aid} reinforces the behaviour that it applies to both
  altered, inflated and deflated cases).
  The first argument is recycled if necessary to become
  \code{n x length(a.mlm)}.
  The second argument becomes
  \code{n x length(i.mlm)}.
  The third argument becomes
  \code{n x length(d.mlm)}.
  Thus these arguments are not used unless
  \code{a.mlm}, \code{i.mlm} and \code{d.mlm} are assigned.
  For deflated models, \code{pdip.mix} and \code{pdip.mlm} are
  positive-valued and \pkg{VGAM} will subtract these quantities;
  the argument \code{deflation} has been deprecated.

  

% If \code{deflation = TRUE} then \code{pstr.mlm} may be negative.




%%% This paragraph only holds if inflation is the sole operator:
% One can think of this matrix as comprising of
% \emph{structural} probabilities.
% Then the matrix augmented with one more column on
% the RHS so it has
% dimension \code{n x (length(inflate) + 1)}
% and whose \code{\link[base]{rowSums}} is a vector of 1s.
% Finally,
% for \code{\link{rgaitdpois.mlm}},
% a multinomial sample is taken and if it belongs to the final
% column then Poisson random variates are drawn.



  }
  
  \item{pobs.mix, pstr.mix, pdip.mix}{
    Vectors of probabilities that are recycled if necessary to
    length \eqn{n}.
  The first  argument is used when \code{a.mix} is not \code{NULL}.
  The second argument is used when \code{i.mix} is not \code{NULL}.
  The third  argument is used when \code{d.mix} is not \code{NULL}.
  
  
  }
% \item{deflation}{
% Logical. If \code{TRUE} then \code{pstr.mlm} is allowed to have
% negative values,
% however, not too negative so that the final PMF becomes negative.
% Of course, if the values are negative then they cannot be
% interpreted as probabilities.
% In theory, one could also allow \code{pstr.mix} to be negative,
% but currently this is disallowed.
% }
}
\details{
  These functions allow any combination of 4 operator types:
  truncation, alteration, inflation and deflation.
The precedence is
truncation, then alteration and lastly inflation and deflation.
Informally, deflation can be thought of as the
opposite of inflation.
This order minimizes the potential interference among the operators.
Loosely, a set of probabilities is set to 0 by truncation
and the remaining probabilities are scaled up.
Then a different set of probabilities are set to some
values \code{pobs.mix} and/or \code{pobs.mlm}
and the remaining probabilities are rescaled up.
Then another different set of probabilities is inflated by
an amount \code{pstr.mlm} and/or proportional
to \code{pstr.mix}
so that individual elements in this set have two sources.
Then another different set of probabilities is deflated by
an amount \code{pdip.mlm} and/or proportional
to \code{pdip.mix}.
Then all the probabilities are
rescaled so that they sum to unity.




Both parametric and nonparametric variants are implemented.
They usually have arguments with suffix
\code{.mix} and \code{.mlm} respectively.
The MLM is a loose coupling that effectively separates
the \emph{parent} (or \emph{base}) distribution from
the altered values.
Values inflated nonparametrically effectively have
their spikes shaved off.
The \code{.mix} variant has associated with it
\code{lambda.a} and \code{lambda.i} and \code{lambda.d}
because it is mixture of 4 Poisson distributions with
partitioned or nested support.




Any value of the support of the distribution that is
altered, inflated, truncated or deflated
is called a \emph{special} value.
A special value that is altered may mean that its probability
increases or decreases relative to the parent distribution.
An inflated special value means that its probability has
increased, provided alteration elsewhere has not made it decrease
in the first case.
There are seven types of special values and they are
represented by
\code{a.mix},
\code{a.mlm},
\code{i.mix},
\code{i.mlm},
\code{d.mix},
\code{d.mlm},
\code{truncate}.




  Terminology-wise, \emph{special} values
  are altered or inflated or truncated or deflated, and
  the remaining support points that correspond directly to
  the parent distribution are \emph{nonspecial} or ordinary.
  These functions do what
  \code{\link{Zapois}},
  \code{\link{Zipois}},
  \code{\link[VGAMdata]{Pospois}}
  collectively did plus much more.

  

In the notation of Yee and Ma (2022)
these functions allow for the special cases:
(i) GAIT--Pois(\code{lambda.p})--Pois(\code{lambda.a},
\code{a.mix}, \code{pobs.mix})--Pois(\code{lambda.i},
\code{i.mix}, \code{pstr.mix});
(ii) GAIT--Pois(\code{lambda.p})--MLM(\code{a.mlm},
\code{pobs.mlm})--MLM(\code{i.mlm}, \code{pstr.mlm}).
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
%  As \eqn{\lambda}{lambda} increases,
%  the positive-Poisson and Poisson
%  distributions become more similar.

% Unlike similar functions for the Poisson distribution,
% a zero value of \code{lambda} is not permitted here.


% These functions merge
% \code{dgaitdpois.mix} and \code{dgaitdpois.mlm},
% \code{pgaitdpois.mix} and \code{pgaitdpois.mlm},
% \code{qgaitdpois.mix} and \code{qgaitdpois.mlm},
% \code{rgaitdpois.mix} and \code{rgaitdpois.mlm}.
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



  

% because the arguments \code{alter}, \code{inflate} and
% \code{truncate}
% have been included in one function.

}
%\section{Warning }{
% See \code{\link{rgaitdpois.mlm}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgaitdpois.mlm(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgaitdpois.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitdpois} gives the density,
  \code{pgaitdpois} gives the distribution function,
  \code{qgaitdpois} gives the quantile function, and
  \code{rgaitdpois} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:Poisson]{dpois}},
  \code{\link[stats:Poisson]{ppois}},
  \code{\link[stats:Poisson]{qpois}},
  \code{\link[stats:Poisson]{rpois}}
  respectively.
  


}
\references{

Yee, T. W. and Ma, C. (2022).
Generally--altered, --inflated, --truncated and --deflated
regression, with application to heaped and seeped data.
\emph{In preparation}.


% count data.
%with application to heaped and seeped counts.



}

\author{ T. W. Yee. }
\note{
  Functions \code{\link[VGAMdata]{Pospois}} and those similar
  have been moved to \pkg{VGAMdata}.
  It is better to use
  \code{dgaitdpois(x, lambda, truncate = 0)} instead of
  \code{dposbinom(x, lambda)}, etc.











% This is repeated info:
%  The MLM variant is called nonparametric while the
%  parent mixture variant is called parametric.
%  Values that are altered nonparametrically are effectively
%  removed from the main analysis because the MLM provides a
%  loose coupling.




}






\seealso{
  \code{\link{gaitdpoisson}},
  \code{\link{multinomial}},
  \code{\link{specials}},
  \code{\link{spikeplot}},
  \code{\link{dgaitdplot}},
  \code{\link{Zapois}},
  \code{\link{Zipois}},
  \code{\link[VGAMdata]{Pospois}}
  \code{\link[stats:Poisson]{Poisson}};
  \code{\link{Gaitdbinom}},
  \code{\link{Gaitdnbinom}},
  \code{\link{Gaitdlog}},
  \code{\link{Gaitdzeta}}.


% \code{\link{Gaitdnbinom}},




  
%  \code{\link{Gaitlog.mix}}  and \code{\link{Gaitlog.mlm}},
%  \code{\link{Gaitdpois.mix}} and \code{\link{Gaitdpois.mlm}},
%  \code{\link{Gaitnbinom.mlm}},

% \code{\link{gaitdpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},



}
\examples{
ivec <- c(6, 14); avec <- c(8, 11); lambda <- 10; xgrid <- 0:25
tvec <- 15; max.support <- 20; pobs.mix <- 0.05; pstr.i <- 0.25
dvec <- 13; pdip.mlm <- 0.05; pobs.mlm <- 0.05
(ddd <- dgaitdpois(xgrid, lambda, lambda.a = lambda + 5,
   truncate = tvec, max.support = max.support, pobs.mix = pobs.mix,
   pobs.mlm = pobs.mlm, a.mlm = avec,
   pdip.mlm = pdip.mlm, d.mlm = dvec,
   pstr.mix = pstr.i, i.mix = ivec))
\dontrun{ dgaitdplot(lambda, ylab = "Probability", xlab = "x",
   truncate = tvec, max.support = max.support, pobs.mix = pobs.mix,
   pobs.mlm = pobs.mlm, a.mlm = avec, all.lwd = 3,
   pdip.mlm = pdip.mlm, d.mlm = dvec,
   pstr.mix = pstr.i, i.mix = ivec, deflation = TRUE,
   main = "GAITD Combo PMF---Poisson Parent")   }
}
\keyword{distribution}


% 20200317; checked identical results to
% gaitdpois.mix() & gaitdpois.mlm().



