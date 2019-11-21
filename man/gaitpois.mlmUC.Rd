\name{Gaitpois.mlm}
\alias{Gaitpois.mlm}
\alias{dgaitpois.mlm}
\alias{pgaitpois.mlm}
\alias{qgaitpois.mlm}
\alias{rgaitpois.mlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered, -Inflated and -Truncated
  Poisson Distribution
  (GAIT--Pois--MLM--MLM)


% (multinomial logit model based; GAIT--Pois--MLM--MLM)


}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally-altered, -inflated and -truncated
  Poisson distribution,
  based on the multinomial logit model (MLM).
  This distribution is sometimes abbreviated as
  GAIT--Pois--MLM--MLM.

}
\usage{
dgaitpois.mlm(x, lambda, alter = NULL, inflate = NULL, truncate = NULL,
              max.support = Inf, pobs.a = 0, pstr.i = 0, byrow.arg = FALSE,
              log.arg = FALSE)
pgaitpois.mlm(q, lambda, alter = NULL, inflate = NULL, truncate = NULL,
              max.support = Inf, pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
qgaitpois.mlm(p, lambda, alter = NULL, inflate = NULL, truncate = NULL,
              max.support = Inf, pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
rgaitpois.mlm(n, lambda, alter = NULL, inflate = NULL, truncate = NULL,
              max.support = Inf, pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log.arg}{
  Same meaning as in \code{\link[stats]{rpois}}.

  }

  \item{lambda}{
  Same meaning as in \code{\link[stats]{rpois}},
  i.e., for an ordinary Poisson distribution.
  Short vectors are recycled.

  }
  \item{alter, inflate, truncate}{
  Vectors of nonnegative integers;
  the altered, inflated and truncated values.
  Each argument must have unique values only.
  The default value of \code{NULL} means an empty set.
  Values in the upper tail of the distribution may be
  truncated by using \code{max.support}.


  
  
% Must be sorted and have unique values only.
  }
  \item{pstr.i, byrow.arg}{
    The first argument can be a \code{length(inflate)}-vector of
    probabilities;
  else a matrix of dimension \code{n x length(inflate)} of probabilities.
  If a vector then this matrix is constructed from the vector
  using \code{byrow.arg} to determine the enumeration of the elements
  (similar to \code{\link[base]{matrix}}).
  These arguments are not used unless \code{inflate} is assigned.



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
  \item{pobs.a}{
  This argument is similar to \code{pstr.i} but is used when
  \code{alter} is assigned a vector.
  The argument \code{byrow.arg} is used similarly to construct
  a matrix of dimension \code{n x length(alter)} of probabilities.
  This argument is not used unless \code{alter} is assigned.

  
  
  }
  \item{max.support}{
  numeric; the maximum support value so that any value larger
  has been truncated.
  This argument is necessary because
  \code{truncate = (max.support + 1):Inf} is impractical,
  so this is the way to truncate the upper tail of the distribution.
  Note that \code{max(truncate) < max.support} must be satisfied
  otherwise an error message will be issued.


  
  }
}
\details{
  These functions allow any combination of 3 operator types:
  truncation, alteration and inflation.
The precedence of these are
truncation, then
alteration and lastly
inflation.
This order minimizes the potential interference among the operators.
Loosely, a set of probabilities is set to 0 by truncation
and the remaining probabilities are scaled up.
Then a different set of probabilities are set to some
values \code{pobs.a}
and the remaining probabilities are rescaled up.
Then another different set of probabilities is inflated by
an amount \code{pstr.i} so that individual elements in this set
have two sources.
Then the remainder of all the probabilities are left
unchanged apart from further rescaling.



Any value of the support of the distribution that is
altered, inflated or truncated is called a \emph{special} value.
A special value that is altered may mean that its probability
increases or decreases relative to the parent distribution.
An inflated special value means that its probability has
increased, provided alteration has not made it decrease
in the first case.





% It can be considered a mixture of a multinomial distribution
% and an ordinary Poisson distribution.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-Poisson and Poisson
%  distributions become more similar.

% Unlike similar functions for the Poisson distribution, a zero value
% of \code{lambda} is not permitted here.


  These functions do what
  \code{dgapois},
  \code{dgipois},
  \code{dgtpois},
  \code{pgapois},
  \code{pgipois},
  \code{pgtpois},
  \code{qgapois},
  \code{qgipois},
  \code{qgtpois},
  \code{rgapois},
  \code{rgipois},
  \code{rgtpois}
  collectively did because the
  arguments \code{alter}, \code{inflate} and \code{truncate}
  have been included in one function.
  The old functions may be relegated to \pkg{VGAMdata} in the
  future.



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
  \code{dgaitpois.mlm} gives the density,
  \code{pgaitpois.mlm} gives the distribution function,
  \code{qgaitpois.mlm} gives the quantile function, and
  \code{rgaitpois.mlm} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:Poisson]{dpois}},
  \code{\link[stats:Poisson]{ppois}},
  \code{\link[stats:Poisson]{qpois}},
  \code{\link[stats:Poisson]{rpois}}
  respectively.
  


}
\references{

Yee, T. W. and Ma, C. C. (2019)
Generally-altered, -inflated and -truncated count regression,
with application to heaped and seeped data.
\emph{In preparation}.


}

\author{ T. W. Yee. }
\note{
  These functions comprise the MLM variant and therefore
  are more nonparametric and unstructured in nature.
  There is another variant based on parent distribution
  mixtures---see \code{\link{Gaitpois.mix}} for information.


}




% \code{\link{gaitpoisson.mlm}},
% \code{\link{gaitpoisson.mix}},
\seealso{
  \code{\link{gatpoisson.mlm}},
  \code{\link{Gaitpois.mix}},
  \code{\link{Gaitnbinom.mlm}},
  \code{\link{Pospois}},
  \code{\link[stats:Poisson]{Poisson}}.



% \code{\link{gaitpoisson.mlm}},
% \code{\link{Gtpois}},
% \code{\link{Gapois.mix}},

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},


}
\examples{
avec <- 5:6    # Alter these (special) values
ivec <- 9:10   # Inflate these (special) values
tvec <- 12:14  # Truncate these (special) values
pobs.a <- c(0.1, 0.2)
pstr.i <- (1:2) / 10
lambda <- 15; xx <- 0:22
y <- rgaitpois.mlm(1000, lambda, alter = avec, inflate = ivec,
                     truncate = tvec, max.support = 20,
                     pstr.i = pstr.i, pobs.a = pobs.a, byrow = TRUE)
table(y)
(ii <- dgaitpois.mlm(xx, lambda, alter = avec, inflate = ivec,
                       truncate = tvec, max.support = 20,
                       pstr.i = pstr.i, pobs.a = pobs.a, byrow = TRUE))
}
\keyword{distribution}

% 20190927; checked identical results to gapoisson() and gtpoisson().



