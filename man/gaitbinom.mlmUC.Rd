\name{Gaitbinom.mlm}
\alias{Gaitbinom.mlm}
\alias{dgaitbinom.mlm}
\alias{pgaitbinom.mlm}
\alias{qgaitbinom.mlm}
\alias{rgaitbinom.mlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered, -Inflated and -Truncated
  Binomial Distribution
  (GAIT--binom--MLM--MLM)

% (multinomial logit model based; GAIT--binom--MLM--MLM)
}
\description{
  Density,
  distribution function,
  quantile function
  and random generation
  for the generally-altered, -inflated and -truncated
  binomial distribution
  based on the multinomial logit model (MLM).
  This distribution is sometimes abbreviated as
  GAIT--Binom--MLM--MLM.


}
\usage{
dgaitbinom.mlm(x, size, prob, alter = NULL, inflate = NULL,
                 truncate = NULL, pobs.a = 0, pstr.i = 0,
                 byrow.arg = FALSE, log.arg = FALSE, .errorcheck = TRUE)
pgaitbinom.mlm(q, size, prob, alter = NULL, inflate = NULL,
                 truncate = NULL, pobs.a = 0, pstr.i = 0,
                 byrow.arg = FALSE, .errorcheck = TRUE)
qgaitbinom.mlm(p, size, prob, alter = NULL, inflate = NULL,
             truncate = NULL, pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
rgaitbinom.mlm(n, size, prob, alter = NULL, inflate = NULL,
             truncate = NULL, pobs.a = 0, pstr.i = 0, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, log.arg}{
  Same meaning as in \code{\link[stats]{rbinom}}.

  }

  \item{size, prob}{
  Same meaning as in \code{\link[stats]{rbinom}},
  i.e., for an ordinary binomial distribution.
  Short vectors are recycled.

  }
  \item{alter, inflate, truncate}{
  Vectors of nonnegative integers;
  the altered, inflated and truncated values.
  Each argument must have unique values only.
  The default value of \code{NULL} means an empty set.


% Value in the upper tail of the distribution may be
% truncated by using \code{max.support}.


  
  
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
% for \code{\link{rgaitbinom.mlm}},
% a multinomial sample is taken and if it belongs to the final
% column then binomial random variates are drawn.



  }
  \item{pobs.a}{
  This argument is similar to \code{pstr.i} but is used when
  \code{alter} is assigned a vector.
  The argument \code{byrow.arg} is used similarly to construct
  a matrix of dimension \code{n x length(alter)} of probabilities.
  This argument is not used unless \code{alter} is assigned.

  
  
  }
  \item{.errorcheck}{
  Logical. Should be ignored by the user.

  }
}
\details{
  See \code{\link{Gaitpois.mlm}} for generic details applying
  to both distributions.


  
  These functions do what
  \code{dgabinom},
  \code{dgibinom},
  \code{dgtbinom},
  \code{pgabinom},
  \code{pgibinom},
  \code{pgtbinom},
  \code{qgabinom},
  \code{qgibinom},
  \code{qgtbinom},
  \code{rgabinom},
  \code{rgibinom},
  \code{rgtbinom}
  collectively did because the
  arguments \code{alter}, \code{inflate} and \code{truncate}
  have been combined.




}
%\section{Warning }{
% See \code{\link{rgaitbinom.mlm}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgaitbinom.mlm(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgaitbinom.mlm(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgaitbinom.mlm} gives the density,
  \code{pgaitbinom.mlm} gives the distribution function,
  \code{qgaitbinom.mlm} gives the quantile function, and
  \code{rgaitbinom.mlm} generates random deviates.
  The default values of the arguments correspond to ordinary
  \code{\link[stats:Binomial]{dbinom}},
  \code{\link[stats:Binomial]{pbinom}},
  \code{\link[stats:Binomial]{qbinom}},
  \code{\link[stats:Binomial]{rbinom}}
  respectively.
  


}
%\references{
%None.
%}

\author{ T. W. Yee. }
%\note{
% 20120405; no longer true to a superior method:
% For \code{rposbinom}, the arguments of the function are fed
% into \code{\link[stats:Binomial]{rbinom}} until \eqn{n} positive
% values are obtained.  This may take a long time if \code{lambda}
% has values close to 0.


%  The family function \code{\link{posbinomson}} estimates
%  \eqn{\lambda}{lambda} by maximum likelihood estimation.


%}




% \code{\link{gaitbinomson}},
\seealso{
  \code{\link{Gaitnbinom.mlm}},
  \code{\link{Posbinom}},
  \code{\link[stats:Binomial]{Binomial}}.



% \code{\link{gabinomial}},
% \code{\link{Gaitbinom.mlm}},

%  \code{\link{zabinomial}},
%  \code{\link{zibinomial}},


}
\examples{
avec <- 5:6   # Alter these (special) values
ivec <- 9:10  # Inflate these (special) values
tvec <- 2:4   # Truncate these (special) values
pobs.a <- c(0.1, 0.2)
pstr.i <- (1:2)/10
size <- 10; prob <- 0.5; xx <- 0:11
y <- rgaitbinom.mlm(1000, size, prob, alter = avec, inflate = ivec,
       truncate = tvec, pstr.i = pstr.i, pobs.a = pobs.a, byrow = TRUE)
table(y)
(ii <- dgaitbinom.mlm(xx, size, prob, alter = avec, inflate = ivec,
       truncate = tvec, pstr.i = pstr.i, pobs.a = pobs.a, byrow = TRUE))
}
\keyword{distribution}

