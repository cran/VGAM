\name{sloglink}
\alias{sloglink}
\alias{lcsloglink}
%\alias{ewsloglink}
%- Also NEED an '\alias' for EACH other
% topic documented here.
\title{ Square root--Log Link Mixtures}
\description{
  Computes some square root--log mixture link
  transformations, including their inverse and
  the first few derivatives.

}
\usage{
  sloglink(theta, bvalue = NULL, taumix.log = 1, tol = 1e-13,
    nmax = 99, inverse = FALSE, deriv = 0, short = TRUE,
    tag = FALSE, c10 = c(2, -2))
lcsloglink(theta, bvalue = NULL, pmix.log = 0.01, tol = 1e-13,
    nmax = 99, inverse = FALSE, deriv = 0, short = TRUE,
    tag = FALSE, c10 = c(2, -2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.
  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }
  \item{taumix.log}{
  Numeric, of length 1.
  Mixing parameter directed at
  \code{\link{loglink}}.  Then
  \code{1 - exp(-taumix.log * theta)} is
  used to weight
  \code{\link{sqrtlink}}.
  Thus a 0 value will result in
  \code{\link{loglink}},
  and a very large numeric such as \code{1e4}
  should be roughly equivalent to
  \code{\link{sqrtlink}} over almost
  all of the parameter space.



  }
  \item{pmix.log}{
  Numeric, of length 1.
  Mixing probability assigned
  to \code{\link{loglink}}.  Then
  \code{1 - pmix.log} is used to weight
  \code{\link{sqrtlink}}.
  Thus a 0 value will result in
  \code{\link{sqrtlink}}
  and 1 is equivalent to \code{\link{loglink}}.



  }
  \item{tol, nmax}{
    Arguments fed into a function implementing a
    vectorized bisection method.


  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
  \item{c10}{
    See \code{\link{sqrtlink}}
    and
    \code{\link{loglink}}.
    
  }

}
\details{
% Poisson regression with the
% log link
% suffers from a fundamental flaw:
% it suffers from the Hauck--Donner effect
% (HDE; Hauck and Donner, 1977)
% whereby, in Wald tests,
% the power function is not a monotonic
% increasing
% function moving away from the null value
% (typically 0).



% These link functions remedy this problem.
% Two variants are implemented:
% an exponential-weighted (EW) sum and a
% simple \emph{linear combination} (LC).
% The latter
% is a mixture of a \emph{s}quare root link
% (with a very heavy weighting)
% and a log link.
% However, the problem with
% \code{\link{sqrtpoislink}}
% is that it is semi-bounded on
% \code{(-2, Inf)}
% so it can be unsuitable for regression.
% By taking a LC of this link with an ordinary
% log link,
% the boundedness problem can be ameliorated.
% It is suggested that a smaller positive values of
% \code{pmix.log} are preferred, since this
% means that the resulting Wald statistics are
% far less likely to suffer from the HDE.


% The LC function operates similarly to
% \code{\link{alogitlink}} in binary regression.
% It is an example of a
% constant information augmented (CIA) link.
% It combats the Hauck--Donner effect (HDE) so
% that the resulting Wald statistic p-values
% are HDE-free.
  For general information see
  \code{\link{alogitlink}}.



}
\value{
  The following holds for the
  \emph{linear combination} (LC)
  variant.
  For \code{deriv = 0},
  \code{(1 - pmix.log) * sqrtlink(mu, c10 = c10)
    + pmix.log * loglink(mu)}
  when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then a nonlinear
  equation is solved for \code{mu},
  given
  \code{eta} passed in as \code{theta}.
  For \code{deriv = 1}, then the function
  returns \emph{d} \code{eta} / \emph{d}
  \code{theta} as a function of \code{theta} if
  \code{inverse = FALSE}, else if \code{inverse
  = TRUE} then it returns the reciprocal.


}
%\references{


%Hauck, J. W. W. and A. Donner (1977).
%Wald's test as applied to hypotheses in
%log analysis.
%\emph{Journal of the American Statistical
%Association},
%\bold{72}, 851--853.


% Corrigenda: JASA, \bold{75}, 482.

% \textit{JASA 72(360): 851--3}] 75 (370), 482

% Yee, T. W. (2023).
% \emph{Constant information augmented link
% functions impervious
% to the Hauck--Donner effect in
% vector generalized linear models}.
% Under review.



%}
\author{ Thomas W. Yee }

\section{Warning }{
  The default values for \code{taumix.log}
  and \code{pmix.log}
  may change in the future.
  The name and order of the arguments
  may change too.

}

%\note{
%  Numerical instability may occur when
%  \code{theta} is close to 1 or 0.  One way of
%  overcoming this is to use \code{bvalue}.

%}

\seealso{
  \code{\link{alogitlink}},
  \code{\link{sqrtlink}},
  \code{\link{loglink}},
  \code{\link{Links}},
  \code{\link{poissonff}},
  \code{\link{hdeff}}.



% \url{https://www.cia.gov}.



}
\examples{
mu  <- seq(0.01, 3, length = 10)
sloglink(mu)
max(abs(sloglink(sloglink(mu), inv = TRUE) - mu))  # 0?
}
\keyword{math}
\keyword{models}
\keyword{regression}

%plot(y, loglink(y, inv = T), type = "l", col = "limegreen",
%     xlab = "transformation", ylab = "p", lwd=2, las=1,
%     main = "Some inverse probability link functions")
%lines(y, probitlink(y, inv = T), col = "purple", lwd=2)
%lines(y, clogloglink(y, inv = T), col = "chocolate", lwd=2)
%abline(h=0.5, v = 0, lty = "dashed")




