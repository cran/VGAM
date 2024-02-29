\name{alogitlink}
\alias{alogitlink}
\alias{lcalogitlink}
%\alias{ewalogitlink}
%- Also NEED an '\alias' for EACH other
% topic documented here.
\title{ Arcsine--Logit Link Mixtures}
\description{
  Computes some arcsine--logit mixture link
  transformations, including their inverse and
  the first few derivatives.

}
\usage{
  alogitlink(theta, bvalue = NULL, taumix.logit = 1,
     tol = 1e-13, nmax = 99, inverse = FALSE, deriv = 0,
     short = TRUE, tag = FALSE, c10 = c(4, -pi))
lcalogitlink(theta, bvalue = NULL, pmix.logit = 0.01,
     tol = 1e-13, nmax = 99, inverse = FALSE, deriv = 0,
     short = TRUE, tag = FALSE, c10 = c(4, -pi))
}
%- maybe also 'usage' for other objs doced here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.
  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }
  \item{taumix.logit}{
  Numeric, of length 1.
  Mixing parameter assigned
  to \code{\link{logitlink}}.  Then
  \code{1 - exp(-taumix.log * theta)} is used
  to weight
  \code{\link{asinlink}}.
  Thus a 0 value will result in
  \code{\link{logitlink}}
  and a very large numeric such as \code{1e4}
  should be roughly equivalent to
  \code{\link{asinlink}} over almost all
  of the parameter space.



  }
  \item{pmix.logit}{
  Numeric, of length 1.
  Mixing probability assigned
  to \code{\link{logitlink}}.  Then
  \code{1 - pmix.logit} is used to weight
  \code{\link{asinlink}}.
  Thus a 0 value will result in
  \code{\link{asinlink}}.
  and 1 is equivalent to \code{\link{logitlink}}.



  }
  \item{tol, nmax}{
    Arguments fed into a function implementing a
    vectorized bisection method.


  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
  \item{c10}{
    See \code{\link{asinlink}}
    and
    \code{\link{logitlink}}.
    
  }

}
\details{
%  Binary regression with the
%  four popular choices of links
%  (\code{\link{logitlink}},
%  \code{\link{probitlink}},
%  \code{\link{clogloglink}},
%  \code{\link{cauchitlink}})
%suffer from a fundamental flaw:
%they suffer from the Hauck--Donner effect
%(HDE; Hauck and Donner, 1977)
%whereby, in Wald tests,
%the power function is not a monotonic
%increasing
%function moving away from the null value
%(typically 0).



%These link functions remedy this problem.
%Two variants are implemented:
%an exponential-weighted (EW) sum and a
%simple \emph{linear combination} (LC).
%The latter
%is a mixture of an \emph{a}rcsine link
%(with a very heavy weighting)
%and a logit link.
%However, the problem with \code{\link{asinlink}}
%is that it is bounded between
%\code{(-pi, pi)}
%so it can be unsuitable for regression.
%By mixing this link with an ordinary
%binary regression link such as the
%logit,
%the boundedness problem can be ameliorated.



%It is suggested that a smaller positive values
%of \code{pmix.logit} are preferred, since this
%means that the resulting Wald statistics are
%far less likely to suffer from the HDE.



%The name of this function was chosen from:
%\emph{ar}csine \emph{m}ixtur\emph{e} with
%the \emph{logit}.  It is referred to as the
%\emph{armed logit}, where the \emph{d} is added
%to give it English sense.



%Numerical values of \code{theta} close to 0
%or 1 or out of range result in \code{Inf},
%\code{-Inf}, \code{NA} or \code{NaN}.



%  More generally, this function is a
%  constant information augmented (CIA) link
%  that allows a form of binary regression
%  to be performed that combats the
%  Hauck--Donner effect (HDE) so that
%  the resulting Wald statistic p-values are
%  HDE-free.
%  Because \code{\link{asinlink}} has a
%  bounded range,
%  it can be unsuitable for regression,
%  hence it is augmented with
%  \code{\link{logitlink}} in the form of a
%  mixture, and this is implemented
%  in \code{\link{alogitlink}}.
%  The tolerance for the HDE can be controlled by
%  \code{pmix.logit}.
%  In general, a suitably crafted CIA link enables
%  parameters of a VGLM to be modelled in a
%  HDE-free manner.




  \code{\link{lcalogitlink}} is a
  \emph{linear combination} (LC) of
  \code{\link{asinlink}} and
  \code{\link{logitlink}}.


}
\value{
  The following holds for the LC variant.
  For \code{deriv >= 0},
\code{(1 - pmix.logit) * asinlink(p, deriv = deriv)
    + pmix.logit * logitlink(p, deriv = deriv)}
  when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then a nonlinear
  equation is solved for the probability,
  given
  \code{eta}.
  For \code{deriv = 1}, then the function
  returns \emph{d} \code{eta} / \emph{d}
  \code{theta} as a function of \code{theta} if
  \code{inverse = FALSE}, else if \code{inverse
  = TRUE} then it returns the reciprocal.


}
\references{


Hauck, J. W. W. and A. Donner (1977).
Wald's test as applied to hypotheses in
logit analysis.
\emph{Journal of the American Statistical
  Association},
\bold{72}, 851--853.


% Corrigenda: JASA, \bold{75}, 482.



% \textit{JASA 72(360): 851--3}] 75 (370), 482


  
% Yee, T. W. (2023).
% \emph{Constant information augmented link
% functions impervious
% to the Hauck--Donner effect in
% vector generalized linear models}.
% Under review.




}
\author{ Thomas W. Yee }

\section{Warning }{
  The default values for \code{taumix.logit}
  and \code{pmix.logit}
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
  \code{\link{asinlink}},
  \code{\link{logitlink}},
  \code{\link{Links}},
  \code{\link{probitlink}},
  \code{\link{clogloglink}},
  \code{\link{cauchitlink}},
  \code{\link{binomialff}},
  \code{\link{sloglink}},
  \code{\link{hdeff}},
  \url{https://www.cia.gov/index.html}.



 }
\examples{
p <- seq(0.01, 0.99, length= 10)
alogitlink(p)
max(abs(alogitlink(alogitlink(p), inv = TRUE) - p))  # 0?

\dontrun{
par(mfrow = c(2, 2), lwd = (mylwd <- 2))
y <- seq(-4, 4, length = 100)
p <- seq(0.01, 0.99, by = 0.01)

for (d in 0:1) {
  matplot(p, cbind(logitlink(p, deriv = d), probitlink(p, deriv = d)),
          type = "n", col = "blue", ylab = "transformation",
          las = 1, main = if (d == 0) "Some probability link functions"
          else "First derivative")
  lines(p,   logitlink(p, deriv = d), col = "green")
  lines(p,  probitlink(p, deriv = d), col = "blue")
  lines(p, clogloglink(p, deriv = d), col = "tan")
  lines(p,  alogitlink(p, deriv = d), col = "red3")
  if (d == 0) {
    abline(v = 0.5, h = 0, lty = "dashed")
    legend(0, 4.5, c("logitlink", "probitlink", "clogloglink",
           "alogitlink"), lwd = mylwd,
           col = c("green", "blue", "tan", "red3"))
  } else
    abline(v = 0.5, lwd = 0.5, col = "gray")
}

for (d in 0) {
  matplot(y, cbind( logitlink(y, deriv = d, inverse = TRUE),
                   probitlink(y, deriv = d, inverse = TRUE)),
          type  = "n", col = "blue", xlab = "transformation", ylab = "p",
          main = if (d == 0) "Some inverse probability link functions"
          else "First derivative", las=1)
  lines(y,   logitlink(y, deriv = d, inverse = TRUE), col = "green")
  lines(y,  probitlink(y, deriv = d, inverse = TRUE), col = "blue")
  lines(y, clogloglink(y, deriv = d, inverse = TRUE), col = "tan")
  lines(y,  alogitlink(y, deriv = d, inverse = TRUE), col = "red3")
  if (d == 0) {
      abline(h = 0.5, v = 0, lwd = 0.5, col = "gray")
      legend(-4, 1, c("logitlink", "probitlink", "clogloglink",
             "alogitlink"), lwd = mylwd,
             col = c("green", "blue", "tan", "red3"))
  }
}
par(lwd = 1)
}
}
\keyword{math}
\keyword{models}
\keyword{regression}

%plot(y, logitlink(y, inverse = TRUE), type = "l", col = "green",
%     xlab = "transformation", ylab = "p",
%     lwd=2, las=1, main = "Some inverse probability link functions")
%lines(y, probitlink(y, inverse = TRUE), col = "blue", lwd=2)
%lines(y, clogloglink(y, inverse = TRUE), col = "tan", lwd=2)
%abline(h=0.5, v = 0, lty = "dashed")




