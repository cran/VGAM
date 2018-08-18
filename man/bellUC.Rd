\name{Bell}
\alias{Bell}
\alias{dbell}
%\alias{pbell}
%\alias{qbell}
\alias{rbell}
\title{The Bell Distribution}
\description{
  Density,
  and random generation
  for the Topp-Leone distribution.

% distribution function,
% quantile function
}
\usage{
dbell(x, shape, log = FALSE)
rbell(n, shape)
}
%pbell(q, shape, lower.tail = TRUE, log.p = FALSE)
%qbell(p, shape)
\arguments{
  \item{x, n}{  % q, p, 
  Same as \code{\link[stats:Uniform]{Uniform}}.


  }
  \item{shape}{the (shape) parameter, which is positive.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
% \item{lower.tail, log.p}{
% Same meaning as in \code{\link[stats:Normal]{pnorm}}
% or \code{\link[stats:Normal]{qnorm}}.


% }

}
\value{
  \code{dbell} gives the density,
  \code{rbell} generates random deviates.
  If \code{shape} is large then \code{rbell} will become
  computationally expensive.


% \code{pbell} gives the distribution function,
% \code{qbell} gives the quantile function, and



}
%\references{
%}


% \author{ T. W. Yee }
\details{
  See \code{\link{bell}}, the \pkg{VGAM} family function for
  estimating the parameter \eqn{s} by maximum likelihood
  estimation.





% , for the formula of the probability density function.


}
%\note{
%  The Bell distribution is related to the Poisson distribution.


%}
\seealso{
  \code{\link{bell}}.
% \code{\link{Triangle}}.



}
\examples{
\dontrun{
plot(0:15, dbell(0:15, shape = 1.5), type = "h", las = 1, col = "blue")
}
}
\keyword{distribution}


