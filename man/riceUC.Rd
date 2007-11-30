\name{Rice}
\alias{Rice}
\alias{drice}
%\alias{price}
%\alias{qrice}
\alias{rrice}
\title{The Rice Distribution}
\description{
  Density
% distribution function, quantile function
  and random generation for the
  Rice distribution.

}
\usage{
drice(x, vee, sigma)
%price(q, vee, sigma)
%qrice(p, vee, sigma)
rrice(n, vee, sigma)
}
\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
 \item{n}{number of observations.
   Must be a positive integer of length 1.}
  \item{vee, sigma}{ See \code{\link{riceff}}.
    }
}
\value{
  \code{drice} gives the density,
% \code{price} gives the distribution function,
% \code{qrice} gives the quantile function, and
  \code{rrice} generates random deviates.

}
\author{ T. W. Yee }
\details{
  See \code{\link{riceff}}, the \pkg{VGAM} family function
  for estimating the two parameters,
  for the formula of the probability density function and other details.

}
%\section{Warning }{
%
%}
\seealso{
  \code{\link{riceff}}.

}
\examples{
\dontrun{
x = seq(0.01, 7, len=201)
plot(x, drice(x, vee=0, sigma=1), type="n", las=1,, ylab="",
     main="Density of Rice distribution for various values of v")
sigma = 1; vee = c(0,0.5,1,2,4)
for(ii in 1:length(vee)) lines(x, drice(x, vee[ii], sigma), col=ii)
legend(x=5, y=0.6, legend=as.character(vee),
       col=1:length(vee), lty=1)
}
}
\keyword{distribution}


