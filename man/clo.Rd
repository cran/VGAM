\name{clo}
\alias{clo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Redirects the User to rrvglm() }
\description{
  Redirects the user to the function \code{\link{rrvglm}}.
}
\usage{
clo(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ Ignored. }
}
\details{
CLO stands for \emph{constrained linear ordination}, and
is fitted with a statistical class of models called
\emph{reduced-rank vector generalized linear models}
(RR-VGLMs). It allows for generalized reduced-rank regression
in that response types such as Poisson counts and presence/absence
data can be handled.

Currently in the \pkg{VGAM} package, \code{\link{rrvglm}} is
used to fit RR-VGLMs. However, the Author's opinion is that
linear responses to a latent variable (composite environmental
gradient) is not as common as unimodal responses, therefore
\code{\link{cqo}} is often more appropriate.

The new CLO/CQO/CAO nomenclature described in Yee (2006).
}
\value{
  Nothing is returned; an error message is issued.



}
\references{

Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

}
\author{Thomas W. Yee}

\seealso{
  \code{\link{rrvglm}},
  \code{\link{cqo}}.
}
\examples{

\dontrun{
clo()
}

}
\keyword{models}
\keyword{regression}

