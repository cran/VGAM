\name{cgo}
\alias{cgo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Redirects the user to cqo }
\description{
  Redirects the user to the function \code{\link{cqo}}.
}
\usage{
cgo(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ Ignored. }
}
\details{
The former function \code{cgo} has been renamed \code{\link{cqo}}
because CGO (for \emph{canonical Gaussian ordination}) is a confusing
and inaccurate name.
CQO (for \emph{constrained quadratic ordination}) is better.
This new nomenclature described in Yee (2006).


}
\value{
  Nothing is returned; an error message is issued.


}
\references{
Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.

Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.

}
\author{Thomas W. Yee}

\section{Warning }{

The code, therefore, in Yee (2004) will not run without changing the
\code{"g"} to a \code{"q"}.

}

\seealso{
  \code{\link{cqo}}.

}
\examples{

\dontrun{
cgo()
}

}
\keyword{models}
\keyword{regression}

