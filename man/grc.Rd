\name{grc}
\alias{grc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fitting Goodman's RC Association Model }
\description{
  Fits a Goodman's RC Association Model to a matrix of counts 

}
\usage{
grc(y, Rank = 1, Index.corner = 2:(1 + Rank),
    Structural.zero = 1, summary.arg = FALSE, h.step = 1e-04, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{
  A matrix of counts. Output from \code{table()} is acceptable;
  it is converted into a matrix.
  Note that \code{y} must be at least 3 by 3. 
  }
  \item{Rank}{
  An integer in the range 
  {1,\ldots,\code{min(nrow(y), ncol(y))}}. 
  This is the dimension of the fit.

  }
  \item{Index.corner}{
  A vector of \code{Rank} integers.
  These are used to store the \code{Rank} by \code{Rank}
  identity matrix in the
  \code{A} matrix; corner constraints are used.

  }
  \item{Structural.zero}{ 
  An integer in the range {1,\ldots,\code{min(nrow(y), ncol(y))}},
  specifying the row that is used as the structural zero.
  }
  \item{summary.arg}{
  Logical. If \code{TRUE}, a summary is returned.
  If \code{TRUE}, \code{y} may be the output (fitted
  object) of \code{grc()}.

  }
  \item{h.step}{
  A small positive value that is passed into
  \code{summary.rrvglm()}. Only used when \code{summary.arg=TRUE}. }
  \item{\dots}{ Arguments that are passed into \code{rrvglm.control()}.

  }
}
\details{
  Goodman's RC association model can fit a reduced-rank approximation
  to a table of counts.  The log of each cell mean is decomposed as an
  intercept plus a row effect plus a column effect plus a reduced-rank
  part. The latter can be collectively written \code{A \%*\% t(C)},
  the product of two `thin' matrices.
  Indeed, \code{A} and \code{C} have \code{Rank} columns.
By default, the first column and row of the interaction matrix
\code{A \%*\% t(C)} is chosen 
to be structural zeros, because \code{Structural.zero=1}.
This means the first row of \code{A} are all zeros.   

This function uses \code{options()$contrasts} to set up the row and 
column indicator variables.

}
\value{
  An object of class \code{"grc"}, which currently is the same as
  an \code{"rrvglm"} object.
}
\references{
Goodman, L. A. (1981)
Association models and canonical correlation in the analysis
of cross-classifications having ordered categories.
\emph{Journal of the American Statistical Association},
\bold{76}, 320--334.

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

Documentation accompanying the \pkg{VGAM} package at 
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information about the setting up of the
indicator variables.
 
}
\author{ Thomas W. Yee }
\note{
  This function sets up variables etc. before calling \code{rrvglm()}.
  The \code{...} is passed into \code{rrvglm.control()}, meaning, e.g.,
  \code{Rank=1} is default. Seting \code{trace=TRUE} may be useful for
  monitoring convergence.

  Using \code{criterion="coefficients"} can result in slow convergence.

  If \code{summary=TRUE}, then \code{y} can be a \code{"grc"} object,
  in which case a summary can be returned. That is,
  \code{grc(y, summary=TRUE)} is equivalent to
  \code{summary(grc(y))}.

}

\section{Warning}{
  This function temporarily creates a permanent data frame called
  \code{.grc.df}, which used to be needed by \code{summary.rrvglm()}.
  Then \code{.grc.df} is deleted before exiting the function.  If an
  error occurs, then \code{.grc.df} may be present in the workspace.

}

\seealso{
\code{\link{rrvglm}},
\code{\link{rrvglm.control}},
\code{\link{rrvglm-class}},
\code{summary.grc},
\code{\link{auuc}}.
}

\examples{
# Some undergraduate student enrolments at the University of Auckland in 1990
data(auuc) 
g1 = grc(auuc, Rank=1)
fitted(g1)
summary(g1)

g2 = grc(auuc, Rank=2, Index.corner=c(2,5))
fitted(g2)
summary(g2)
}
\keyword{models}
\keyword{regression}
