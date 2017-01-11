\name{aux.posbernoulli.t}
\alias{aux.posbernoulli.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary Function for the
  Positive Bernoulli Family Function with Time Effects }
\description{
  Returns behavioural effects indicator variables from a
  capture history matrix.


}
\usage{
aux.posbernoulli.t(y, check.y = FALSE, rename = TRUE, name = "bei")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  Capture history matrix.
  Rows are animals, columns are sampling occasions, and
  values should be 0s and 1s only.

  }
  \item{check.y}{
  Logical, if \code{TRUE} then some basic checking is performed.

  }
  \item{rename, name}{
  If \code{rename = TRUE} then the behavioural effects indicator
  are named using the value of \code{name} as the prefix.
  If \code{FALSE} then use the same column names as \code{y}.

  }
}
\details{
  This function can help fit certain capture--recapture models
  (commonly known as \eqn{M_{tb}} or \eqn{M_{tbh}}
  (no prefix \eqn{h} means it is an intercept-only model)
  in the literature).
  See \code{\link{posbernoulli.t}} for details.


}
\value{
  A list with the following components.

  \describe{
    \item{cap.hist1}{
      A matrix the same dimension as \code{y}.
      In any particular row there are 0s up to
      the first capture. Then there are 1s thereafter.


    }
    \item{cap1}{
      A vector specifying which time occasion the animal
      was first captured.


    }
    \item{y0i}{
      Number of noncaptures before the first capture.

    }
    \item{yr0i}{
      Number of noncaptures after the first capture.


    }
    \item{yr1i}{
      Number of recaptures after the first capture.


    }
  }
}

% \author{ Thomas W. Yee. }

%\note{
%  Models \eqn{M_{tbh}}{M_tbh} can be fitted using the
%  \code{xij} argument (see \code{\link{vglm.control}})
%  to input the behavioural effect indicator variables.
%  Rather than manually setting these up, they may be more conveniently
%  obtained by \code{\link{aux.posbernoulli.t}}. See
%  the example below.
%
%
%}

%\section{Warning }{
%
%  See \code{\link{posbernoulli.tb}}.
%
%
%}

\seealso{
  \code{\link{posbernoulli.t}},
  \code{\link{deermice}}.

}

\examples{
# Fit a M_tbh model to the deermice data:
(pdata <- aux.posbernoulli.t(with(deermice, cbind(y1, y2, y3, y4, y5, y6))))

deermice <- data.frame(deermice,
                    bei = 0,  # Add this
                    pdata$cap.hist1)  # Incorporate these
head(deermice)  # Augmented with behavioural effect indicator variables
tail(deermice)
}
\keyword{models}
\keyword{regression}

