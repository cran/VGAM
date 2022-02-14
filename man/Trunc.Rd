\name{Trunc}
\alias{Trunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Truncated Values for the GT-Expansion Method

}
\description{
  Given the minimum and maximum values in a response variable,
  and a positive multiplier,
  returns the truncated values 
  for generally-truncated regression

  


}
\usage{
Trunc(Range, mux = 2, location = 0, omits = TRUE)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{Range}{
    Numeric, of length 2 containing the minimum and maximum
    (in that order) of the untransformed data.
    Alternatively, if \code{length(Range) > 2} then it is assumed
    that the entire untransformed data is passed in so that
    \code{\link[base]{range}} is applied.


  }
  \item{mux}{
    Numeric,
    the multiplier.
    A positive integer.

  }
  \item{location}{
    Numeric,
    the location parameter, allows a shift to the right.

  }
  \item{omits}{
    Logical.
    The default is to return the truncated values (those being
    omitted).
    If \code{FALSE} then the multiples are returned.



  }
}
\details{
  Generally-truncated regression can handle underdispersion
  with respect to some parent or base distribution such as the
  Poisson.
  Yee and Ma (2022) call this the
  \emph{GT-Expansion} (GTE) method,
  which is a special case of the GT-location-scale (GT-LS) method.
  This is a utility function to help make life easier.
  It is assumed that the response is a count variable.
  
  

}
\value{
  A vector of values to be fed into the \code{truncate} argument
  of a \pkg{VGAM} family function such as \code{\link{gaitdpoisson}}.
  If \code{mux = 1} then the function will return a \code{NULL}
  rather than \code{integer(0)}.



%---in general setting mux = 1 should be avoided.
  


}
% \references{
% Yee, T. W. and Ma, C. (2022).
%Generally--altered, --inflated, --truncated and --altered
%count regression.
% \emph{In preparation}.


%% , \bold{3}, 15--41.


% }




%\section{Warning }{
%    See \code{\link{gaitdpoisson}}



%}



\author{ T. W. Yee}
%\note{
%  See \code{\link{gaitdpoisson}}


% See \code{\link{gaitlog.mlm}} for other general details.


%}

\seealso{
  \code{\link{gaitdpoisson}},
  \code{\link{gaitdlog}},
  \code{\link{gaitdzeta}},
  \code{\link[base]{range}},
  \code{\link[base]{setdiff}},
  \code{\link{goffset}}.


% \code{\link{gaitlog.mlm}},





}

\examples{
Trunc(c(1, 8), 2)

set.seed(1)  # The following example is based on the normal
mymean <- 20; m.truth <- 3  # approximation to the Poisson.
gdata <- data.frame(y1 = round(rnorm((nn <- 1000), mymean,
                                     sd = sqrt(mymean / m.truth))))
org1 <- with(gdata, range(y1))  # Original range of the raw data
m.max <- 5  # Try multipliers 1:m.max
logliks <- numeric(m.max)
names(logliks) <- as.character(1:m.max)
for (i in 1:m.max) {
  logliks[i] <- logLik(vglm(i * y1 ~ offset(rep(log(i), nn)),
    gaitdpoisson(truncate = Trunc(org1, i)), data = gdata))
}
sort(logliks, decreasing = TRUE)  # Best to worst
\dontrun{ par(mfrow = c(1, 2))
plot(with(gdata, table(y1)))  # Underdispersed wrt Poisson
plot(logliks, col = "blue", type = "b", xlab = "Multiplier") }
}
\keyword{models}
\keyword{regression}



%                 truncate = tvec, max.support = max.support
%                           eq.ip = TRUE, max.support = max.support


