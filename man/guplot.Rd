\name{guplot}
\alias{guplot}
\alias{guplot.default}
\alias{guplot.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gumbel Plot }
\description{
  Produces a Gumbel plot,
  a diagnostic plot for checking whether the data appears to be
  from a Gumbel distribution.

}
\usage{
guplot(object, ...)
guplot.default(y, main = "Gumbel Plot",
    xlab = "Reduced data", ylab = "Observed data", type = "p", ...)
guplot.vlm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ A numerical vector. \code{NA}s etc. are not allowed.}
  \item{main}{Character. Overall title for the plot. }
  \item{xlab}{Character. Title for the x axis. }
  \item{ylab}{Character. Title for the y axis. }
  \item{type}{Type of plot. The default means points are
    plotted. }
  \item{object}{ An object that inherits class \code{"vlm"},
    usually of class \code{\link{vglm-class}} or
    \code{\link{vgam-class}}.

  }
  \item{\dots}{ Graphical argument passed into
    \code{\link[graphics]{plot}}. See \code{\link[graphics]{par}}
    for an exhaustive list. The arguments \code{xlim} and
    \code{ylim} are particularly useful.

  }
}
\details{
  If \eqn{Y} has a Gumbel distribution then plotting the sorted
  values \eqn{y_i} versus the \emph{reduced values} \eqn{r_i} should
  appear linear. The reduced values are given by
  \deqn{r_i = -\log(-\log(p_i)) }{%
        r_i = - log(- log(p_i)) }
  where \eqn{p_i} is the \eqn{i}th plotting position, taken
  here to be \eqn{(i-0.5)/n}.
  Here, \eqn{n} is the number of observations.
  Curvature upwards/downwards may indicate a Frechet/Weibull
  distribution, respectively. Outliers may also be detected
  using this plot.


  The function \code{guplot} is generic, and
  \code{guplot.default} and \code{guplot.vlm} are some
  methods functions for Gumbel plots.


}
\value{
  A list is returned invisibly with the following components.
  \item{x }{The reduced data. }
  \item{y }{The sorted y data. }
}

%% zz not sure about the reference
\references{

  Coles, S. (2001)
  \emph{An Introduction to Statistical Modeling of Extreme Values}.
  London: Springer-Verlag.


  Gumbel, E. J. (1958)
  \emph{Statistics of Extremes}.
  New York, USA: Columbia University Press.


}
\author{ T. W. Yee }
\note{
  The Gumbel distribution is a special case of the
  GEV distribution with shape parameter equal to zero.


}
\seealso{
  \code{\link{gumbel}},
  \code{\link{gumbelff}},
  \code{\link{gev}},
  \code{\link{venice}}.


}
\examples{\dontrun{guplot(rnorm(500), las = 1) -> ii
names(ii)

guplot(with(venice, r1), col = "blue")  # Venice sea levels data
}}
\keyword{models}
\keyword{regression}

