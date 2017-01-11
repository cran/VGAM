\name{gew}
\alias{gew}
\docType{data}
\title{ General Electric and Westinghouse Data }
\description{
  General Electric and Westinghouse capital data.

}
\usage{data(gew)}
\format{
  A data frame with 20 observations on the following 7 variables.
  All variables are numeric vectors.
  Variables ending in \code{.g} correspond to General Electric and
  those ending in \code{.w} are Westinghouse.
  \describe{
    \item{year}{The observations are the years from 1934 to 1953}
    \item{invest.g, invest.w}{investment figures.
These are \eqn{I=} Gross investment =
additions to plant and equipment plus maintenance and repairs
in millions of dollars deflated by \eqn{P_1}.


  }
    \item{capital.g, capital.w}{capital stocks.
These are \eqn{C=} The stock of plant and equipment =
accumulated sum of net additions to plant and equipment deflated
by \eqn{P_1} minus depreciation allowance deflated by \eqn{P_3}.


  }
    \item{value.g, value.w}{market values.
These are \eqn{F=} Value of the firm =
price of common and preferred shares at December 31
(or average price of December 31 and January 31 of the following year)
times number of common and preferred shares outstanding plus
total book value of debt at December 31 in millions of
dollars deflated by \eqn{P_2}.


}
  }
}
\details{
  These data are a subset of a table in Boot and de Wit (1960),
  also known as the Grunfeld data.
  It is used a lot in econometrics,
  e.g., for seemingly unrelated regressions
  (see \code{\link[VGAM:SURff]{SURff}}).


  Here,
  \eqn{P_1 =} Implicit price deflator of producers durable
  equipment (base 1947),
  \eqn{P_2 =} Implicit price deflator of G.N.P.
  (base 1947),
  \eqn{P_3 =} Depreciation expense deflator = ten years
  moving average of wholesale price index of metals and metal
  products (base 1947).





}
\source{

  Table 10 of:
  Boot, J. C. G. and de Wit, G. M. (1960)
  Investment Demand: An Empirical Contribution to the Aggregation Problem.
  \emph{International Economic Review},
  \bold{1}, 3--30.


  Grunfeld, Y. (1958)
  The Determinants of Corporate Investment.
  Unpublished PhD Thesis (Chicago).


}
\seealso{
  \code{\link[VGAM:SURff]{SURff}},
  \url{http://statmath.wu.ac.at/~zeileis/grunfeld}.


}

\references{

  Zellner, A. (1962)
  An efficient method of estimating seemingly unrelated regressions
  and tests for aggregation bias.
  \emph{Journal of the American Statistical Association},
  \bold{57}, 348--368.


}
\examples{
str(gew)
}
\keyword{datasets}
