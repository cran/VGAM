\name{beggs}
\alias{beggs}
\docType{data}
\title{Bacon and Eggs Data}
\description{
  Purchasing of bacon and eggs.

}
\usage{
data(beggs)
}
\format{
  Data frame of a two way table.

  \describe{

    \item{b0, b1, b2, b3, b4}{
    The \code{b} refers to bacon.
    The number of times bacon was purchased was 0, 1, 2, 3, or 4.


    }
    \item{e0, e1, e2, e3, e4}{
    The \code{e} refers to eggs.
    The number of times eggs was purchased was 0, 1, 2, 3, or 4.


    }

  }
}
\details{

The data is from Information Resources, Inc., a consumer panel
based in a large US city [see Bell and Lattin (1998) for further
details]. Starting in June 1991, the purchases in the bacon and
fresh eggs product categories for a sample of 548 households over
four consecutive store trips was tracked.  Only those grocery
shopping trips with a total basket value of at least five dollars
was considered.  For each household, the total number of bacon
purchases in their four eligible shopping trips and the total
number of egg purchases (usually a package of eggs) for the same
trips, were counted.


%    Data from Bell and Latin (1998).
%    Also see Danaher and Hardie (2005).


}
\source{

  Bell, D. R. and Lattin, J. M. (1998)
  Shopping Behavior and Consumer Preference
  for Store Price Format: Why `Large Basket' Shoppers Prefer EDLP.
  \emph{Marketing Science},
  \bold{17}, 66--88.


}
\references{

  Danaher, P. J. and Hardie, B. G. S. (2005)
  Bacon with Your Eggs?
  Applications of a New Bivariate Beta-Binomial Distribution.
  \emph{American Statistician},
  \bold{59}(4), 282--286.



}
\seealso{
  \code{\link[VGAM]{rrvglm}},
  \code{\link[VGAM]{rcim}},
  \code{\link[VGAM]{grc}}.


}
\examples{
beggs
colSums(beggs)
rowSums(beggs)
}
\keyword{datasets}


%
%
