\name{olympics}
\alias{olym08}
\alias{olym12}
\docType{data}
\title{ 2008 and 2012 Summer Olympic Final Medal Count Data}
\description{
  Final medal count, by country, for the Summer
  2008 and 2012 Olympic Games.


}
\usage{
data(olym08)
data(olym12)
}
\format{
  A data frame with 87 or 85 observations on the following 6 variables.
  \describe{
    \item{\code{rank}}{a numeric vector, overall ranking of the countries. }
    \item{\code{country}}{a factor. }
    \item{\code{gold}}{a numeric vector, number of gold medals. }
    \item{\code{silver}}{a numeric vector, number of silver medals. }
    \item{\code{bronze}}{a numeric vector, number of bronze medals. }
    \item{\code{totalmedal}}{a numeric vector, total number of medals. }

%   \item{\code{country}}{a factor. character vector. }

  }
}
\details{
  The events were held during
  (i) August 8--24, 2008, in Beijing; and
  (ii) 27 July--12 August, 2012, in London.


% This is a simple two-way contingency table of counts.

}
% \source{


% url{http://www.associatedcontent.com/article/979484/2008_summer_olympic_medal_count_total.html},
% url{http://www.london2012.com/medals/medal-count/}.


% }
\references{
  The official English website was/is \code{http://en.beijing2008.cn}
  and \code{http://www.london2012.com}.
  Help from Viet Hoang Quoc is gratefully acknowledged.


}
\seealso{
  \code{\link[VGAM]{grc}}.


}

\examples{
summary(olym08)
summary(olym12)
## maybe str(olym08) ; plot(olym08) ...
\dontrun{ par(mfrow = c(1, 2))
myylim <- c(0, 55)
with(head(olym08, n = 8),
barplot(rbind(gold, silver, bronze),
        col = c("gold", "grey", "brown"),  # No "silver" or "bronze"!
#               "gold", "grey71", "chocolate4",
        names.arg = country, cex.names = 0.5, ylim = myylim,
        beside = TRUE, main = "2008 Summer Olympic Final Medal Count",
        ylab = "Medal count", las = 1,
        sub = "Top 8 countries; 'gold'=gold, 'grey'=silver, 'brown'=bronze"))
with(head(olym12, n = 8),
barplot(rbind(gold, silver, bronze),
        col = c("gold", "grey", "brown"),  # No "silver" or "bronze"!
        names.arg = country, cex.names = 0.5, ylim = myylim,
        beside = TRUE, main = "2012 Summer Olympic Final Medal Count",
        ylab = "Medal count", las = 1,
        sub = "Top 8 countries; 'gold'=gold, 'grey'=silver, 'brown'=bronze")) }
}
\keyword{datasets}
