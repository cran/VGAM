\name{olympic}
\alias{olympic}
\docType{data}
\title{ 2008 Summer Olympic Final Medal Count Data}
\description{
  Final count of medal winners by country for the 2008 Summer Olympic
  games in Beijing.


}
\usage{data(olympic)}
\format{
  A data frame with 87 observations on the following 6 variables.
  \describe{
    \item{\code{rank}}{a numeric vector, overall ranking of the countries. }
    \item{\code{country}}{a character vector. }
    \item{\code{gold}}{a numeric vector, number of gold medals. }
    \item{\code{silver}}{a numeric vector, number of silver medals. }
    \item{\code{bronze}}{a numeric vector, number of bronze medals. }
    \item{\code{totalmedal}}{a numeric vector, total number of medals. }
  }
}
\details{
  The event was held during August 8--24, 2008, in Beijing.


% This is a simple two-way contingency table of counts.

}
\source{
\url{http://www.associatedcontent.com/article/979484/2008_summer_olympic_medal_count_total.html}.

}
\references{
  The official English webite was
  \url{http://en.beijing2008.cn}. 

}
\seealso{
  \code{\link{grc}}.


}

\examples{
summary(olympic)
## maybe str(olympic) ; plot(olympic) ...
\dontrun{ with(head(olympic, n = 8),
barplot(rbind(gold,silver,bronze),
        col = c("gold","grey","brown"),  # No "silver" or "bronze"!
        names.arg = country, cex.names = 0.5,
        beside = TRUE, main = "2008 Summer Olympic Final Medal Count",
        ylab = "Medal count", las = 1,
        sub = "Top 8 countries; 'gold'=gold, 'grey'=silver, 'brown'=bronze")) }
}
\keyword{datasets}
