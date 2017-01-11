\name{waitakere}
\alias{waitakere}
\docType{data}
\title{Waitakere Ranges Data}
\description{
  The \code{waitakere} data frame has 579 rows and 18 columns.
  Altitude is explanatory, and there are binary responses
  (presence/absence = 1/0 respectively) for 17 plant species.

}
\usage{data(waitakere)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{agaaus}{Agathis australis, or Kauri}
    \item{beitaw}{Beilschmiedia tawa, or Tawa}
    \item{corlae}{Corynocarpus laevigatus}
    \item{cyadea}{Cyathea dealbata}
    \item{cyamed}{Cyathea medullaris}
    \item{daccup}{Dacrydium cupressinum}
    \item{dacdac}{Dacrycarpus dacrydioides}
    \item{eladen}{Elaecarpus dentatus}
    \item{hedarb}{Hedycarya arborea}
    \item{hohpop}{Species name unknown}
    \item{kniexc}{Knightia excelsa, or Rewarewa}
    \item{kuneri}{Kunzea ericoides}
    \item{lepsco}{Leptospermum scoparium}
    \item{metrob}{Metrosideros robusta}
    \item{neslan}{Nestegis lanceolata}
    \item{rhosap}{Rhopalostylis sapida}
    \item{vitluc}{Vitex lucens, or Puriri}
    \item{altitude}{meters above sea level}
  }
}
\details{
  These were collected from the Waitakere Ranges, a small forest in northern
  Auckland, New Zealand. At 579 sites in the forest, the presence/absence
  of 17 plant species was recorded, as well as the altitude.
  Each site was of area size 200\eqn{m^2}{m^2}.


}
\source{
  Dr Neil Mitchell, University of Auckland.

}
%\references{
%None.
%}
\seealso{
    \code{\link{hunua}}.

}
\examples{
fit <- vgam(agaaus ~ s(altitude, df = 2), binomialff, waitakere)
head(predict(fit, waitakere, type = "response"))
\dontrun{ plot(fit, se = TRUE, lcol = "orange", scol = "blue") }
}
\keyword{datasets}
