\name{hunua}
\alias{hunua}
\docType{data}
\title{Hunua Ranges Data}
\description{
  The \code{hunua} data frame has 392 rows and 18 columns.
  Altitude is explanatory, and there are binary responses 
  (presence/absence = 1/0 respectively) for 17 plant species.

}
\usage{data(hunua)}
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
  These were collected from the Hunua Ranges, a small forest in southern
  Auckland, New Zealand. At 392 sites in the forest, the presence/absence
  of 17 plant species was recorded, as well as the altitude. 
  Each site was of area size 200\eqn{m^2}{m^2}.

}
\source{
  Dr Neil Mitchell, University of Auckland. 
}
%\references{
%  None. 
%}
\seealso{
    \code{\link{waitakere}}.
}
\examples{
# Fit a GAM using vgam() and compare it with the Waitakere Ranges one
data(hunua)
fit.h = vgam(agaaus ~ s(altitude, df=2), binomialff, hunua)
\dontrun{
plot(fit.h, se=TRUE, lcol="red", scol="red",
     main="Red is Hunua, Blue is Waitakere") }
predict(fit.h, hunua, type="response")[1:3]

data(waitakere)
fit.w = vgam(agaaus ~ s(altitude, df=2), binomialff, waitakere)
\dontrun{
plot(fit.w, se=TRUE, lcol="blue", scol="blue", add=TRUE) }
predict(fit.w, hunua, type="response")[1:3]        # Same as above? 
}
\keyword{datasets}
