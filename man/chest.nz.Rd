\name{chest.nz}
\alias{chest.nz}
\docType{data}
\title{ Chest Pain in NZ Adults Data}
\description{
  Presence/absence of chest pain in 10186 New Zealand adults.
}
\usage{data(chest.nz)}
\format{
  A data frame with 73 rows and the following 5 variables.
  \describe{
    \item{age}{a numeric vector; age (years).}
    \item{nolnor}{a numeric vector of counts; no pain on LHS or RHS.}
    \item{nolr}{a numeric vector of counts; no pain on LHS but pain on RHS.}
    \item{lnor}{a numeric vector of counts; no pain on RHS but pain on LHS.}
    \item{lr}{a numeric vector of counts; pain on LHS and RHS of chest.}
  }
}
\details{
  Each adult was asked their age and whether they experienced any
  pain or discomfort in their chest over the last six months.
  If yes, they indicated whether it
  was on their LHS and/or RHS of their chest.


}
\source{
  MacMahon, S., Norton, R., Jackson, R., Mackie, M. J.,
  Cheng, A., Vander Hoorn, S., Milne, A., McCulloch, A. (1995)
  Fletcher Challenge-University of Auckland Heart &
  Health Study: design and baseline findings.
  \emph{New Zealand Medical Journal},
  \bold{108}, 499--502.


}
\examples{
\dontrun{
fit <- vgam(cbind(nolnor, nolr, lnor, lr) ~ s(age, c(4, 3)),
            binom2.or(exchan = TRUE, zero = NULL), data = chest.nz)
coef(fit, matrix = TRUE)
}
\dontrun{ plot(fit, which.cf = 2, se = TRUE) }
}
\keyword{datasets}


