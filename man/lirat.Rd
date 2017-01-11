\name{lirat}
\alias{lirat}
\docType{data}
\title{ Low-iron Rat Teratology Data }
\description{
  Low-iron rat teratology data.


}
\usage{data(lirat)}
\format{
  A data frame with 58 observations on the following 4 variables.
  \describe{
    \item{\code{N}}{Litter size.}
    \item{\code{R}}{Number of dead fetuses.}
    \item{\code{hb}}{Hemoglobin level.}
    \item{\code{grp}}{Group number.
          Group 1 is the untreated (low-iron) group,
          group 2 received injections on day 7 or day 10 only,
          group 3 received injections on days 0 and 7, and
          group 4 received injections weekly.}
  }
}
\details{
The following description comes from Moore and Tsiatis (1991).  The data
comes from the experimental setup from Shepard et al. (1980), which is
typical of studies of the effects of chemical agents or dietary regimens
on fetal development in laboratory rats.


Female rats were put in iron-deficient diets and divided into 4
groups. One group of controls was given weekly injections of iron
supplement to bring their iron intake to normal levels, while another
group was given only placebo injections.  Two other groups were given
fewer iron-supplement injections than the controls.  The rats were made
pregnant, sacrificed 3 weeks later, and the total number of fetuses and
the number of dead fetuses in each litter were counted.


For each litter the number of dead fetuses may be considered to be
Binomial(\eqn{N,p}) where \eqn{N} is the litter size and \eqn{p}
is the probability of a fetus dying. The parameter \eqn{p} is expected
to vary from litter to litter, therefore the total variance of the
proportions will be greater than that predicted by a binomial model,
even when the covariates for hemoglobin level and experimental group
are accounted for.


}
\source{
  Moore, D. F. and Tsiatis, A. (1991)
  Robust Estimation of the Variance in Moment Methods for
  Extra-binomial and Extra-Poisson Variation.
  \emph{Biometrics},
  \bold{47}, 383--401.


}
\references{
   Shepard, T. H., Mackler, B. and Finch, C. A. (1980)
   Reproductive studies in the iron-deficient rat.
   \emph{Teratology}, \bold{22}, 329--334.


}
\examples{
\dontrun{
# cf. Figure 3 of Moore and Tsiatis (1991)
plot(R / N ~ hb, data = lirat, pch = as.character(grp), col = grp,
     las = 1, xlab = "Hemoglobin level", ylab = "Proportion Dead") }
}
\keyword{datasets}
