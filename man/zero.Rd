\name{zero}
% \alias{zeroarg}
\alias{zero}
\title{ The zero Argument in VGAM Family Functions }
\description{
  The \code{zero} argument allows users to conveniently
  model certain linear/additive predictors as intercept-only.


}
% \usage{
% VGAMfamilyFunction(zero = 3)
% }
\value{
  Nothing is returned.
  It is simply a convenient argument for constraining
  certain linear/additive predictors to be an intercept only.



}
\section{Warning }{
  The use of other arguments may conflict with the \code{zero}
  argument. For example, using \code{constraints} to input constraint
  matrices may conflict with the \code{zero} argument.
  Another example is the argument \code{parallel}.
  In general users
  should not assume any particular order of precedence when
  there is potential conflict of definition.
  Currently no checking for consistency is made.


  The argument \code{zero} may be renamed in the future to
  something better.


}
\section{Side Effects}{
  The argument creates the appropriate constraint matrices
  internally.


}
\details{
  Often a certain parameter needs to be modelled simply while other
  parameters in the model may be more complex, for example, the
  \eqn{\lambda}{lambda} parameter in LMS-Box-Cox quantile regression
  should be modelled more simply compared to its \eqn{\mu}{mu} parameter.
  Another example is the \eqn{\xi}{xi} parameter in a GEV distribution
  which is should be modelled simpler than its \eqn{\mu}{mu} parameter.
  Using the \code{zero} argument allows this to be fitted conveniently
  without having to input all the constraint matrices explicitly.


  The \code{zero} argument can be assigned an integer vector from the
  set \{\code{1:M}\} where \code{M} is the number of linear/additive
  predictors. Full details about constraint matrices can be found in
  the references.
  See \code{\link{CommonVGAMffArguments}} for more information.




}

\author{T. W. Yee }
\note{
  In all \pkg{VGAM} family functions \code{zero = NULL} means
  none of the linear/additive predictors are modelled as
  intercepts-only.
  Almost all \pkg{VGAM} family function have \code{zero = NULL}
  as the default, but there are some exceptions, e.g.,
  \code{\link{binom2.or}}.


  Typing something like \code{coef(fit, matrix = TRUE)} is a useful
  way to ensure that the \code{zero} argument has worked as expected.


}
\references{

Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


%\url{http://www.stat.auckland.ac.nz/~yee}


}

\seealso{
  \code{\link{CommonVGAMffArguments}},
  \code{\link{constraints}}.


}

\examples{
args(multinomial)
args(binom2.or)
args(gpd)

#LMS quantile regression example
fit <- vglm(BMI ~ sm.bs(age, df = 4), lms.bcg(zero = c(1, 3)),
            data = bmi.nz, trace = TRUE)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
\keyword{programming}


%zz Here is a conflict which is not picked up (no internal checking done)


