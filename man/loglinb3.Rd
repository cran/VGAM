\name{loglinb3}
\alias{loglinb3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Loglinear Model for Three Binary Responses }
\description{
  Fits a loglinear model to three binary responses.

}
\usage{
loglinb3(exchangeable = FALSE, zero = c("u12", "u13", "u23",
         if (u123.arg) "u123" else NULL),
         u123.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exchangeable}{ Logical.
    If \code{TRUE}, the three marginal probabilities
    are constrained to be equal,
    as well as their second-order interaction terms.

  }
  \item{zero}{ Which linear/additive predictors
    are modelled as intercept-only?
    A \code{NULL} means none.
    See \code{\link{CommonVGAMffArguments}} for
    further information.


  }
  \item{u123.arg}{ Logical.
    Include the 3rd-order interaction \code{u123}?
    Note that \code{u123.arg = TRUE} might
    become the default some time in the future.


  }

  
}
\details{
  The full model is
\eqn{P(Y_1=y_1,Y_2=y_2,Y_3=y_3) =}{P(Y1=y1,Y2=y2,Y3=y3) =}
\deqn{\exp(u_0+u_1 y_1+u_2 y_2+u_3 y_3+u_{12} y_1 y_2+
  u_{13} y_1 y_3+u_{23} y_2 y_3 +
  u_{123} y_1 y_2 y_3 )}{%
  exp(u0 + u1*y1 + u2*y2 + u3*y3 + u12*y1*y2 + u13*y1*y3+
  u23*y2*y3+ u123*y1*y2*y3)}
where \eqn{y_1}{y1}, \eqn{y_2}{y2} and
\eqn{y_3}{y3} are 0
or 1, and the parameters are \eqn{u_1}{u1},
\eqn{u_2}{u2},
\eqn{u_3}{u3}, \eqn{u_{12}}{u12},
\eqn{u_{13}}{u13},
  \eqn{u_{23}}{u23}, and if
  \code{u123.arg} then \eqn{u_{123}}{u123} too.
  The normalizing parameter \eqn{u_0}{u0} can
  be expressed as a
  function of the other parameters.
  The the parameters are estimated by
  \code{\link{identitylink}}.
  Note that a third-order association parameter,
  \eqn{u_{123}}{u123} for the product
  \eqn{y_1 y_2 y_3}{y1*y2*y3},
  is assumed to be zero for this family function
  by default;
  it is estimated if \code{u123.arg} is \code{TRUE}.
  Note the default for this argument might change
  in the future.


  

  The linear/additive predictors are,
  for the full model,
  \eqn{(\eta_1,\eta_2,\ldots,\eta_6,\eta_7)^T =
(u_1,u_2,u_3,u_{12},u_{13},u_{23},u_{123})^T}{(eta1,eta2,...,eta6,eta7) =
(u1,u2,u3,u12,u13,u23,u123)}.
  By default, the last element is not there since
  \code{u123.arg = FALSE}.



}

\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.



  When fitted,
  the \code{fitted.values} slot of the object
  contains the eight joint probabilities,
  labelled as
  \eqn{(Y_1,Y_2,Y_3)}{(Y1,Y2,Y3)} = (0,0,0),
  (0,0,1), (0,1,0),
  (0,1,1), (1,0,0), (1,0,1), (1,1,0),
  (1,1,1), respectively.


}

\references{

Yee, T. W. and Wild, C. J. (2001).
Discussion to:
``Smoothing spline ANOVA for multivariate
Bernoulli observations, with application
to ophthalmology data
(with discussion)''
by Gao, F., Wahba, G., Klein, R., Klein, B.
\emph{Journal of the American Statistical
  Association},
\bold{96}, 127--160.


McCullagh, P. and Nelder, J. A. (1989).
\emph{Generalized Linear Models}, 2nd ed.
London: Chapman & Hall.


%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}
\author{ Thomas W. Yee and
  Yunhao (Harry) Han who added \code{u123}.

}
\note{
  The response must be a 3-column matrix
  of ones and zeros only.
  Note that each of the 8 combinations of
  the multivariate
  response need to appear in the data set,
  therefore data sets
  will need to be large in order for this
  family function to work.
  After estimation, the response attached
  to the object is also
  a 3-column matrix; possibly in the future
  it might change into
  a 8-column matrix.




}

\seealso{
  \code{\link{binom3.or}},
  \code{\link{loglinb2}},
  \code{\link{loglinb4}},
  \code{\link{binom2.or}},
  \code{\link{hunua}}.


}
\examples{
lfit1 <- vglm(cbind(cyadea, beitaw, kniexc) ~ altitude,
             loglinb3, data = hunua, trace = TRUE)
coef(lfit1, matrix = TRUE)
lfit2 <- vglm(cbind(cyadea, beitaw, kniexc) ~ altitude,
             loglinb3(u123 = TRUE), hunua, trace = TRUE)
coef(lfit2, matrix = TRUE)
head(fitted(lfit2))
summary(lfit2)
}
\keyword{models}
\keyword{regression}
