\name{loglinb4}
\alias{loglinb4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Loglinear Model for Four Binary Responses }
\description{
  Fits a loglinear model to four binary responses.

}
\usage{
loglinb4(order4 = 4, zero = c("u12", "u13", "u14", "u23",
         "u24", "u34", if (order4 > 2) c("u123", "u124",
         "u134", "u234") else NULL, if (order4 > 3) "u1234"
         else NULL), exchangeable = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exchangeable}{ Logical.
    If \code{TRUE}, the four marginal probabilities
    are constrained to be equal,
    as well as their higher-order interaction terms.

  }
  \item{zero}{ Which linear/additive predictors
    are modelled as intercept-only?
    A \code{NULL} means none.
    See \code{\link{CommonVGAMffArguments}} for
    further information.


  }
  \item{order4}{ Logical or either 2 or 3 or 4.
    If logical and \code{TRUE} then 4, else 3.
    It ends up an integer that is either 2 or 3 or 4.
    Any higher-order term is 0, e.g.,
    \code{u1234 = 0} if \code{order4 = 3}.
  }

  
}
\details{
  The full model is
\eqn{P(Y_1=y_1,Y_2=y_2,Y_3=y_3,Y_4=y_4) =}{P(Y1=y1,Y2=y2,Y3=y3,Y4=y4) =}
\deqn{\exp(u_0+u_1 y_1+u_2 y_2+u_3 y_3+u_4 y_4+
  u_{12} y_1 y_2+
  u_{13} y_1 y_3+
  u_{14} y_1 y_4 +
  u_{23} y_2 y_3 +
  u_{24} y_2 y_4 +
  u_{34} y_3 y_4 +
  u_{123} y_1 y_2 y_3 +
  u_{124} y_1 y_2 y_4 +
  u_{134} y_1 y_3 y_4 +
  u_{234} y_2 y_3 y_4 +
  u_{1234} y_1 y_2 y_3 y_4)}{%
exp(u0 + u1*y1 + u2*y2 + u3*y3 + u4*y4 +
u12*y1*y2 + u13*y1*y3+ u14*y1*y4+
u23*y2*y3 + u24*y2*y4 + u34*y3*y4 +
u123*y1*y2*y3 +
u124*y1*y2*y4 +
u134*y1*y3*y4 +
u234*y2*y3*y4 +
u123*y1*y2*y3 +
u1234*y1*y2*y3*y4)}
where \eqn{y_1}{y1}, \eqn{y_2}{y2} and
\eqn{y_3}{y3} and
\eqn{y_4}{y4} are 0
or 1, and the parameters are \eqn{u_1}{u1},
\eqn{u_2}{u2},
\eqn{u_3}{u3},
\eqn{u_4}{u4},
\eqn{u_{12}}{u12},
\eqn{u_{13}}{u13},
\eqn{u_{14}}{u14},
\eqn{u_{23}}{u23},
\eqn{u_{24}}{u24},
\eqn{u_{34}}{u34}, and if
\code{order4 >= 3} then
\eqn{u_{123}}{u123},
\eqn{u_{124}}{u124},
\eqn{u_{134}}{u134},
\eqn{u_{234}}{u234}, too, and if
\code{order4 == 4} then
\eqn{u_{1234}}{u1234} too.
  The normalizing parameter \eqn{u_0}{u0} can
  be expressed as a function of the others.
  The the parameters are estimated by
  \code{\link{identitylink}}.
  Unlike \code{\link{loglinb3}},
  a fourth-order (full) association parameter,
  \eqn{u_{1234}}{u1234} for the product
  \eqn{y_1 y_2 y_3 y_4}{y1*y2*y3*y4},
  is \emph{not} assumed to be zero for this
  family function by default.
  Note the default for this argument might change
  in the future.
  If the data cannot support such a high
  order interaction term then reduce \code{order4}.

  
  The linear/additive predictors are,
  for the full model,
  \eqn{(\eta_1,\eta_2,\ldots,\eta_{15})^T =
(u_1,u_2,u_3,u_4,u_{12},u_{13},\ldots,u_{1234})^T}{(eta1,eta2,...,eta15) =
(u1,u2,u3,u4,u12,u13,u14,u23,u24,u34,u123,u124,u134,u234,u1234)}.
The ordering agrees with \code{\link[utils]{combn}}
piecemeal.



% By default, the last element is not there since
% \code{u123.arg = FALSE}.


}

\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.



  When fitted,
  the \code{fitted.values} slot of the object
  contains the joint probabilities,
  labelled as
  \eqn{(Y_1,Y_2,Y_3,Y_4)}{(Y1,Y2,Y3,Y4)} =
  (0,0,0,0), (0,0,0,1), (0,0,1,0),
  (0,0,1,1), (0,1,0,0), (0,1,0,1), (0,1,1,0),
  (0,1,1,1),
  (1,0,0,0), (1,0,0,1), (1,0,1,0),
  (1,0,1,1), (1,1,0,0), (1,1,0,1), (1,1,1,0),
  (1,1,1,1),
  respectively.


}

%references{
%}
\author{ 
  Yunhao (Harry) Han wrote
  \code{@deriv} and \code{@weight},
  and Thomas Yee wrote the rest.

}
\note{
  The response must be a 4-column matrix
  of ones and zeros only.
  By default, each of the 16 combinations of
  the multivariate
  response need to appear in the data set,
  therefore data sets will usually
  need to be large in order for this
  family function to work most simply.
  As stated above, reduce \code{order4}
  if there are problems.
  After estimation, the response attached
  to the object is also
  a 4-column matrix; possibly in the future
  it might change into
  a 16-column matrix.




}

\seealso{
  \code{\link{loglinb3}},
  \code{\link{loglinb2}},
  \code{\link[utils]{combn}},
  \code{\link{hunua}}.


}
\examples{
lfit4 <- vglm(cbind(cyadea, beitaw, kniexc, vitluc) ~ altitude,
              loglinb4, hunua, trace = TRUE)
coef(lfit4, matrix = TRUE)
head(fitted(lfit4))
head(predict(lfit4))
summary(lfit4, HDEtest = FALSE)
}
\keyword{models}
\keyword{regression}
