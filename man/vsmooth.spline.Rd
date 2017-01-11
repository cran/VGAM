\name{vsmooth.spline}
\alias{vsmooth.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Vector cubic smoothing spline }
\description{
  Fits a vector cubic smoothing spline.
}
\usage{
vsmooth.spline(x, y, w = NULL, df = rep(5, M), spar = NULL,
               i.constraint = diag(M),
               x.constraint = diag(M),
               constraints = list("(Intercepts)" = i.constraint,
                                  x = x.constraint),
               all.knots = FALSE, var.arg = FALSE, scale.w = TRUE,
               nk = NULL, control.spar = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A vector, matrix or a list.
  If a list, the  \code{x} component is used.
  If a matrix, the first column is used.
  \code{x} may also be a complex vector, in which case
  the real part is used, and the imaginary part is used for the response.
  In this help file, \code{n} is the number of unique values of \code{x}.


}
  \item{y}{
  A vector, matrix or a list.
  If a list, the  \code{y} component is used.
  If a matrix, all but the first column is used.
  In this help file, \code{M} is the number of columns of \code{y} if
  there are no constraints on the functions.


}
  \item{w}{
  The weight matrices or the number of observations.
  If the  weight matrices, then this must be a \code{n}-row matrix
  with the elements in matrix-band form (see \code{iam}).
  If a vector, then these are the number of observations.
  By default, \code{w} is the \code{M} by \code{M} identity
  matrix, denoted by  \code{matrix(1, n, M)}.


}
  \item{df}{
Numerical vector containing the degrees of
freedom for each component function (smooth).
If necessary, the vector is recycled to have length equal
to the number of component functions to be estimated
(\code{M} if there are no constraints), which
equals the number of columns of the \code{x}-constraint matrix.
A value of 2 means a linear fit, and each element of
\code{df} should lie between 2 and \code{n}.
The larger the values of \code{df} the more wiggly the
smooths.


}
  \item{spar}{
Numerical vector containing the non-negative smoothing
parameters for each component function (smooth).
If necessary, the vector is recycled to have length equal
to the number of component functions to be estimated
(\code{M} if there are no constraints), which
equals the number of columns of the \code{x}-constraint matrix.
A value of zero means the smooth goes through the data and hence
is wiggly.
A value of \code{Inf} may be assigned, meaning the smooth will
be linear.
By default, the \code{NULL} value of \code{spar} means
\code{df} is used to determine the smoothing
parameters.


}
  \item{all.knots}{ Logical. If \code{TRUE} then each distinct
value of \code{x} will be a knot. By default, only a subset of
the unique values  of \code{x} are used; typically, the number
of knots is \code{O(n^0.25)} for \code{n} large,
but if \code{n <= 40} then all the unique values of \code{x} are used.


}
  \item{i.constraint}{ A \code{M}-row constraint matrix for the
intercepts. It must be of full column rank.
By default, the constraint matrix for the intercepts is the
\code{M} by \code{M} identity matrix, meaning no constraints.


}
  \item{x.constraint}{ A \code{M}-row constraint matrix for \code{x}.
It must be of full column rank.
By default, the constraint matrix for the intercepts is the
\code{M} by \code{M} identity matrix, meaning no constraints.


}
  \item{constraints}{
An alternative to specifying \code{i.constraint} and \code{x.constraint},
this is a list with two components corresponding to the
intercept and \code{x} respectively. They must both be a
\code{M}-row constraint matrix with full column rank.


}
  \item{var.arg}{ Logical: return the pointwise variances
of the fit?
Currently, this corresponds only to the nonlinear part of the
fit, and may be wrong.


}
  \item{scale.w}{
Logical.
By default, the weights \code{w} are scaled so that the
diagonal elements have mean 1.


}
  \item{nk}{ Number of knots.
If used, this argument overrides \code{all.knots}, and
must lie between 6 and \code{n}+2 inclusive.


}
  \item{control.spar}{
See \code{\link[stats]{smooth.spline}}.


}
}
\details{
  The algorithm implemented is detailed in Yee (2000).
  It involves decomposing the component functions into a linear and
  nonlinear part, and using B-splines.
  The cost of the computation is \code{O(n M^3)}.


  The argument \code{spar} contains \emph{scaled} smoothing parameters.


}
\value{
  An object of class \code{"vsmooth.spline"}
  (see \code{vsmooth.spline-class}).


}
\references{
Yee, T. W. (2000)
Vector Splines and Other Vector Smoothers.
Pages 529--534.
In: Bethlehem, J. G. and van der Heijde, P. G. M.
\emph{Proceedings in Computational Statistics COMPSTAT 2000}.
Heidelberg: Physica-Verlag.


}
\author{ Thomas W. Yee }
\note{
  This function is quite similar to \code{\link[stats]{smooth.spline}}
  but offers less functionality.
  For example, cross validation is not implemented here.
  For \code{M = 1}, the results will be generally different,
  mainly due to the different way the knots are selected.


  The vector cubic smoothing spline which \code{s()} represents is
  computationally demanding for large \eqn{M}.
  The cost is approximately \eqn{O(n M^3)} where
  \eqn{n} is the number of unique abscissae.


  Yet to be done: return the \emph{unscaled} smoothing parameters.


}

%~Make other sections like WARNING with \section{WARNING }{....} ~
\section{WARNING}{
  See \code{\link{vgam}} for information about an important bug.


}



\seealso{
  \code{vsmooth.spline-class},
  \code{plot.vsmooth.spline},
  \code{predict.vsmooth.spline},
  \code{iam},
  \code{\link{sm.os}},
  \code{\link[VGAM]{s}},
  \code{\link[stats]{smooth.spline}}.


}
\examples{
nn <- 20; x <- 2 + 5*(nn:1)/nn
x[2:4] <- x[5:7]  # Allow duplication
y1 <- sin(x) + rnorm(nn, sd = 0.13)
y2 <- cos(x) + rnorm(nn, sd = 0.13)
y3 <- 1 + sin(x) + rnorm(nn, sd = 0.13)  # Run this for constraints
y <- cbind(y1, y2, y3)
ww <- cbind(rep(3, nn), 4, (1:nn)/nn)

(fit <- vsmooth.spline(x, y, w = ww, df = 5))
\dontrun{
plot(fit)  # The 1st and 3rd functions do not differ by a constant
}

mat <- matrix(c(1,0,1, 0,1,0), 3, 2)
(fit2 <- vsmooth.spline(x, y, w = ww, df = 5, i.constr = mat,
                        x.constr = mat))
# The 1st and 3rd functions do differ by a constant:
mycols <- c("orange", "blue", "orange")
\dontrun{ plot(fit2, lcol = mycols, pcol = mycols, las = 1) }

p <- predict(fit, x = model.matrix(fit, type = "lm"), deriv = 0)
max(abs(depvar(fit) - with(p, y)))  # Should be 0; and fit@y is not good

par(mfrow = c(3, 1))
ux <- seq(1, 8, len = 100)
for (dd in 1:3) {
  pp <- predict(fit, x = ux, deriv = dd)
\dontrun{with(pp, matplot(x, y, type = "l", main = paste("deriv =", dd),
                          lwd = 2, ylab = "", cex.axis = 1.5,
                          cex.lab = 1.5, cex.main = 1.5)) }
}
}
\keyword{regression}
\keyword{smooth}
