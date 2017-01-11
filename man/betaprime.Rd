\name{betaprime}
\alias{betaprime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Beta-Prime Distribution }
\description{
  Estimation of the two shape parameters of the beta-prime distribution
  by maximum likelihood estimation.


}
\usage{
betaprime(lshape = "loge", ishape1 = 2, ishape2 = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape}{
  Parameter link function applied to the two (positive) shape parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape1, ishape2, zero}{
  See
  \code{\link{CommonVGAMffArguments}}.


% Initial values for the first and second shape parameters.
% A \code{NULL} value means it is obtained in the \code{initialize} slot.
% Note that \code{ishape2} is obtained using \code{ishape1}.


  }
% \item{zero}{
% An integer-valued vector specifying which linear/additive predictors
% are modelled as intercepts only.  The value must be from the set
% \{1,2\} corresponding respectively to \code{shape1} and \code{shape2}
% respectively.  If \code{zero=NULL} then both parameters are modelled
% with the explanatory variables.
% }
}

%% what is the mean if shape2 < 1?

\details{
  The beta-prime distribution is given by
  \deqn{f(y) = y^{shape1-1}   (1+y)^{-shape1-shape2} / B(shape1,shape2)}{%
        f(y) = y^(shape1-1) * (1+y)^(-shape1-shape2) / B(shape1,shape2) }
  for \eqn{y > 0}.
  The shape parameters are positive, and
  here, \eqn{B} is the beta function.
  The mean of \eqn{Y} is \eqn{shape1 / (shape2-1)} provided \eqn{shape2>1};
  these are returned as the fitted values.



  If \eqn{Y} has a \eqn{Beta(shape1,shape2)} distribution then
  \eqn{Y/(1-Y)} and \eqn{(1-Y)/Y} have a \eqn{Betaprime(shape1,shape2)}
  and \eqn{Betaprime(shape2,shape1)} distribution respectively.
  Also, if \eqn{Y_1}{Y1} has a \eqn{gamma(shape1)} distribution
  and \eqn{Y_2}{Y2} has a \eqn{gamma(shape2)} distribution
  then \eqn{Y_1/Y_2}{Y1/Y2} has a \eqn{Betaprime(shape1,shape2)}
  distribution.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}

%% zz not sure about the JKB reference.
\references{

Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
Chapter 25 of:
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 2,
New York: Wiley.


%Documentation accompanying the \pkg{VGAM} package at
%\url{https://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}
\author{ Thomas W. Yee }
\note{
  The response must have positive values only.


  The beta-prime distribution is also known as the
  \emph{beta distribution of the second kind} or the
  \emph{inverted beta distribution}.


}

\seealso{
  \code{\link{betaff}},
  \code{\link[stats]{Beta}}.


}
\examples{
nn <- 1000
bdata <- data.frame(shape1 = exp(1), shape2 = exp(3))
bdata <- transform(bdata, yb = rbeta(nn, shape1, shape2))
bdata <- transform(bdata, y1 = (1-yb) /    yb,
                          y2 =    yb  / (1-yb),
                          y3 = rgamma(nn, exp(3)) / rgamma(nn, exp(2)))

fit1 <- vglm(y1 ~ 1, betaprime, data = bdata, trace = TRUE)
coef(fit1, matrix = TRUE)

fit2 <- vglm(y2 ~ 1, betaprime, data = bdata, trace = TRUE)
coef(fit2, matrix = TRUE)

fit3 <- vglm(y3 ~ 1, betaprime, data = bdata, trace = TRUE)
coef(fit3, matrix = TRUE)

# Compare the fitted values
with(bdata, mean(y3))
head(fitted(fit3))
Coef(fit3)  # Useful for intercept-only models
}
\keyword{models}
\keyword{regression}

