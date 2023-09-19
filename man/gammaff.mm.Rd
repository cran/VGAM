\name{gammaff.mm}
\alias{gammaff.mm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate Gamma Family Function:
  Mathai and Moschopoulos (1992) }
\description{
  Estimate the scale parameter and
  shape parameters
  of the Mathai and Moschopoulos (1992)
  multivariate gamma
  distribution by maximum likelihood estimation.

}
\usage{
gammaff.mm(lscale = "loglink", lshape = "loglink",
           iscale = NULL, ishape = NULL, imethod = 1,
           eq.shapes = FALSE, sh.byrow = TRUE, zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape}{
  Link functions applied to the (positive)
  parameters \eqn{b},
  and \eqn{s_1}{s1}, \ldots, \eqn{s_Q}{sQ}
  respectively.
  See \code{\link{Links}} for more choices.
  In the future, \code{lshapes} might be used instead;
  of course, this link applies to all the shape parameters.


%  where there are \eqn{Q} shape parameters.

  

  }
  \item{iscale, ishape, sh.byrow}{
    Optional initial values.
    The default is to compute them internally.
    Argument \code{sh.byrow} is fed into
    \code{byrow} in \code{\link[base]{matrix}}
    and concerns the ordering of the initial shape
    parameters;
    a matrix of dimension \eqn{n} by \eqn{Q} is
    ultimately constructed.
  See also \code{\link{CommonVGAMffArguments}}.


  }
  \item{eq.shapes}{
    Logical.
  Constrain the shape parameters to be equal?
  See also \code{\link{CommonVGAMffArguments}}.


  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  This distribution has the
  bivariate gamma distribution
  \code{\link[VGAMdata]{bigamma.mckay}}
  as a special case.
  Let \eqn{Q  > 1} be the number of columns of the
  response matrix \code{y}.
  Then the
  joint probability density function is given by
  \deqn{f(y_1,\ldots,y_Q; b, s_1, \ldots, s_Q) =
    y_1^{s_1} (y_2 - y_1)^{s_2}
    \cdots (y_Q - y_{Q-1})^{s_Q}
    \exp(-y_Q / b) / [b^{s_Q^*}
    \Gamma(s_1) \cdots \Gamma(s_Q)]}{%
    f(y1,...,yQ; b, s1, ..., sQ) =
    y1^{s1} (y2 - y1)^{s2}
    ... (yQ - y(Q-1))^{sQ}
    \exp(-yQ / b) / [b^(sQ*)
    Gamma(s1) ... Gamma(sQ)]}
  for \eqn{b > 0},
  \eqn{s_1 > 0}{s1 > 0}, \ldots,
  \eqn{s_Q > 0}{sQ > 0} and
  \eqn{0<y_1< y_2<\cdots<y_Q<\infty}{0<y1<y2<...<yQ<Inf}.
  Also,
  \eqn{s_Q^* = s_1+\cdots+s_Q}{sQ* = s1+...+sQ}.
  Here, \eqn{\Gamma}{Gamma} is
  the \code{\link[base:Special]{gamma}} function,
  By default, the linear/additive predictors are
  \eqn{\eta_1=\log(b)}{eta1=log(b)},
  \eqn{\eta_2=\log(s_1)}{eta2=log(s1)},
  \ldots,
  \eqn{\eta_M=\log(s_Q)}{etaM=log(sQ)}.
  Hence \eqn{Q = M - 1}.
  The marginal distributions are gamma,
  with shape parameters
  \eqn{s_1}{s1} up to \eqn{s_Q}{sQ}, but they have a
  common scale parameter \eqn{b}.



  The fitted value returned
  is a matrix with columns equalling
  their respective means;
  for column \eqn{j} it is
  \code{sum(shape[1:j]) * scale}.



  The correlations are always positive;
  for columns \eqn{j} and \eqn{k}
  with \eqn{j < k},
  the correlation is
  \code{sqrt(sum(shape[1:j]) /sum(shape[1:k]))}.
  Hence the variance of column \eqn{j}
  is \code{sum(shape[1:j]) * scale^2}.



}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}

%% improve the references
\references{


Mathai, A. M. and Moschopoulos, P. G. (1992).
A form of multivariate gamma distribution.
\emph{Ann. Inst. Statist. Math.},
\bold{44}, 97--106.



}
\author{ T. W. Yee }
\note{
  The response must be a matrix with at least
  two columns.
  Apart from the first column,
  the differences between a column and
  its LHS adjacent column
  must all be positive.
  That is, each row must be strictly increasing.



% The data are sorted internally and the user need not
% input the data presorted.



}

\seealso{
  \code{\link[VGAMdata]{bigamma.mckay}},
  \code{\link{gammaff}}.


}
\examples{
\dontrun{
data("mbflood", package = "VGAMdata")
mbflood <- transform(mbflood, VdivD = V / D)
fit <- vglm(cbind(Q, y2 = Q + VdivD) ~ 1,
            gammaff.mm, trace = TRUE, data = mbflood)
coef(fit, matrix = TRUE)
Coef(fit)
vcov(fit)
colMeans(depvar(fit))  # Check moments
head(fitted(fit), 1)
}}
\keyword{models}
\keyword{regression}

