\name{dirmultinomial}
\alias{dirmultinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a Dirichlet-Multinomial Distribution }
\description{
  Fits a Dirichlet-multinomial distribution to a matrix response.

}
\usage{
dirmultinomial(lphi = "logit", iphi = 0.10, parallel = FALSE, zero = "M")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi}{
  Link function applied to the \eqn{\phi}{phi}
  parameter, which lies in the open unit interval \eqn{(0,1)}.
  See \code{\link{Links}} for more choices.


  }
  \item{iphi}{
  Numeric. Initial value for \eqn{\phi}{phi}.
  Must be in the open unit interval \eqn{(0,1)}.
  If a failure to converge occurs then try assigning this argument
  a different value.


  }
  \item{parallel}{
  A logical (formula not allowed here) indicating whether the
  probabilities \eqn{\pi_1,\ldots,\pi_{M-1}}{pi_1,\ldots,pi_{M-1}}
  are to be equal via equal coefficients.
  Note \eqn{\pi_M}{pi_M} will generally be different from the
  other probabilities.
  Setting \code{parallel = TRUE} will only work if you also set
  \code{zero = NULL} because of interference between these arguments
  (with respect to the intercept term).

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \eqn{\{1,2,\ldots,M\}}.
  If the character \code{"M"} then this means the numerical value
  \eqn{M}, which corresponds to linear/additive predictor associated
  with \eqn{\phi}{phi}.
  Setting \code{zero = NULL} means none of the values from
  the set \eqn{\{1,2,\ldots,M\}}.

  }
}
\details{
  The Dirichlet-multinomial distribution
  arises from a multinomial distribution where
  the probability parameters are not constant but are generated from a
  multivariate distribution called the Dirichlet distribution.
  The Dirichlet-multinomial distribution has probability function
  \deqn{P(Y_1=y_1,\ldots,Y_M=y_M) =
     {N_{*} \choose {y_1,\ldots,y_M}}
        \frac{
        \prod_{j=1}^{M}
        \prod_{r=1}^{y_{j}}
        (\pi_j (1-\phi) + (r-1)\phi)}{
        \prod_{r=1}^{N_{*}}
        (1-\phi + (r-1)\phi)}}{%
        P(Y_1=y_1,\ldots,Y_M=y_M) = C_{y_1,\ldots,y_M}^{N_{*}}
        prod_{j=1}^{M}
        prod_{r=1}^{y_{j}}
        (pi_j (1-phi) + (r-1)phi) /
        prod_{r=1}^{N_{*}}
        (1-phi + (r-1)phi)}
  where \eqn{\phi}{phi} is the \emph{over-dispersion} parameter
  and \eqn{N_{*} = y_1+\cdots+y_M}{N_* = y_1+\cdots+y_M}.
  Here, \eqn{a \choose b}{C_b^a} means ``\eqn{a} choose \eqn{b}'' and
  refers to combinations (see \code{\link[base]{choose}}).
  The above formula applies to each row of the matrix response.
  In this \pkg{VGAM} family function the first \eqn{M-1} linear/additive
  predictors correspond to the first \eqn{M-1} probabilities via
  \deqn{\eta_j = \log(P[Y=j]/ P[Y=M]) = \log(\pi_j/\pi_M)}{%
         eta_j =  log(P[Y=j]/ P[Y=M]) = log(pi_j/pi_M)}
  where \eqn{\eta_j}{eta_j} is the \eqn{j}th linear/additive predictor
  (\eqn{\eta_M=0}{eta_M=0} by definition for \eqn{P[Y=M]} but not for
   \eqn{\phi}{phi})
  and
  \eqn{j=1,\ldots,M-1}.
  The \eqn{M}th linear/additive predictor corresponds to \code{lphi}
  applied to \eqn{\phi}{phi}.


  Note that \eqn{E(Y_j) = N_* \pi_j}{E(Y_j) = N_* pi_j} but
  the probabilities (returned as the fitted values)
  \eqn{\pi_j}{pi_j} are bundled together as a \eqn{M}-column matrix.
  The quantities \eqn{N_*} are returned as the prior weights.


  The beta-binomial distribution is a special case of
  the Dirichlet-multinomial distribution when \eqn{M=2};
  see \code{\link{betabinomial}}.  It is easy to show that
  the first shape parameter of the beta distribution is
  \eqn{shape1=\pi(1/\phi-1)}{shape1=pi*(1/phi-1)} and the second shape
  parameter is \eqn{shape2=(1-\pi)(1/\phi-1)}{shape2=(1-pi)*(1/phi-1)}.
  Also, \eqn{\phi=1/(1+shape1+shape2)}{phi=1/(1+shape1+shape2)}, which
  is known as the \emph{intra-cluster correlation} coefficient.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


  If the model is an intercept-only model then \code{@misc} (which is a
  list) has a component called \code{shape} which is a vector with the
  \eqn{M} values \eqn{\pi_j(1/\phi-1)}{pi_j * (1/phi-1)}.


% zz not sure: These are the shape parameters of the underlying
% Dirichlet distribution.


}
\references{


Paul, S. R., Balasooriya, U. and Banerjee, T. (2005)
Fisher information matrix of the Dirichlet-multinomial distribution.
\emph{Biometrical Journal}, \bold{47}, 230--236.



Tvedebrink, T. (2010)
Overdispersion in allelic counts and \eqn{\theta}-correction in
forensic genetics.
\emph{Theoretical Population Biology}, \bold{78}, 200--210.



Yu, P. and Shaw, C. A. (2014).
An Efficient Algorithm for Accurate Computation of
the Dirichlet-Multinomial Log-Likelihood Function.
\emph{Bioinformatics},
\bold{30}, 1547--54.



% url {doi:10.1093/bioinformatics/btu079}.
% number =      {11},



% first published online February 11, 2014





}
\author{ Thomas W. Yee }

\section{Warning }{
  This \pkg{VGAM} family function is prone to numerical problems,
  especially when there are covariates.


}
\note{
  The response can be a matrix of non-negative integers, or else
  a matrix of sample proportions and the total number of counts in
  each row specified using the \code{weights} argument.
  This dual input option is similar to \code{\link{multinomial}}.


  To fit a `parallel' model with the \eqn{\phi}{phi} parameter
  being an intercept-only you will need to use the \code{constraints}
  argument.


  Currently, Fisher scoring is implemented. To compute the expected
  information matrix a \code{for} loop is used; this may be very slow
  when the counts are large.
  Additionally, convergence may be slower than usual due to round-off
  error when computing the expected information matrices.


}

\seealso{
  \code{\link{dirmul.old}},
  \code{\link{betabinomial}},
  \code{\link{betabinomialff}},
  \code{\link{dirichlet}},
  \code{\link{multinomial}}.


}

\examples{
nn <- 5; M <- 4; set.seed(1)
ydata <- data.frame(round(matrix(runif(nn * M, max = 100), nn, M)))  # Integer counts
colnames(ydata) <- paste("y", 1:M, sep = "")

fit <- vglm(cbind(y1, y2, y3, y4) ~ 1, dirmultinomial,
            data = ydata, trace = TRUE)
head(fitted(fit))
depvar(fit)  # Sample proportions
weights(fit, type = "prior", matrix = FALSE)  # Total counts per row

\dontrun{
ydata <- transform(ydata, x2 = runif(nn))
fit <- vglm(cbind(y1, y2, y3, y4) ~ x2, dirmultinomial,
            data = ydata, trace = TRUE)
Coef(fit)
coef(fit, matrix = TRUE)
(sfit <- summary(fit))
vcov(sfit)
}
}
\keyword{models}
\keyword{regression}


% zz \eqn{\alpha_j = P[Y=j] \times (1/\phi - 1)}{alpha_j = P[Y=j] *
% (1/phi - 1)} are the shape parameters,
% for \eqn{j=1,\ldots,M}.

% Currently, initial values can be improved upon.


% \dontrun{ # This does not work:





