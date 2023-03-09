\name{genpoisson0}
\alias{genpoisson0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Poisson Regression
  (Original Parameterization) }
\description{
  Estimation of the two-parameter generalized Poisson distribution
  (original parameterization).

}
\usage{
genpoisson0(ltheta = "loglink", llambda = "logitlink",
            itheta = NULL, ilambda = NULL, imethod = c(1, 1),
            ishrinkage = 0.95, glambda = ppoints(5),
            parallel = FALSE, zero = "lambda")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ltheta, llambda}{
  Parameter link functions for \eqn{\theta} and \eqn{\lambda}.
  See \code{\link{Links}} for more choices.
  In theory the \eqn{\lambda} parameter is allowed to be negative to
  handle underdispersion, however this is no longer supported,
  hence \eqn{0 < \lambda < 1}.
  The \eqn{\theta} parameter is positive, therefore the default is the
  log link.


  % lies at least within the interval \eqn{[-1,1]};
  % see below for more details,
  % and an alternative link is \code{\link{rhobitlink}}.


  }
  \item{itheta, ilambda}{
  Optional initial values for \eqn{\lambda} and \eqn{\theta}.
  The default is to choose values internally.


  }
% \item{use.approx}{
% Logical. If \code{TRUE} then an approximation to the expected
% information matrix is used, otherwise Newton-Raphson is used.


% }
  \item{imethod}{
    See \code{\link{CommonVGAMffArguments}} for information.
    Each value is an integer \code{1} or \code{2} or \code{3} which
  specifies the initialization method for each of the parameters.
  If failure to converge occurs try another value
  and/or else specify a value for \code{ilambda} and/or \code{itheta}.
  The argument is recycled to length 2, and the first value
  corresponds to \code{theta}, etc.
  

  }
  \item{ishrinkage, zero}{
    See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{glambda, parallel}{
    See \code{\link{CommonVGAMffArguments}} for information.
    Argument \code{glambda} is similar to \code{gsigma}
    there and is currently used only if \code{imethod[2] = 1}.


  }
% \item{zero}{
% An integer vector, containing the value 1 or 2.
% If so, \eqn{\lambda} or \eqn{\theta} respectively
% are modelled as an intercept only.
% If set to \code{NULL} then both linear/additive predictors are
% modelled as functions of the explanatory variables.


% }


}
\details{
  The generalized Poisson distribution (GPD) was proposed by
  Consul and Jain (1973),
  and it has PMF
\deqn{f(y)=\theta(\theta+\lambda y)^{y-1} \exp(-\theta-\lambda y) / y!}{
  f(y)=\theta(\theta+\lambda * y)^(y-1) * exp(-\theta-\lambda * y) / y!}
for \eqn{0 < \theta} and \eqn{y = 0,1,2,\ldots}.
Theoretically,
    \eqn{\max(-1,-\theta/m) \leq \lambda \leq 1}{
          max(-1,-\theta/m) \le   lambda \le  1}
where \eqn{m}  \eqn{(\geq 4)}{(\ge 4)} is the greatest positive
integer satisfying \eqn{\theta + m\lambda > 0}
when \eqn{\lambda < 0}
[and then \eqn{Pr(Y=y) = 0} for \eqn{y > m}].
However, there are problems with a negative \eqn{\lambda} such as
it not being normalized, so this family function restricts
\eqn{\lambda} to \eqn{(0, 1)}.



%Note the complicated support for this distribution means,
%for some data sets,
%the default link for \code{llambda} will not always work, and
%some tinkering may be required to get it running.



This original parameterization is called the GP-0 by \pkg{VGAM},
partly because there are two other common parameterizations
called the GP-1 and GP-2 (see Yang et al. (2009)),
\code{\link{genpoisson1}}
and \code{\link{genpoisson2}})
that are more suitable for regression.
However, \code{genpoisson()} has been simplified to
\code{\link{genpoisson0}} by only handling positive parameters,
hence only overdispersion relative to the Poisson is accommodated.
Some of the reasons for this are described in
Scollnik (1998), e.g., the probabilities do not
sum to unity when \code{lambda} is negative.
To simply things, \pkg{VGAM} 1.1-4 and later will only
handle positive \code{lambda}.




%Note that the GP-1 and GP-2 will handle underdispersion
%relative to the Poisson.




%As Consul and Famoye (2006) state on p.165, the lower limits
%on \eqn{\lambda} and \eqn{m \ge 4}{m >= 4} are imposed
%to ensure that there are at least 5 classes with nonzero
%probability when \eqn{\lambda} is negative.



An ordinary Poisson distribution corresponds
to \eqn{\lambda = 0}{lambda = 0}.
The mean (returned as the fitted values) is
\eqn{E(Y) = \theta / (1 - \lambda)}
and the variance is \eqn{\theta / (1 - \lambda)^3}
so that the variance is proportional to the mean,
just like the NB-1 and quasi-Poisson.





For more information see Consul and Famoye (2006) for a summary
and Consul (1989) for more details.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Consul, P. C. and Jain, G. C. (1973).
A generalization of the Poisson distribution.
\emph{Technometrics},
\bold{15}, 791--799.


Consul, P. C. and Famoye, F. (2006).
\emph{Lagrangian Probability Distributions},
Boston, USA: Birkhauser.


Jorgensen, B. (1997).
\emph{The Theory of Dispersion Models}.
London: Chapman & Hall.


Consul, P. C. (1989).
\emph{Generalized Poisson Distributions: Properties and Applications}.
New York, USA: Marcel Dekker.



Yang, Z., Hardin, J. W., Addy, C. L. (2009).
A score test for overdispersion in Poisson regression based
on the generalized Poisson-2 model.
\emph{J. Statist. Plann. Infer.},
\bold{139}, 1514--1521.



Yee, T. W. (2020).
  On generalized Poisson regression.
  \emph{In preparation}.



}
\section{Warning }{
  Although this family function is far less fragile compared to
  what used to be called \code{genpoisson()} it is still a
  good idea to monitor convergence because
  equidispersion may result in numerical problems;
  try \code{\link{poissonff}} instead.
  And underdispersed data will definitely result in
  numerical problems and warnings;
  try \code{\link[stats]{quasipoisson}} instead.




}

\author{ T. W. Yee.
  Easton Huch derived the EIM and it has been implemented
  in the \code{weights} slot.


}
\note{
  This family function replaces \code{genpoisson()}, and some of the
  major changes are:
  (i) the swapping of the linear predictors;
  (ii) the change from \code{\link{rhobitlink}} to
  \code{\link{logitlink}} in \code{llambda}
  to reflect the no longer handling of underdispersion;
  (iii) proper Fisher scoring is implemented to give improved
  convergence.

  

  Notationally, and in the literature too,
  don't get confused because \code{theta}
  (and not \code{lambda}) here really
  matches more closely with \code{lambda} of
  \code{\link[stats:Poisson]{dpois}}.




  This family function handles multiple responses.
  This distribution is potentially useful for dispersion modelling.
  Convergence and numerical problems may occur when \code{lambda}
  becomes very close to 0 or 1.

  

% If a failure occurs then you might want to try something like
% \code{llambda = extlogitlink(min = -0.9, max = 1)}
% to handle the LHS complicated constraint,
% and if that doesn't work, try
% \code{llambda = extlogitlink(min = -0.8, max = 1)}, etc.


}

\seealso{
  \code{\link{Genpois0}},
  \code{\link{genpoisson1}},
  \code{\link{genpoisson2}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}},
  \code{\link[stats]{Poisson}},
  \code{\link[stats]{quasipoisson}}.



% \code{\link{rhobitlink}},
% \code{\link{rgenpois0}},
% \code{\link{extlogitlink}}.



}
\examples{
gdata <- data.frame(x2 = runif(nn <- 500))
gdata <- transform(gdata, y1 = rgenpois0(nn, theta = exp(2 + x2),
                                         logitlink(1, inverse = TRUE)))
gfit0 <- vglm(y1 ~ x2, genpoisson0, data = gdata, trace = TRUE)
coef(gfit0, matrix = TRUE)
summary(gfit0)
}
\keyword{models}
\keyword{regression}

% yettodo: see csda 2009, 53(9): 3478--3489.

%{% f(y) = theta*(theta+lambda*y)^(y-1) exp(-theta-lambda*y) / y!}
%# gdata <- transform(gdata, y1=rnbinom(nn, exp(1), mu=exp(2 - x2)))




