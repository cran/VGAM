\name{dexpbinomial}
\alias{dexpbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Double Exponential Binomial Distribution Family Function }
\description{
  Fits a double exponential binomial distribution by
  maximum likelihood estimation.
  The two parameters here are the mean and dispersion parameter.

}
\usage{
dexpbinomial(lmean="logit", ldispersion="logit",
             emean=list(),  edispersion=list(),
             idispersion=0.25, zero=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean, ldispersion}{ 
  Link functions applied to the two parameters, called
  \eqn{\mu}{mu} and \eqn{\theta}{theta} respectively below.
  See \code{\link{Links}} for more choices.
  The defaults cause the parameters to be restricted to \eqn{(0,1)}. 

  }
  \item{emean, edispersion}{ 
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{idispersion}{ 
  Initial value for the dispersion parameter.
  If given, it must be in range, and is recyled to the necessary length.
  Use this argument if convergence failure occurs.

  }
  \item{zero}{ 
  An integer specifying which
  linear/additive predictor is to be modelled as an intercept only.
  If assigned, the single value should be either \code{1} or \code{2}.
  The default is to have a single dispersion parameter value.
  To model both parameters as functions of the covariates assign
  \code{zero=NULL}.

  }
}
\details{
  This distribution provides a way for handling overdispersion in
  a binary response.
  The double exponential binomial distribution belongs the family of
  double exponential distributions proposed by Efron (1986).
  Below, equation numbers refer to that original article.
  Briefly, the idea is that an ordinary one-parameter exponential family
  allows the addition of a second parameter \eqn{\theta}{theta}
  which varies the dispersion
  of the family without changing the mean.
  The extended family behaves like the original family with
  sample size changed from \eqn{n} to \eqn{n\theta}{n*theta}.
  The extended family is an exponential family in \eqn{\mu}{mu} when
  \eqn{n} and \eqn{\theta}{theta} are fixed, and
  an exponential family in \eqn{\theta}{theta} when
  \eqn{n} and \eqn{\mu}{mu} are fixed.
  Having \eqn{0 < \theta < 1}{0 < theta < 1} corresponds to overdispersion
  with respect to the binomial distribution.
  See Efron (1986) for full details.

  This \pkg{VGAM} family function implements an \emph{approximation}
  (2.10) to the exact density (2.4). It replaces the normalizing
  constant by unity since the true value nearly equals 1.
  The default model fitted is \eqn{\eta_1 = logit(\mu)}{eta1 =logit(mu)}
  and \eqn{\eta_2 = logit(\theta)}{eta2 = logit(theta)}.
  This restricts both parameters to lie between 0 and 1, although
  the dispersion parameter can be modelled over a larger parameter space by
  assigning the arguments \code{ldispersion} and \code{edispersion}.

  Approximately, the mean (of \eqn{Y}) is \eqn{\mu}{mu}.
  The \emph{effective sample size} is the dispersion parameter multiplied
  by the original sample size,
  i.e., \eqn{n\theta}{n*theta}.
  This family function uses Fisher scoring, and the two estimates are
  asymptotically independent because the expected information matrix
  is diagonal.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}.

}
\references{

  Efron, B. (1986)
  Double exponential families and their use in generalized linear regression.
  \emph{Journal of the American Statistical Association},
  \bold{81}, 709--721.

}

\author{ T. W. Yee }
\note{
  This function processes the input in the same way
  as \code{\link{binomialff}}, however multivariate responses are
  not allowed (\code{binomialff(mv=FALSE)}).

}
\section{Warning }{
  Numerical difficulties can occur; if so, try using \code{idispersion}.

}
\seealso{
  \code{\link{binomialff}},
  \code{\link{toxop}}.

}
\examples{
# This example mimics the example in Efron (1986). The results here
# differ slightly.

# Scale the variables
toxop = transform(toxop,
                  phat = positive / ssize,
                  srainfall = scale(rainfall),  # (6.1)
                  sN = scale(ssize))            # (6.2)

# A fit similar (should be identical) to Section 6 of Efron (1986).
# But does not use poly(), and M=1.25 here, as in (5.3)
cmlist = list("(Intercept)"=diag(2),
              "I(srainfall)"=rbind(1,0),
              "I(srainfall^2)"=rbind(1,0),
              "I(srainfall^3)"=rbind(1,0),
              "I(sN)"=rbind(0,1),
              "I(sN^2)"=rbind(0,1))
dlist = list(min=0, max=1.25)
fit = vglm(phat ~ I(srainfall) + I(srainfall^2) + I(srainfall^3) +
                  I(sN) + I(sN^2),
           fam = dexpbinomial(ldisp="elogit", idisp=0.2,
                              edisp=dlist, zero=NULL),
           data=toxop, weight=ssize, trace=TRUE, constraints=cmlist)

# Now look at the results
coef(fit)
coef(fit, matrix=TRUE)
head(fitted(fit))
summary(fit)
vcov(fit)
sqrt(diag(vcov(fit)))   # Standard errors

# Effective sample size (not quite the last column of Table 1)
head(predict(fit))
Dispersion = elogit(predict(fit)[,2], earg=dlist, inverse=TRUE)
c(round(weights(fit, type="prior") * Dispersion, dig=1))


# Ordinary logistic regression (gives same results as (6.5))
ofit = vglm(phat ~ I(srainfall) + I(srainfall^2) + I(srainfall^3),
            fam = binomialff, data=toxop, weight=ssize, trace=TRUE)


# Same as fit but it uses poly(), and can be plotted (cf. Figure 1)
cmlist2 = list("(Intercept)"=diag(2),
               "poly(srainfall, 3)"=rbind(1,0),
               "poly(sN, 2)"=rbind(0,1))
fit2 = vglm(phat ~ poly(srainfall, 3) + poly(sN, 2),
            fam = dexpbinomial(ldisp="elogit", idisp=0.2,
                               edisp=dlist, zero=NULL),
            data=toxop, weight=ssize, trace=TRUE, constraints=cmlist2)
\dontrun{
par(mfrow=c(1,2))
plotvgam(fit2, se=TRUE, lcol="blue", scol="red")  # Cf. Figure 1


# Cf. Figure 1(a)
par(mfrow=c(1,2))
o = with(toxop, sort.list(rainfall))
with(toxop, plot(rainfall[o], fitted(fit2)[o], type="l", col="blue",
                 las=1, ylim=c(0.3, 0.65)))
with(toxop, points(rainfall[o], fitted(ofit)[o], col="red", type="b",
                   pch=19))

# Cf. Figure 1(b)
o = with(toxop, sort.list(ssize))
with(toxop, plot(ssize[o], Dispersion[o], type="l", col="blue", las=1,
                 xlim=c(0, 100)))
}
}
\keyword{models}
\keyword{regression}
