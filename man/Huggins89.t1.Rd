\name{Huggins89.t1}
\alias{Huggins89.t1}
\docType{data}
\title{
  Table 1 of Huggins (1989)
}
\description{
  Simulated capture data set for the linear logistic model
  depending on an occasion covariate and an individual
  covariate for 10 trapping occasions and 20 individuals.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Huggins89.t1)}
\format{
  The format is a data frame.


%chr "Huggins89.t1"


}
\details{
  Table 1 of Huggins (1989) gives this toy data set.
  Note that variables \code{z1},\ldots,\code{z10} are
  occasion-specific variables. They correspond to the
  response variables \code{y1},\ldots,\code{y10} which
  have values 1 for capture and 0 for not captured.


%%  ~~ If necessary, more details than the __description__ above ~~
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{

Huggins, R. M. (1989)
On the statistical analysis of capture experiments.
\emph{Biometrika},
\bold{76}, 133--140.


%%  ~~ possibly secondary sources and usages ~~
}
\examples{
\dontrun{
small.Huggins89.t1 <- transform(Huggins89.t1, Zedd = z1, Z2 = z2, Z3 = z3)
small.Huggins89.t1 <- subset(small.Huggins89.t1, y1 + y2 + y3 > 0)
# fit1 is the bottom equation on p.133, but this is only for the 1st 3 responses.
# Currently posbernoulli.tb() cannot handle more than 3 Bernoulli variates.
# The fit is not very good.
fit1 <-
  vglm(cbind(y1, y2, y3) ~  x2 + Zedd,
       xij = list(Zedd ~ z1 + z2 + z3 + Z2 + Z3 - 1),
       posbernoulli.tb(parallel.t = TRUE), maxit = 155,
       data = small.Huggins89.t1, trace = TRUE,
       form2 = ~ x2 + Zedd + z1 + z2 + z3 + Z2 + Z3)
coef(fit1)
coef(fit1, matrix = TRUE)  # M_t model
constraints(fit1)
summary(fit1)
fit1@extra$N.hat     # Estimate of the population size N
fit1@extra$SE.N.hat  # Its standard error


fit.t <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
              posbernoulli.t, data = Huggins89.t1, trace = TRUE)
coef(fit.t)
coef(fit.t, matrix = TRUE)  # M_t model
summary(fit.t)
fit.t@extra$N.hat     # Estimate of the population size N
fit.t@extra$SE.N.hat  # Its standard error


fit.b <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
              posbernoulli.b, data = Huggins89.t1, trace = TRUE)
coef(fit.b)
coef(fit.b, matrix = TRUE)  # M_b model
summary(fit.b)
fit.b@extra$N.hat
fit.b@extra$SE.N.hat


fit.0 <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
             posbernoulli.b(parallel.b = TRUE), data = Huggins89.t1,
             trace = TRUE)
coef(fit.0, matrix = TRUE)  # M_0 model (version 1)
coef(fit.0)
summary(fit.0)
fit.0@extra$N.hat
fit.0@extra$SE.N.hat


Fit.0 <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
              posbernoulli.t(parallel.t = TRUE), data = Huggins89.t1,
              trace = TRUE)
coef(Fit.0)
coef(Fit.0, matrix = TRUE)  # M_0 model (version 2)
summary(Fit.0)
Fit.0@extra$N.hat
Fit.0@extra$SE.N.hat
}
}
\keyword{datasets}
% data(Huggins89.t1)
%## maybe str(Huggins89.t1) ; plot(Huggins89.t1) ...


