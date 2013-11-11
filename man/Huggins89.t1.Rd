\name{Huggins89.t1}
\alias{Huggins89.t1}
\alias{Huggins89table1}
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
\usage{
data(Huggins89table1)
data(Huggins89.t1)
}
\format{
  The format is a data frame.


%chr "Huggins89.t1"


}
\details{
  Table 1 of Huggins (1989) gives this toy data set.
  Note that variables \code{t1},\ldots,\code{t10} are
  occasion-specific variables. They correspond to the
  response variables \code{y1},\ldots,\code{y10} which
  have values 1 for capture and 0 for not captured.



  Both \code{Huggins89table1} and \code{Huggins89.t1} are identical.
  The latter used variables beginning with \code{z}, 
  not \code{t}, and may be withdrawn very soon.



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
Huggins89table1 <- transform(Huggins89table1, x3.tij = t1,
                             T2 = t2, T3 = t3, T4 = t4,  T5 =  t5, T6 = t6,
                             T7 = t7, T8 = t8, T9 = t9, T10 = t10)
small.table1 <- subset(Huggins89table1,
                       y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 > 0)
# fit.tbh is the bottom equation on p.133.
# It is a M_tbh model.
fit.tbh <-
  vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~  x2 + x3.tij,
       xij = list(x3.tij ~ t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 +
                           T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9 + T10 - 1),
       posbernoulli.tb(parallel.t = TRUE ~ x2 + x3.tij),
       data = small.table1, trace = TRUE,
       form2 = ~  x2 + x3.tij +
                  t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 +
                       T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9 + T10)

# These results differ a bit from Huggins (1989), probably because
# two animals had to be removed here (they were never caught):
coef(fit.tbh)  # First element is the behavioural effect
sqrt(diag(vcov(fit.tbh)))  # SEs
constraints(fit.tbh, matrix = TRUE)
summary(fit.tbh, presid = FALSE)
fit.tbh@extra$N.hat     # Estimate of the population site N; cf. 20.86
fit.tbh@extra$SE.N.hat  # Its standard error; cf. 1.87 or 4.51

fit.th <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
               posbernoulli.t, data = small.table1, trace = TRUE)
coef(fit.th)
constraints(fit.th)
coef(fit.th, matrix = TRUE)  # M_th model
summary(fit.th, presid = FALSE)
fit.th@extra$N.hat     # Estimate of the population size N
fit.th@extra$SE.N.hat  # Its standard error

fit.bh <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
               posbernoulli.b(I2 = FALSE), data = small.table1, trace = TRUE)
coef(fit.bh)
constraints(fit.bh)
coef(fit.bh, matrix = TRUE)  # M_bh model
summary(fit.bh, presid = FALSE)
fit.bh@extra$N.hat
fit.bh@extra$SE.N.hat

fit.h <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
             posbernoulli.b, data = small.table1, trace = TRUE)
coef(fit.h, matrix = TRUE)  # M_h model (version 1)
coef(fit.h)
summary(fit.h, presid = FALSE)
fit.h@extra$N.hat
fit.h@extra$SE.N.hat

Fit.h <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ x2,
              posbernoulli.t(parallel.t = TRUE ~ x2),
              data = small.table1, trace = TRUE)
coef(Fit.h)
coef(Fit.h, matrix = TRUE)  # M_h model (version 2)
summary(Fit.h, presid = FALSE)
Fit.h@extra$N.hat
Fit.h@extra$SE.N.hat
}
\keyword{datasets}


%\dontrun{
%}
% data(Huggins89table1)
%## maybe str(Huggins89table1) ; plot(Huggins89table1) ...


% coef(fit1, matrix = TRUE)  # M_t model
% Huggins89.t1 <- transform(Huggins89.t1, xx2 = c(matrix(x2, 2, 10, byrow = TRUE)))




%This code below is equivalent to the above fit.tbh (same name).
%But this version uses manual construction of the constraint matrices:
%tau <- 10
%Hlist <-list("(Intercept)" = cbind(bhvr.effect = c(rep(0, len = tau),
%                                                   rep(1, len = tau-1)),
%                                   overall.intercept = 1),
%             x2            = cbind(rep(1, len = 2*tau-1)),
%             Zedd          = cbind(rep(1, len = 2*tau-1)))
%fit.tbh <-
%  vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~  x2 + Zedd,
%       xij = list(Zedd ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 +
%                              Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10 - 1),
%       posbernoulli.tb, data = small.t1, trace = TRUE,
%       constraints = Hlist,
%       form2 = ~  x2 + Zedd +
%                 z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 +
%                      Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10)


