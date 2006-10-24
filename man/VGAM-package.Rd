\name{VGAM-package}
\alias{VGAM-package}
\alias{VGAM}
\docType{package}
\title{
Vector Generalized Linear and Additive Models
}
\description{
Vector generalized linear and additive models, and associated models
(Reduced-Rank VGLMs, Quadratic RR-VGLMs, Reduced-Rank VGAMs). This package
fits many models and distribution by maximum likelihood estimation (MLE)
or penalized MLE. Also fits constrained ordination models in ecology.

}
\details{
\tabular{ll}{
Package: \tab VGAM\cr
Version: \tab 0.7-1\cr
Date: \tab 2006-10-24\cr
Depends: \tab R (>= 2.4.0), splines, methods, stats, stats4\cr
License: \tab GPL version 2\cr
URL: \tab http://www.stat.auckland.ac.nz/~yee/VGAM\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Packaged: \tab Tue Oct 24 15:36:03 2006; yee\cr
Built: \tab R 2.4.0; i686-pc-linux-gnu; 2006-10-24 15:44:11; unix\cr
}

Index:
\preformatted{
AA.Aa.aa                The AA-Aa-aa Blood Group System
AB.Ab.aB.ab             The AB-Ab-aB-ab Blood Group System
AB.Ab.aB.ab2            The AB-Ab-aB-ab2 Blood Group System
ABO                     The ABO Blood Group System
Benini                  The Benini Distribution
Betabin                 The Beta-Binomial Distribution
Betageom                The Beta-Geometric Distribution
Betanorm                The Beta-Normal Distribution
Bisa                    The Birnbaum-Saunders Distribution
Brat                    Inputting Data to fit a Bradley Terry Model
Coef                    Computes Model Coefficients and Quantities
Coef.qrrvglm            Returns Important Matrices etc. of a QO Object
Coef.qrrvglm-class      Class "Coef.qrrvglm"
Coef.rrvglm             Returns Important Matrices etc. of a RR-VGLM
                        Object
Coef.rrvglm-class       Class "Coef.rrvglm"
Coef.vlm                Extract Model Coefficients for VLM Objects
Dagum                   The Dagum Distribution
Fisk                    The Fisk Distribution
Frank                   Frank's Bivariate Distribution
Frechet                 The Frechet Distribution
G1G2G3                  The G1G2G3 Blood Group System
Hzeta                   Haight's Zeta Function
Inv.gaussian            The Inverse Gaussian Distribution
Invlomax                The Inverse Lomax Distribution
Invparalogistic         The Inverse Paralogistic Distribution
Links                   Link functions for VGLM/VGAM/etc. families
Lino                    The Generalized Beta Distribution (Libby and
                        Novick, 1982)
Log                     Logarithmic Distribution
Lomax                   The Lomax Distribution
MNSs                    The MNSs Blood Group System
Max                     Maxima
Maxwell                 The Maxwell Distribution
Nakagami                Nakagami Distribution
Opt                     Maxima
Paralogistic            The Paralogistic Distribution
Pareto                  The Pareto Distribution
ParetoIV                The Pareto(IV/III/II) Distributions
Polono                  The Poisson Lognormal Distribution
Posbinom                Positive-Binomial Distribution
Posnorm                 The Positive-Normal Distribution
Pospois                 Positive-Poisson Distribution
Rayleigh                The Rayleigh Distribution
Sinmad                  The Singh-Maddala Distribution
Tikuv                   A Short-tailed Symmetric Distribution
Tol                     Tolerances
Tpareto                 The Truncated Pareto Distribution
Zeta                    The Zeta Distribution
Zibinom                 Zero-Inflated Binomial Distribution
Zipf                    The Zipf Distribution
Zipois                  Zero-Inflated Poisson Distribution
acat                    Ordinal Regression with Adjacent Categories
                        Probabilities
auuc                    Auckland University Undergraduate Counts
benini                  Benini Distribution Family Function
betaII                  Beta Distribution of the Second Kind
betabin.ab              Beta-binomial Distribution Family Function
betabinomial            Beta-binomial Distribution Family Function
betaff                  The Two-parameter Beta Distribution Family
                        Function
betageometric           Beta-geometric Distribution Family Function
betaprime               The Beta-Prime Distribution
bilogis4                Bivariate Logistic Distribution
bilogistic4             Bivariate Logistic Distribution Family Function
binom2.or               Bivariate Logistic Regression
binom2.rho              Bivariate Probit Model
binomialff              Binomial Family Function
bisa                    Birnbaum-Saunders Distribution Family Function
bminz                   Body Mass Index of New Zealand Adults
brat                    Bradley Terry Model
bratt                   Bradley Terry Model With Ties
calibrate               Model Calibrations
calibrate.qrrvglm       Calibration for CQO, UQO and CAO models
calibrate.qrrvglm.control
                        Control function for CQO/UQO/CAO calibration
cao                     Fitting Constrained Additive Ordination (CAO)
cao.control             Control Function for RR-VGAMs (CAO)
cauchit                 Cauchit Link Function
cauchy1                 Cauchy Distribution Family Function
ccoef                   Extract Model Constrained/Canonical
                        Coefficients
cdf.lmscreg             Cumulative Distribution Function for LMS
                        Quantile Regression
cexpon                  Censored Exponential Distribution
cgo                     Redirects the user to cqo
cgumbel                 Censored Gumbel Distribution
chest                   Chest Pain in NZ Adults
chisq                   Chi-squared Distribution
clo                     Redirects the user to rrvglm
cloglog                 Complementary Log-log Link Function
cnormal1                Censored Normal Distribution
coalminers              Breathlessness and Wheeze Amongst Coalminers
constraints             Constraint Matrices
cqo                     Fitting Constrained Quadratic Ordination (CQO)
cratio                  Ordinal Regression with Continuation Ratios
cumulative              Ordinal Regression with Cumulative
                        Probabilities
dagum                   Dagum Distribution Family Function
dcnormal1               Univariate Normal Distribution with Double
                        Censoring
deplot.lmscreg          Density Plot for LMS Quantile Regression
dgumbel                 The Gumbel Distribution
dirichlet               Fitting a Dirichlet Distribution
dirmul.old              Fitting a Dirichlet-Multinomial Distribution
dirmultinomial          Fitting a Dirichlet-Multinomial Distribution
dlaplace                The Laplace Distribution
enzyme                  Enzyme data
erf                     Error Function
erlang                  Erlang Distribution
expexp                  Exponentiated Exponential Distribution
expexp1                 Exponentiated Exponential Distribution
exponential             Exponential Distribution
fff                     F Distribution Family Function
fill                    Creates a Matrix of Appropriate Dimension
fisherz                 Fisher's Z Link Function
fisk                    Fisk Distribution family function
fitted.vlm              Fitted Values of a VLM object
frank                   Frank's Bivariate Distribution Family Function
frechet2                Frechet Distribution Family Function
freund61                Freund's (1961) Bivariate Extension of the
                        Exponential Distribution
gamma1                  1-parameter Gamma Distribution
gamma2                  2-parameter Gamma Distribution
gamma2.ab               2-parameter Gamma Distribution
gammahyp                Gamma Hyperbola Bivariate Distribution
garma                   GARMA (Generalized Autoregressive
                        Moving-Average) Models
gaussianff              Gaussian (normal) Family Function
genbetaII               Generalized Beta Distribution of the Second
                        Kind
genpoisson              Generalized Poisson distribution
geometric               Geometric Distribution
get.smart               Retrieve One Component of ".smart.prediction"
get.smart.prediction    Retrieves ".smart.prediction"
gev                     Generalized Extreme Value Distribution Family
                        Function
gevUC                   The Generalized Extreme Value Distribution
gew                     General Electric and Westinghouse Data
ggamma                  Generalized Gamma distribution family function
ggammaUC                The Generalized Gamma Distribution
golf                    Gamma-Ordinal Link Function
gpd                     Generalized Pareto Distribution Family Function
gpdUC                   The Generalized Pareto Distribution
grc                     Fitting Goodman's RC Association Model
gumbel                  Gumbel Distribution Family Function
guplot                  Gumbel Plot
hspider                 Hunting Spider Data
hunua                   Hunua Ranges data
hyper                   Hypergeometric Family Function
hzeta                   Haight's Zeta Family Function
iam                     Index from Array to Matrix
identity                Identity Link Function
inv.gaussianff          Inverse Gaussian Distribution Family Function
invlomax                Inverse Lomax Distribution Family Function
invparalogistic         Inverse Paralogistic Distribution Family
                        Function
is.smart                Test For a Smart Object
leipnik                 Leipnik Distribution Family Function
lerch                   Lerch Phi Function
levy                    Levy Distribution Family Function
lgammaUC                The Log-Gamma Distribution
lgammaff                Log-gamma Distribution Family Function
lino                    Generalized Beta Distribution Family Function
lirat                   Low-iron Rat Teratology Data
lms.bcg                 LMS Quantile Regression with a Box-Cox
                        transformation to a Gamma Distribution
lms.bcn                 LMS Quantile Regression with a Box-Cox
                        Transformation to Normality
lms.yjn                 LMS Quantile Regression with a Yeo-Johnson
                        Transformation to Normality
logc                    Complementary-log Link Function
loge                    Log link function
logff                   Logarithmic Distribution
logistic                Logistic Distribution Family Function
logit                   Logit Link Function
loglinb2                Loglinear Model for Two Binary Responses
loglinb3                Loglinear Model for Three Binary Responses
loglog                  Log-log Link Function
lognormal               Lognormal Distribution
logoff                  Log link function with an offset
lomax                   Lomax Distribution Family Function
lv                      Latent Variables
lvplot                  Latent Variable Plot
lvplot.qrrvglm          Latent Variable Plot for QO models
lvplot.rrvglm           Latent Variable Plot for RR-VGLMs
maxwell                 Maxwell Distribution Family Function
mccullagh89             McCullagh (1989) Distribution Family Function
mckaygamma2             McKay's Bivariate Gamma Distribution
meplot                  Mean Excess Plot
micmen                  Michaelis-Menten Model
mix2normal1             Mixture of Two Univariate Normal Distributions
mix2poisson             Mixture of Two Poisson Distributions
model.framevlm          Construct the Model Frame of a VLM Object
model.matrixvlm         Construct the Design Matrix of a VLM Object
multinomial             Multinomial Logit Model
nakagami                Nakagami Distribution Family Function
nbolf                   Negative Binomial-Ordinal Link Function
negbinomial             Negative Binomial Distribution Family Function
normal1                 Univariate normal distribution
notdocumentedyet        Undocumented and Internally Used Functions and
                        Classes
nzc                     Chinese Population in New Zealand 1867-2001
oxtemp                  Oxford Temperature Data
paralogistic            Paralogistic Distribution Family Function
pareto1                 Pareto and Truncated Pareto Distribution Family
                        Functions
paretoIV                Pareto(IV/III/II) Distribution Family Functions
persp.qrrvglm           Perspective plot for QRR-VGLMs
plotdeplot.lmscreg      Density Plot for LMS Quantile Regression
plotqrrvglm             Model Diagnostic Plots for QRR-VGLMs
plotqtplot.lmscreg      Quantile Plot for LMS Quantile Regression
plotvgam                Default VGAM Plotting
plotvgam.control        Control Function for plotvgam()
pneumo                  Pneumoconiosis amongst a group of coalminers
poissonff               Poisson Family Function
polf                    Poisson-Ordinal Link Function
posbinomial             Positive Binomial Distribution Family Function
posnegbinomial          Positive Negative Binomial Distribution Family
                        Function
posnormal1              Positive Normal Distribution Family Function
pospoisson              Positive Poisson Distribution Family Function
predict.vglm            Predict Method for a VGLM fit
prentice74              Prentice (1974) Log-gamma Distribution
probit                  Probit Link Function
put.smart               Adds a List to the End of the List
                        ".smart.prediction"
qrrvglm.control         Control function for QRR-VGLMs (CQO)
qtplot.gumbel           Quantile Plot for Gumbel Regression
qtplot.lmscreg          Quantile Plot for LMS Quantile Regression
quasibinomialff         Quasi-Binomial Family Function
quasipoissonff          Quasi-Poisson Family Function
rayleigh                Rayleigh Distribution Family Function
rcqo                    Constrained Quadratic Ordination
rdiric                  The Dirichlet distribution
recexp1                 Upper Record Values from a 1-parameter
                        Exponential Distribution
reciprocal              Reciprocal link function
recnormal1              Upper Record Values from a Univariate Normal
                        Distribution
rhobit                  Rhobit Link Function
rig                     Reciprocal Inverse Gaussian distribution
rlplot.egev             Return Level Plot for GEV Fits
rposnegbin              Positive-negative binomial distribution random
                        variates
rrar                    Nested reduced-rank autoregressive models for
                        multiple time series
rrvglm                  Fitting Reduced-Rank Vector Generalized Linear
                        Models (RR-VGLMs)
rrvglm-class            Class "rrvglm"
rrvglm.control          Control function for rrvglm
rrvglm.optim.control    Control function for rrvglm() calling optim()
s                       Defining smooths in VGAM formulae
setup.smart             Smart Prediction Setup
simplex                 Simplex distribution
sinmad                  Singh-Maddala Distribution Family Function
skewnormal1             Univariate Skew-Normal Distribution Family
                        Function
smart.expression        S Expression for Smart Functions
smart.mode.is           Determine What Mode the Smart Prediction is In
smartpred               Smart Prediction
snorm                   Skew-Normal Distribution
sratio                  Ordinal Regression with Stopping Ratios
studentt                Student t Distribution
tikuv                   Short-tailed Symmetric Distribution Family
                        Function
tobit                   Tobit Model
trplot                  Trajectory Plot
trplot.qrrvglm          Trajectory plot for QRR-VGLMs
uqo                     Fitting Unconstrained Quadratic Ordination
                        (UQO)
uqo.control             Control Function for UQO models
usagrain                USA grain prices
venice                  Venice Maximum Sea Levels
vgam                    Fitting Vector Generalized Additive Models
vgam-class              Class "vgam"
vgam.control            Control function for vgam
vglm                    Fitting Vector Generalized Linear Models
vglm-class              Class "vglm"
vglm.control            Control function for vglm
vglmff-class            Class "vglmff"
vonmises                von Mises Distribution Family Function
vsmooth.spline          Vector cubic smoothing spline
waitakere               Waitakere Ranges data
wald                    Wald Distribution Family Function
weibull                 Weibull Distribution Family Function
weightsvglm             Prior and Working Weights of a VGLM fit
wrapup.smart            Cleans Up After Smart Prediction
yeo.johnson             Yeo-Johnson Transformation
yip88                   Zero-Inflated Poisson Distribution (Yip (1988)
                        algorithm)
zanegbinomial           Zero-Altered Negative Binomial Distribution
zapoisson               Zero-Altered Poisson Distribution
zero                    The zero Argument in VGAM Family Functions
zeta                    Riemann's Zeta Function
zetaff                  Zeta Distribution Family Function
zibinomial              Zero-Inflated Binomial Distribution Family
                        Function
zipf                    Zipf Distribution Family Function
zipoisson               Zero-Inflated Poisson Distribution Family
                        Function
}

%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~


}
\author{
Thomas W. Yee <t.yee@auckland.ac.nz>

Maintainer: Thomas Yee <t.yee@auckland.ac.nz>
}
\references{

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.

Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.

}

\keyword{ package }
\keyword{models}
\keyword{regression}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# Fit the proportional odds model, p.179, in McCullagh and Nelder (1989)
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
(fit = vglm(cbind(normal, mild, severe) ~ let,
            cumulative(parallel=TRUE, reverse=TRUE), pneumo))
fit@y   # Sample proportions
weights(fit, type="prior")   # Number of observations
coef(fit, matrix=TRUE)
constraints(fit)   # Constraint matrices


# Fit a two species GAM simultaneously
data(hunua)
fit2 = vgam(cbind(agaaus, kniexc) ~ s(altitude),
            fam = binomialff(mv=TRUE), hunua)
coef(fit2, mat=TRUE)   # Not really interpretable
\dontrun{
plot(fit2, se=TRUE, overlay=TRUE, lcol=1:2, scol=1:2)
attach(hunua)
o = order(altitude)
matplot(altitude[o], fitted(fit2)[o,], type="l", lwd=2, las=1,
    xlab="Altitude (m)", ylab="Probability of presence",
    main="Two plant species' response curves", ylim=c(0,.8))
rug(altitude)
detach(hunua)
}



# LMS quantile regression
data(bminz)
fit = vgam(BMI ~ s(age, df=c(4,2)), fam=lms.bcn(zero=1), data=bminz, tr=TRUE)
predict(fit)[1:3,]
fitted(fit)[1:3,]
bminz[1:3,]
# Person 1 is near the lower quartile of BMI amongst people his age
cdf(fit)[1:3]

\dontrun{
# Quantile plot
par(bty="l", mar=c(5,4,4,3)+0.1, xpd=TRUE)
qtplot(fit, percentiles=c(5,50,90,99), main="Quantiles",
       xlim=c(15,90), las=1, ylab="BMI", lwd=2, lcol=4)

# Density plot
ygrid = seq(15, 43, len=100)  # BMI ranges
par(mfrow=c(1,1), lwd=2)
a = deplot(fit, x0=20, y=ygrid, xlab="BMI", col="black",
    main="Density functions at Age = 20 (black), 42 (red) and 55 (blue)")
a
a = deplot(fit, x0=42, y=ygrid, add=TRUE, llty=2, col="red")
a = deplot(fit, x0=55, y=ygrid, add=TRUE, llty=4, col="blue", Attach=TRUE)
a@post$deplot  # Contains density function values
}


# GEV distribution for extremes
data(oxtemp)
(fit = vglm(maxtemp ~ 1, egev, data=oxtemp, trace=TRUE))
fitted(fit)[1:3,]
coef(fit, mat=TRUE)
Coef(fit)
vcov(fit)
vcov(fit, untransform=TRUE)
sqrt(diag(vcov(fit)))   # Approximate standard errors
\dontrun{ rlplot(fit) }
}


