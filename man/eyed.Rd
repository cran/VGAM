\name{eyed}
\alias{eyed}
\docType{data}
\title{
  Eye Disease Data:
  Diabetic Retinopathy}
\description{
  From Agresti's 2010 ordinal book, an
  eye disease data set involving ordinal
  measurements of retinopathy in left and
  right eyes of insulin-taking young diabetics.



  %%  ~~ A concise (1-5 lines) description
  %% of the dataset. ~~
}
\usage{data("eyed")}
\format{
  A data frame with 720 observations on the
  following 23 variables.
  \describe{
    \item{\code{rrl}}{a numeric vector}
    \item{\code{lrl}}{a numeric vector}
    \item{\code{rl}}{a numeric vector}
    \item{\code{rme}}{Right eye macular oedema
     (absent = 0, present = 1) }
    \item{\code{lme}}{Left eye macular oedema
     (absent = 0, present = 1) }
    \item{\code{rre}}{Right eye refraction index}
    \item{\code{lre}}{Left eye refraction index}
    \item{\code{riop}}{Right eye intraocular pressure}
    \item{\code{liop}}{Left eye intraocular pressure}
    \item{\code{age}}{Age at diagnosis
      (values are less than 30 years)  }
    \item{\code{diab}}{Duration of diabetes (years)}
    \item{\code{gh}}{Glycosylated haemoglobin level}
    \item{\code{sbp}}{Systolic blood pressure (mmHg)}
    \item{\code{dbp}}{Diastolic blood pressure (mmHg)}
    \item{\code{bmi}}{Body mass index, in kg / m2,
    rounded to an integer }
    \item{\code{pr}}{a numeric vector}
    \item{\code{sex}}{male = 0, female = 1}
    \item{\code{prot}}{Proteinuria
      (absent = 0, present = 1)}
    \item{\code{dose}}{a numeric vector}
    \item{\code{resid}}{a numeric vector}
    \item{\code{rerl}}{Right eye severity of retinopathy,
      an ordered factor with levels
      \code{None}, \code{Mild}, \code{Moderate},
      \code{Proliferative}}
    \item{\code{lerl}}{Same as \code{rerl} but for
    left eye}
    \item{\code{rl2}}{a numeric vector}
  }
}
\details{
  The data is used as the example in
  Section 9.1.3 of Agresti (2010),
  entitled
  \emph{Eye Disease Risk Factors}.
  Williamson and Kim (1996) describe a
  Wisconsin epidemiological study
  in which insulin-taking young diabetics
  were examined for diabetic retinopathy
  (both prevalence and severity)
  in both eyes.
  One could fit a proportional odds model
  to each eye.
  The retinopathy severity is ordinal:
  \code{None}, \code{Mild}, \code{Moderate},
  \code{Proliferative}.
  Variables \code{lerl} and \code{rerl} can
  be considered the primary responses.



 }
\source{
  \url{https://users.stat.ufl.edu/~aa/ordinal/data.html}
  retrieved in Nov 2025.


  
  %%  ~~ reference to a publication or URL
  %%from which the data were obtained ~~
}
\references{


Agresti, A. (2010).
\emph{Analysis of Ordered Categorical Data},
Second Edition. Wiley: Hoboken, NJ, USA. 


Williamson, J. and Kim, K. (1996).
A global odds ratio regression model for
bivariate ordered categorical data
from opthalmologic studies.
\emph{Statistics in Medicine} \bold{15}, 1507--1518.



}
\examples{
summary(eyed)
}
\keyword{datasets}
