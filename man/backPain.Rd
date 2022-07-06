\name{backPain}
\alias{backPain}
\alias{backPain2}
\docType{data}
\title{ Data on Back Pain Prognosis, from Anderson (1984) }
\description{
  Data from a study of patients suffering from back
  pain. Prognostic variables were recorded at presentation and
  progress was categorised three weeks after treatment.

}
\usage{data(backPain)}
\format{
  A data frame with 101 observations on the following 4 variables.
  \describe{
    \item{x2}{length of previous attack.}
    \item{x3}{pain change.}
    \item{x4}{lordosis.}
    \item{pain}{an ordered factor describing the progress of each
      patient with levels \code{worse} < \code{same} <
      \code{slight.improvement} < \code{moderate.improvement} <
      \code{marked.improvement} < \code{complete.relief}. }
  }
}
\source{
  \code{http://ideas.repec.org/c/boc/bocode/s419001.html}


% \url{http://ideas.repec.org/c/boc/bocode/s419001.html}


  The data set and this help file was copied from \pkg{gnm}
  so that a vignette in \pkg{VGAM} could be run; the analysis is
  described in Yee (2010).
  

  The data frame \code{backPain2} is a modification of
  \code{backPain} where the variables have been renamed
  (\code{x1} becomes \code{x2},
   \code{x2} becomes \code{x3},
   \code{x3} becomes \code{x4})
  and
  converted into factors.

  


}
\references{
  Anderson, J. A. (1984).
  Regression and Ordered Categorical Variables.
  \emph{J. R. Statist. Soc. B}, \bold{46(1)}, 1-30.


Yee, T. W. (2010).
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\doi{10.18637/jss.v032.i10}.



% \url{https://www.jstatsoft.org/article/view/v032i10/}.
% \url{https://www.jstatsoft.org/v32/i10/}.



}
\examples{
summary(backPain)
summary(backPain2)
}
\keyword{datasets}
% set.seed(1)
% data(backPain)




