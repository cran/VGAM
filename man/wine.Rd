\name{wine}
\alias{wine}
\docType{data}
\title{ Bitterness in Wine Data
%%   ~~ data name/kind ... ~~

}
\description{
 This oenological data frame concerns the amount of bitterness
 in 78 bottles of white wine.


}
\usage{
data(wine)
}

\format{
  A data frame with 4 rows on the following 7 variables.

\describe{
  \item{temp}{temperature, with levels cold and warm.


  }
  \item{contact}{whether contact of the juice with the
    skin was allowed or avoided, for a specified period.
    Two levels: no or yes.

  }
  \item{bitter1, bitter2, bitter3, bitter4, bitter5}{
  numeric vectors, the counts.
  The order is none to most intense.


  }
}

}
\details{

The data set comes from Randall (1989) and concerns a factorial
experiment for investigating factors that affect the bitterness of white
wines. There are two factors in the experiment: temperature at the time
of crushing the grapes and contact of the juice with the skin.
Two bottles of wine were fermented for each of the treatment combinations.
A panel of 9 judges were selected and trained for the ability to
detect bitterness.
Thus there were 72 bottles in total.
Originally, the bitterness of the wine were taken on a continuous
scale in the interval from 0 (none) to 100 (intense)
but later they were
grouped using equal lengths into five ordered categories 1, 2, 3, 4 and 5.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{

% Further information is at:


% September 30, 2013
Christensen, R. H. B. (2013)
Analysis of ordinal data with cumulative link models---estimation
with the R-package \pkg{ordinal}.
R Package Version 2013.9-30.
\url{https://CRAN.R-project.org/package=ordinal}.
%\url{https://www.R-project.org/package=ordinal}.
%\url{https://www.CRAN.R-project.org/package=ordinal}.  % Prior to 20150728




  Randall, J. H. (1989)
  The analysis of sensory data by generalized linear model.
  \emph{Biometrical Journal} \bold{31}(7), 781--793.


  Kosmidis, I. (2014)
  Improved estimation in cumulative link models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{76}, in press.


}
\examples{
wine
summary(wine)
}
\keyword{datasets}



