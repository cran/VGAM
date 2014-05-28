\name{prinia}
\alias{prinia}
\docType{data}
\title{Yellow-bellied Prinia
%%   ~~ data name/kind ... ~~

}
\description{
 A data frame with yellow-bellied Prinia.

}
\usage{
data(prinia)
}

\format{
  A data frame with 151 observations on the following 23 variables.

\describe{
  \item{length}{a numeric vector,
  the scaled wing length (zero mean and unit variance).


  }
  \item{fat}{a numeric vector, fat index;
  originally 1 (no fat) to 4 (very fat) but
  converted to 0 (no fat) versus 1 otherwise.


  }
  \item{cap}{a numeric vector,
  number of times the bird was captured or recaptured.


  }
  \item{noncap}{a numeric vector,
  number of times the bird was not captured.


  }
  \item{y01, y02, y03, y04, y05, y06}{
  a numeric vector of 0s and 1s; for noncapture and capture resp.


  }
  \item{y07, y08, y09, y10, y11, y12}{
  same as above.


  }
  \item{y13, y14, y15, y16, y17, y18, y19}{
  same as above.


  }
}

}
\details{
  The yellow-bellied Prinia \emph{Prinia flaviventris}
  is a common bird species located in Southeast Asia. A
  capture--recapture experiment was conducted at the Mai Po
  Nature Reserve in Hong Kong during 1991, where captured
  individuals had their wing lengths measured and fat index
  recorded. A total of 19 weekly capture occasions were
  considered, where 151 distinct birds were captured.


  More generally, the prinias are a genus of small
  insectivorous birds, and are sometimes referred to as
  \emph{wren-warblers}.  They are a little-known group of the
  tropical and subtropical Old World, the roughly 30 species
  being divided fairly equally between Africa and Asia.



% 20131030; this is old:
% The example below illustrates the necessity of creating
% variables \code{y1}, \code{y2}, \ldots in order for
% \code{\link{posbernoulli.b}},
% \code{\link{posbernoulli.t}} and
% \code{\link{posbernoulli.tb}} to work.
% In contrast, \code{\link{posbinomial}} may have a simple 2-column
% matrix as the response.


% \emph{Prinia inornate} is from the SS paper, not exactly this bird.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{

  Thanks to Paul Yip for permission to make this data available.


% Further information is at:

% Huggins, R. M. and Yip, P. S. F. (1997).
% Statistical analysis of removal experiments with the use of auxillary variables.
% \emph{Statistica Sinica} \bold{7}, 705--712.


  Hwang, W.-H. and Huggins, R. M. (2007)
  Application of semiparametric regression models in the
  analysis of capture--recapture experiments.
  \emph{Australian and New Zealand Journal of Statistics} \bold{49}, 191--202.


}
\examples{
head(prinia)
summary(prinia)
rowSums(prinia[, c("cap", "noncap")])  # 19s

#  Fit a positive-binomial distribution (M.h) to the data:
fit1 <- vglm(cbind(cap, noncap) ~ length + fat, posbinomial, data = prinia)

#  Fit another positive-binomial distribution (M.h) to the data:
#  The response input is suitable for posbernoulli.*-type functions.
fit2 <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08, y09, y10,
                   y11, y12, y13, y14, y15, y16, y17, y18, y19) ~
             length + fat, posbernoulli.b(drop.b = FALSE ~ 0), data = prinia)
}
\keyword{datasets}



