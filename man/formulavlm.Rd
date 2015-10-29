\name{formulavlm}
%\name{confint}
\alias{formula.vlm}
\alias{formulavlm}
\alias{term.names}
\alias{term.namesvlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model Formulae and Term Names for VGLMs }
\description{
     The methods function for \code{formula} to
     extract the formula from a fitted object,
     as well as a methods function to return the names
     of the terms in the formula.

}
\usage{
\method{formula}{vlm}(x, \dots)
formulavlm(x, form.number = 1, \dots)
term.names(model, \dots)
term.namesvlm(model, form.number = 1, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, model}{ A fitted model object.


  }
 \item{form.number}{Formula number, is 1 or 2.
       which correspond to the arguments \code{formula}
       and \code{form2} respectively.


  }
 \item{\dots}{Same as \code{\link[stats]{formula}}.
 }
}
\details{
  The \code{formula} methods function is
  based on \code{\link[stats]{formula}}.


}
\value{
  The \code{formula} methods function should return something similar to
  \code{\link[stats]{formula}}.
  The \code{term.names} methods function should return a character string
  with the terms in the formula; this includes any intercept (which
  is denoted by \code{"(Intercept)"} as the first element.)



}
%\references{
%}
\author{ Thomas W. Yee }



%\note{
%}
%\section{Warning }{
%}



\seealso{
  \code{\link{has.interceptvlm}}.
% \code{termsvlm}.


}

\examples{
# Example: this is based on a glm example
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
vglm.D93 <- vglm(counts ~ outcome + treatment, family = poissonff)
formula(vglm.D93)
pdata <- data.frame(counts, outcome, treatment)  # Better style
vglm.D93 <- vglm(counts ~ outcome + treatment, poissonff, data = pdata)
formula(vglm.D93)
term.names(vglm.D93)
responseName(vglm.D93)
has.intercept(vglm.D93)
}
\keyword{models}
\keyword{regression}

