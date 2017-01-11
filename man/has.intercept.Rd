\name{has.interceptvlm}
%\name{confint}
\alias{has.intercept}
%\alias{has.intercept.vlm}
\alias{has.interceptvlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Has a Fitted VGLM Got an Intercept Term? }
\description{
     Looks at the \code{formula} to
     see if it has an intercept term.

}
\usage{
has.intercept(object, \dots)
has.interceptvlm(object, form.number = 1, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A fitted model object.


  }
 \item{form.number}{Formula number, is 1 or 2.
       which correspond to the arguments \code{formula}
       and \code{form2} respectively.


  }
 \item{\dots}{Arguments that are might be passed from
  one function to another.

 }
}

\details{
  This methods function is a simple way to determine whether a
  fitted \code{\link{vglm}} object etc. has an intercept term
  or not.
  It is not entirely foolproof because one might suppress the
  intercept from the formula and then add in a variable in the
  formula that has a constant value.


}
\value{
  Returns a single logical.



}
%\references{
%}
\author{ Thomas W. Yee }

%\note{
%}
%\section{Warning }{
%}

\seealso{
  \code{\link{formulavlm}},
  \code{termsvlm}.


}

\examples{
# Example: this is based on a glm example
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
pdata <- data.frame(counts, outcome, treatment)  # Better style
vglm.D93 <- vglm(counts ~ outcome + treatment, poissonff, data = pdata)
formula(vglm.D93)
term.names(vglm.D93)
responseName(vglm.D93)
has.intercept(vglm.D93)
}
\keyword{models}
\keyword{regression}



% \method{has.intercept}{vlm}(object, \dots)



