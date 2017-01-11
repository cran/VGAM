\name{SurvS4-class}
\docType{class}
\alias{SurvS4-class}
%%%% 20120216 \alias{print,SurvS4-method}
\alias{show,SurvS4-method}

\title{Class "SurvS4" }
\description{ S4 version of the Surv class. }
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{
%Class \code{"\linkS4class{Surv}"}, directly.
Class \code{"Surv"}, directly.
Class \code{"\linkS4class{matrix}"}, directly.
Class \code{"\linkS4class{oldClass}"}, by class "Surv", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 3, with explicit coerce.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.

}
\section{Methods}{
  \describe{
%   \item{print}{\code{signature(x = "SurvS4")}: ... }
    \item{show}{\code{signature(object = "SurvS4")}: ... }
   }
}
\references{
  See \pkg{survival}.


}
\author{
  T. W. Yee.


}
\note{
  The purpose of having \code{\link{SurvS4}} in \pkg{VGAM} is so that
  the same input can be fed into \code{\link{vglm}} as functions in
  \pkg{survival} such as \code{\link[survival]{survreg}}.


}

\section{Warning }{
  This code has not been thoroughly tested.

}

\seealso{
  \code{\link{SurvS4}}.


% or \code{\linkS4class{CLASSNAME}} for links to other classes



}
\examples{
showClass("SurvS4")
}
\keyword{classes}
