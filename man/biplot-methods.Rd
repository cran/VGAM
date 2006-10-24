\name{biplot-methods}
\docType{methods}
\alias{biplot,rrvglm-method}
\alias{biplot,qrrvglm-method}
\title{ Biplot of Constrained Regression Models }
\description{
  \code{biplot} is a generic function applied to RR-VGLMs and QRR-VGLMs etc.
  These apply to rank-1 and rank-2 models of these only.
  For RR-VGLMs these plot the second latent variable scores against the
  first latent variable scores.

}
%\usage{
%  \S4method{biplot}{cao,Coef.cao}(object, ...)
%}

\section{Methods}{
\describe{

\item{x}{
  The object from which the latent variables are
  extracted and/or plotted.

}
}
}
\note{
  See \code{\link{lvplot}} which is very much related to biplots.

}

\keyword{methods}
\keyword{classes}
%\keyword{ ~~ other possible keyword(s)}
\keyword{models}
\keyword{regression}



