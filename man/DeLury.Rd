\name{DeLury}
\alias{DeLury}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ DeLury's Method for Population Size Estimation }
\description{
  Computes DeLury's method or Leslie's method for estimating a
  biological population size.

}
\usage{
DeLury(catch, effort, type = c("DeLury","Leslie"), ricker = FALSE)
       
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{catch, effort}{
  Catch and effort. These should be numeric vectors of equal length.

  }
  \item{type}{
  Character specifying which of the DeLury or Leslie models is to be fitted.
  The default is the first value.

  }
  \item{ricker}{
  Logical. If \code{TRUE} then the Ricker (1975) modification is computed.

  }
}
\details{
  This simple function implements the methods of DeLury (1947).
  These are called the DeLury and Leslie models.
  Note that there are many assumptions.
  These include: (i) Catch and effort records are available for a series
  of consecutive time intervals. The catch for a given time
  interval, specified by \eqn{t}, is \eqn{c(t)}, and the
  corresponding effort by \eqn{e(t)}.
  The \emph{catch per unit effort} (CPUE) for the time interval \eqn{t} is
  \eqn{C(t) = c(t)/e(t)}.
  Let \eqn{d(t)} represent the proportion of the population
  captured during the time interval \eqn{t}.
  Then \eqn{d(t) = k(t) e(t)} so that \eqn{k(t)} is the 
  proportion of the population captured during interval \eqn{t} by one unit
  of effort. Then \eqn{k(t)} is called the \emph{catchability},
  and the \emph{intensity} of effort is \eqn{e(t)}.
  Let \eqn{E(t)} and \eqn{K(t)} be the total effort and total catch
  up to interval \eqn{t}, and \eqn{N(t)} be the number of individuals
  in the population at time \eqn{t}.
  It is good idea to plot
  \eqn{\log(C(t))} against \eqn{E(t)} for \code{type = "DeLury"} and
  \eqn{C(t)} versus \eqn{K(t)} for \code{type = "Leslie"}.

The other assumptions are as follows.
%
   (ii) The population is closed---the population must be closed to sources
   of animals such as recruitment and immigration and losses of animals
   due to natural mortality and emigration.
%
   (iii) Catchability is constant over the period of removals.
%
   (iv) The units of effort are independent, i.e., the individual units
   of the method of capture (i.e., nets, traps, etc) do not compete with
   each other.
%
   (v) All fish are equally vulnerable to the method of capture---source
   of error may include gear saturation and trap-happy or trap-shy
   individuals.
%
   (vi) Enough fish must be removed to substantially reduce the CPUE.
%
   (vii) The catches may remove less than 2\% of the population.
%
   Also, the usual assumptions of simple regression such as
%
   (viii) random sampling,
%
   (ix) the independent variable(s) are measured without error---both
   catches and effort should be known, not estimated,
%
   (x) a line describes the data,
%
  (xi) the errors are independent and normally distributed.



}
\value{
  A list with the following components.


  \item{catch, effort }{
  Catch and effort. Same as the original vectors.
  These correspond to \eqn{c(t)} and \eqn{e(t)} respectively.


  }
  \item{type, ricker}{Same as input.


  }
  \item{N0}{an estimate of the population size at time 0.
  Only valid if the assumptions are satisfied.


  }
  \item{CPUE}{Catch Per Unit Effort \eqn{=C(t)}. }
  \item{K, E}{\eqn{K(t)}, \eqn{E(t)}. Only one is computed
  depending on \code{type}. }
  \item{lmfit}{
  the \code{\link[stats:lm]{lm}} object from the
  fit of \code{log(CPUE)} on \code{K} (when \code{type = "Leslie"}).
  Note that the \code{x} component of the object is the model matrix.


  }

}
\references{


DeLury, D. B. (1947)
On the estimation of biological populations.
\emph{Biometrics},
\bold{3}, 145--167.


Ricker, W. E. (1975)
Computation and interpretation of biological
statistics of fish populations.
\emph{Bull. Fish. Res. Bd. Can.},
\bold{191}, 382--


Yee, T. W. (2010)
VGLMs and VGAMs: an overview for applications in fisheries research.
\emph{Fisheries Research},
\bold{101}, 116--126.




}
\author{ T. W. Yee. }
\note{
The data in the example below comes from DeLury (1947), and
some plots of his are reproduced.
Note that he used log to base 10 whereas natural logs are used here.
His plots had some observations obscured by the y-axis!


The DeLury method is not applicable to the data frame
\code{\link{wffc.nc}} since the 2008 World Fly Fishing Competition was
strictly catch-and-release.


}
\seealso{ \code{\link{wffc.nc}}. }
\examples{
pounds <- c(  147, 2796, 6888, 7723, 5330, 8839, 6324, 3569, 8120, 8084,
            8252, 8411, 6757, 1152, 1500, 11945, 6995, 5851, 3221, 6345,
            3035, 6271, 5567, 3017, 4559, 4721, 3613,  473,  928, 2784,
            2375, 2640, 3569)
traps  <- c(  200, 3780, 7174, 8850, 5793, 9504, 6655, 3685, 8202, 8585,
            9105, 9069, 7920, 1215, 1471, 11597, 8470, 7770, 3430, 7970,
            4740, 8144, 7965, 5198, 7115, 8585, 6935, 1060, 2070, 5725,
            5235, 5480, 8300)
table1 <- DeLury(pounds/1000, traps/1000)

\dontrun{
with(table1, plot(1+log(CPUE) ~ E, las = 1, pch = 19, main = "DeLury method",
     xlab = "E(t)", ylab = "1 + log(C(t))", col = "blue"))
}
omitIndices <- -(1:16)
table1b <- DeLury(pounds[omitIndices]/1000, traps[omitIndices]/1000)
\dontrun{
with(table1b, plot(1+log(CPUE) ~ E, las = 1, pch = 19, main = "DeLury method",
     xlab = "E(t)", ylab = "1 + log(C(t))", col = "blue"))
mylmfit <- with(table1b, lmfit)
lines(mylmfit$x[, 2], 1 + predict.lm(mylmfit), col = "red", lty = "dashed")
}


omitIndices <- -(1:16)
table2 <- DeLury(pounds[omitIndices]/1000, traps[omitIndices]/1000, type = "L")
\dontrun{
with(table2, plot(CPUE ~ K, las = 1, pch = 19,
     main = "Leslie method; Fig. III",
     xlab = "K(t)", ylab = "C(t)", col = "blue"))
mylmfit <- with(table2, lmfit)
abline(a = coef(mylmfit)[1], b = coef(mylmfit)[2],
       col = "orange", lty = "dashed")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }


