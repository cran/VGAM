\name{lms.yjn}
\alias{lms.yjn}
\alias{lms.yjn2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMS Quantile Regression with a Yeo-Johnson Transformation to Normality }
\description{
  LMS quantile regression with the Yeo-Johnson transformation
  to normality.
}
\usage{
lms.yjn(percentiles = c(25, 50, 75), zero = c("lambda", "sigma"),
        llambda = "identitylink", lsigma = "loge",
        idf.mu = 4, idf.sigma = 2,
        ilambda = 1, isigma = NULL, rule = c(10, 5),
        yoffset = NULL, diagW = FALSE, iters.diagW = 6)
lms.yjn2(percentiles=c(25,50,75), zero = c("lambda", "sigma"),
         llambda = "identitylink", lmu = "identitylink", lsigma = "loge",
         idf.mu = 4, idf.sigma = 2, ilambda = 1.0,
         isigma = NULL, yoffset = NULL, nsimEIM = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{percentiles}{
  A numerical vector containing values between 0 and 100,
  which are the quantiles. They will be returned as `fitted values'.

  }
  \item{zero}{
  See \code{\link{lms.bcn}}.

  }
  \item{llambda, lmu, lsigma}{
  See \code{\link{lms.bcn}}.

  }
  \item{idf.mu, idf.sigma}{
  See \code{\link{lms.bcn}}.

  }

  \item{ilambda, isigma}{
  See \code{\link{lms.bcn}}.

  }

  \item{rule}{
  Number of abscissae used in the Gaussian integration
  scheme to work out elements of the weight matrices.
  The values given are the possible choices, with the first value
  being the default.
  The larger the value, the more accurate the approximation is likely
  to be but involving more computational expense.

  }
  \item{yoffset}{
  A value to be added to the response y, for the purpose
  of centering the response before fitting the model to the data.
  The default value, \code{NULL}, means \code{-median(y)} is used, so that
  the response actually used has median zero. The \code{yoffset} is
  saved on the object and used during prediction.

  }
  \item{diagW}{
  Logical.
  This argument is offered because the expected information matrix may not
  be positive-definite. Using the diagonal elements of this matrix results
  in a higher chance of it being positive-definite, however convergence will
  be very slow.
  If \code{TRUE}, then the first \code{iters.diagW} iterations will
  use the diagonal of the expected information matrix.
  The default is \code{FALSE}, meaning faster convergence.

  }
  \item{iters.diagW}{
  Integer. Number of iterations in which the
  diagonal elements of the expected information matrix are used.
  Only used if \code{diagW = TRUE}.

  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
 Given a value of the covariate, this function applies a Yeo-Johnson
 transformation to the response to best obtain normality. The parameters
 chosen to do this are estimated by maximum likelihood or penalized
 maximum likelihood.
 The function \code{lms.yjn2()} estimates the expected information
 matrices using simulation (and is consequently slower) while
 \code{lms.yjn()} uses numerical integration.
  Try the other if one function fails.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.



}
\references{
Yeo, I.-K. and Johnson, R. A. (2000)
A new family of power transformations to improve normality or symmetry.
\emph{Biometrika},
\bold{87}, 954--959.

Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.

Yee, T. W. (2002)
An Implementation for Regression Quantile Estimation.
Pages 3--14.
In: Haerdle, W. and Ronz, B.,
\emph{Proceedings in Computational Statistics COMPSTAT 2002}.
Heidelberg: Physica-Verlag.


%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}
\author{ Thomas W. Yee }
\note{
The response may contain both positive and negative values.
In contrast, the LMS-Box-Cox-normal and LMS-Box-Cox-gamma
methods only handle a positive response because the
Box-Cox transformation cannot handle negative values.


Some other notes can be found at \code{\link{lms.bcn}}.


}

\section{Warning }{
The computations are not simple, therefore convergence may fail.
In that case, try different starting values.


The generic function \code{predict}, when applied to a
\code{lms.yjn} fit, does not add back the \code{yoffset} value.


}
\seealso{
  \code{\link{lms.bcn}},
  \code{\link{lms.bcg}},
  \code{\link{qtplot.lmscreg}},
  \code{\link{deplot.lmscreg}},
  \code{\link{cdf.lmscreg}},
  \code{\link{bmi.nz}},
  \code{\link{amlnormal}}.


}
\examples{
fit <- vgam(BMI ~ s(age, df = 4), lms.yjn, bmi.nz, trace = TRUE)
head(predict(fit))
head(fitted(fit))
head(bmi.nz)
# Person 1 is near the lower quartile of BMI amongst people his age
head(cdf(fit))

\dontrun{
# Quantile plot
par(bty = "l", mar = c(5, 4, 4, 3) + 0.1, xpd = TRUE)
qtplot(fit, percentiles = c(5, 50, 90, 99), main = "Quantiles",
       xlim = c(15, 90), las = 1, ylab = "BMI", lwd = 2, lcol = 4)

# Density plot
ygrid <- seq(15, 43, len = 100)  # BMI ranges
par(mfrow = c(1, 1), lwd = 2)
(aa <- deplot(fit, x0 = 20, y = ygrid, xlab = "BMI", col = "black",
    main = "Density functions at Age = 20 (black), 42 (red) and 55 (blue)"))
aa <- deplot(fit, x0 = 42, y = ygrid, add = TRUE, llty = 2, col = "red")
aa <- deplot(fit, x0 = 55, y = ygrid, add = TRUE, llty = 4, col = "blue",
             Attach = TRUE)
with(aa@post, deplot)  # Contains density function values; == a@post$deplot
}
}
\keyword{models}
\keyword{regression}

