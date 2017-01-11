\name{rcqo}
\alias{rcqo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constrained Quadratic Ordination }
\description{
  Random generation for constrained quadratic ordination (CQO).
}
\usage{
rcqo(n, p, S, Rank = 1,
     family = c("poisson", "negbinomial", "binomial-poisson",
                "Binomial-negbinomial", "ordinal-poisson",
                "Ordinal-negbinomial", "gamma2"),
     eq.maximums = FALSE, eq.tolerances = TRUE, es.optimums = FALSE,
     lo.abundance = if (eq.maximums) hi.abundance else 10,
     hi.abundance = 100, sd.latvar = head(1.5/2^(0:3), Rank),
     sd.optimums = ifelse(es.optimums, 1.5/Rank, 1) *
                       ifelse(scale.latvar, sd.latvar, 1),
     sd.tolerances = 0.25, Kvector = 1, Shape = 1,
     sqrt.arg = FALSE, log.arg = FALSE, rhox = 0.5, breaks = 4,
     seed = NULL, optimums1.arg = NULL, Crow1positive = TRUE,
     xmat = NULL, scale.latvar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Number of sites. It is denoted by \eqn{n} below.

  }
  \item{p}{
    Number of environmental variables, including an intercept term.
    It is denoted by \eqn{p} below.
    Must be no less than \eqn{1+R} in value.

  }
  \item{S}{
    Number of species.
    It is denoted by \eqn{S} below.

  }
  \item{Rank}{
    The rank or the number of latent variables or true dimension
    of the data on the reduced space.
    This must be either 1, 2, 3 or 4.
    It is denoted by \eqn{R}.

  }
  \item{family}{
    What type of species data is to be returned.
    The first choice is the default.
    If binomial then a 0 means absence and 1 means presence.
    If ordinal then the \code{breaks} argument is passed into
    the \code{breaks} argument of \code{\link[base]{cut}}.
    Note that either the Poisson or negative binomial distributions
    are used to generate binomial and ordinal data, and that
    an upper-case choice is used for the negative binomial distribution
    (this makes it easier for the user).
    If \code{"gamma2"} then this is the 2-parameter gamma distribution.

%    ,
%    and the resulting values are
%    1,2,\ldots,\code{breaks} if \code{breaks} is a single integer zz
%    else zz.

  }
  \item{eq.maximums}{
    Logical. Does each species have the same maximum?
    See arguments \code{lo.abundance} and \code{hi.abundance}.


  }
  \item{eq.tolerances}{
    Logical. Does each species have the
    same tolerance? If \code{TRUE} then the common value is 1 along
    every latent variable, i.e., all species' tolerance matrices are the
    order-\eqn{R} identity matrix.


  }
  \item{es.optimums}{
    Logical. Do the species have equally spaced optimums?
    If \code{TRUE} then the quantity
    \eqn{S^{1/R}}{S^(1/R)} must be an
    integer with value 2 or more. That is, there has to be an
    appropriate number of species in total. This is so that a grid
    of optimum values is possible in \eqn{R}-dimensional
    latent variable space
    in order to place the species' optimums.
    Also see the argument \code{sd.tolerances}.


  }
  \item{lo.abundance, hi.abundance}{
    Numeric. These are recycled to a vector of length \eqn{S}.
    The species have a maximum
    between \code{lo.abundance} and \code{hi.abundance}. That is,
    at their optimal environment, the mean abundance of each
    species is between the two componentwise values. If \code{eq.maximums}
    is \code{TRUE} then \code{lo.abundance} and \code{hi.abundance}
    must have the same values.
    If \code{eq.maximums} is \code{FALSE} then the
    logarithm of the maximums are uniformly distributed between
    \code{log(lo.abundance)} and \code{log(hi.abundance)}.

  }
  \item{sd.latvar}{
    Numeric, of length \eqn{R}
    (recycled if necessary). Site scores along
    each latent variable have these standard deviation values.
    This must be a decreasing sequence of values because the first
    ordination axis contains the greatest spread of the species'
    site scores, followed by the second axis, followed by the third
    axis, etc.

  }
  \item{sd.optimums}{
    Numeric, of length \eqn{R} (recycled if necessary).
    If \code{es.optimums = FALSE} then,
    for the \eqn{r}th latent variable axis,
    the optimums of the species are generated from a
    normal distribution centered about 0.
    If \code{es.optimums = TRUE} then the \eqn{S} optimums
    are equally spaced about 0 along every latent variable axis.
    Regardless of the value of \code{es.optimums}, the optimums
    are then scaled to give standard deviation \code{sd.optimums[r]}.

  }
  \item{sd.tolerances}{
    Logical. If \code{eq.tolerances = FALSE} then, for the
    \eqn{r}th latent variable, the
    species' tolerances are
    chosen from a normal distribution with mean 1 and
    standard deviation
    \code{sd.tolerances[r]}.
    However, the first species \code{y1} has its tolerance matrix
    set equal to the order-\eqn{R} identity matrix.
    All tolerance matrices for all species are diagonal in this function.
    This argument is ignored if \code{eq.tolerances} is \code{TRUE},
    otherwise it is recycled to length \eqn{R} if necessary.

  }
  \item{Kvector}{
    A vector of positive \eqn{k} values (recycled to length \eqn{S}
    if necessary) for the negative binomial distribution
    (see \code{\link{negbinomial}} for details).
    Note that a natural default value does not exist, however the default
    value here is probably a realistic one, and that for large values
    of \eqn{\mu} one has \eqn{Var(Y) = \mu^2 / k}{Var(Y) = mu^2 / k}
    approximately.

  }
  \item{Shape}{
    A vector of positive \eqn{\lambda}{lambda} values (recycled to length
    \eqn{S} if necessary) for the 2-parameter gamma distribution (see
    \code{\link{gamma2}} for details).  Note that a natural default value
    does not exist, however the default value here is probably a realistic
    one, and that \eqn{Var(Y) = \mu^2 / \lambda}{Var(Y) = mu^2 / lambda}.

  }
  \item{sqrt.arg}{
    Logical. Take the square-root of the negative binomial counts?
    Assigning \code{sqrt.arg = TRUE} when \code{family="negbinomial"} means
    that the resulting species data can be considered very crudely to be
    approximately Poisson distributed.
    They will not integers in general but much easier (less numerical
    problems) to estimate using something like \code{cqo(..., family="poissonff")}.

  }
  \item{log.arg}{
    Logical. Take the logarithm of the gamma random variates?
    Assigning \code{log.arg = TRUE} when \code{family="gamma2"} means
    that the resulting species data can be considered very crudely to be
    approximately Gaussian distributed about its (quadratic) mean.
    The result is that it is much easier (less numerical
    problems) to estimate using something like \code{cqo(..., family="gaussianff")}.

  }
  \item{rhox}{
    Numeric, less than 1 in absolute value.
    The correlation between the environmental variables.
    The correlation matrix is a matrix of 1's along the diagonal
    and \code{rhox} in the off-diagonals.
    Note that each environmental variable is normally distributed
    with mean 0. The standard deviation of each environmental variable
    is chosen so that the site scores have the determined standard
    deviation, as given by argument \code{sd.latvar}.

  }
  \item{breaks}{
    If \code{family} is assigned an ordinal value then this argument
    is used to define the cutpoints. It is fed into the
    \code{breaks} argument of \code{\link[base]{cut}}.

  }
  \item{seed}{
    If given, it is passed into \code{\link[base:Random]{set.seed}}.
    This argument can be used to obtain reproducible results.
    If set, the value is saved as the \code{"seed"}
    attribute of the returned value. The default will
    not change the random generator state, and return
    \code{\link[base:Random]{.Random.seed}} as \code{"seed"} attribute.

  }
  \item{optimums1.arg}{
    If assigned and \code{Rank = 1} then these are the explicity optimums.
    Recycled to length \code{S}.


  }
  \item{Crow1positive}{
    See \code{\link{qrrvglm.control}} for details.

  }
  \item{xmat}{
   The
   \eqn{n} by  \eqn{p-1}
   environmental matrix can be inputted.

  }
  \item{scale.latvar}{
   Logical. If \code{FALSE} the argument
   \code{sd.latvar} is ignored and no scaling of the latent variable
   values is performed.

  }
}
\details{
  This function generates data coming from a constrained quadratic
  ordination (CQO) model. In particular,
  data coming from a \emph{species packing model} can be generated
  with this function.
  The species packing model states that species have equal tolerances,
  equal maximums, and optimums which are uniformly distributed over
  the latent variable space. This can be achieved by assigning
  the arguments \code{es.optimums = TRUE}, \code{eq.maximums = TRUE},
  \code{eq.tolerances = TRUE}.

  At present, the Poisson and negative binomial abundances are
  generated first using \code{lo.abundance} and \code{hi.abundance},
  and if \code{family} is binomial or ordinal then it is converted into
  these forms.

  In CQO theory the \eqn{n} by \eqn{p} matrix \eqn{X} is partitioned
  into two parts \eqn{X_1} and \eqn{X_2}. The matrix
  \eqn{X_2} contains the `real' environmental variables whereas
  the variables in \eqn{X_1} are just for adjustment purposes;
  they contain the intercept terms and other variables that one
  wants to adjust for when (primarily) looking at the
  variables in \eqn{X_2}.
  This function has \eqn{X_1} only being a matrix of ones,
  i.e., containing an intercept only.

}
\value{
  A \eqn{n} by \eqn{p-1+S} data frame with
  components and attributes.
  In the following the attributes are labelled with double
  quotes.
  \item{x2, x3, x4, \ldots, xp}{
    The environmental variables. This makes up the
    \eqn{n} by \eqn{p-1} \eqn{X_2} matrix.
    Note that \code{x1} is not present; it is effectively a vector
    of ones since it corresponds to an intercept term when
    \code{\link{cqo}} is applied to the data.

  }
  \item{y1, y2, x3, \ldots, yS}{
    The species data. This makes up the
    \eqn{n} by \eqn{S} matrix \eqn{Y}.
    This will be of the form described by the argument
    \code{family}.

  }
  \item{"concoefficients"}{
    The \eqn{p-1} by \eqn{R} matrix of
    constrained coefficients
    (or canonical coefficients).
    These are also known as weights or loadings.

  }
  \item{"formula"}{
    The formula involving the species and environmental variable names.
    This can be used directly in the \code{formula} argument of
    \code{\link{cqo}}.

  }
  \item{"log.maximums"}{
    The \eqn{S}-vector of species' maximums, on a log scale.
    These are uniformly distributed between
    \code{log(lo.abundance)} and \code{log(hi.abundance)}.


  }
  \item{"latvar"}{
    The \eqn{n} by \eqn{R} matrix of site scores.
    Each successive column (latent variable) has
    sample standard deviation
    equal to successive values of \code{sd.latvar}.

  }
  \item{"eta"}{
    The linear/additive predictor value.

  }
  \item{"optimums"}{
    The \eqn{S} by \eqn{R} matrix of species' optimums.

  }
  \item{"tolerances"}{
    The \eqn{S} by \eqn{R} matrix of species' tolerances.
    These are the
    square root of the diagonal elements of the tolerance matrices
    (recall that all tolerance matrices are restricted to being
    diagonal in this function).

  }
  Other attributes are \code{"break"},
  \code{"family"}, \code{"Rank"},
  \code{"lo.abundance"}, \code{"hi.abundance"},
  \code{"eq.tolerances"}, \code{"eq.maximums"},
  \code{"seed"} as used.

}
\references{

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


ter Braak, C. J. F. and Prentice, I. C. (1988)
A theory of gradient analysis.
\emph{Advances in Ecological Research},
\bold{18}, 271--317.

}
\author{ T. W. Yee }
\note{
  This function is under development and is not finished yet.
  There may be a few bugs.


  Yet to do: add an argument that allows absences to be equal
  to the first level if ordinal data is requested.


}
\seealso{
  \code{\link{cqo}},
  \code{\link{qrrvglm.control}},
  \code{\link[base]{cut}},
  \code{\link{binomialff}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}},
  \code{\link{gamma2}},
  \code{\link{gaussianff}}.


}

\examples{
\dontrun{
# Example 1: Species packing model:
n <- 100; p <- 5; S <- 5
mydata <- rcqo(n, p, S, es.opt = TRUE, eq.max = TRUE)
names(mydata)
(myform <- attr(mydata, "formula"))
fit <- cqo(myform, poissonff, mydata, Bestof = 3)  # eq.tol = TRUE
matplot(attr(mydata, "latvar"), mydata[,-(1:(p-1))], col = 1:S)
persp(fit, col = 1:S, add = TRUE)
lvplot(fit, lcol = 1:S, y = TRUE, pcol = 1:S)  # The same plot as above

# Compare the fitted model with the 'truth'
concoef(fit)  # The fitted model
attr(mydata, "concoefficients")  # The 'truth'

c(apply(attr(mydata, "latvar"), 2, sd),
  apply(latvar(fit), 2, sd))  # Both values should be approx equal


# Example 2: negative binomial data fitted using a Poisson model:
n <- 200; p <- 5; S <- 5
mydata <- rcqo(n, p, S, fam = "negbin", sqrt = TRUE)
myform <- attr(mydata, "formula")
fit <- cqo(myform, fam = poissonff, dat = mydata)  # I.tol = TRUE,
lvplot(fit, lcol = 1:S, y = TRUE, pcol = 1:S)
# Compare the fitted model with the 'truth'
concoef(fit)  # The fitted model
attr(mydata, "concoefficients")  # The 'truth'


# Example 3: gamma2 data fitted using a Gaussian model:
n <- 200; p <- 5; S <- 3
mydata <- rcqo(n, p, S, fam = "gamma2", log.arg = TRUE)
fit <- cqo(attr(mydata, "formula"),
           fam = gaussianff, data = mydata)  # I.tol = TRUE,
matplot(attr(mydata, "latvar"),
        exp(mydata[, -(1:(p-1))]), col = 1:S)  # 'raw' data
# Fitted model to transformed data:
lvplot(fit, lcol = 1:S, y = TRUE, pcol = 1:S)
# Compare the fitted model with the 'truth'
concoef(fit)  # The fitted model
attr(mydata, "concoefficients")  # The 'truth'
}
}
\keyword{distribution}
\keyword{datagen}


