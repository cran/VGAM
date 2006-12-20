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
     EqualMaxima = FALSE, EqualTolerances = TRUE, ESOptima = FALSE,
     loabundance = if (EqualMaxima) hiabundance else 10,
     hiabundance = 100, sdlv = c(1.5/2^(0:3))[1:Rank],
     sdOptima = ifelse(ESOptima, 1.5/Rank, 1) * sdlv,
     sdTolerances = 0.25, Kvector = 1, Shape = 1,
     sqrt = FALSE, Log = FALSE, rhox = 0.5, breaks = 4,
     seed = NULL, Crow1positive=TRUE)
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
  \item{EqualMaxima}{
    Logical. Does each species have the same maxima?
    See arguments \code{loabundance} and \code{hiabundance}.
    
  }
  \item{EqualTolerances}{
    Logical. Does each species have the
    same tolerance? If \code{TRUE} then the common value is 1 along
    every latent variable, i.e., all species' tolerance matrices are the
    order-\eqn{R} identity matrix.
    
  }
  \item{ESOptima}{
    Logical. Do the species have equally spaced optima?
    If \code{TRUE} then the quantity
    \eqn{S^{1/R}}{S^(1/R)} must be an
    integer with value 2 or more. That is, there has to be an
    appropriate number of species in total. This is so that a grid
    of optimum values is possible in \eqn{R}-dimensional
    latent variable space
    in order to place the species' optima.
    Also see the argument \code{sdTolerances}.
    
  }
  \item{loabundance, hiabundance}{
    Numeric. These are recycled to a vector of length \eqn{S}.
    The species have a maximum
    between \code{loabundance} and \code{hiabundance}. That is,
    at their optimal environment, the mean abundance of each
    species is between the two componentwise values. If \code{EqualMaxima}
    is \code{TRUE} then \code{loabundance} and \code{hiabundance}
    must have the same values.
    If \code{EqualMaxima} is \code{FALSE} then the
    logarithm of the maxima are uniformly distributed between
    \code{log(loabundance)} and \code{log(hiabundance)}.
    
  }
  \item{sdlv}{
    Numeric, of length \eqn{R}
    (recycled if necessary). Site scores along
    each latent variable have these standard deviation values.
    This must be a decreasing sequence of values because the first
    ordination axis contains the greatest spread of the species'
    site scores, followed by the second axis, followed by the third
    axis, etc.
        
  }
  \item{sdOptima}{
    Numeric, of length \eqn{R} (recycled if necessary).
    If \code{ESOptima=FALSE} then,
    for the \eqn{r}th latent variable axis,
    the optima of the species are generated from a
    normal distribution centered about 0.
    If \code{ESOptima=TRUE} then the \eqn{S} optima
    are equally spaced about 0 along every latent variable axis.
    Regardless of the value of \code{ESOptima}, the optima
    are then scaled to give standard deviation \code{sdOptima[r]}.
    
  }
  \item{sdTolerances}{
    Logical. If \code{EqualTolerances=FALSE} then, for the
    \eqn{r}th latent variable, the
    species' tolerances are
    chosen from a normal distribution with mean 1 and
    standard deviation
    \code{sdTolerances[r]}.
    However, the first species \code{y1} has its tolerance matrix
    set equal to the order-\eqn{R} identity matrix.
    All tolerance matrices for all species are diagonal in this function.
    This argument is ignored if \code{EqualTolerances} is \code{TRUE},
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
  \item{sqrt}{
    Logical. Take the square-root of the negative binomial counts?
    Assigning \code{sqrt=TRUE} when \code{family="negbinomial"} means
    that the resulting species data can be considered very crudely to be
    approximately Poisson distributed.
    They will not integers in general but much easier (less numerical
    problems) to estimate using something like \code{cqo(..., family="poissonff")}.

  }
  \item{Log}{
    Logical. Take the logarithm of the gamma random variates?
    Assigning \code{Log=TRUE} when \code{family="gamma2"} means
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
    deviation, as given by argument \code{sdlv}.
    
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
  \item{Crow1positive}{
    See \code{\link{qrrvglm.control}} for details.
    
  }
}
\details{
  This function generates data coming from a constrained quadratic
  ordination (CQO) model. In particular,
  data coming from a \emph{species packing model} can be generated
  with this function.
  The species packing model states that species have equal tolerances,
  equal maxima, and optima which are uniformly distributed over
  the latent variable space. This can be achieved by assigning
  the arguments \code{ESOptima = TRUE}, \code{EqualMaxima = TRUE},
  \code{EqualTolerances = TRUE}.

  At present, the Poisson and negative binomial abundances are
  generated first using \code{loabundance} and \code{hiabundance},
  and if \code{family} is binomial or ordinal then it is converted into
  these forms.

  In CQO theory the \eqn{n \times p}{n * p}
  matrix \eqn{X} is partitioned
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
  A \eqn{n \times (p-1+S)}{n * (p-1+S)} data frame with
  components and attributes.
  In the following the attributes are labelled with double
  quotes.
  \item{x2, x3, x4, \ldots, xp}{
    The environmental variables. This makes up the
    \eqn{n \times (p-1)}{n * (p-1)} \eqn{X_2} matrix.
    Note that \code{x1} is not present; it is effectively a vector
    of ones since it corresponds to an intercept term when
    \code{\link{cqo}} is applied to the data.
    
  }
  \item{y1, y2, x3, \ldots, yS}{
    The species data. This makes up the
    \eqn{n \times S}{n * S} matrix \eqn{Y}.
    This will be of the form described by the argument
    \code{family}.
    
  }
  \item{"ccoefficients"}{
    The \eqn{(p-1) \times R}{(p-1) * R} matrix of
    constrained coefficients
    (or canonical coefficients).
    These are also known as weights or loadings.
    
  }
  \item{"formula"}{
    The formula involving the species and environmental variable names.
    This can be used directly in the \code{formula} argument of
    \code{\link{cqo}}.
    
  }
  \item{"logmaxima"}{
    The \eqn{S}-vector of species' maxima, on a log scale.
    These are uniformly distributed between
    \code{log(loabundance)} and \code{log(hiabundance)}.
    

  }
  \item{"lv"}{
    The \eqn{n \times R}{n * R} matrix of site scores.
    Each successive column (latent variable) has
    sample standard deviation
    equal to successive values of \code{sdlv}.
    
  }
  \item{"optima"}{
    The \eqn{S \times R}{S * R} matrix of species' optima.

  }
  \item{"tolerances"}{
    The \eqn{S \times R}{S * R} matrix of species' tolerances.
    These are the
    square root of the diagonal elements of the tolerance matrices
    (recall that all tolerance matrices are restricted to being
    diagonal in this function).
    
  }
  Other attributes are \code{"break"},
  \code{"family"}, \code{"Rank"},
  \code{"loabundance"}, \code{"hiabundance"},
  \code{"EqualTolerances"}, \code{"EqualMaxima"},
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
  \code{gaussianff}.
}
\examples{

# Example 1: Species packing model:
n = 100; p = 5; S = 5
mydata = rcqo(n, p, S, ESOpt=TRUE, EqualMax=TRUE)
names(mydata)
myform = attr(mydata, "formula")
fit = cqo(myform, fam=poissonff, ITol=TRUE, data=mydata,
          Bestof=3)  # Fit a CQO model to the data
\dontrun{
matplot(attr(mydata, "lv"), mydata[,-(1:(p-1))], col=1:S)
persp(fit, col=1:S, add=TRUE)
lvplot(fit, lcol=1:S, y=TRUE, pcol=1:S)  # The same plot as above
}

# Compare the fitted model with the 'truth'
ccoef(fit)  # The fitted model
attr(mydata, "ccoefficients") # The 'truth'

c(sd(attr(mydata, "lv")), sd(lv(fit))) # Both values should be approx equal


# Example 2: negative binomial data fitted using a Poisson model:
n = 200; p = 5; S = 5
mydata = rcqo(n, p, S, fam="negbin", sqrt=TRUE)
myform = attr(mydata, "formula")
fit = cqo(myform, fam=poissonff, ITol=TRUE, dat=mydata)
\dontrun{
lvplot(fit, lcol=1:S, y=TRUE, pcol=1:S)
}
# Compare the fitted model with the 'truth'
ccoef(fit)  # The fitted model
attr(mydata, "ccoefficients") # The 'truth'


# Example 3: gamma2 data fitted using a Gaussian model:
n = 200; p = 5; S = 3
mydata = rcqo(n, p, S, fam="gamma2", Log=TRUE)
fit = cqo(attr(mydata, "formula"), fam=gaussianff, ITol=TRUE, dat=mydata)
\dontrun{
matplot(attr(mydata, "lv"), exp(mydata[,-(1:(p-1))]), col=1:S) # 'raw' data
lvplot(fit, lcol=1:S, y=TRUE, pcol=1:S)  # Fitted model to transformed data
}
# Compare the fitted model with the 'truth'
ccoef(fit)  # The fitted model
attr(mydata, "ccoefficients") # The 'truth'
}
\keyword{distribution}


