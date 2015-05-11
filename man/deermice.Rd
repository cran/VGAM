\name{deermice}
%\alias{Perom}
\alias{deermice}
\docType{data}
\title{
  Captures of Peromyscus maniculatus, also known as deer mice.


%%   ~~ data name/kind ... ~~
}
\description{

  Captures of \emph{Peromyscus maniculatus} collected at East
  Stuart Gulch, Colorado, USA.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
% data(Perom)
\usage{
  data(deermice)
}
\format{
  The format is a data frame.


}
\details{

  \emph{Peromyscus maniculatus} is a rodent native to North America.
  The deer mouse is small in size, only about 8 to 10 cm long,
  not counting the length of the tail.


  Originally,
  the columns of this data frame
  represent the sex (\code{m} or \code{f}),
  the ages (\code{y}: young, \code{sa}: semi-adult, \code{a}: adult),
  the weights in grams, and the
  capture histories of 38 individuals over 6 trapping
  occasions (1: captured, 0: not captured).


  The data set was collected by V. Reid and distributed
  with the \pkg{CAPTURE} program of Otis et al. (1978).


  \code{deermice} has 38 deermice whereas
  \code{Perom} had 36 deermice
  (\code{Perom} has been withdrawn.)
  In \code{deermice} the two semi-adults have been classified as adults.
  The \code{sex} variable has 1 for female, and 0 for male.


%%  ~~ If necessary, more details than the __description__ above ~~
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{


Huggins, R. M. (1991)
Some practical aspects of a conditional likelihood
approach to capture experiments.
\emph{Biometrics},
\bold{47}, 725--732.


  Otis, D. L. et al. (1978)
  Statistical inference from capture data on closed animal populations,
  \emph{Wildlife Monographs},
  \bold{62}, 3--135.


%%  ~~ possibly secondary sources and usages ~~
}

\seealso{
    \code{\link[VGAM:posbernoulli.b]{posbernoulli.b}},
    \code{\link[VGAM:posbernoulli.t]{posbernoulli.t}},
    \code{\link{fill1}}.


}

\examples{
head(deermice)
\dontrun{
fit1 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + age,
             posbernoulli.t(parallel.t = TRUE), data = deermice, trace = TRUE)
coef(fit1)
coef(fit1, matrix = TRUE)
}
}
\keyword{datasets}
