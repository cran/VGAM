\name{ducklings}
\alias{ducklings}
\docType{data}
\title{
Relative Frequencies of Serum Proteins in White Pekin Ducklings 

%%   ~~ data name/kind ... ~~
}
\description{
  Relative frequencies of serum proteins in white Pekin ducklings 
  as determined by electrophoresis.


}
\usage{data(ducklings)}
\format{
  The format is:
 chr "ducklings"

}
\details{
  Columns \code{p1}, \code{p2}, \code{p3}
  stand for pre-albumin, albumin, globulins respectively.
  These were collected from 3-week old white Pekin ducklings.
  Let \eqn{Y_1}{Y1} be proportional to the total milligrams of
  pre-albumin in the blood serum of a duckling.
  Similarly,
  let \eqn{Y_2}{Y2} and \eqn{Y_3}{Y3} be directly proportional
  to the same factor as \eqn{Y_1}{Y1} to the total milligrams 
  respectively of albumin and globulins in its blood serum.
  The proportion of pre-albumin is given by
  \eqn{Y_1/(Y_1 + Y_2 + Y_3)}{Y1/(Y1 + Y2 + Y3)},
  and similarly for the others.





%  Each set of 3 measurements is based on from 7 to 12 individual ducklings.






%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Mosimann, J. E.  (1962)
  On the compound multinomial distribution,
  the multivariate \eqn{\beta}{beta}-distribution,
  and correlations among proportions,
  {Biometrika},
  \bold{49}, 65--82.
}


\seealso{
  \code{\link{dirichlet}}.


}

%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
print(ducklings)
}
\keyword{datasets}
