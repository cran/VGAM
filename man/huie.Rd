\name{huie}
\alias{huie}
\docType{data}
\title{
Harvard University International Enrollments

%%   ~~ data name/kind ... ~~
}
\description{
  A two-way table of counts; there are 12
  degrees and 8 areas of the world.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(huie)}
\format{
  The format is:
 chr "huie"
}
\details{
  The rownames and colnames have been edited.
  The full colnames are:
Africa,
Asia,
Europe,
Caribbean and Central and and South America,
Middle East,
North America,
Oceania,
Stateless.

The data was for the autumn (Fall) of 2010.
GSAS stands for Graduate School of Arts and Sciences.
See the URL below for more technical details supporting the data.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  \url{http://www.provost.harvard.edu/institutional_research/factbook.php}


%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\seealso{
    \code{\link{hued}},
    \code{\link{huse}}.

}

%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
print(huie)
## maybe str(huie) ; plot(huie) ...
}
\keyword{datasets}
