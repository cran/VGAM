\name{huse}
\alias{huse}
\docType{data}
\title{
Harvard University Numbers of Ladder Faculty by School and Ethnicity

%%   ~~ data name/kind ... ~~
}
\description{
  A two-way table of counts; there are 14
  schools and 5 race/ethnicities.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(huse)}
\format{
  The format is:
 chr "huse"
}
\details{
  Ladder faculty members of Harvard University are cross-classified by
  their school and their race/ethnicity.
  This was for the period 2010--1.
  Ladder Faculty are defined as Assistant Professors or Convertible
  Instructors, Associate Professors, and Professors that have
  been appointed in certain Schools.


  Abbreviations:
  FAS = Faculty of Arts and Sciences = Humanities + Social Sciences +
  Natural Sciences + SEAS,
  Natural Sciences = Life Sciences + Physical Sciences,
  SEAS = School of Engineering and Applied Sciences,
  HBS = Harvard Business School,
  HMS = Harvard Medical School,
  HSPH = Harvard School of Public Health,
  HLS = Harvard Law School,
  HKS = Harvard Kennedy School,
  HGSE = Harvard Graduate School of Education,
  GSD = Graduate School of Design ,
  HDS = Harvard Divinity School,
  HSDM = Harvard School of Dental Medicine.


  See the URL below for many technical details supporting the data.
  The table was constructed from pp.31--2 from the source.


%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
 \url{http://www.provost.harvard.edu/institutional_research/factbook.php}

%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  \emph{Harvard University Office of the Senior Vice Provost Faculty
  Development \& Diversity: 2010 Annual Report}.


%%  ~~ possibly secondary sources and usages ~~
}
\seealso{
    \code{\link{hued}},
    \code{\link{huie}}.

}

\examples{
print(huse)
## maybe str(huse) ; plot(huse) ...
}
\keyword{datasets}
