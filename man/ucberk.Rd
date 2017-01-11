\name{ucberk}
\alias{ucberk}
\docType{data}
\title{ University California Berkeley Graduate Admissions }
\description{
  University California Berkeley Graduate Admissions: counts
  cross-classified by acceptance/rejection and gender, for
  the six largest departments.

}
\usage{data(ucberk)}
\format{
  A data frame with 6 departmental groups with the following 5 columns.
  \describe{
    \item{m.deny}{Counts of men denied admission. }
    \item{m.admit}{Counts of men admitted. }
    \item{w.deny}{Counts of women denied admission. }
    \item{w.admit}{Counts of women admitted. }
    \item{dept}{Department (the six largest),
                called \code{A}, code{B}, \dots, code{F}.
    }
  }
}
\details{
  From Bickel et al. (1975),
  the data consists of applications for admission to graduate
  study at the University of California, Berkeley, for the
  fall 1973 quarter.
  In the admissions cycle for that quarter,
  the Graduate Division at Berkeley received approximately
  15,000 applications, some of which were later withdrawn or
  transferred to a different proposed entry quarter by the
  applicants.  Of the applications finally remaining for the
  fall 1973 cycle 12,763 were sufficiently complete to permit
  a decision.
  There were about 101 graduate department and
  interdepartmental graduate majors.  There were 8442 male
  applicants and 4321 female applicants. About 44 percent of
  the males and about 35 percent of the females were admitted.
  The data are well-known for illustrating Simpson's paradox.




}
%\source{
%
%
%}
\references{
  Bickel, P. J., Hammel, E. A. and O'Connell, J. W. (1975)
  Sex bias in graduate admissions: data from Berkeley.
  \emph{Science}, \bold{187}(4175): 398--404.


  Freedman, D., Pisani, R. and Purves, R. (1998)
  Chapter 2 of \emph{Statistics}, 3rd. ed.,
  W. W. Norton & Company.



}
\examples{
summary(ucberk)
}
\keyword{datasets}
% 7 February 1975
% Bickel, et al., 187 (4175): 398-404

