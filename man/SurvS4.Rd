\name{SurvS4}
\alias{SurvS4}
\alias{is.SurvS4}
%%%% 20120216 \alias{print.SurvS4}
\alias{show.SurvS4}
\alias{Math.SurvS4}
\alias{Summary.SurvS4}
\alias{[.SurvS4}
\alias{format.SurvS4}
\alias{as.data.frame.SurvS4}
\alias{as.character.SurvS4}
\alias{is.na.SurvS4}
\alias{Ops.SurvS4}
\title{
Create a Survival Object
}
\description{
  Create a survival object, usually used as a response
  variable in a model formula.

}
\usage{
SurvS4(time, time2, event, type =, origin = 0)
is.SurvS4(x)
}
\arguments{
  \item{time}{
    for right censored data, this is the follow up time.  For interval
    data, the first argument is the starting time for the interval.


  }
  \item{x}{
    any R object.


  }
  \item{event}{
    The status indicator, normally 0=alive, 1=dead.  Other choices are
    \code{TRUE}/\code{FALSE} (\code{TRUE} = death) or 1/2 (2=death). For
    interval censored data, the status indicator is 0=right censored,
    1=event at \code{time}, 2=left censored, 3=interval censored.
    Although unusual, the event indicator can be omitted, in which case
    all subjects are assumed to have an event.


  }
  \item{time2}{
    ending time of the interval for interval censored  or counting
    process data only.  Intervals are assumed to be open on the left and
    closed on the right, \code{(start, end]}.  For counting process
    data, \code{event} indicates whether an event occurred at the end of
    the interval.


  }
  \item{type}{
    character string specifying the type of censoring. Possible values
    are \code{"right"}, \code{"left"}, \code{"counting"},
    \code{"interval"}, or \code{"interval2"}.  The default is
    \code{"right"} or \code{"counting"} depending on whether the
    \code{time2} argument is absent or present, respectively.


  }
  \item{origin}{
    for counting process data, the hazard function origin.  This is most
    often used in conjunction with a model containing time dependent
    strata in order to align the subjects properly when they cross over
    from one strata to another.


  }
}
\value{
  An object of class \code{SurvS4} (formerly \code{Surv}).
  There are methods for \code{print}, \code{is.na}, and
  subscripting survival objects. \code{SurvS4} objects are
  implemented as a matrix of 2 or 3 columns.

  In the case of \code{is.SurvS4}, a logical value
  \code{TRUE} if \code{x} inherits from class
  \code{"SurvS4"}, otherwise a \code{FALSE}.


}
\details{
  Typical usages are
\preformatted{
SurvS4(time, event)
SurvS4(time, time2, event, type=, origin=0)
}


In theory it is possible to represent interval censored
data without a third column containing the explicit status.
Exact, right censored, left censored and interval censored
observation would be represented as intervals of (a,a),
(a, infinity), (-infinity,b), and (a,b) respectively; each
specifying the interval within which the event is known to
have occurred.


If \code{type = "interval2"} then the representation given
above is assumed, with NA taking the place of infinity.
If `type="interval" \code{event} must be given.
If \code{event} is \code{0}, \code{1}, or \code{2},
the relevant information is assumed to be contained in
\code{time}, the value in \code{time2} is ignored, and the
second column of the result will contain a placeholder.


Presently, the only methods allowing interval
censored data are the parametric models computed by
\code{\link[survival]{survreg}}, so the distinction between
open and closed intervals is unimportant.  The distinction
is important for counting process data and the Cox model.


The function tries to distinguish between the use of 0/1
and 1/2 coding for left and right censored data using
\code{if (max(status)==2)}.  If 1/2 coding is used and all
the subjects are censored, it will guess wrong.  Use 0/1
coding in this case.


}


\author{
  The code and documentation comes from \pkg{survival}.
  Slight modifications have been made for conversion to S4 by T. W. Yee.
  Also, for \code{"interval"} data, \code{as.character.SurvS4()} has
  been modified to print intervals of the form
  \code{(start, end]} and not
  \code{[start, end]} as previously.
  (This makes a difference for discrete data, such as for
  \code{\link{cens.poisson}}).
  All \pkg{VGAM} family functions beginning with \code{"cen"} require
  the packaging function \code{Surv} to format the input.


}
\note{
  The purpose of having \code{SurvS4} in \pkg{VGAM} is so that
  the same input can be fed into \code{\link{vglm}} as functions in
  \pkg{survival} such as \code{\link[survival]{survreg}}.  The class
  name has been changed from \code{"Surv"} to \code{"SurvS4"}; see
  \code{\link{SurvS4-class}}.


  The format \code{J+} is interpreted in \pkg{VGAM} as \eqn{\ge J}.
  If \code{type="interval"} then these should not be used in \pkg{VGAM}:
  \code{(L,U-]} or  \code{(L,U+]}.
% zz is this for type="count" only?


}

\seealso{
  \code{\link{SurvS4-class}},
  \code{\link{cens.poisson}},
  \code{\link[survival]{survreg}},
  \code{\link{leukemia}}.


% \code{\link[survival]{coxph}},
% \code{\link[survival]{survfit}},


}
\examples{
with(leukemia, SurvS4(time, status))
class(with(leukemia, SurvS4(time, status)))
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
% with(heart, SurvS4(start,stop,event))
