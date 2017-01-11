\name{lakeO}
\alias{lakeO}
\docType{data}
\title{
  Annual catches on Lake Otamangakau from October 1974 to October 1989

%%   ~~ data name/kind ... ~~
}
\description{
  Rainbow and brown trout catches by a Mr Swainson at
  Lake Otamangakau in the central North Island of New Zealand
  during the 1970s and 1980s.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(lakeO)}
\format{
  A data frame with 15 observations on the following 5 variables.
  \describe{
    \item{\code{year}}{a numeric vector,
  the season began on 1 October of the year and ended 12 months later.
% Hence the fishing ended around October 1989.


  }
  \item{\code{total.fish}}{a numeric vector,
  the total number of fish caught during the season.
  Simply the sum of brown and rainbow trout.


  }
  \item{\code{brown}}{a numeric vector,
  the number of brown trout
  (\emph{Salmo trutta})
  caught.


  }
  \item{\code{rainbow}}{a numeric vector,
  the number of rainbow trout
  (\emph{Oncorhynchus mykiss})
  caught.


  }
  \item{\code{visits}}{a numeric vector,
  the number of visits during the season that the angler made to
  the lake.
  It is necessary to assume that the visits were of an equal time length
  in order to interpret the usual Poisson regressions.

  }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~


  The data was extracted from the season summaries at Lake Otamangakau
  by Anthony Swainson
  for the seasons 1974--75 to 1988--89.
% Note however that the final year's data
% was cut off from the scanned version.
  Mr Swainson was one of a small group of regular fly fishing anglers
  and kept a diary of his catches.
  Lake Otamangakau is a lake of area 1.8 squared km and has a
  maximum depth of about 12m, and is located
  in the central North Island of New Zealand.
  It is trout-infested and known for its trophy-sized fish.


  See also \code{\link[VGAMdata]{trapO}}.


}
\source{

  Table 7.2 of the reference below.
  Thanks to Dr Michel Dedual for a copy of the report and for help
  reading the final year's data.
  The report is available from TWY on request.


% p.43


%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{

  {Dedual, M.  and MacLean, G. and Rowe, D. and Cudby, E.},
  \emph{The Trout Population and Fishery of {L}ake {O}tamangakau---Interim Report}.
  {National Institute of Water and Atmospheric Research},
  {Hamilton, New Zealand}.
  Consultancy Report Project No. {ELE70207},
  (Dec 1996).



%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(lakeO)
lakeO
summary(lakeO)
}
\keyword{datasets}


