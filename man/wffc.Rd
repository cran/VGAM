\name{wffc}
\alias{wffc}
\docType{data}
\title{ 2008 World Fly Fishing Championships }
\description{
  Capture records of the 2008 FIPS-MOUCHE
  World Fly Fishing Championships held in Rotorua, New Zealand during
  22--30 March 2008.

}
\usage{data(wffc)}
\format{
  A data frame with 4267 observations on the following 8 variables.
  Each row is a recorded capture.
  \describe{
    \item{\code{length}}{a numeric vector; length of fish in mm.}
    \item{\code{water}}{a factor with levels \code{Waihou},
    \code{Waimakariri}, \code{Whanganui}, \code{Otamangakau}, \code{Rotoaira}.
    These are known as Sectors IV, V, I, II, III respectively, and
    are also represented by the variable \code{sector}.
    }
    \item{\code{session}}{a numeric vector; a value from the set 1,2,\ldots,6.
    These are time ordered, and there were two sessions per competition day.}
    \item{\code{sector}}{a numeric vector; a value from the set 1,2,\ldots,5.}
    \item{\code{beatboat}}{a numeric vector; beat or boat number,
    a value from the set 1,2,\ldots,19.}
    \item{\code{comid}}{a numeric vector; the competitor's ID number. Uniquely
    identifies each competitor.
    These ID numbers actually correspond to their rankings
    in the individual level. }
    \item{\code{iname}}{a character vector; the individual competitor's name. }
    \item{\code{country}}{a character vector;
    what country the competitors represented.
The countries represented were
Australia (AUS),
Canada (CAN),
Croatia (CRO),
Czech Republic (CZE),
England (ENG),
Finland (FIN),
France (FRA),
Holland (NED),
Ireland (IRE),
Italy (ITA),
Japan (JPN),
Malta (MAL),
New Zealand (NZL),
Poland (POL),
Portugal (POR),
South Africa (RSA),
Slovakia (SVK),
USA (USA),
Wales (WAL).
  }
  }
}
\details{
  Details may be obtained at the official website and Yee (2008).
  Here is a brief summary.
  The three competition days were 28--30 March.
  Each session was fixed at 9.00am--12.00pm and 2.30--5.30pm daily.
  One of the sessions was a rest session.
  Each of 19 teams had 5 members, called A, B, C, D and E
  (there was a composite team, actually).
  The scoring system allocated 100 points to each eligible fish
  (minimum length was 18 cm)
  and 20 points for each cm of its length
  (rounded up to the nearest centimeter).
  Thus a 181mm or 190mm fish was worth 480 points.
  Each river was divided into 19 contiguous downstream beats labelled
  1,2,\ldots,19.
  Each lake was fished by 9 boats, each with two competitors
  except for one boat which only had one.
  Each competitor was randomly assigned to a beat/boat.

Competitors were ranked according to their placings at each sector-session
combination, and then these placings were summed. Those with the minimum
total placings were the winners, thus it was not necessarily those who had
the maximum points who won. For example, in Session 1 at the Waihou River,
each of the 19 competitors was ranked 1 (best) to 19 (worst) according
to the point system. This is the ``placing'' for that session. These
placings were added up over the 5 sessions to give the ``total placings''.

All sectors have naturally wild Rainbow trout (\emph{Oncorhynchus mykiss})
while Lake Otamangakau and the Whanganui River also holds Brown trout
(\emph{Salmo trutta}). Only these two species were targetted.
The species was not recorded electronically, however a post-analysis
of the paper score sheets from the two lakes showed that, approximately,
less than 5 percent were Brown trout.
It may be safely assumed that all the Waihou and Waimakariri
fish were Rainbow trout.
The gender of the fish were also not recorded electronically, and
anyway, distinguishing between male and female was very difficult
for small fish.

Although species and gender data were supposed to have been
collected at the time of capture the quality of these variables 
is rather poor and furthermore they were not recorded electronically.

% 11 out of (11 + 210) were brown trout, in Otamangakau.
%  52 were NAs.
%
%  3 out of ( 3 + 179) were brown trout, in Rotoaira.
%  19 were NAs.

  Note that some fish may have been caught more than once, hence
  these data do not represent individual fish but rather recorded captures.

}

\seealso{
  \code{\link{wffc.indiv}},
  \code{\link{wffc.teams}},
  \code{\link{wffc.nc}},
  \code{\link{wffc.P1}}.

}
\source{
  This data frame was adapted from the WFFC's spreadsheet.
  Special thanks goes to 
  Paul Dewar,
  Jill Mandeno,
  Ilkka Pirinen,
  and the other members of the Organising Committee of the 28th FIPS-Mouche
  World Fly Fishing Championships for access to the data.
  The assistance and feedback of Colin Shepherd is gratefully
  acknowledged.

}
\references{
  \url{http://www.2008worldflyfishingchamps.com}
  is the official website.

  Yee, T. W. (2008)
  Vector generalized linear and additive models, with
  applications to the 2008 World Fly Fishing Championships.
  In preparation.

}
\examples{
data(wffc)
summary(wffc)
with(wffc, table(water,session))

# Obtain some simple plots
waihou = with(wffc, wffc[water=="Waihou",])
waimak = with(wffc, wffc[water=="Waimakariri",])
whang = with(wffc, wffc[water=="Whanganui",])
otam = with(wffc, wffc[water=="Otamangakau",])
roto = with(wffc, wffc[water=="Rotoaira",])
minlength = min(wffc[,"length"])
maxlength = max(wffc[,"length"])
nwater = c("Waihou"=nrow(waihou), "Waimakariri"=nrow(waimak),
           "Whanganui"=nrow(whang), "Otamangakau"=nrow(otam),
           "Rotoaira"=nrow(roto))
\dontrun{
par(mfrow=c(2,3), las=1)
# Overall distribution of length
with(wffc, boxplot(length/10 ~ water, ylim=c(minlength, maxlength)/10,
                   border="blue", main="Length (cm)"))

# Overall distribution of LOG length
with(wffc, boxplot(length/10 ~ water, ylim=c(minlength, maxlength)/10,
                   border="blue", log="y",
                   main="Length (cm) on a log scale"))

# Overall distribution of number of captures
pie(nwater, border="blue", main="Proportion of captures",
    labels=names(nwater), density=10, col=1:length(nwater),
    angle=85+30* 1:length(nwater))

# Overall distribution of number of captures
with(wffc, barplot(nwater, main="Number of captures", cex.names=0.5))

# Overall distribution of proportion of number of captures
with(wffc, barplot(nwater / sum(nwater), cex.names=0.5,
                   main="Proportion of captures"))
# An interesting lake
with(roto, hist(length/10, xlab="Fish lengths (cm)", breaks=seq(18, 70, by=3),
                prob=TRUE, ylim=c(0, 0.08), border="blue", ylab="",
                main="Lake Rotoaira", lwd=2))
}
}
\keyword{datasets}
