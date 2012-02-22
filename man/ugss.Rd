\name{ugss}
\alias{ugss}
\docType{data}
\title{
  Undergraduate Statistics Students Lifestyle Questionnaire

%%   ~~ data name/kind ... ~~
}
\description{
  About 800 students studying undergraduate statistics
  were asked many lifestyle questions.
  
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(ugss)}
\format{
  A data frame with 804 observations on the following 29 variables.
\describe{
    \item{\code{sex}}{Gender, a factor, (female or male) }
    \item{\code{age}}{age in years, a numeric vector}
    \item{\code{eyes}}{eye colour, a factor,
                       (blue, brown, green, hazel or other)}
    \item{\code{piercings}}{Number of body piercings, a numeric vector}
    \item{\code{pierced}}{Any body piercings? a factor, (Yes or No)}
    \item{\code{tattoos}}{Number of tattoos, a numeric vector}
    \item{\code{tattooed}}{Any tattoos? a factor, (Yes or No) }
    \item{\code{glasses}}{Wears glasses etc.? a factor, (Yes or No)}
    \item{\code{sleep}}{Average number of hours of sleep per night,
                        a numeric vector}
    \item{\code{study}}{Average number of hours of study per week,
                        a numeric vector}
    \item{\code{tv}}{Average number of hours watching TV per week,
                        a numeric vector}
    \item{\code{movies}}{
  Number of movies seen at a cinema during the last 3 months,
                         a numeric vector
  }
    \item{\code{movies3m}}{Seen movies in last 3 months? 
                           a factor, (Yes or No)}
    \item{\code{sport}}{Favourite sport, a factor,
                        about 19 of them }
    \item{\code{entertainment}}{Favourite entertainment, a factor,
                        about 15 of them }
    \item{\code{fruit}}{Favourite fruit a factor,
                        about 13 of them }
    \item{\code{income}}{Average income during semester per week,
                         a numeric vector }
    \item{\code{rent}}{Amount spent on rent or room and board per week,
                         a numeric vector  }
    \item{\code{clothes}}{Average amount spent on clothes per month,
                         a numeric vector }
    \item{\code{hair}}{Average cost to get a hair-cut,
                         a numeric vector }
    \item{\code{tobacco}}{Average amount spent on tobacco per week,
                          a numeric vector}
    \item{\code{smokes}}{Smokes? a factor, (Yes or No) }
    \item{\code{alcohol}}{Average amount spent on alcohol per week,
                          a numeric vector }
    \item{\code{buy.alcohol}}{Buys (purchases) alcohol? a factor, (Yes or No) }
    \item{\code{sendtxt}}{Average number text messages sent per day,
                          a numeric vector.}
    \item{\code{receivetxt}}{Average number text messages received per day,
                          a numeric vector.}
    \item{\code{txts}}{Uses text messaging? a factor, (Yes or No) }
    \item{\code{country}}{Country of birth,  a factor,
                          about 54 of them  }
    \item{\code{status}}{Student status, a factor,
                         (International, NZ.Citizen, NZ.Resident) }
  }
}
%%\format{
%%  The format is:
%% chr "ugss"
%%}
\details{
  This data was collected online and anonymously in 2010.
  The respondents were students studying an undergraduate statistics
  course at a New Zealand university.
  Possibly there are duplicate students (due to failing and
  re-enrolling).
  All monies are in NZD.
  Note the data has had minimal checking.
  Most numerical variables tend to have measurement error, and all of
  them happen to be all integer-valued.


%%  ~~ If necessary, more details than the __description__ above
}


%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
summary(ugss)
}
\keyword{datasets}
