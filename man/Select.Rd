\name{Select}
\alias{Select}
\alias{subsetcol}
% \alias{subsetc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Select Variables for a Formula Response or the RHS of a Formula


%%  ~~function to do ... ~~
}
\description{
  Select variables from a data frame whose names
  begin with a certain character string.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Select(data = list(), prefix = "y",
       lhs = NULL, rhs = NULL, rhs2 = NULL, rhs3 = NULL,
       as.character = FALSE, as.formula.arg = FALSE, tilde = TRUE,
       exclude = NULL, sort.arg = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A data frame or a matrix.


%%     ~~Describe \code{data} here~~
}
\item{prefix}{
  A vector of character strings, or a logical.
  If a character then
  the variables chosen from \code{data} begin with the
  value of \code{prefix}.
  If a logical then
  only \code{TRUE} is accepted and all the variables
  in \code{data} are chosen.



%%     ~~Describe \code{prefix} here~~
}
\item{lhs}{
  A character string.
  The response of a formula.


%%     ~~Describe \code{lhs} here~~
}
\item{rhs}{
  A character string.
  Included as part of the RHS a formula.
  Set \code{rhs = "0"} to suppress the intercept.


%%     ~~Describe \code{rhs} here~~
}
\item{rhs2, rhs3}{
  Same as \code{rhs} but appended to its RHS,
  i.e., \code{paste0(rhs, " + ", rhs2, " + ", rhs3)}.
  If used, \code{rhs} should be used first,
  and then possibly \code{rhs2}
  and then possibly \code{rhs3}.




%%     ~~Describe \code{rhs} here~~
}
\item{as.character}{
  Logical.
  Return the answer as a character string?


%%     ~~Describe \code{as.character} here~~
}
\item{as.formula.arg}{
  Logical.
  Is the answer a formula?


%%     ~~Describe \code{as.formula.arg} here~~
}
\item{tilde}{
  Logical.
  If \code{as.character} and \code{as.formula.arg}
  are both \code{TRUE}
  then include the tilde in the formula?


}
\item{exclude}{
  Vector of character strings.
  Exclude these variables explicitly.


%%     ~~Describe \code{exclude} here~~
}
\item{sort.arg}{
  Logical.
  Sort the variables?


%%     ~~Describe \code{sort.arg} here~~
}
}
\details{
  This is meant as a utility function to avoid manually:
  (i) making a \code{\link[base]{cbind}} call to construct
      a big matrix response,
  and
  (ii) constructing a formula involving a lot of terms.
  The savings can be made because the variables of interest
  begin with some prefix, e.g., with the character \code{"y"}.



}
\value{
  If \code{as.character = FALSE} and
  \code{as.formula.arg = FALSE} then a matrix such
  as \code{cbind(y1, y2, y3)}.
  If \code{as.character = TRUE} and
  \code{as.formula.arg = FALSE} then a character string such
  as \code{"cbind(y1, y2, y3)"}.



  If \code{as.character = FALSE} and
  \code{as.formula.arg = TRUE} then a \code{\link[stats]{formula}} such
  as \code{lhs ~ y1 + y2 + y3}.
  If \code{as.character = TRUE} and
  \code{as.formula.arg = TRUE} then a character string such
  as \code{"lhs ~ y1 + y2 + y3"}.
  See the examples below.
  By default, if no variables beginning the the value of \code{prefix}
  is found then a \code{NULL} is returned.
  Setting \code{prefix = " "} is a way of selecting no variables.




%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
  T. W. Yee.


%%  ~~who you are~~
}
\note{
  This function is a bit experimental at this stage and
  may change in the short future.
  Some of its utility may be better achieved using
  \code{\link[base]{subset}} and its \code{select} argument,
  e.g., \code{subset(pdata, TRUE, select = y01:y10)}.



  For some models such as \code{\link{posbernoulli.t}} the
  order of the variables in the \code{xij} argument is
  crucial, therefore care must be taken with the
  argument \code{sort.arg}.
  In some instances, it may be good to rename variables
  \code{y1} to \code{y01},
  \code{y2} to \code{y02}, etc.
  when there are variables such as
  \code{y14}.



  Currently \code{subsetcol()} and \code{Select()} are identical.
  One of these functions might be withdrawn in the future.



%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{vglm}},
  \code{\link[base]{cbind}},
  \code{\link[base]{subset}},
  \code{\link[stats]{formula}},
  \code{\link{fill}}.




%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Pneumo <- pneumo
colnames(Pneumo) <- c("y1", "y2", "y3", "x2")  # The "y" variables are response
Pneumo$x1 <- 1; Pneumo$x3 <- 3; Pneumo$x <- 0; Pneumo$x4 <- 4  # Add these

Select(data = Pneumo)  # Same as with(Pneumo, cbind(y1, y2, y3))
Select(Pneumo, "x")
Select(Pneumo, "x", sort = FALSE, as.char = TRUE)
Select(Pneumo, "x", exclude = "x1")
Select(Pneumo, "x", exclude = "x1", as.char = TRUE)
Select(Pneumo, c("x", "y"))
Select(Pneumo, "z")  # Now returns a NULL
Select(Pneumo, " ")  # Now returns a NULL
Select(Pneumo, prefix = TRUE, as.formula = TRUE)
Select(Pneumo, "x", exclude = c("x3", "x1"), as.formula = TRUE,
       lhs = "cbind(y1, y2, y3)", rhs = "0")
Select(Pneumo, "x", exclude = "x1", as.formula = TRUE, as.char = TRUE,
       lhs = "cbind(y1, y2, y3)", rhs = "0")

# Now a 'real' example:
Huggins89table1 <- transform(Huggins89table1, x3.tij = t01)
tab1 <- subset(Huggins89table1,
               rowSums(Select(Huggins89table1, "y")) > 0)
# Same as
# subset(Huggins89table1, y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 > 0)

# Long way to do it:
fit.th <-
   vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08, y09, y10) ~ x2 + x3.tij,
        xij = list(x3.tij ~ t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 +
                            t09 + t10 - 1),
        posbernoulli.t(parallel.t = TRUE ~ x2 + x3.tij),
        data = tab1, trace = TRUE,
        form2 = ~ x2 + x3.tij + t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 +
                                t09 + t10)
# Short way to do it:
Fit.th <- vglm(Select(tab1, "y") ~ x2 + x3.tij,
               xij = list(Select(tab1, "t", as.formula = TRUE,
                                 sort = FALSE, lhs = "x3.tij", rhs = "0")),
               posbernoulli.t(parallel.t = TRUE ~ x2 + x3.tij),
               data = tab1, trace = TRUE,
               form2 = Select(tab1, prefix = TRUE, as.formula = TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}



% 20140524; For Fit.th before prefix = TRUE was allowed:
%              form2 = Select(tab1, "t", as.formula = TRUE,
%                             rhs = "x2 + x3.tij"))


% dim(subset(prinia, TRUE, select = grepl("^y", colnames(prinia))))

