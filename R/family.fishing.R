# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.





DeLury <- function(catch, effort,
                   type = c("DeLury", "Leslie"),
                   ricker = FALSE) {
  type <- match.arg(type, c("DeLury", "Leslie"))[1]
  if (!is.logical(ricker))
    stop("bad input for argument 'ricker'")
  if ((LLL <- Lcatch <- length(catch)) != (Leffort <- length(effort)))
    stop("length(catch) != length(effort)")

  CPUE <- catch / effort
  if (type == "DeLury") {
    Et <- cumsum(effort) - ifelse(ricker, 0.5, 1) * effort
    logCPUE <- log(CPUE)
    lmfit <- lm(logCPUE ~ Et, x = TRUE)
    myq <- catchabilityCoefficient <- -coef(lmfit)[2]
    N0 <- exp(coef(lmfit)["(Intercept)"]) / myq
  } else {
    Kt <- cumsum(catch) - ifelse(ricker, 0.5, 1) * catch
    lmfit <- lm(CPUE ~ Kt, x = TRUE)
    myq <- catchabilityCoefficient <- -coef(lmfit)[2]
    N0 <- coef(lmfit)["(Intercept)"] / myq
  }

  rlist <-
  list(catch = catch,
       effort = effort,
       type = type,
       N0 = N0,
       CPUE = CPUE,
       lmfit = lmfit)
  if (type == "DeLury") {
    rlist$E <- Et
  } else {
    rlist$K <- Kt
  }
  rlist
}






wffc.P1     <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
  ifelse(length >= min.eligible, c1 + (ppm/100) *
         ceiling(  signif(100 * length, digits = 8)  ), 0)


wffc.P1star <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
  ifelse(length >= min.eligible, c1 + ppm * length, 0)















wffc.P2     <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
  wffc.P1(length, c1 = c1, min.eligible = min.eligible, ppm = ppm) +
  ifelse(length >= min.eligible,
           ceiling(100*(length-min.eligible))^2, 0)

wffc.P2star <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
  wffc.P1star(length, c1 = c1, min.eligible = min.eligible, ppm = ppm) +
  ifelse(length >= min.eligible, 10000 * (length-min.eligible)^2, 0)





wffc.P3     <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000) {

  temp1 <- floor((ceiling(100*length)/100) / min.eligible) # zz not sure
  temp1 <- floor(length / min.eligible)
  ans <- ifelse(temp1 >= 1, c1, length * 0) # Handles NAs
  ans <- ans + ifelse(temp1 >= 1, ppm * (ceiling(100*length)/100), 0)
  maxtemp1 <- max(temp1, na.rm = TRUE)
  if (maxtemp1 > 1)
    for (ii in 2:maxtemp1) {
      ans <- ans +
            ifelse(ii <  temp1,         min.eligible  * (ii-1) * ppm, 0) +
            ifelse(ii == temp1, (ceiling(100*length)/100 -
                   ii*min.eligible) * (ii-1) * ppm, 0)
    }
  ans
}



wffc.P3star <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000) {
  temp1 <- floor(length / min.eligible)
  ans <- ifelse(temp1 >= 1, c1, length * 0) # Handles NAs
  ans <- ans + ifelse(temp1 >= 1, length * ppm, 0)
  maxtemp1 <- max(temp1, na.rm = TRUE)
  if (maxtemp1 > 1)
    for (ii in 2:maxtemp1) {
      ans <- ans + ifelse(ii <  temp1,  min.eligible  * (ii-1) * ppm, 0) +
                   ifelse(ii == temp1, (length - ii*min.eligible) *
                                       (ii-1) * ppm, 0)
    }
  ans
}







