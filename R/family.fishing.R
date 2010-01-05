# "family.fishing.q"
# Last modified: 01/12/08, 02/12/08

# These functions are Copyright (C) 2009-2010 T. W. Yee   All rights reserved.


# ====================================================================
# 20081201

DeLury = function(catch, effort,
                  type=c("DeLury","Leslie"),
                  ricker=FALSE) {
# 20081202; this function has been checked not ok
    type = match.arg(type, c("DeLury","Leslie"))[1]
    if (!is.logical(ricker)) stop("bad input for 'ricker'")
    if ((LLL <- Lcatch <- length(catch)) != (Leffort <- length(effort)))
        stop("length(catch) != length(effort)")

    CPUE = catch / effort
    if (type == "DeLury") {
        Et = cumsum(effort) - ifelse(ricker, 0.5, 1) * effort
        logCPUE = log(CPUE)
        lmfit = lm(logCPUE ~ Et, x=TRUE)
        myq = catchabilityCoefficient = -coef(lmfit)[2]
        N0 = exp(coef(lmfit)["(Intercept)"]) / myq
    } else {
        Kt = cumsum(catch) - ifelse(ricker, 0.5, 1) * catch
        lmfit = lm(CPUE ~ Kt, x=TRUE)
        myq = catchabilityCoefficient = -coef(lmfit)[2]
        N0 = coef(lmfit)["(Intercept)"] / myq
    }

    rlist =
    list(catch=catch,
         effort=effort,
         type=type,
         N0 = N0,
         CPUE = CPUE,
         lmfit=lmfit)
    if (type == "DeLury") {
        rlist$E = Et
    } else {
        rlist$K = Kt
    }
    rlist
}



# ======================================================================
# 20081201
# Transferred over from my own files and then modified here.

# length is in metres
wffc.P1     = function(length, min.eligible=0.18)
    ifelse(length >= min.eligible, 100 + 20 * ceiling(100*length), 0)
wffc.P1star = function(length, min.eligible=0.18)
    ifelse(length >= min.eligible, 100 + 2000 * length, 0)

# This was in the original mss. but problem is P2 does not return an integer
#wffc.P2     = function(y, min.eligible=0.18)
#    P1(y) + ifelse(y >= min.eligible, 0.7*ceiling(100*(y-min.eligible))^2, 0)
#wffc.P2star = function(y, min.eligible=0.18)
#    P1star(y) + ifelse(y >= min.eligible, 7000 * (y-min.eligible)^2, 0)

# 7/6/08; This returns an integer
wffc.P2     = function(length, min.eligible=0.18)
    wffc.P1(length) +
    ifelse(length >= min.eligible,
           ceiling(100*(length-min.eligible))^2, 0)
wffc.P2star = function(length, min.eligible=0.18)
    wffc.P1star(length) +
    ifelse(length >= min.eligible, 10000 * (length-min.eligible)^2, 0)


# ======================================================================



