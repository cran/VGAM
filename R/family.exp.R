# These functions are Copyright (C) 1998-2010 T. W. Yee  All rights reserved.

# Families for expectile regression are put in this file 
# 20100324;
# Last modified: 20100324, 20100326, 20100329, 20100331,

# Yet to do:
# 1. lms.bcn(expectiles = FALSE). If lms.bcn(expectiles = TRUE) then
#    expectiles, and not quantiles, are the fitted values.
#    This is LMS-BCN expectile regression, a new method.
# 2. Improve the approximations (initial values) for each of the
#    three distributions. See the zzs below.
# 3. For peunif(q) etc.: use e or q as first argument??
#    For deunif(x) etc.: use e or x as first argument??
#    For qeunif(x) etc.: rename to eeunif(x)?

# Done:
# 1. For norm, exp and unif distributions:
#    qenorm(0.25) returns the 0.25-expectile of a standard normal,
#    penorm(1.25) returns the tau (in (0,1)) for an expectile of 1.25.
#    This is based on the paper by M C Jones (1994) in Stat Prob Letters.

# Notes:
# 1. 

# ======================================================================
# Expectiles for uniform distribution ----------------------------------
# 20100324
# The [et]norm() here were adapted from MC Jones paper.

qeunif <- function(p, min = 0, max = 1, Maxit_nr = 10, Tol_nr = 1.0e-6) {
# Using Newton-Raphson may be a problem at the boundaries.
# The secant method may be better.

  ppp = p
  vsmallno = sqrt(.Machine$double.eps)
   smallno = 0.10
  if (any(min >= max))
    stop("argument 'min' has values greater or equal to argument 'max'")
  if (!is.Numeric( Tol_nr, allow = 1, posit = TRUE) || Tol_nr > 0.10)
    stop("argument 'Tol_nr' is not a single positive value, or is too large")
  nrok = ppp >= vsmallno & ppp <= 1.0 - vsmallno & is.finite(ppp)

# A beta function seems to approximate it ok near the middle.
# This can be improved zz.
  eee = qbeta(ppp, shape1 = 3, shape2 = 3)
# A different quadratic fits each boundary well (asymptotic expansion).
  eee[ppp <        smallno] = sqrt(ppp[ppp <  smallno])
  eee[ppp > 1.0 -  smallno] = 1.0 - sqrt(1.0 - ppp[ppp > 1.0 -  smallno])

#lines(ppp, eee, col="purple", type="b")
#print("initial eee")
#isample = sample(length(eee))
#isample = 1:length(eee)
#print( head(eee[isample]) )
#print(     (eee[isample]) )
#cat("\n")

  for(iii in 1:Maxit_nr) {
    realdiff <- (peunif(eee[nrok]) - ppp[nrok]) / deunif(eee[nrok])
#  #print("max(abs(realdiff))")
#  #print( max(abs(realdiff)) )
    eee[nrok] = eee[nrok] - realdiff
#   cat("Iteration ", iii, "\n")
#  #print( head(eee[isample]) )
#  #print(     (eee[isample]) )
#   cat("\n")
    if (all(abs(realdiff) / (1.0 + abs(realdiff)) < Tol_nr )) break
    if (iii == Maxit_nr) warning("did not converge")
  }

# Check again (on the standard uniform distribution);
  if (max(abs(peunif(eee[nrok]) - ppp[nrok])) > Tol_nr)
    warning("did not converge on the second check")

# zz; Needs checking, esp. near the boundary of 1.0:
  eee[ppp <       vsmallno] =       sqrt(      ppp[ppp <       vsmallno])
  eee[ppp > 1.0 - vsmallno] = 1.0 - sqrt(1.0 - ppp[ppp > 1.0 - vsmallno])
  eee[ppp == 0] = 0
  eee[ppp == 1] = 1
  eee[ppp <  0] = NA
  eee[ppp >  1] = NA
  min + eee * (max - min)
}


peunif <- function(q, min = 0, max = 1, log = FALSE) {
# zz use e or x ??
# This is G(y).
  if (!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
  rm(log)
  if (any(min >= max))
    stop("argument 'min' has values greater or equal to argument 'max'")

  eee = (q - min) / (max - min)
  if (log.arg) {
    logGofy = 2 * log(eee) - log1p(2 * eee * (eee - 1))
    logGofy[eee < 0] = -Inf
    logGofy[eee > 1] = 0.0
    logGofy
  } else {
    Gofy = eee^2 / (1 + 2 * eee * (eee - 1))
    Gofy[eee < 0] = 0.0
    Gofy[eee > 1] = 1.0
    Gofy
  }
}



deunif <- function(x, min = 0, max = 1, log = FALSE) {
# This is g(x).
  if (!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
  rm(log)
  if (any(min >= max))
    stop("argument 'min' has values greater or equal to argument 'max'")

  eee = (x - min) / (max - min)

  if (log.arg) {
    ans = log(2) + log(eee) + log1p(-eee) -
          2.0 * log(2*eee*(1-eee) - 1) - log(max - min)
    ans[eee <= 0.0] = log(0.0)
    ans[eee >= 1.0] = log(0.0)
  } else {
    gunif <- function(y)
        as.numeric(y >= 0 & y <= 1) * 2*y*(1-y) / (2*y*(1-y) - 1)^2
    ans = gunif(eee) / (max - min)
#   ans[eee <  0.0] = 0.0
#   ans[eee >  1.0] = 0.0
  }
  ans
}




reunif <- function(n, min = 0, max = 1) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ=TRUE, allow=1, posit=TRUE))
                stop("bad input for argument 'n'") else n
    qeunif(runif(use.n), min = min, max = max)
}




# ======================================================================
# Expectiles for normal distribution -----------------------------------
# 20100324
# The [et]norm() here were adapted from MC Jones paper.

qenorm <- function(p, mean = 0, sd = 1, Maxit_nr = 10, Tol_nr = 1.0e-6) {
  ppp = p
  if (!is.Numeric( Tol_nr, allow = 1, posit = TRUE) || Tol_nr > 0.10)
    stop("argument 'Tol_nr' is not a single positive value, or is too large")
# if (!is.Numeric( sd, posit = TRUE))
#   stop("argument 'sd' must contain positive values")
  nrok = is.finite(ppp)

# A N(0, sd = 2/3) approximation is good according to the paper.
  eee =  qnorm(ppp, sd = 2/3)

# lines(ppp, eee, col="purple", type="b")
##print("initial eee")
#isample = sample(length(eee))
#isample = 1:length(eee)
##print( head(eee[isample]) )
##print(     (eee[isample]) )
# cat("\n")

  gnorm = function(y) dnorm(y) / (y * (1-2*pnorm(y)) - 2*dnorm(y))^2

  for(iii in 1:Maxit_nr) {
    realdiff <- (penorm(eee[nrok]) - ppp[nrok]) / gnorm(eee[nrok])
#  #print("max(abs(realdiff))")
#  #print( max(abs(realdiff)) )
    eee[nrok] = eee[nrok] - realdiff
#   cat("Iteration ", iii, "\n")
#  #print( head(eee[isample]) )
#  #print(     (eee[isample]) )
#   cat("\n")
    if (all(abs(realdiff) / (1.0 + abs(realdiff)) < Tol_nr )) break
    if (iii == Maxit_nr) warning("did not converge")
  }

# Check again (on the standard normal distribution);
  if (max(abs(penorm(eee[nrok]) - ppp[nrok])) > Tol_nr)
    warning("did not converge on the second check")

# zz; Needs checking, esp. near the boundary of 1.0:
  eee[ppp == 0] = -Inf
  eee[ppp == 1] =  Inf
  eee[ppp <  0] = NA
  eee[ppp >  1] = NA
  eee * ifelse(sd >= 0, sd, NaN) + mean
}


penorm <- function(q, mean = 0, sd = 1, log = FALSE) {
  if (!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
  rm(log)

  eee = (q - mean) / sd
  tmp1 = -dnorm(eee) - eee * pnorm(eee)
  if (log.arg) {
    logGofy = log(tmp1) - log(2 * tmp1 + eee)
    logGofy[eee <= -Inf] = -Inf
    logGofy[eee >=  Inf] = 0.0
    logGofy
  } else {
    Gofy = tmp1 / (2 * tmp1 + eee)
    Gofy[eee <= -Inf] = 0.0
    Gofy[eee >=  Inf] = 1.0
    Gofy
  }
}


denorm <- function(x, mean = 0, sd = 1, log = FALSE) {
# This is g(x).
  if (!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
  rm(log)

  eee = (x - mean) / sd
  if (log.arg) {
    ans = dnorm(eee, log = TRUE) -
          2.0 * log(eee * (1-2*pnorm(eee)) - 2*dnorm(eee)) - log(sd)
  } else {
    gnorm = function(y) dnorm(y) / (y * (1-2*pnorm(y)) - 2*dnorm(y))^2
    ans = gnorm(eee) / sd
    ans[sd  <=  0.0] = NaN
  }
  ans
}




renorm <- function(n, mean = 0, sd = 1) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ=TRUE, allow=1, posit=TRUE))
                stop("bad input for argument 'n'") else n
    qenorm(runif(use.n), mean = mean, sd = sd)
}





# ======================================================================
# Expectiles for exponential distribution ------------------------------
# 20100324
# The [et]exp() here were adapted from MC Jones paper.


qeexp <- function(p, rate = 1, Maxit_nr = 10, Tol_nr = 1.0e-6) {
  ppp = p
  vsmallno = sqrt(.Machine$double.eps)
  if (!is.Numeric( Tol_nr, allow = 1, posit = TRUE) || Tol_nr > 0.10)
    stop("argument 'Tol_nr' is not a single positive value, or is too large")
  nrok = ppp >= vsmallno & is.finite(ppp)

# 20100401; An approximation: (zz improve this!!)
# eee =  qf(0.8 * ppp, df1 =  4.0, df2 = 44) * 1.5

# 20100408; This is a piecewise approximation, and looks ok.
  eee = qf(1.0 * ppp, df1 =  4.0, df2 = 44)
  if ( any(rangex <- ppp < 0.8) )
      eee[rangex] = qrayleigh(ppp[rangex], a =  0.8)


# A different quadratic fits each boundary well (asymptotic expansion). zz
  eee[ppp <       vsmallno] = sqrt(ppp[ppp < vsmallno])

#lines(ppp,eee,col="purple",type="b") # See what the initial values were like
##print("initial eee")
#isample = sample(length(eee))
#isample = 1:length(eee)
##print( head(eee[isample]) )
##print(     (eee[isample]) )
##cat("\n")

  for(iii in 1:Maxit_nr) {
    realdiff <- (peexp(eee[nrok]) - ppp[nrok]) / deexp(eee[nrok])
#  #print("max(abs(realdiff))")
#  #print( max(abs(realdiff)) )
    eee[nrok] = eee[nrok] - realdiff
#   cat("Iteration ", iii, "\n")
#  #print( head(eee[isample]) )
#  #print(     (eee[isample]) )
#   cat("\n")
    if (all(abs(realdiff) / (1.0 + abs(realdiff)) < Tol_nr )) break
    if (iii == Maxit_nr) warning("did not converge")
  }

# Check again (on the standard exponential distribution);
  if (max(abs(peexp(eee[nrok]) - ppp[nrok])) > Tol_nr)
    warning("did not converge on the second check")

# zz; Needs checking, esp. near the boundary of 1.0:
  eee[ppp < vsmallno] = sqrt(ppp[ppp < vsmallno])
  eee[ppp == 0] = 0
  eee[ppp == 1] = Inf
  eee[ppp <  0] = NaN
  eee[ppp >  1] = NaN
  eee / rate
}


peexp <- function(q, rate = 1, log = FALSE) {
  if (!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
  rm(log)

  eee = q * rate
  if (log.arg) {
    tmp1 = -expm1(-eee) - eee
#   logGofy = log(tmp1) - log(2 * tmp1 + eee - 1.0)
    logGofy = log1p(- eee - exp(-eee)) - log(2 * tmp1 + eee - 1.0)
    logGofy[eee <    0] = log(0.0)
    logGofy[eee >= Inf] = log(1.0)
    logGofy
  } else {
#   tmp1 = 1 - eee - exp(-eee)
    tmp1 = -expm1(-eee) - eee
    Gofy = tmp1 / (2 * tmp1 + eee - 1.0)
    Gofy[eee <    0] = 0.0
    Gofy[eee >= Inf] = 1.0
    Gofy
  }
}



deexp <- function(x, rate = 1, log = FALSE) {
# This is g(x).
  if (!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
  rm(log)
  if (any(rate <= 0))
    stop("argument 'rate' must have positive values")

  eee = x * rate

  if (log.arg) {
    ans = log(eee) - eee + 2.0 * log((1-y) - 2 * exp(-y)) + log(rate)
  } else {
    gexp = function(y)
      as.numeric(y >= 0) * y * exp(-y) / ((1-y) - 2 * exp(-y))^2
    ans = gexp(eee) * rate
    ans[rate <=  0.0] = NaN
  }
  ans
}



reexp <- function(n, rate = 1) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ=TRUE, allow=1, posit=TRUE))
                stop("bad input for argument 'n'") else n
    qeexp(runif(use.n), rate = rate)
}


# ======================================================================

# ======================================================================


