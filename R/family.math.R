# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.











if (FALSE)
log1pexp <- function(x) {

  ans <- log1p(exp(x))
  big <- (x > 10)
  ans[big] <- x[big] + log1p(exp(-x[big]))
  ans
}







erf <- function(x, inverse = FALSE) {
  if (inverse) {
    ans <- qnorm((x+1)/2) / sqrt(2)
    ans[x <  -1] <- NA
    ans[x >  +1] <- NA
    ans[x == -1] <- -Inf
    ans[x == +1] <-  Inf
    ans
  } else {
    2 * pnorm(x * sqrt(2)) - 1
  }
}



erfc <- function(x, inverse = FALSE) {
  if (inverse) {
    ans <- qnorm(x/2, lower.tail = FALSE) / sqrt(2)
    ans[x <  0] <- NA
    ans[x >  2] <- NA
    ans[x == 0] <-  Inf
    ans[x == 2] <- -Inf
    ans
  } else {
    2 * pnorm(x * sqrt(2), lower.tail = FALSE)
  }
}







lambertW <- function(x, tolerance = 1.0e-10, maxit = 50) {
  if (any(Im(x) != 0.0))
    stop("argument 'x' must be real, not complex!")

  ans <- x
  ans[!is.na(x) & x <  -exp(-1)] <- NA
  ans[!is.na(x) & x >= -exp(-1)] <- log1p(x[!is.na(x) & x >= -exp(-1)])
  ans[!is.na(x) & x >= 0       ] <-  sqrt(x[!is.na(x) & x >= 0       ]) / 2

  cutpt <- 3.0
  if (any(myTF <- !is.na(x) & x > cutpt)) {
    L1 <- log(x[!is.na(x) & x > cutpt])  # log(as.complex(x))
    L2 <- log(L1)  # log(as.complex(L1))
    wzinit <- L1 - L2 +
          (L2 +
          (L2*( -2 + L2)/(2) +
          (L2*(  6 + L2*(-9 + L2*   2)) / (6) +
           L2*(-12 + L2*(36 + L2*(-22 + L2*3))) / (12*L1)) / L1) / L1) / L1

    ans[myTF] <- wzinit
  }

  for (ii in 1:maxit) {
    exp1 <- exp(ans)
    exp2 <- ans * exp1
    delta <- (exp2 - x) / (exp2 + exp1 -
                ((ans + 2) * (exp2 - x) / (2 * (ans + 1.0))))
    ans <- ans - delta
    if (all(is.na(delta) ||
        max(abs(delta), na.rm = TRUE) < tolerance)) break
    if (ii == maxit)
      warning("did not converge")
  }
  ans[x == Inf] <- Inf
  ans
}






 pgamma.deriv <- function(q, shape, tmax = 100) {

  nnn <- max(length(q), length(shape))
  if (length(q)     != nnn) q     <- rep_len(q,     nnn)
  if (length(shape) != nnn) shape <- rep_len(shape, nnn)

  if (!is.Numeric(q, positive = TRUE))
    stop("bad input for argument 'q'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")

  if (!is.Numeric(tmax, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'tmax'")
  if (tmax < 10)
    warning("probably argument 'tmax' is too small")


  gplog  <- lgamma(shape)
  gp1log <- gplog + log(shape)
  psip   <- digamma(shape)
  psip1  <- psip + 1 / shape
  psidp  <- trigamma(shape)
  psidp1 <- psidp - 1 / shape^2

  fred <-
    .C("VGAM_C_vdigami",
         d = as.double(matrix(0, 6, nnn)),
         x = as.double(q), p = as.double(shape),
         as.double(gplog), as.double(gp1log), as.double(psip),
         as.double(psip1), as.double(psidp), as.double(psidp1),
         ifault = integer(nnn),
         tmax = as.double(tmax),
         as.integer(nnn))
  answer <- matrix(fred$d, nnn, 6, byrow = TRUE)
  dimnames(answer) <- list(names(q),
                           c("q", "q^2", "shape", "shape^2",
                             "q.shape", "pgamma(q, shape)"))

  if (any(fred$ifault != 0)) {
    indices <- which(fred$ifault != 0)
    warning("convergence problems with elements ",
             indices)
  }

  answer
}








expint <- function (x, deriv = 0) {
  if (deriv == 0) {
    LLL <- length(x)
    answer <- .C("sf_C_expint", x = as.double(x), size = as.integer(LLL),
                 ans = double(LLL))$ans
    answer[x < 0] <- NA
    answer[x == 0] <- NA
    answer
  } else {
    if (!is.Numeric(deriv, integer.valued = TRUE, positive = TRUE) ||
        deriv > 3)
      stop("Bad input for argument 'deriv'")
    answer <- rep_len(0, length(x))
    if (deriv == 1) {
      answer <- exp(x) / x
    }
    if (deriv == 2) {
      answer <- exp(x) / x - exp(x) / x^2
    }
    if (deriv == 3) {
      answer <- exp(x) / x - 2 * exp(x) / x^2 +
        2 * exp(x) / x^3
    }
    answer
  }
}


expexpint <- function (x, deriv = 0) {
  LLL <- length(x)
  answer <- .C("sf_C_expexpint", x = as.double(x), size = as.integer(LLL),
               ans = double(LLL))$ans
  answer[x <  0] <- NA
  answer[x == 0] <- NA
  if (deriv > 0) {
    if (!is.Numeric(deriv, integer.valued = TRUE, positive = TRUE) ||
        deriv > 3)
      stop("Bad input for argument 'deriv'")
    if (deriv >= 1) {
      answer <- -answer + 1 / x
    }
    if (deriv >= 2) {
      answer <- -answer - 1 / x^2
    }
    if (deriv == 3) {
      answer <- -answer + 2 / x^3
    }
  }
  answer
}


expint.E1 <- function (x, deriv = 0) {
  if (deriv == 0) {
    LLL <- length(x)
    answer <- .C("sf_C_expint_e1", x = as.double(x), size = as.integer(LLL),
                 ans = double(LLL))$ans
    answer[x < 0] <- NA
    answer[x == 0] <- NA
  } else {
    if (!is.Numeric(deriv, integer.valued = TRUE, positive = TRUE) ||
        deriv > 3)
      stop("Bad input for argument 'deriv'")
    answer <- rep_len(0, length(x))
    if (deriv == 1) {
      answer <- exp(-x) / x
    }
    if (deriv == 2) {
      answer <- exp(-x) / x + exp(-x) / x^2
    }
    if (deriv == 3) {
      answer <- exp(-x) / x + 2 * exp(-x) / x^2 +
        2 * exp(-x) / x^3
    }
    answer <- (-1)^deriv * answer
  }
  answer
}











if (FALSE)
expint <- function(x) {


  LLL <- length(x)
  answer <- .C("sf_C_expint",
                 x = as.double(x),
                 size = as.integer(LLL),
                 ans = double(LLL))$ans

  answer[x  < 0] <- NA
  answer[x == 0] <- NA

  answer
}



if (FALSE)
expexpint <- function(x) {




  LLL <- length(x)
  answer <- .C("sf_C_expexpint",
                 x = as.double(x),
                 size = as.integer(LLL),
                 ans = double(LLL))$ans

  answer[x  < 0] <- NA
  answer[x == 0] <- NA

  answer
}


if (FALSE)
pochhammer <- function (x, n) {
  exp(lgamma(x+n) - lgamma(x))
}







if (FALSE)
expint.E1 <- function(x) {




  LLL <- length(x)
  answer <- .C("sf_C_expint_e1",
                 x = as.double(x),
                 size = as.integer(LLL),
                 ans = double(LLL))$ans

  answer[x  < 0] <- NA
  answer[x == 0] <- NA

  answer
}




 Zeta.aux <- function(shape, qq, shift = 1) {



  LLL <- max(length(shape), length(qq))
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(qq   ) != LLL) qq    <- rep_len(qq,    LLL)

  if (any(qq < 12-1))
    warning("all values of argument 'q' should be 12 or more")
  aa <- qq


  B2 <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  kk <- length(B2)  # 8
  ans <- 1 / ((shape-1) * (shift + aa)^(shape-1)) +
         0.5 / (shift + aa)^shape

  term <- (shape/2) / (shift + aa)^(shape+1)
  ans <- ans + term * B2[1]

  for (mm in 2:kk) {
    term <- term * (shape+2*mm-3) *
            (shape+2*mm-2) / ((2*mm-1) * 2 * mm * (shift + aa)^2)
    ans <- ans + term * B2[mm]
  }
  ifelse(aa - 1 <= qq, ans, rep(0, length(ans)))  # Handled above
}







 zeta <- function(x, deriv = 0, shift = 1) {




  deriv.arg <- deriv
  rm(deriv)
  if (!is.Numeric(deriv.arg, length.arg = 1,
                  integer.valued = TRUE))
    stop("'deriv' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv' must be 0, 1, or 2")


  if (deriv.arg > 0)
    return(Zeta.derivative(x, deriv.arg = deriv.arg, shift = shift))



  if (any(special <- Re(x) <= 1)) {
    ans <- x
    ans[special] <- Inf  # For Re(x) == 1

    special3 <- Re(x) < 1
    ans[special3] <- NA  # For 0 < Re(x) < 1

    special4 <- (0 < Re(x)) & (Re(x) < 1) & (Im(x) == 0)
    ans[special4] <-
      Zeta.derivative(x[special4], deriv.arg = deriv.arg, shift = shift)


    special2 <- Re(x) < 0
    if (any(special2)) {
      x2 <- x[special2]
      cx <- 1 - x2
      ans[special2] <- 2^(x2) * pi^(x2-1) * sin(pi*x2/2) *
                       gamma(cx) * Recall(cx)
    }  # special2

    if (any(!special)) {
      ans[!special] <- Recall(x[!special])
    }
    return(ans)
  }  # special

  aa <- 12
  ans <- 0
  for (ii in 0:(aa-1))
    ans <- ans + 1 / (shift + ii)^x

  ans <- ans + Zeta.aux(shape = x, aa, shift = shift)
  ans[shift <= 0] <- NaN
  ans
}  # zeta



 Zeta.derivative <- function(x, deriv.arg = 0, shift = 1) {



  if (!all(shift == 1))
    stop("currently 'shift' must all be 1")


  if (!is.Numeric(deriv.arg, length.arg = 1,
                  integer.valued = TRUE))
    stop("'deriv.arg' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv.arg' must be 0, 1, or 2")

  if (any(Im(x) != 0))
    stop("Sorry, currently can only handle x real, not complex")
  if (any(x < 0))
    stop("Sorry, currently cannot handle x < 0")

  ok <- is.finite(x) & x > 0 & x != 1  # Handles NAs
  ans <- rep_len(NA_real_, length(x))
  nn <- sum(ok)  # Effective length (excludes x < 0 and x = 1 values)
  if (nn)
    ans[ok] <- .C("vzetawr", as.double(x[ok]), ans = double(nn),
                  as.integer(deriv.arg), as.integer(nn))$ans



  if (deriv.arg == 0)
    ans[is.finite(x) & abs(x) < 1.0e-12] <- -0.5

  ans
}






