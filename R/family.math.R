# These functions are
# Copyright (C) 1998-2014 T.W. Yee, University of Auckland.
# All rights reserved.











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
  if (length(q)     != nnn) q     <- rep(q,     length = nnn)
  if (length(shape) != nnn) shape <- rep(shape, length = nnn)

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





