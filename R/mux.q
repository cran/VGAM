# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.




mux34 <- function(xmat, cc, symmetric = FALSE) {


  if (!is.matrix(xmat))
    xmat <- as.matrix(xmat)
  d <- dim(xmat)
  nnn <- d[1]
  RRR <- d[2]
  if (length(cc) == 1)
    cc <- matrix(cc, 1, 1)
  if (!is.matrix(cc))
    stop("'cc' is not a matrix")
  c( .C("VGAM_C_mux34", as.double(xmat), as.double(cc),
         as.integer(nnn), as.integer(RRR),
         as.integer(symmetric), ans = as.double(rep_len(0.0, nnn)),
         NAOK = TRUE)$ans)
}




if (FALSE)
mux34 <- function(xmat, cc, symmetric = FALSE) {

  if (!is.matrix(xmat))
    xmat <- as.matrix(xmat)
  d <- dim(xmat)
  n <- d[1]
  R <- d[2]
  if (length(cc) == 1)
    cc <- matrix(cc, 1, 1)
  if (!is.matrix(cc))
    stop("'cc' is not a matrix")
  c( .Fortran("vgamf90mux34", as.double(xmat), as.double(cc),
               as.integer(n), as.integer(R),
               as.integer(symmetric), ans = as.double(rep_len(0.0, n)),
               NAOK = TRUE)$ans)
}





mux2 <- function(cc, xmat) {


  if (!is.matrix(xmat))
    xmat <- as.matrix(xmat)
  d <- dim(xmat)
  n <- d[1]
  p <- d[2]
  if (is.matrix(cc))
    cc <- array(cc, c(dim(cc), n))
  d <- dim(cc)
  M <- d[1]
  if (d[2] != p || d[3] != n)
    stop("dimension size inconformable")
  ans <- rep_len(NA_real_, n*M)
  fred <- .C("mux2", as.double(cc), as.double(t(xmat)),
               ans = as.double(ans), as.integer(p), as.integer(n),
               as.integer(M), NAOK = TRUE)
  matrix(fred$ans, n, M, byrow = TRUE)
}






mux22 <- function(cc, xmat, M, upper = FALSE, as.matrix = FALSE) {

  n <- ncol(cc)

  index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
  dimm.value <- nrow(cc)  # Usually M or M(M+1)/2

  ans <- rep_len(NA_real_, n*M)
  fred <- .C("mux22", as.double(cc), as.double(t(xmat)),
               ans = as.double(ans), as.integer(dimm.value),
               as.integer(index$row), as.integer(index$col),
               as.integer(n), as.integer(M), wk = double(M*M),
               as.integer(as.numeric(upper)), NAOK = TRUE)
  if (!as.matrix) fred$ans else {
    dim(fred$ans) <- c(M, n)
    t(fred$ans)
  }
}





mux5 <- function(cc, x, M, matrix.arg = FALSE) {



  dimx <- dim(x)
  dimcc <- dim(cc)
  r <- dimx[2]

  if (matrix.arg) {
    n <- dimcc[1]
    neltscci <- ncol(cc)
    cc <- t(cc)
  } else {
    n <- dimcc[3]
    if (dimcc[1] != dimcc[2] ||
        dimx[1]  != dimcc[1] ||
        (length(dimx) == 3 && dimx[3] != dimcc[3]))
      stop('input nonconformable')
    neltscci <- M*(M+1)/2
  }

  if (is.matrix(x))
    x <- array(x,c(M, r, n))
  index.M <- iam(NA, NA, M, both = TRUE, diag = TRUE)
  index.r <- iam(NA, NA, r, both = TRUE, diag = TRUE)

  size <- if (matrix.arg) dimm(r)*n else r*r*n
  fred <- .C("mux5", as.double(cc), as.double(x),
               ans = double(size),
               as.integer(M), as.integer(n), as.integer(r),
               as.integer(neltscci),
               as.integer(dimm(r)),
               as.integer(as.numeric(matrix.arg)),
               double(M*M), double(r*r),
               as.integer(index.M$row), as.integer(index.M$col),
               as.integer(index.r$row), as.integer(index.r$col),
               ok3 = as.integer(1), NAOK = TRUE)
  if (fred$ok3 == 0)
    stop("can only handle matrix.arg == 1")


  if (matrix.arg) {
    ans <- fred$ans
    dim(ans) <- c(dimm(r), n)
    t(ans)
  } else {
    array(fred$ans, c(r, r, n))
  }
}




mux55 <- function(evects, evals, M) {

  d <- dim(evects)
  n <- ncol(evals)
  if (d[1] != M || d[2] != M || d[3] != n ||
      nrow(evals)!= M || ncol(evals) != n)
    stop("input nonconformable")
  MMp1d2 <- M*(M+1)/2  # The answer is a full-matrix
  index <- iam(NA, NA, M, both = TRUE, diag = TRUE)

  fred <- .C("mux55", as.double(evects), as.double(evals),
             ans = double(MMp1d2 * n),
             double(M*M), double(M*M),
             as.integer(index$row), as.integer(index$col),
             as.integer(M), as.integer(n), NAOK = TRUE)
  dim(fred$ans) <- c(MMp1d2, n)
  fred$ans
}




mux7 <- function(cc, x) {

  dimx <- dim(x)
  dimcc <- dim(cc)
  if (dimx[1]!= dimcc[2] ||
     (length(dimx) == 3 && dimx[3]!= dimcc[3]))
    stop('input nonconformable')
  M  <- dimcc[1]
  qq <- dimcc[2]
  n  <- dimcc[3]
  r <- dimx[2]
  if (is.matrix(x))
    x <- array(x, c(qq, r, n))

  ans <- array(NA, c(M, r, n))
  fred <- .C("mux7", as.double(cc), as.double(x),
               ans = as.double(ans),
               as.integer(M), as.integer(qq), as.integer(n),
               as.integer(r), NAOK = TRUE)
  array(fred$ans, c(M, r, n))
}





mux9 <- function(cc, xmat) {

  if (is.vector(xmat))
    xmat <- cbind(xmat)
  dimxmat <- dim(xmat)
  dimcc <- dim(cc)

  if (dimcc[1]   != dimcc[2] ||
      dimxmat[1] != dimcc[3] ||
      dimxmat[2] != dimcc[1])
    stop('input nonconformable')
  M <- dimcc[1]
  n <- dimcc[3]

  ans <-  matrix(NA_real_, n, M)
  fred <- .C("mux9", as.double(cc), as.double(xmat),
               ans = as.double(ans),
               as.integer(M), as.integer(n), NAOK = TRUE)
  matrix(fred$ans, n, M)
}




mux11 <- function(cc, xmat) {


  dcc <- dim(cc)
  d <- dim(xmat)
  M <- dcc[1]
  R <- d[2]
  n <- dcc[3]
  if (M != dcc[2] || d[1] != n*M)
    stop("input inconformable")

  Xmat <- array(c(t(xmat)), c(R, M, n))
  Xmat <- aperm(Xmat, c(2, 1, 3))   # Xmat becomes M x R x n
  mat <- mux7(cc, Xmat)             # mat is M x R x n
  mat <- aperm(mat, c(2, 1, 3))     # mat becomes R x M x n
  mat <- matrix(c(mat), n*M, R, byrow = TRUE)
  mat
}



mux111 <- function(cc, xmat, M, upper = TRUE) {


  R <- ncol(xmat)
  n <- nrow(xmat) / M
  index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
  dimm.value <- nrow(cc)  # M or M(M+1)/2

  fred <- .C("mux111", as.double(cc),
               b = as.double(t(xmat)),
               as.integer(M),
               as.integer(R), as.integer(n), wk = double(M * M),
               wk2 = double(M * R), as.integer(index$row),
               as.integer(index$col), as.integer(dimm.value),
               as.integer(as.numeric(upper)), NAOK = TRUE)

  ans <- fred$b
  dim(ans) <- c(R, nrow(xmat))
  d <- dimnames(xmat)
  dimnames(ans) <- list(d[[2]], d[[1]])
  t(ans)
}






mux15 <- function(cc, xmat) {

  n <- nrow(xmat)
  M <- ncol(xmat)
  if (nrow(cc) != M || ncol(cc) != M)
    stop("input inconformable")
  if (max(abs(t(cc)-cc))>0.000001)
    stop("argument 'cc' is not symmetric")

  ans <- rep_len(NA_real_, n*M*M)
  fred <- .C("mux15", as.double(cc), as.double(t(xmat)),
               ans = as.double(ans), as.integer(M),
               as.integer(n), NAOK = TRUE)
  array(fred$ans, c(M, M, n))
}





vforsub <- function(cc, b, M, n) {



    index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
    dimm.value <- nrow(cc)  # M or M(M+1)/2


  fred <- .C("vforsub", as.double(cc), b = as.double(t(b)),
               as.integer(M), as.integer(n), wk = double(M*M),
               as.integer(index$row), as.integer(index$col),
               as.integer(dimm.value), NAOK = TRUE)

  dim(fred$b) <- c(M, n)
  fred$b
}




vbacksub <- function(cc, b, M, n) {
  index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
  dimm.value <- nrow(cc)
  if (nrow(b) != M || ncol(b) != n)
    stop("dimension size inconformable")

  fred <- .C("vbacksub", as.double(cc), b = as.double(b),
               as.integer(M), as.integer(n), wk = double(M*M),
               as.integer(index$row), as.integer(index$col),
               as.integer(dimm.value), NAOK = TRUE)

  if (M == 1) {
    fred$b
  } else {
    dim(fred$b) <- c(M, n)
    t(fred$b)
  }
}



vchol <- function(cc, M, n, silent = FALSE, callno = 0) {





  index <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
  cc <- t(cc)
  MM <- nrow(cc)    # cc is big enough to hold its Cholesky decom.

  fred <- .C("vchol", cc = as.double(cc), as.integer(M),
               as.integer(n), ok = integer(n),
               wk = double(M*M), as.integer(index$row),
               as.integer(index$col),
               as.integer(MM),
               NAOK = TRUE)

  failed <- (fred$ok != 1)
  if ((correction.needed <- any(failed))) {
    index <- (1:n)[failed]
    if (!silent) {
      if (length(index) < 11)
        warning("weight matri",
                ifelse(length(index) > 1, "ces ","x "),
                paste(index, collapse = ", "),
                " not positive-definite")
    }
  }

  ans <- fred$cc
  dim(ans) <- c(MM, n)

  if (correction.needed) {
    temp <- cc[, index, drop = FALSE]
    tmp777 <- vchol.greenstadt(temp, M = M, silent = silent,
                               callno = callno + 1)


    if (length(index) == n) {
      ans <- tmp777[1:nrow(ans), , drop = FALSE]
    } else {


      ans[, index] <- tmp777  # restored 20031016
    }
  }
  dim(ans) <- c(MM, n)  # Make sure

  ans
}



vchol.greenstadt <- function(cc, M, silent = FALSE,
                             callno = 0) {




  MM <- dim(cc)[1]
  n <- dim(cc)[2]

  if (!silent)
    cat(paste("Applying Greenstadt modification to ", n, " matri",
              ifelse(n > 1, "ces", "x"), "\n", sep = ""))





  temp <- veigen(cc, M = M)  # , mat = TRUE)
  dim(temp$vectors) <- c(M, M, n)  # Make sure (when M = 1) for mux5
  dim(temp$values)  <- c(M, n)    # Make sure (when M = 1) for mux5

  is.neg <- (temp$values < .Machine$double.eps)
  is.pos <- (temp$values > .Machine$double.eps)
  zilch  <- (!is.pos & !is.neg)

  temp$values <- abs(temp$values)

  temp.small.value <- quantile(temp$values[!zilch], prob = 0.15)
  if (callno > 2) {
    temp.small.value <- abs(temp.small.value) * 1.50^callno


      small.value <- temp.small.value


      temp$values[zilch] <- small.value
  }

  if (callno > 9) {
    warning("taking drastic action; setting all wz to ",
            "scaled versions of the order-M identity matrix")

    cc2mean <- abs(colMeans(cc[1:M, , drop = FALSE]))
    temp$values  <- matrix(cc2mean, M, n, byrow = TRUE)
    temp$vectors <- array(c(diag(M)), c(M, M, n))
  }



  temp3 <- mux55(temp$vectors, temp$values, M = M)  #, matrix.arg = TRUE)
  ans <- vchol(t(temp3), M = M, n = n, silent = silent,
               callno = callno + 1)  #, matrix.arg = TRUE)



  if (nrow(ans) == MM) ans else ans[1:MM, , drop = FALSE]
}




if (FALSE)
myf <- function(x) {
    .Fortran("VGAM_F90_fill9",
               x = as.double(x), lenx = as.integer(length(x)),
               answer = as.double(x),
               NAOK = TRUE)$answer
}




