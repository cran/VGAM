# These functions are
# Copyright (C) 1998-2016 T.W. Yee, University of Auckland.
# All rights reserved.







ps <-
  function(x,
           ...,
           ps.intervals = NULL,
           lambda = 0, degree = 2, order = 2,
           ridge.adj = 1e-5, ridge.inv = 0.0001,
           spillover = 0.01, maxlambda = 1e4) {


  xs <- substitute(x)
  ans <- as.character(xs)
  x.index <- as.vector(x)




  x.orig <- x.index
  xdots <- list(...)
  uses.xij <- length(xdots) > 0
  if (uses.xij)
    x.index <- as.vector(c(x.index, unlist(xdots)))
  if (is.null(ps.intervals))
    ps.intervals <- ceiling(1.5 * log(length(unique(x.index))))



  number.knots <- ps.intervals + 2 * degree + 1
  xl <- min(x.index)
  xr <- max(x.index)



 if (smart.mode.is("read")) {
    smartlist  <- get.smart()
    xl <- smartlist$xl  # Overwrite its value
    xr <- smartlist$xr  # Overwrite its value
    ps.intervals <- smartlist$ps.intervals  # Ditto
    number.knots <- ps.intervals + 2 * degree + 1  # Redo
    lambda    <- smartlist$lambda
    degree    <- smartlist$degree
    order     <- smartlist$order
    ridge.adj <- smartlist$ridge.adj
    ridge.inv <- smartlist$ridge.inv
    spillover <- smartlist$spillover
    maxlambda <- smartlist$maxlambda
    maXX      <- smartlist$maXX
    Cmat      <- smartlist$Cmat
  } else {
    maXX      <- NULL
    Cmat      <- NULL
  }

  xmax <- xr + spillover * (xr - xl)
  xmin <- xl - spillover * (xr - xl)
  dx <- (xmax - xmin) / ps.intervals
  nx <- names(x.index)
  nax <- is.na(x.index)
  if (nas <- any(nax))
    x.index <- x[!nax]
  sorder <- degree + 1
  if (length(ps.intervals)) {
    nAknots <- ps.intervals - 1
    if (nAknots < 1) {
      nAknots <- 1
      warning("ps.intervals was too small; have used 2")
    }






    if (nAknots > 0) {
      Aknots <- seq(from = xmin - degree * dx,
                    to   = xmax + degree * dx, by = dx)
    } else {
      knots <- NULL
    }
  }
  basis <- splineDesign(Aknots, x.index, sorder, 0 * x.index)
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA_real_, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(1:nrow(basis), 1:n.col)
  if ((order - n.col + 1) > 0) {
    order <- n.col - 1
    warning("order was too large; have used ", n.col - 1)
  }
  if (any(lambda < 0)) {
    lambda[lambda < 0] <- 0
    warning("some lambda values are negative : have used lambda = ",
            paste(lambda, collapse = ", "))
  }

  if (any(lambda > maxlambda)) {
    lambda[lambda > maxlambda] <- maxlambda
    warning("some lambda values are > ", maxlambda, ": ",
                  "for stability have used lambda = ",
            paste(lambda, collapse = ", "))
  }
  aug <- if (order > 0) diff(diag(n.col), diff = order) else diag(n.col)

  
  Pen <- t(aug) %*% aug
  pen.aug <- (Pen + t(Pen))/2



  if (is.null(maXX))
    maXX <- mean(abs(t(basis) %*% basis))
 maS <- mean(abs(pen.aug))/maXX



  pen.aug <- pen.aug / maS
  kk <- ncol(basis)
  if (is.null(Cmat))
    Cmat <- matrix(colSums(basis), 1, kk)


  qrCt <- qr(t(Cmat))
  jay <- nrow(Cmat)  # 1
  XZ <- t(qr.qty(qrCt, t(basis))[(jay+1):kk, ])


  ZtSZ <- qr.qty(qrCt, t(qr.qty(qrCt, t(pen.aug))))[(jay+1):kk,
                                                    (jay+1):kk]


  basis <- XZ


  if (smart.mode.is("write"))
    put.smart(list(xl        = xl,
                   xr        = xr,
                   ps.intervals = ps.intervals,
                   lambda    = lambda,
                   degree    = degree,
                   order     = order,
                   ridge.adj = ridge.adj,
                   ridge.inv = ridge.inv,
                   spillover = spillover,
                   maxlambda = maxlambda,
                   maXX      = maXX,
                   Cmat      = Cmat))




  basis <- basis[seq_along(x.orig), , drop = FALSE]




  attr(basis, "S.arg") <- ZtSZ

  attr(basis, "degree") <- degree
  attr(basis, "knots") <- Aknots
  attr(basis, "lambda") <- lambda  # Vector
  attr(basis, "order") <- order
  attr(basis, "ps.intervals") <- ps.intervals
  attr(basis, "ps.xargument") <- ans
  attr(basis, "ridge.adj") <- ridge.adj
  attr(basis, "ridge.inv") <- ridge.inv
  basis
}




