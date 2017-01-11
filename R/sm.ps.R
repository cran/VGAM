# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.







 sm.ps <-
  function(x,
           ...,
           ps.int = NULL,
           spar = -1,  # was 0 prior to 20160810
           degree = 3, p.order = 2,
           ridge.adj = 1e-5,  # ridge.inv = 0.0001,
           spillover = 0.01, maxspar = 1e12,
           outer.ok = FALSE,
           mux = NULL,  # 1.25,
           fixspar = FALSE) {


  xs <- substitute(x)
  ans <- as.character(xs)
  x.index <- as.vector(x)




  x.orig <- x.index
  xdots <- list(...)
  uses.xij <- length(xdots) > 0
  if (uses.xij)
    x.index <- as.vector(c(x.index, unlist(xdots)))
  if (is.null(ps.int)) {
    ps.int <- if (length(mux)) {
      nux <- length(unique(x.index))
      ceiling(mux * log(nux))
    } else {
      min(max(degree, 7), length(x.index) - 2)
    }
  }  # if (is.null(ps.int))
  if (length(x.index) - 1 <= ps.int)
    stop("argument 'ps.int' is too large")


  xl <- min(x.index)
  xr <- max(x.index)



 if (smart.mode.is("read")) {
    smartlist  <- get.smart()
    xl <- smartlist$xl  # Overwrite its value
    xr <- smartlist$xr  # Overwrite its value
    ps.int    <- smartlist$ps.int  # Ditto
    spar      <- smartlist$spar
    degree    <- smartlist$degree
    p.order   <- smartlist$p.order
    ridge.adj <- smartlist$ridge.adj
    spillover <- smartlist$spillover
    maxspar   <- smartlist$maxspar
    maXX      <- smartlist$maXX
    Cmat      <- smartlist$Cmat
    outer.ok  <- smartlist$outer.ok
    mux       <- smartlist$mux
    fixspar   <- smartlist$fixspar
  } else {
    maXX      <- NULL
    Cmat      <- NULL
  }

  xmax <- xr + spillover * (xr - xl)
  xmin <- xl - spillover * (xr - xl)
  dx <- (xmax - xmin) / ps.int
  nx <- names(x.index)
  nax <- is.na(x.index)
  if (nas <- any(nax))
    x.index <- x[!nax]
  s.order <- degree + 1
  if (length(ps.int)) {
    nAknots <- ps.int - 1
    if (nAknots < 1) {
      nAknots <- 1
      warning("'ps.int' was too small; have used 2")
    }


  if (FALSE && # nux < 6 &&
      smart.mode.is("write"))
    warning("smoothing when there are less than 6 distinct 'x' values",
            " is not advised")



    if (nAknots > 0) {
      Aknots <- seq(from = xmin - degree * dx,
                    to   = xmax + degree * dx, by = dx)
    } else {
      knots <- NULL
    }
  }  # length(ps.int)





  basis <- splineDesign(Aknots, x.index, s.order, 0 * x.index,
                        outer.ok = outer.ok)
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA_real_, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(1:nrow(basis), 1:n.col)
  if ((p.order - n.col + 1) > 0) {
    p.order <- n.col - 1
    warning("argument 'p.order' was too large; have used ", n.col - 1)
  }
  fixspar <- rep_len(fixspar, max(length(fixspar), length(spar)))
     spar <- rep_len(   spar, max(length(fixspar), length(spar)))

  if (any(spar < 0 & fixspar)) {
    spar[spar < 0 & fixspar] <- 0
    warning("some 'spar' values are negative : have used 'spar' = ",
            paste(spar, collapse = ", "))
  }

  if (any(spar > maxspar)) {
    spar[spar > maxspar] <- maxspar
    warning("some 'spar' values are > ", maxspar, ": ",
                  "for stability have used 'spar' = ",
            paste(spar, collapse = ", "))
  }
  aug <- if (p.order > 0) diff(diag(n.col), diff = p.order) else diag(n.col)

  
  pen.aug <- crossprod(aug)



  if (is.null(maXX))
    maXX <- mean(abs(crossprod(basis)))
 maS <- mean(abs(pen.aug)) / maXX

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
                   ps.int    = ps.int,
                   spar      = spar,
                   degree    = degree,
                   p.order   = p.order,
                   ridge.adj = ridge.adj,
                   spillover = spillover,
                   maxspar   = maxspar,
                   maXX      = maXX,
                   Cmat      = Cmat,
                   outer.ok  = outer.ok,
                   mux       = mux,
                   fixspar   = fixspar))




  basis <- basis[seq_along(x.orig), , drop = FALSE]




  attr(basis, "S.arg") <- ZtSZ

  attr(basis, "degree") <- degree
  attr(basis, "knots") <- Aknots
  attr(basis, "spar") <- spar  # Vector
  attr(basis, "p.order") <- p.order
  attr(basis, "ps.int") <- ps.int
  attr(basis, "ridge.adj") <- ridge.adj
  attr(basis, "outer.ok") <- outer.ok
  attr(basis, "mux") <- mux
  attr(basis, "fixspar") <- fixspar
  basis
}




