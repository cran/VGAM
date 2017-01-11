# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





 sm.os <-
  function(x, ...,
           niknots = 6,  # NULL if 'alg.niknots' is to be used.
           spar = -1,  # was 0 prior to 20160810
           o.order = 2,
           alg.niknots = c("s", ".nknots.smspl")[1],
           all.knots = FALSE,  # 20161013
           ridge.adj = 1e-5,
           spillover = 0.01, maxspar = 1e12,
           outer.ok = FALSE,
           fixspar = FALSE) {






  niknots.orig <- niknots
  if (all.knots && length(niknots))
    warning("ignoring 'niknots' because 'all.knots = TRUE'")






Penalty.os <- function(a, b, intKnots, o.order = 2) {



  if (any(diff(intKnots) <= 0))
    stop("argument 'intKnots' must be sorted in increasing order")
  if (length(unique(intKnots)) != (KK <- length(intKnots)))
    stop("argument 'intKnots' must have unique values")
  if (KK == 0)
    stop("no interior knots (intKnots == 0)")


  allKnots <- c(rep(a, 2 * o.order), intKnots, rep(b, 2 * o.order))
  mkmat <- matrix(c(1, rep(NA, 6),
           1/3, 4/3, 1/3, rep(NA, 4),
           14/45, 64/45, 8/15, 64/45, 14/45, NA, NA,
           41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140),
           4, 7, byrow = TRUE)

  vec.ell  <- 1:(KK + 4 * o.order - 1)  # length(allKnots) - 1
  vec.ellp <- 0:(2 * o.order - 2)

  hmell <- if (o.order == 1)
    diff(allKnots) else diff(allKnots) /  (2 * o.order - 2)

  xtilde <- wts <- numeric((2*o.order - 1) * (KK * 4*o.order - 1))
  index1 <- (2*o.order - 1) * (vec.ell - 1) + 1

  for (ellp in vec.ellp) {
    xtilde[index1 + ellp] <- hmell * ellp + allKnots[vec.ell]
       wts[index1 + ellp] <- hmell * mkmat[o.order, ellp + 1]
  }

  Bdd <- splineDesign(allKnots, xtilde,
                      ord = 2 * o.order,
                      derivs = rep(o.order, length(xtilde)),
                      outer.ok = TRUE)
  Omega <- crossprod(Bdd * wts, Bdd)
  attr(Omega, "allKnots") <- allKnots
  Omega
}






  xs <- substitute(x)
  ans <- as.character(xs)
  x.index <- as.vector(x)


  alg.niknots <- match.arg(alg.niknots,
                          c("s", ".nknots.smspl"))[1]
  if (!is.Numeric(o.order, length.arg = 1, integer.valued = TRUE,
                  positive = TRUE) ||
      o.order > 4)
    stop("argument 'o.order' must be one value from the set 1:4")



  x.orig <- x.index
  xdots <- list(...)
  uses.xij <- length(xdots) > 0
  if (uses.xij)
    x.index <- as.vector(c(x.index, unlist(xdots)))


  xl <- min(x.index)
  xr <- max(x.index)


 if (smart.mode.is("read")) {
    smartlist  <- get.smart()
    xl <- smartlist$xl  # Overwrite its value
    xr <- smartlist$xr  # Overwrite its value
    alg.niknots <- smartlist$alg.niknots  # Ditto
    spar       <- smartlist$spar
    o.order     <- smartlist$o.order
    all.knots  <- smartlist$all.knots
    ridge.adj  <- smartlist$ridge.adj
    spillover  <- smartlist$spillover
    maxspar    <- smartlist$maxspar
    maXX       <- smartlist$maXX
    Cmat       <- smartlist$Cmat
    intKnots   <- smartlist$intKnots
    outer.ok   <- smartlist$outer.ok
    fixspar    <- smartlist$fixspar
  } else {
    intKnots   <- NULL
    maXX       <- NULL
    Cmat       <- NULL
  }




  xmax <- xr + spillover * (xr - xl)
  xmin <- xl - spillover * (xr - xl)
  nx <- names(x.index)
  nax <- is.na(x.index)
  if (nas <- any(nax))
    x.index <- x[!nax]




  usortx <- unique(sort(as.vector(x.index)))
  neff <- length(usortx)
  if (neff < 2) {
    stop("not enough unique 'x' values (need 2 or more)")
  }

  noround <- TRUE   # Improvement 20020803


  if (all.knots) {
    xbar <- (usortx - usortx[1]) / (usortx[neff] - usortx[1])
    knot <- if (noround) {
      valid.vknotl2(c(rep_len(xbar[   1], 2 * o.order - 1),  # 3
                      xbar,
                      rep_len(xbar[neff], 2 * o.order - 1)))  # 3
    } else {
      c(rep_len(xbar[   1], 2 * o.order - 1),
        xbar,
        rep_len(xbar[neff], 2 * o.order - 1))
    }
    if (length(niknots.orig)) {
      warning("overriding 'niknots' by 'all.knots = TRUE'")
    }
    niknots <- length(knot) - 2 * o.order  # TWYee
  } else if (is.null(niknots.orig)) {
    xbar <- (usortx - usortx[1]) / (usortx[neff] - usortx[1])
    if (alg.niknots == "s") {
      chosen <- length(niknots)
      if (chosen && (niknots > neff + 2 || niknots <= 5)) {
        stop("bad value for 'niknots'")
      }
      if (!chosen) {
        niknots <- 0
      }
      knot.list <-
        .C("vknootl2", as.double(xbar),
           as.integer(neff),
           knot = double(neff + 4 * o.order - 2),  # (neff+6), zz unsure
           k = as.integer(niknots + 2 * o.order),  # (niknots+4), zz unsure
           chosen = as.integer(chosen))
      if (noround) {
        knot <- valid.vknotl2(knot.list$knot[1:(knot.list$k)])
        knot.list$k <- length(knot)
      } else {
        knot <- knot.list$knot[1:(knot.list$k)]
      }
      niknots <- knot.list$k - 2 * o.order  # TWYee
    } else {
      niknots <- .nknots.smspl(neff)
    }

  }  # !all.knots



  if (!is.Numeric(niknots, positive = TRUE, integer.valued = TRUE,
                  length.arg = 1)) {
    stop("bad value of 'niknots'")
  }






  numIntKnots <- niknots


  if (is.null(intKnots))
    intKnots <- quantile(usortx,  # unique(x.index),
                         probs = seq(0, 1, length = numIntKnots + 2)[
                                 -c(1, numIntKnots + 2)])





  Basis <- bs(x, knots = intKnots,
              degree = 2 * o.order - 1,  # 3 by default
              Boundary.knots = c(a = xmin, b = xmax),  # zz not sure
              intercept = TRUE)




  n.col <- ncol(Basis)
  if (nas) {
    nmat <- matrix(NA_real_, length(nax), n.col)
    nmat[!nax, ] <- Basis
    Basis <- nmat
  }
  dimnames(Basis) <- list(1:nrow(Basis), 1:n.col)


  fixspar <- rep_len(fixspar, max(length(fixspar), length(spar)))
     spar <- rep_len(   spar, max(length(fixspar), length(spar)))

  if (any(spar < 0 & fixspar)) {
     spar[spar < 0 & fixspar] <- 0
     warning("some 'spar' values are negative : have used 'spar' = ",
             paste(spar, collapse = ", "))
  }

  if (any(maxspar < spar)) {
     spar[maxspar < spar] <- maxspar
     warning("some 'spar' values are > ", maxspar, ": ",
                   "for stability have used 'spar' = ",
             paste(spar, collapse = ", "))
  }






  pen.aug <- Penalty.os(a = xmin, b = xmax, intKnots, o.order = o.order)
  allKnots <- attr(pen.aug, "allKnots")  # Retrieved





  if (is.null(maXX))
    maXX <- mean(abs(crossprod(Basis)))
 maS <- mean(abs(pen.aug)) / maXX

  pen.aug <- pen.aug / maS




  kk <- ncol(Basis)
  if (is.null(Cmat))
    Cmat <- matrix(colSums(Basis), 1, kk)

  qrCt <- qr(t(Cmat))
  jay <- nrow(Cmat)  # 1
  XZ <- t(qr.qty(qrCt, t(Basis))[(jay+1):kk, ])

  Basis <- XZ


  ZtSZ <- qr.qty(qrCt, t(qr.qty(qrCt, t(pen.aug))))[(jay+1):kk,
                                                    (jay+1):kk]





  if (smart.mode.is("write"))
    put.smart(list(xl          = xl,
                   xr          = xr,
                   alg.niknots = alg.niknots,
                   spar        = spar,
                   o.order     = o.order,
                   all.knots   = all.knots,
                   ridge.adj   = ridge.adj,
                   spillover   = spillover,
                   maxspar     = maxspar,
                   maXX        = maXX,
                   Cmat        = Cmat,
                   intKnots    = intKnots,
                   outer.ok    = outer.ok,
                   fixspar     = fixspar))




  Basis <- Basis[seq_along(x.orig), , drop = FALSE]





  attr(Basis, "S.arg") <- ZtSZ

  attr(Basis, "knots") <- allKnots  # zz might be intKnots
  attr(Basis, "intKnots") <- intKnots
  attr(Basis, "spar") <- spar  # Vector
  attr(Basis, "o.order") <- o.order  # Save argument
  attr(Basis, "ps.int") <- NA_real_  # For the psint() methods function
  attr(Basis, "all.knots") <- all.knots  # Save logical argument
  attr(Basis, "alg.niknots") <- alg.niknots  # Save argument
  attr(Basis, "ridge.adj") <- ridge.adj  # Save argument
  attr(Basis, "outer.ok") <- outer.ok  # Save argument
  attr(Basis, "fixspar") <- fixspar  # Save argument

  Basis
}



