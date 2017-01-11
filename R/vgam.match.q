# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.




vgam.match <- function(x, all.knots = FALSE, nk = NULL) {

  if (is.list(x)) {
    nvar <- length(x)
    if (length(nk))
      nk <- rep_len(nk, nvar)
    temp <- vgam.match(x[[1]], all.knots = all.knots, nk = nk[1])

    ooo <- matrix(temp$matcho, length(temp$matcho), nvar)
    neffec <- rep(temp$neffec, nvar)
    xmin <- rep(temp$xmin, nvar)
    xmax <- rep(temp$xmax, nvar)
    nknots <- rep(temp$nknots, nvar)
    knots <- vector("list", nvar)
    knots[[1]] <- temp$knots

    if (nvar > 1)
      for (ii in 2:nvar) {
        temp <- vgam.match(x[[ii]], all.knots = all.knots, nk = nk[ii])
        ooo[, ii] <- temp$matcho
        neffec[ii] <- temp$neffec
        nknots[ii] <- temp$nknots
        knots[[ii]] <- temp$knots
        xmin[ii] <- temp$xmin
        xmax[ii] <- temp$xmax
      }
    names(nknots) <- names(knots) <-
    names(neffec) <- names(xmin) <- names(xmax) <- names(x)
    dimnames(ooo) <- list(NULL, names(x))

    return(list(matcho = ooo, neffec = neffec, nknots = nknots, knots = knots,
                xmin = xmin, xmax = xmax))
  }

  if (!is.null(attributes(x)$NAs) || anyNA(x))
    stop("cannot smooth on variables with NAs")

  sx <- unique(sort(as.vector(x)))  # "as.vector()" strips off attributes
  ooo <- match(x, sx)  # as.integer(match(x, sx))      # sx[o]==x
  neffec <- length(sx)  # as.integer(length(sx))

  if (neffec < 7)
    stop("smoothing variables must have at least 7 unique values")

  xmin <- sx[1]     # Don't use rounded value
  xmax <- sx[neffec]
  xbar <- (sx - xmin) / (xmax - xmin)

    noround <- TRUE   # Improvement 20020803
  if (all.knots) {
    knot <- if (noround) {
      valid.vknotl2(c(rep(xbar[1], 3), xbar, rep(xbar[neffec], 3)))
    } else {
      c(rep(xbar[1], 3), xbar, rep(xbar[neffec], 3))
    }
    if (length(nk))
      warning("overriding nk by all.knots = TRUE")
    nk <- length(knot) - 4 # No longer: neffec + 2
  } else {
    chosen <- length(nk)
    if (chosen && (nk > neffec+2 || nk <= 5))
      stop("bad value for 'nk'")
    if (!chosen)
      nk <- 0
    knot.list <- .C("vknootl2",
                    as.double(xbar),
                    as.integer(neffec), knot = double(neffec+6),
                    k = as.integer(nk+4), chosen = as.integer(chosen))
    if (noround) {
      knot <- valid.vknotl2(knot.list$knot[1:(knot.list$k)])
      knot.list$k <- length(knot)
    } else {
      knot <- knot.list$knot[1:(knot$k)]
    }
    nk <- knot.list$k - 4
  }
  if (nk <= 5)
    stop("not enough distinct knots found")

  return(list(matcho = ooo, neffec = neffec, nknots = nk, knots = knot,
              xmin = xmin, xmax = xmax))
}



