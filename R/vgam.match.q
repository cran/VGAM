# These functions are
# Copyright (C) 1998-2010 T.W. Yee, University of Auckland. All rights reserved.



vgam.match <- function(x, all.knots=FALSE, nk=NULL) {

  if (is.list(x)) {
      nvar <- length(x)
      if (length(nk))
          nk = rep(nk, length=nvar)
      temp <- vgam.match(x[[1]], all.knots=all.knots, nk=nk[1])

      ooo <- matrix(temp$o, length(temp$o), nvar)
      nef <- rep(temp$nef, nvar)
      xmin <- rep(temp$xmin, nvar)
      xmax <- rep(temp$xmax, nvar)
      nknots <- rep(temp$nknots, nvar)
      knots <- vector("list", nvar)
      knots[[1]] <- temp$knots

      if (nvar > 1) 
          for (ii in 2:nvar) {
              temp = vgam.match(x[[ii]], all.knots = all.knots, nk = nk[ii])
              ooo[, ii] <- temp$o
              nef[ii] <- temp$nef
              nknots[ii] <- temp$nknots
              knots[[ii]] <- temp$knots
              xmin[ii] <- temp$xmin
              xmax[ii] <- temp$xmax
          }
      names(nknots) <- names(knots) <- 
      names(nef) <- names(xmin) <- names(xmax) <- names(x)
      dimnames(ooo) <- list(NULL, names(x))

      return(list(o=ooo, nef=nef, nknots=nknots, knots=knots,
                  xmin=xmin, xmax=xmax))
  }

  if (!is.null(attributes(x)$NAs) || any(is.na(x)))
      stop("cannot smooth on variables with NAs") 

  sx <- unique(sort(as.vector(x))) # "as.vector()" strips off attributes
  ooo <- match(x, sx)  # as.integer(match(x, sx))      # sx[o]==x
  nef <- length(sx)  # as.integer(length(sx))

  if (nef < 7)
      stop("smoothing variables must have at least 7 unique values")

  xmin <- sx[1]     # Don't use rounded value 
  xmax <- sx[nef]
  xbar <- (sx - xmin) / (xmax - xmin)

    noround = TRUE   # Improvement 3/8/02
  if (all.knots) {
      if (noround) {
          knot = valid.vknotl2(c(rep(xbar[1], 3), xbar, rep(xbar[nef], 3)))
      } else {
          knot <- c(rep(xbar[1], 3), xbar, rep(xbar[nef], 3))
      }
      if (length(nk)) warning("overriding nk by all.knots = TRUE")
      nk <- length(knot) - 4    # No longer: nef + 2
  } else {
      chosen = length(nk)
      if (chosen && (nk > nef+2 || nk <= 5))
          stop("bad value for 'nk'")
      if (!chosen) nk = 0
      knot.list <- dotC(name="vknootl2", as.double(xbar),
                        as.integer(nef), knot=double(nef+6),
                        k=as.integer(nk+4), chosen=as.integer(chosen))
      if (noround) {
          knot = valid.vknotl2(knot.list$knot[1:(knot.list$k)])
          knot.list$k = length(knot)
      } else {
          knot <- knot.list$knot[1:(knot$k)]
      }
      nk <- knot.list$k - 4
  }
  if (nk <= 5) stop("not enough distinct knots found")

  return(list(o=ooo, nef=nef, nknots=nk, knots=knot, xmin=xmin, xmax=xmax))
}



