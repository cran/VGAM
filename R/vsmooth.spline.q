# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.






setClass("vsmooth.spline.fit", representation(
         "Bcoefficients" = "matrix",
         "knots"         = "numeric",
         "xmin"          = "numeric",
         "xmax"          = "numeric"))



setClass("vsmooth.spline", representation(
         "call"         = "call",
         "constraints"  = "list",
         "df"           = "numeric",
         "nlfit"        = "vsmooth.spline.fit",  # is the nonlinear component
         "lev"          = "matrix",
         "lfit" = "vlm",  # 20020606 was "vlm.wfit"; is the linear component
         "spar"         = "numeric",
         "lambda"       = "numeric",
         "var"          = "matrix",
         "w"            = "matrix",
         "x"            = "numeric",
         "y"            = "matrix",
         "yin"          = "matrix"))


setMethod("coefficients", signature(object = "vsmooth.spline"),
          function(object, ...)
          coefvsmooth.spline(object, ...))
setMethod("coef", signature(object = "vsmooth.spline"),
          function(object, ...)
          coefvsmooth.spline(object, ...))

setMethod("coefficients", signature(object = "vsmooth.spline.fit"),
          function(object, ...)
          coefvsmooth.spline.fit(object, ...))
setMethod("coef", signature(object = "vsmooth.spline.fit"),
          function(object, ...)
          coefvsmooth.spline.fit(object, ...))

setMethod("fitted.values", signature(object = "vsmooth.spline"),
          function(object, ...)
          fittedvsmooth.spline(object, ...))
setMethod("fitted", signature(object = "vsmooth.spline"),
          function(object, ...)
          fittedvsmooth.spline(object, ...))

setMethod("residuals", signature(object = "vsmooth.spline"),
          function(object, ...)
          residvsmooth.spline(object, ...))
setMethod("resid", signature(object = "vsmooth.spline"),
          function(object, ...)
          residvsmooth.spline(object, ...))

setMethod("predict", signature(object="vsmooth.spline"),
          function(object, ...)
          predictvsmooth.spline(object, ...))



setMethod("show",  "vsmooth.spline",
          function(object)
          show.vsmooth.spline(object))




setMethod("plot", "vsmooth.spline",
          function(x, y, ...) {
          if (!missing(y)) stop("cannot process the 'y' argument")
          invisible(plotvsmooth.spline(x, ...))})
setMethod("predict",  "vsmooth.spline.fit",
          function(object, ...)
          predictvsmooth.spline.fit(object, ...))



setMethod("model.matrix",  "vsmooth.spline",
          function(object, ...)
          model.matrixvlm(object, ...))




depvar.vsmooth.spline <- function(object, ...) {
  object@y
}


if (!isGeneric("depvar"))
    setGeneric("depvar", function(object, ...) standardGeneric("depvar"),
               package = "VGAM")


setMethod("depvar",  "vsmooth.spline", function(object, ...)
           depvar.vsmooth.spline(object, ...))




vsmooth.spline <-
  function(x, y, w = NULL, df = rep(5, M),
           spar = NULL,  #rep(0,M),
           i.constraint = diag(M),
           x.constraint = diag(M),
           constraints = list("(Intercepts)" = i.constraint,
                              x = x.constraint),
           all.knots = FALSE,
           var.arg = FALSE,
           scale.w = TRUE,
           nk = NULL,
           control.spar = list()) {



  if (var.arg) {
    warning("@var will be returned, but no use will be made of it")
  }


  missing.constraints <- missing(constraints)
  if (!(missing.spar <- missing(spar)) && !missing(df)) {
    stop("cannot specify both 'spar' and 'df'")
  }



  contr.sp <- list(low = -1.5,## low = 0.      was default till R 1.3.x
                   high = 1.5,
                   tol = 1e-4,## tol = 0.001   was default till R 1.3.x
                   eps = 2e-8,## eps = 0.00244 was default till R 1.3.x
                   maxit = 500 )




  contr.sp[(names(control.spar))] <- control.spar
  if (!all(sapply(contr.sp[1:4], is.numeric)) ||
      contr.sp$tol < 0 || contr.sp$eps <= 0 || contr.sp$maxit <= 0)
    stop("invalid 'control.spar'")


  my.call <- match.call()
  if (missing(y)) {
    if (is.list(x)) {
      if (anyNA(match(c("x", "y"), names(x))))
        stop("cannot find 'x' and 'y' in list")
      y <- x$y
      x <- x$x
    } else if (is.complex(x)) {
      y <- Im(x)
      x <- Re(x)
    } else if (is.matrix(x)) {
      y <- x[,-1]
      x <- x[,1]
    } else {
      y <- x
      x <- time(x)
    }
  }

  xvector <- x
  n_lm <- length(xvector)
  ymat <- as.matrix(y)
  ny2 <- dimnames(ymat)[[2]]  # NULL if vector
  M <- ncol(ymat)
  if (n_lm != nrow(ymat)) {
    stop("lengths of arguments 'x' and 'y' must match")
  }

  if (anyNA(xvector) || anyNA(ymat)) {
    stop("NAs not allowed in arguments 'x' or 'y'")
  }

  if (is.null(w)) {
    wzmat <- matrix(1, n_lm, M)
  } else {
    if (anyNA(w)) {
      stop("NAs not allowed in argument 'w'")
    }
    wzmat <- as.matrix(w)

    if (nrow(ymat) != nrow(wzmat) || ncol(wzmat) > M * (M+1) / 2) {
      stop("arguments 'w' and 'y' don't match")
    }

    if (scale.w) {
      wzmat <- wzmat / mean(wzmat[,1:M])    # 'Average' value is 1
    }
  }
  dim2wz <- ncol(wzmat)


  if (missing.constraints) {
    constraints <- list("(Intercepts)" = eval(i.constraint),
                        "x"            = eval(x.constraint))
  }
  constraints <- eval(constraints)
  if (is.matrix(constraints)) {
    constraints <- list("(Intercepts)" = constraints,
                        "x"            = constraints)
  }
  if (!is.list(constraints) || length(constraints) != 2) {
    stop("'constraints' must equal a list (of length 2) or a matrix")
  }
  for (ii in 1:2)
    if (!is.numeric(constraints[[ii]]) ||
        !is.matrix (constraints[[ii]]) ||
        nrow(constraints[[ii]]) != M   ||
        ncol(constraints[[ii]]) >  M)
      stop("something wrong with argument 'constraints'")
  names(constraints) <- c("(Intercepts)", "x")


    usortx <- unique(sort(as.vector(xvector)))
    ooo <- match(xvector, usortx)  # usortx[ooo] == x
    neff <- length(usortx)
    if (neff < 7) {
      stop("not enough unique 'x' values (need 7 or more)")
    }

    dim1U <- dim2wz  # 20000110; was M * (M+1) / 2

    collaps <- .C("vsuff9",
      as.integer(n_lm), as.integer(neff), as.integer(ooo),
      as.double(xvector), as.double(ymat), as.double(wzmat),

      xbar = double(neff), ybar = double(neff * M),
          wzbar = double(neff * dim2wz),
      uwzbar = double(1), wzybar = double(neff * M), okint = as.integer(0),
      as.integer(M), dim2wz = as.integer(dim2wz), dim1U = as.integer(dim1U),
      Hlist1 = as.double(diag(M)), ncolb = as.integer(M),
      trivc = as.integer(1), wuwzbar = as.integer(0),
      dim1Uwzbar = as.integer(dim1U), dim2wzbar = as.integer(dim2wz))

    if (collaps$okint != 1) {
      stop("some non-positive-definite weight matrices ",
           "detected in 'vsuff9'")
    }
    dim(collaps$ybar)   <- c(neff, M)


    if (FALSE) {
    } else {
      yinyin <- collaps$ybar  # Includes both linear and nonlinear parts
      x <- collaps$xbar  # Could call this xxx for location finder

      lfit <- vlm(yinyin ~ 1 + x,  # xxx
                 constraints = constraints,
                 save.weights = FALSE,
                 qr.arg = FALSE, x.arg = FALSE, y.arg = FALSE,
                 smart = FALSE,
                 weights = matrix(collaps$wzbar, neff, dim2wz))
    }

    ncb0  <- ncol(constraints[[2]])  # Of xxx and not of the intercept
    spar  <- rep_len(if (length(spar)) spar else 0, ncb0)
    dfvec <- rep_len(df, ncb0)

    if (!missing.spar) {
      ispar <- 1
      if (any(spar <= 0) || !is.numeric(spar)) {
        stop("not allowed non-positive or non-numeric ",
             "smoothing parameters")
      }
      nonlin <- (spar != Inf)
    } else {
      ispar <- 0
      if (!is.numeric(dfvec) || any(dfvec < 2 | dfvec > neff)) {
        stop("you must supply '2 <= df <= ", neff, "'")
      }
      nonlin <- (abs(dfvec - 2) > contr.sp$tol)
    }


    if (all(!nonlin)) {

      junk.fill <- new("vsmooth.spline.fit",
                       "Bcoefficients" = matrix(NA_real_, 1, 1),
                       "knots"         = numeric(0),
                       "xmin"          = numeric(0),
                       "xmax"          = numeric(0))  # 20031108

      dratio <- NA_real_

      object <-
      new("vsmooth.spline",
          "call"         = my.call,
          "constraints"  = constraints,
          "df"     = if (ispar == 0) dfvec else rep_len(2, length(spar)),
          "lfit"         = lfit,
          "nlfit"        = junk.fill,
          "spar"   = if (ispar == 1) spar   else rep_len(Inf, length(dfvec)),
          "lambda" = if (ispar == 1) dratio * 16.0^(spar * 6.0 - 2.0) else
                                     rep_len(Inf, length(dfvec)),
          "w"            = matrix(collaps$wzbar, neff, dim2wz),
          "x"            = usortx,
          "y"            = lfit@fitted.values,
          "yin"          = yinyin)


      return(object)
  }


  xbar <- (usortx - usortx[1]) / (usortx[neff] - usortx[1])
  noround <- TRUE   # Improvement 20020803
  nknots <- nk
  if (all.knots) {
    knot <- if (noround) {
      valid.vknotl2(c(rep_len(xbar[1], 3), xbar, rep_len(xbar[neff], 3)))
    } else {
      c(rep_len(xbar[1], 3), xbar, rep_len(xbar[neff], 3))
    }
    if (length(nknots)) {
      warning("overriding 'nk' by 'all.knots = TRUE'")
    }
    nknots <- length(knot) - 4  # No longer neff + 2
  } else {
    chosen <- length(nknots)
    if (chosen && (nknots > neff+2 || nknots <= 5)) {
      stop("bad value for 'nk'")
    }
    if (!chosen) {
      nknots <- 0
    }
    knot.list <- .C("vknootl2", as.double(xbar),
                      as.integer(neff), knot = double(neff+6),
                      k = as.integer(nknots+4),
                      chosen = as.integer(chosen))
    if (noround) {
      knot <- valid.vknotl2(knot.list$knot[1:(knot.list$k)])
      knot.list$k <- length(knot)
    } else {
      knot <- knot.list$knot[1:(knot.list$k)]
    }
    nknots <- knot.list$k - 4
  }
  if (nknots <= 5) {
    stop("not enough distinct knots found")
  }


  conmat <- (constraints[[2]])[, nonlin, drop = FALSE]
  ncb <- sum(nonlin)
  trivc <- trivial.constraints(conmat)
  resmat <- collaps$ybar - lfit@fitted.values     # neff by M
  spar.nl <-  spar[nonlin]
  dofr.nl <- dfvec[nonlin]

   dim1Uwzbar <- if (trivc) dim1U  else ncb * (ncb+1) / 2
   dim2wzbar  <- if (trivc) dim2wz else ncb * (ncb+1) / 2
   ooo <- 1:neff # Already sorted



  collaps <- .C("vsuff9",
      as.integer(neff), as.integer(neff), as.integer(ooo),
      as.double(collaps$xbar), as.double(resmat), as.double(collaps$wzbar),

      xbar = double(neff), ybar = double(neff * ncb),
          wzbar = double(neff * dim2wzbar),
      uwzbar = double(1), wzybar = double(neff * ncb), okint = as.integer(0),
      as.integer(M), as.integer(dim2wz), as.integer(dim1U),
      Hlist1 = as.double(conmat), ncolb = as.integer(ncb),
      as.integer(trivc), wuwzbar = as.integer(0),
      as.integer(dim1Uwzbar), as.integer(dim2wzbar))

  if (collaps$okint != 1) {
   stop("some non-positive-definite weight matrices ",
        "detected in 'vsuff9' during the second call.")
  }

  dim(collaps$ybar) <- dim(collaps$wzybar) <- c(neff, ncb)
  dim(collaps$wzbar) <- c(neff, dim2wzbar)







  wzyb.c <-
  zedd.c <- matrix(0, neff, ncb)
  Wmat.c <- array(0, c(ncb, ncb, neff))
 if (FALSE)
  for (ii in 1:neff) {
    Wi.indiv <- m2a(wzmat[ii, , drop = FALSE], M = ncb)
    Wi.indiv <- Wi.indiv[,, 1]  # Drop the 3rd dimension
    Wmat.c[,, ii] <- t(conmat) %*% Wi.indiv %*% conmat
    one.Wmat.c <- matrix(Wmat.c[,, ii], ncb, ncb)
    zedd.c[ii, ] <- solve(Wmat.c[,, ii],
                          t(conmat) %*% Wi.indiv %*% cbind(resmat[ii, ]))
    wzyb.c[ii, ] <- one.Wmat.c %*% zedd.c[ii, ]
  }










  ldk <- 3 * ncb + 1  # 20020710; Previously 4 * ncb
  varmat <- if (var.arg) matrix(0, neff, ncb) else double(1)





  vsplin <- .C("Yee_spline",
     xs = as.double(xbar),
     yyy = as.double(collaps$wzybar),  # zz

         as.double(collaps$wzbar), xknot = as.double(knot),
     n = as.integer(neff), nknots = as.integer(nknots), as.integer(ldk),
         M = as.integer(ncb), dim2wz = as.integer(dim2wzbar),

     spar.nl = as.double(spar.nl), lamvec = as.double(spar.nl),

         iinfo = integer(1), fv = double(neff * ncb),
     Bcoef = double(nknots * ncb), varmat = as.double(varmat),

     levmat = double(neff * ncb), as.double(dofr.nl),

     ifvar = as.integer(var.arg), ierror = as.integer(0),
     n_lm = as.integer(neff),
     double(nknots), double(nknots), double(nknots), double(nknots),
     double(1), as.integer(0),

     icontrsp = as.integer(contr.sp$maxit),
      contrsp = as.double(unlist(contr.sp[1:4])))







  if (vsplin$ierror != 0) {
    stop("vsplin$ierror == ", vsplin$ierror,
         ". Something gone wrong in 'vsplin'")
  }
  if (vsplin$iinfo != 0) {
    stop("leading minor of order ", vsplin$iinfo,
         " is not positive-definite")
  }

  dim(vsplin$levmat) <- c(neff, ncb)   # A matrix even when ncb == 1
  if (ncb > 1) {
    dim(vsplin$fv) <- c(neff, ncb)
    if (var.arg)
      dim(vsplin$varmat) <- c(neff, ncb)
  }

  dofr.nl <- colSums(vsplin$levmat)  # Actual EDF used




  fv <- lfit@fitted.values + vsplin$fv %*% t(conmat)
  if (M > 1) {
    dimnames(fv) <- list(NULL, ny2)
  }

  dfvec[!nonlin] <- 2.0
  dfvec[ nonlin] <- dofr.nl
  if (ispar == 0) {
    spar[!nonlin] <- Inf
    spar[ nonlin] <- vsplin$spar.nl   # Actually used
  }

  fit.object <- new("vsmooth.spline.fit",
                   "Bcoefficients" = matrix(vsplin$Bcoef, nknots, ncb),
                   "knots"         = knot,
                   "xmax"          = usortx[neff],
                   "xmin"          = usortx[1])

  object <-
  new("vsmooth.spline",
      "call"         = my.call,
      "constraints"  = constraints,
      "df"           = dfvec,
      "nlfit"        = fit.object,
      "lev"          = vsplin$levmat,
      "lfit"         = lfit,
      "spar"         = spar,   # if (ispar == 1) spar else vsplin$spar,
      "lambda"       = vsplin$lamvec,  #
      "w"            = collaps$wzbar,
      "x"            = usortx,
      "y"            = fv,
      "yin"          = yinyin)

  if (var.arg)
    object@var <- vsplin$varmat

  object
}


show.vsmooth.spline <- function(x, ...) {
  if (!is.null(cl <- x@call)) {
    cat("Call:\n")
    dput(cl)
  }

  ncb <- if (length(x@nlfit)) ncol(x@nlfit@Bcoefficients) else NULL
  cat("\nSmoothing Parameter (Spar):",
    if (length(ncb) && ncb == 1) format(x@spar) else
        paste(format(x@spar), collapse = ", "), "\n")

  cat("\nEquivalent Degrees of Freedom (Df):",
    if (length(ncb) && ncb == 1) format(x@df) else
        paste(format(x@df), collapse = ", "), "\n")

  if (!all(trivial.constraints(x@constraints) == 1)) {
    cat("\nConstraint matrices:\n")
    print(x@constraints)
  }

  invisible(x)
}


coefvsmooth.spline.fit <- function(object, ...) {
  object@Bcoefficients
}


coefvsmooth.spline <- function(object, matrix = FALSE, ...) {

        list(lfit = coefvlm(object@lfit, matrix.out = matrix),
             nlfit = coefvsmooth.spline.fit(object@nlfit))
}


fittedvsmooth.spline <- function(object, ...) {
  object@y
}

residvsmooth.spline <- function(object, ...) {
  as.matrix(object@yin - object@y)
}



plotvsmooth.spline <- function(x, xlab = "x", ylab = "", points = TRUE,
                               pcol = par()$col, pcex = par()$cex,
                               pch = par()$pch, lcol = par()$col,
                               lwd = par()$lwd, lty = par()$lty,
                               add = FALSE, ...) {
  points.arg <- points; rm(points)
  M <- ncol(x@y)
  pcol <- rep_len(pcol, M)
  pcex <- rep_len(pcex, M)
  pch  <- rep_len(pch,  M)
  lcol <- rep_len(lcol, M)
  lwd  <- rep_len(lwd,  M)
  lty  <- rep_len(lty,  M)
  if (!add)
    matplot(x@x, x@yin, type = "n", xlab = xlab, ylab = ylab, ...)
  for (ii in 1:ncol(x@y)) {
    if (points.arg)
      points(x@x, x@yin[,ii], col = pcol[ii], pch = pch[ii], cex = pcex[ii])
    lines(x@x, x@y[,ii], col = lcol[ii], lwd = lwd[ii], lty = lty[ii])
  }
  invisible(x)
}



predictvsmooth.spline <- function(object, x, deriv = 0, se.fit = FALSE) {
  if (se.fit)
    warning("'se.fit = TRUE' is not currently implemented. ",
            "Using 'se.fit = FALSE'")

   lfit <- object@lfit    #    Linear part of the vector spline
  nlfit <- object@nlfit   # Nonlinear part of the vector spline

  if (missing(x)) {
    if (deriv == 0) {
      return(list(x = object@x, y = object@y))
    } else {
      x <- object@x
      return(Recall(object, x, deriv))
    }

  }

  mat.coef <- coefvlm(lfit, matrix.out = TRUE)
  coeflfit <- t(mat.coef)   # M x p now
  M <- nrow(coeflfit)  # if (is.matrix(object@y)) ncol(object@y) else 1

  pred <- if (deriv == 0)
           predict(lfit, data.frame(x = x)) else
          if (deriv == 1)
            matrix(coeflfit[,2], length(x), M, byrow = TRUE) else
            matrix(0, length(x), M)
  if (!length(nlfit@knots)) {
    return(list(x = x, y = pred))
  }

  nonlin <- (object@spar != Inf)

  conmat <- if (!length(lfit@constraints)) diag(M) else
              lfit@constraints[[2]]
  conmat <- conmat[, nonlin, drop = FALSE]  # Of nonlinear functions

  list(x = x, y = pred + predict(nlfit, x, deriv)$y %*% t(conmat))
}




predictvsmooth.spline.fit <- function(object, x, deriv = 0) {
  nknots <- nrow(object@Bcoefficients)
  drangex <- object@xmax - object@xmin
  if (missing(x))
    x <- seq(from = object@xmin, to = object@xmax, length.out = nknots-4)

  xs <- as.double((x - object@xmin) / drangex)

  bad.left  <- (xs <  0)
  bad.right <- (xs >  1)
  good <- !(bad.left | bad.right)

  ncb <- ncol(object@Bcoefficients)
  y <- matrix(NA_real_, length(xs), ncb)
  if (ngood <- sum(good)) {
    junk <- .C("Yee_vbvs", as.integer(ngood),
          as.double(object@knots), as.double(object@Bcoefficients),
          as.double(xs[good]), smomat = double(ngood * ncb),
          as.integer(nknots), as.integer(deriv), as.integer(ncb))
    y[good,] <- junk$smomat

    if (TRUE && deriv > 1) {
      edges <- xs <= 0 | xs >= 1 # Zero the edges & beyond explicitly
      y[edges,] <- 0
    }
  }
  if (any(!good)) {
    xrange <- c(object@xmin, object@xmax)
    if (deriv == 0) {
      end.object <- Recall(object, xrange)$y
      end.slopes <- Recall(object, xrange, 1)$y * drangex

      if (any(bad.left)) {
        y[bad.left,] <-  rep(end.object[1,], rep(sum(bad.left), ncb)) +
                         rep(end.slopes[1,], rep(sum(bad.left), ncb)) *
                         xs[bad.left]
      }
      if (any(bad.right)) {
        y[bad.right,] <- rep(end.object[2,], rep(sum(bad.right), ncb)) +
                         rep(end.slopes[2,], rep(sum(bad.right), ncb)) *
                         (xs[bad.right] - 1)
      }
    } else if (deriv == 1) {
      end.slopes <- Recall(object, xrange, 1)$y * drangex
      y[bad.left,]  <- rep(end.slopes[1,], rep(sum(bad.left),  ncb))
      y[bad.right,] <- rep(end.slopes[2,], rep(sum(bad.right), ncb))
    } else
      y[!good,] <- 0
  }
  if (deriv > 0)
    y <- y / (drangex^deriv)
  list(x = x, y = y)
}



valid.vknotl2 <- function(knot, tol = 1/1024) {

  junk <- .C("Yee_pknootl2", knot = as.double(knot),
               as.integer(length(knot)),
               keep = integer(length(knot)), as.double(tol))
  keep <- as.logical(junk$keep)
  knot <- junk$knot[keep]
  if (length(knot) <= 11) {
    stop("too few (distinct) knots")
  }
  knot
}


