# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.




















smartpredenv <- new.env()


smart.mode.is <- function(mode.arg = NULL) {
  if (!length(mode.arg)) {
    if (exists(".smart.prediction", envir = VGAM:::smartpredenv)) {
      get(".smart.prediction.mode", envir = VGAM:::smartpredenv)
    } else {
      "neutral"
    }
  } else {
    if (mode.arg != "neutral" &&
        mode.arg != "read" &&
        mode.arg != "write")
      stop("argument \"mode.arg\" must be one of",
           " \"neutral\", \"read\" or \"write\"")
    if (exists(".smart.prediction", envir = VGAM:::smartpredenv)) {
      get(".smart.prediction.mode", envir = VGAM:::smartpredenv) ==
        mode.arg
    } else {
      mode.arg == "neutral"
    }
  }
}


setup.smart <- function(mode.arg, smart.prediction = NULL,
                        max.smart = 30) {
  actual <- if (mode.arg == "write") vector("list", max.smart) else 
            if (mode.arg == "read") smart.prediction else
            stop("value of 'mode.arg' unrecognized")

  wrapup.smart()  # make sure

  if (length(actual)) {


    assign(".smart.prediction", actual, envir = VGAM:::smartpredenv)
    assign(".smart.prediction.counter", 0, envir = VGAM:::smartpredenv)
    assign(".smart.prediction.mode", mode.arg, envir = VGAM:::smartpredenv)
    assign(".max.smart", max.smart, envir = VGAM:::smartpredenv)
    assign(".smart.prediction", actual, envir = VGAM:::smartpredenv)
  }
}

wrapup.smart <- function() {
  if (exists(".smart.prediction", envir = VGAM:::smartpredenv))
    rm(".smart.prediction", envir = VGAM:::smartpredenv)
  if (exists(".smart.prediction.counter", envir = VGAM:::smartpredenv))
    rm(".smart.prediction.counter", envir = VGAM:::smartpredenv)
  if (exists(".smart.prediction.mode", envir = VGAM:::smartpredenv))
    rm(".smart.prediction.mode", envir = VGAM:::smartpredenv)
  if (exists(".max.smart", envir = VGAM:::smartpredenv))
    rm(".max.smart", envir = VGAM:::smartpredenv)
}


get.smart.prediction <- function() {

  smart.prediction.counter <- get(".smart.prediction.counter",
                                  envir = VGAM:::smartpredenv)
  max.smart <- get(".max.smart", envir = VGAM:::smartpredenv)

  if (smart.prediction.counter > 0) {
    # Save this on the object for smart prediction later
    smart.prediction <- get(".smart.prediction", envir = VGAM:::smartpredenv)
    if (max.smart >= (smart.prediction.counter + 1))
      for(i in max.smart:(smart.prediction.counter + 1))
        smart.prediction[[i]] <- NULL
    smart.prediction
  } else 
    NULL
}


put.smart <- function(smart) {



  max.smart <- get(".max.smart", envir = VGAM:::smartpredenv)
  smart.prediction.counter <- get(".smart.prediction.counter",
                                  envir = VGAM:::smartpredenv)
  smart.prediction <- get(".smart.prediction", envir = VGAM:::smartpredenv)
  smart.prediction.counter <- smart.prediction.counter + 1

  if (smart.prediction.counter > max.smart) {
    # if list is too small, make it larger
    max.smart <- max.smart + (inc.smart <- 10) # can change inc.smart
    smart.prediction <- c(smart.prediction, vector("list", inc.smart))
    assign(".max.smart", max.smart, envir = VGAM:::smartpredenv)
  }

  smart.prediction[[smart.prediction.counter]] <- smart
  assign(".smart.prediction", smart.prediction, envir = VGAM:::smartpredenv)
  assign(".smart.prediction.counter", smart.prediction.counter,
         envir = VGAM:::smartpredenv)
}


get.smart <- function() {
  # Returns one list component of information
  smart.prediction <- get(".smart.prediction", envir = VGAM:::smartpredenv)
  smart.prediction.counter <- get(".smart.prediction.counter",
                                  envir = VGAM:::smartpredenv)
  smart.prediction.counter <- smart.prediction.counter + 1
  assign(".smart.prediction.counter", smart.prediction.counter,
         envir = VGAM:::smartpredenv)
  smart <- smart.prediction[[smart.prediction.counter]]
  smart
}

smart.expression <- expression({


  smart  <- get.smart()
  assign(".smart.prediction.mode", "neutral", envir = VGAM:::smartpredenv)

  .smart.match.call <- as.character(smart$match.call)
  smart$match.call <- NULL  # Kill it off for the do.call 

  ans.smart <- do.call(.smart.match.call[1], c(list(x=x), smart))
  assign(".smart.prediction.mode", "read", envir = VGAM:::smartpredenv)

  ans.smart
})



is.smart <- function(object) {
  if (is.function(object)) {
    if (is.logical(a <- attr(object, "smart"))) a else FALSE
  } else {
    if (length(slotNames(object))) {
        if (length(object@smart.prediction) == 1 &&
            is.logical(object@smart.prediction$smart.arg))
        object@smart.prediction$smart.arg else
            any(slotNames(object) == "smart.prediction")
    } else {
      if (length(object$smart.prediction) == 1 &&
          is.logical(object$smart.prediction$smart.arg))
        object$smart.prediction$smart.arg else
        any(names(object) == "smart.prediction")
    }
  }
}










library(splines) 



bs <-
function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
    Boundary.knots = range(x)) 
{
    x <- x  # Evaluate x
    if (smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1]) | (or <- x > 
            Boundary.knots[2L])
    } else outside <- FALSE
    ord <- 1 + (degree <- as.integer(degree))
    if (ord <= 1) 
        stop("'degree' must be integer >= 1")
    if (!missing(df) && missing(knots)) {
        nIknots <- df - ord + (1 - intercept)
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used  ", ord - 
                (1 - intercept))
        }
        knots <- if (nIknots > 0) {
            knots <- seq(from = 0, to = 1, length = nIknots + 
                2)[-c(1, nIknots + 2)]
            stats::quantile(x[!outside], knots)
        }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))
    if (any(outside)) {
        warning("some 'x' values beyond boundary knots may ",
                "cause ill-conditioned bases")
        derivs <- 0:degree
        scalef <- gamma(1L:ord)
        basis <- array(0, c(length(x), length(Aknots) - degree - 
            1L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord, 
                derivs)$design
            basis[ol, ] <- xl %*% (tt/scalef)
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord, 
                derivs)$design
            basis[or, ] <- xr %*% (tt/scalef)
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                ord)$design
    } else basis <- spline.des(Aknots, x, ord)$design
    if (!intercept) 
        basis <- basis[, -1L, drop = FALSE]
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    a <- list(degree = degree,
        knots = if (is.null(knots)) numeric(0L) else knots, 
        Boundary.knots = Boundary.knots,
        intercept = intercept,
        Aknots = Aknots)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("bs", "basis", "matrix")

    if (smart.mode.is("write"))
        put.smart(list(df = df,
                       knots = knots,
                       degree = degree,
                       intercept = intercept,
                       Boundary.knots = Boundary.knots,
                       match.call = match.call()))

    basis
}
attr(bs, "smart") <- TRUE






ns <-
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) 
{
    x <- x  # Evaluate x
    if (smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    } else outside <- FALSE
    if (!missing(df) && missing(knots)) {
        nIknots <- df - 1 - intercept
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used ", 1 + intercept)
        }
        knots <- if (nIknots > 0) {
            knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, nIknots +
                2L)]
            stats::quantile(x[!outside], knots)
        }
    } else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, x[ol] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[ol, ] <- xl %*% tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, x[or] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[or, ] <- xr %*% tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 4)$design
    } else basis <- spline.des(Aknots, x, 4)$design
    const <- spline.des(Aknots, Boundary.knots, 4, c(2, 2))$design
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),
        drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    a <- list(degree = 3,
              knots = if (is.null(knots)) numeric(0) else knots, 
              Boundary.knots = Boundary.knots,
              intercept = intercept,
              Aknots = Aknots)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("ns", "basis", "matrix")

    if (smart.mode.is("write"))
        put.smart(list(df = df,
                       knots = knots,
                       intercept = intercept,
                       Boundary.knots = Boundary.knots,
                       match.call = match.call()))

    basis
}
attr(ns, "smart") <- TRUE








poly <-
function (x, ..., degree = 1, coefs = NULL, raw = FALSE) 
{
    x <- x  # Evaluate x
    if (!raw && smart.mode.is("read")) {
        smart <- get.smart()
        degree <- smart$degree
        coefs  <- smart$coefs
        raw  <- smart$raw
    }

    dots <- list(...)
    if (nd <- length(dots)) {
        if (nd == 1 && length(dots[[1]]) == 1L)
            degree <- dots[[1L]] else
        return(polym(x, ..., degree = degree, raw = raw))
    }
    if (is.matrix(x)) {
        m <- unclass(as.data.frame(cbind(x, ...)))
        return(do.call("polym", c(m, degree = degree, raw = raw)))
    }
    if (degree < 1) 
        stop("'degree' must be at least 1")



    # At prediction time x may be less than the degree
    if (smart.mode.is("write") || smart.mode.is("neutral"))
    if (degree >= length(x))
        stop("degree must be less than number of points")




    if (any(is.na(x))) 
        stop("missing values are not allowed in 'poly'")
    n <- degree + 1
    if (raw) {
        if (degree >= length(unique(x)))
            stop("'degree' must be less than number of unique points")
        Z <- outer(x, 1L:degree, "^")
        colnames(Z) <- 1L:degree
        attr(Z, "degree") <- 1L:degree
        class(Z) <- c("poly", "matrix")
        return(Z)
    }
    if (is.null(coefs)) {
        if (degree >= length(unique(x))) 
            stop("'degree' must be less than number of unique points")
        xbar <- mean(x)
        x <- x - xbar
        X <- outer(x, seq_len(n) - 1, "^")
        QR <- qr(X)

        if (QR$rank < degree)
            stop("'degree' must be less than number of unique points")

        z <- QR$qr
        z <- z * (row(z) == col(z))
        raw <- qr.qy(QR, z)
        norm2 <- colSums(raw^2)
        alpha <- (colSums(x * raw^2)/norm2 + xbar)[1L:degree]
        Z <- raw/rep(sqrt(norm2), each = length(x))
        colnames(Z) <- 1L:n - 1L
        Z <- Z[, -1, drop = FALSE]
        attr(Z, "degree") <- 1:degree
        attr(Z, "coefs") <- list(alpha = alpha, norm2 = c(1, 
            norm2))
        class(Z) <- c("poly", "matrix")
    } else {
        alpha <- coefs$alpha
        norm2 <- coefs$norm2
        Z <- matrix(, length(x), n)
        Z[, 1] <- 1
        Z[, 2] <- x - alpha[1L]
        if (degree > 1) 
            for (i in 2:degree) Z[, i + 1] <- (x - alpha[i]) * 
                Z[, i] - (norm2[i + 1]/norm2[i]) * Z[, i - 1]
        Z <- Z/rep(sqrt(norm2[-1L]), each = length(x))
        colnames(Z) <- 0:degree
        Z <- Z[, -1, drop = FALSE]
        attr(Z, "degree") <- 1L:degree
        attr(Z, "coefs") <- list(alpha = alpha, norm2 = norm2)
        class(Z) <- c("poly", "matrix")
    }

    if (smart.mode.is("write"))
        put.smart(list(degree = degree,
                       coefs = attr(Z, "coefs"),
                       raw = FALSE,  # raw is changed above
                       match.call = match.call()))

    Z
}
attr(poly, "smart") <- TRUE






scale.default <-
function (x, center = TRUE, scale = TRUE) 
{
    x <- as.matrix(x)

    if (smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    nc <- ncol(x)
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(x, na.rm = TRUE)
            x <- sweep(x, 2L, center, check.margin = FALSE)
        }
    } else if (is.numeric(center) && (length(center) == nc)) 
        x <- sweep(x, 2L, center, check.margin = FALSE) else
    stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
        if (scale) {
            f <- function(v) {
                v <- v[!is.na(v)]
                sqrt(sum(v^2)/max(1, length(v) - 1L))
            }
            scale <- apply(x, 2L, f)
            x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
        }
    } else if (is.numeric(scale) && length(scale) == nc) 
        x <- sweep(x, 2L, scale, "/", check.margin = FALSE) else 
    stop("length of 'scale' must equal the number of columns of 'x'")
    if (is.numeric(center)) 
        attr(x, "scaled:center") <- center
    if (is.numeric(scale)) 
        attr(x, "scaled:scale") <- scale

    if (smart.mode.is("write")) {
        put.smart(list(center = center, scale = scale,
                       match.call = match.call()))
    }

    x
}
attr(scale.default, "smart") <- TRUE


attr(scale, "smart") <- TRUE





"my1" <- function(x, minx=min(x)) {

    x <- x   # Evaluate x

    if (smart.mode.is("read")) {
        smart  <- get.smart()
        minx <- smart$minx          # Overwrite its value 
    } else 
    if (smart.mode.is("write"))
        put.smart(list(minx=minx))

    (x-minx)^2
}
attr(my1, "smart") <- TRUE




"my2" <- function(x, minx=min(x)) {

    x <- x   # Evaluate x

    if (smart.mode.is("read")) {
        return(eval(smart.expression))
    } else 
    if (smart.mode.is("write"))
        put.smart(list(minx=minx, match.call=match.call()))

    (x-minx)^2
}

attr(my2, "smart") <- TRUE




"stdze1" <- function(x, center=TRUE, scale=TRUE) {

    x <- x  # Evaluate x

    if (!is.vector(x))
        stop("x must be a vector")

    if (smart.mode.is("read")) {
        smart  <- get.smart()
        return((x-smart$center)/smart$scale)
    }

    if (is.logical(center))
        center <- if (center) mean(x) else 0
    if (is.logical(scale))
        scale <- if (scale) sqrt(var(x)) else 1

    if (smart.mode.is("write"))
        put.smart(list(center=center,
                       scale=scale))
    # Normal use
    (x-center)/scale
}
attr(stdze1, "smart") <- TRUE

"stdze2" <- function(x, center=TRUE, scale=TRUE) {

    x <- x  # Evaluate x

    if (!is.vector(x))
        stop("x must be a vector")

    if (smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    if (is.logical(center))
        center <- if (center) mean(x) else 0
    if (is.logical(scale))
        scale <- if (scale) sqrt(var(x)) else 1

    if (smart.mode.is("write"))
        put.smart(list(center=center,
                       scale=scale,
                       match.call=match.call()))

    (x-center)/scale
}
attr(stdze2, "smart") <- TRUE




