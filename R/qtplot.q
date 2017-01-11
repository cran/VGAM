# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.














qtplot.lms.bcn <- function(percentiles = c(25, 50, 75),
                           eta = NULL, yoffset = 0) {

  lp <- length(percentiles)
  answer <- matrix(NA_real_, nrow(eta), lp,
                   dimnames = list(dimnames(eta)[[1]],
                   paste(as.character(percentiles), "%", sep = "")))
  for (ii in 1:lp) {
    answer[, ii] <- qlms.bcn(p      = percentiles[ii]/100,
                             lambda = eta[, 1],
                             mu     = eta[, 2],
                             sigma  = eta[, 3])
  }
  answer
}



qtplot.lms.bcg <- function(percentiles = c(25,50,75),
                           eta = NULL, yoffset = 0) {

  cc <- percentiles
  lp <- length(percentiles)
  answer <- matrix(NA_real_, nrow(eta), lp,
                   dimnames = list(dimnames(eta)[[1]],
                   paste(as.character(percentiles), "%", sep = "")))
  lambda <- eta[, 1]
  sigma <- eta[, 3]
  shape <- 1 / (lambda * sigma)^2
  for (ii in 1:lp) {
    ccc <- rep_len(cc[ii]/100, nrow(eta))
    ccc <- ifelse(lambda > 0, ccc, 1-ccc)
    answer[, ii] <- eta[, 2] *
                    (qgamma(ccc, shape = shape)/shape)^(1/lambda)
  }
  answer
}


qtplot.lms.yjn2 <-
qtplot.lms.yjn <- function(percentiles = c(25,50,75),
                           eta = NULL, yoffset = 0) {

  cc <- percentiles
  lp <- length(percentiles)
  answer <- matrix(NA_real_, nrow(eta), lp,
                   dimnames = list(dimnames(eta)[[1]],
                   paste(as.character(percentiles), "%", sep = "")))
  lambda <- eta[, 1]
  mu <- eta[, 2]
  sigma <- eta[, 3]  # Link function already taken care of above
  for (ii in 1:lp) {
    ccc <- mu + sigma * qnorm(cc[ii]/100)
    answer[, ii] <- yeo.johnson(ccc, lambda, inverse= TRUE) - yoffset
  }
  answer
}


qtplot.default <- function(object, ...) {

  warning("no methods function. Returning the object")
  invisible(object)
}



"qtplot.vglm" <- function(object, Attach= TRUE, ...) {

  LL <- length(object@family@vfamily)
  newcall <- paste("qtplot.", object@family@vfamily[LL],
                   "(object, ...)", sep = "")
  newcall <- parse(text = newcall)[[1]]

  if (Attach) {
    object@post$qtplot <- eval(newcall)
    invisible(object)
  } else
    eval(newcall)
}


qtplot.lmscreg <- function(object,
                           newdata = NULL,
                           percentiles = object@misc$percentiles,
                           show.plot = TRUE, ...) {

  same <- length(percentiles) == length(object@misc$percentiles) &&
          all(percentiles == object@misc$percentiles)

  lp <- length(percentiles)
  if (same) {
    fitted.values <- if (!length(newdata))
      object@fitted.values else {
      predict(object, newdata = newdata, type = "response")
    }
    fitted.values <- as.matrix(fitted.values)
  } else {
    if (!is.numeric(percentiles))
        stop("'percentiles' must be specified")

    eta <- if (length(newdata))
             predict(object, newdata = newdata, type = "link") else
             object@predictors


    if (!length(double.check.earg <- object@misc$earg))
      double.check.earg <- list(theta = NULL)
      eta  <- eta2theta(eta, link = object@misc$link,
                        earg = double.check.earg)  # lambda, mu, sigma



    if (!is.logical(expectiles <- object@misc$expectiles)) {
      expectiles <- FALSE
    }

    newcall <- paste(if (expectiles) "explot." else "qtplot.",
                    object@family@vfamily[1],
                    "(percentiles = percentiles",
                    ", eta = eta, yoffset=object@misc$yoffset)",
                    sep = "")
    newcall <- parse(text = newcall)[[1]]
    fitted.values <- as.matrix( eval(newcall) )
    dimnames(fitted.values) <-
      list(dimnames(eta)[[1]],
           paste(as.character(percentiles), "%", sep = ""))
  }

  if (show.plot) {
    plotqtplot.lmscreg(fitted.values = fitted.values,
                       object = object,
                       newdata = newdata,
                       lp = lp,
                       percentiles = percentiles, ...)
  }

  list(fitted.values = fitted.values, percentiles = percentiles)
}



plotqtplot.lmscreg <-
  function(fitted.values, object,
           newdata = NULL,
           percentiles = object@misc$percentiles,
           lp = NULL,
           add.arg = FALSE,
           y = if (length(newdata)) FALSE else TRUE,
           spline.fit = FALSE,
           label = TRUE,
           size.label = 0.06,
           xlab = NULL, ylab = "",
           pch = par()$pch, pcex = par()$cex,
           pcol.arg = par()$col,
           xlim = NULL, ylim = NULL,
           llty.arg = par()$lty,
           lcol.arg = par()$col, llwd.arg = par()$lwd,
           tcol.arg = par()$col,
           tadj = 1, ...) {



  if (!length(newdata)) {
    X <- model.matrixvlm(object, type = "lm")
    if (is.matrix(X) && length(object@y) && ncol(X)==2 &&
       dimnames(X)[[2]][1] == "(Intercept)") {
      xx <- X[, 2]
      if (is.null(xlab)) {
        xlab <- if (object@misc$nonparametric)
                as.vector(slot(object, "s.xargument")) else
                names(object@assign)[2]
      }

      if (!add.arg) {
        if (!is.numeric(xlim))
          xlim <- if (label)
              c(min(xx), max(xx) + size.label*diff(range(xx))) else
              c(min(xx), max(xx))
        fred <- cbind(object@y, fitted.values)
        if (!is.numeric(ylim))
          ylim <- c(min(fred), max(fred))
        matplot(x = xx, y = fred,
                xlab = xlab, ylab = ylab, type = "n",
                xlim = xlim, ylim = ylim, ...)
      }

      if (y && length(object@y))
        matpoints(x = xx, y = object@y, pch = pch, cex = pcex,
                  col = pcol.arg)
    } else {
      warning("there is not a single covariate. ",
              "Returning the object.")
      return(fitted.values)
    }
  } else {

    firstterm <- attr(terms(object), "term.labels")[1]

    if (object@misc$nonparametric &&
       length(object@s.xargument[firstterm]))
      firstterm <-  object@s.xargument[firstterm]

    xx <- newdata[[firstterm]]
    if (!is.numeric(xx))
      stop("couldn't extract the 'primary' variable from newdata")

    if (!add.arg) {
      if (is.null(xlab))
        xlab <- firstterm
      if (!is.numeric(xlim))
        xlim <- if (label)
            c(min(xx), max(xx)+size.label*diff(range(xx))) else
            c(min(xx), max(xx))
      if (!is.numeric(ylim))
        ylim <- c(min(fitted.values), max(fitted.values))
      matplot(x = xx, y = fitted.values,
              xlab = xlab, ylab = ylab, type = "n",
                  xlim = xlim, ylim = ylim, col = pcol.arg)
    }
    if (y && length(object@y))
      matpoints(x = xx, y = object@y, pch = pch, cex = pcex,
                col = pcol.arg)

  }

    tcol.arg <- rep_len(tcol.arg, lp)
    lcol.arg <- rep_len(lcol.arg, lp)
    llwd.arg <- rep_len(llwd.arg, lp)
    llty.arg <- rep_len(llty.arg, lp)
    for (ii in 1:lp) {
      temp <- cbind(xx, fitted.values[, ii])
      temp <- temp[sort.list(temp[, 1]), ]
      index <- !duplicated(temp[, 1])
      if (spline.fit) {
        lines(spline(temp[index, 1], temp[index, 2]),
              lty = llty.arg[ii], col = lcol.arg[ii], err = -1,
              lwd = llwd.arg[ii])
      } else {
        lines(temp[index, 1], temp[index, 2],
              lty = llty.arg[ii], col = lcol.arg[ii], err = -1,
              lwd = llwd.arg[ii])
      }
      if (label)
        text(par()$usr[2], temp[nrow(temp), 2],
             paste( percentiles[ii], "%", sep = ""),
             adj = tadj, col = tcol.arg[ii], err = -1)
    }

    invisible(fitted.values)
}


if (TRUE) {
  if (!isGeneric("qtplot"))
    setGeneric("qtplot", function(object, ...)
    standardGeneric("qtplot"))


  setMethod("qtplot", signature(object = "vglm"),
            function(object, ...)
            invisible(qtplot.vglm(object, ...)))
  setMethod("qtplot", signature(object = "vgam"),
            function(object, ...)
            invisible(qtplot.vglm(object, ...)))
}







"qtplot.vextremes" <- function(object, ...) {


  newcall <- paste("qtplot.", object@family@vfamily[1],
                  "(object = object, ... )", sep = "")
  newcall <- parse(text = newcall)[[1]]
  eval(newcall)
}


qtplot.gumbelff <-
qtplot.gumbel <-
    function(object, show.plot = TRUE, y.arg = TRUE,
             spline.fit = FALSE, label = TRUE,
             R = object@misc$R,
             percentiles = object@misc$percentiles,
             add.arg = FALSE,
             mpv = object@misc$mpv,
             xlab = NULL, ylab = "", main = "",
             pch = par()$pch, pcol.arg = par()$col,
             llty.arg = par()$lty, lcol.arg = par()$col,
             llwd.arg = par()$lwd,
             tcol.arg = par()$col, tadj = 1, ...) {
  if (!is.logical(mpv) || length(mpv) != 1)
    stop("bad input for 'mpv'")
  if (!length(percentiles) ||
     (!is.Numeric(percentiles, positive = TRUE) ||
      max(percentiles) >= 100))
    stop("bad input for 'percentiles'")



  eta <- predict(object)


  if (is.Numeric(R))
    R <- rep_len(R, nrow(eta))

  if (!is.Numeric(percentiles))
    stop("the 'percentiles' argument needs to be assigned a value")


  extra <- object@extra
  extra$mpv <- mpv  # Overwrite if necessary
  extra$R <- R
  extra$percentiles <- percentiles
  fitted.values <- object@family@linkinv(eta = eta, extra = extra)

  answer <- list(fitted.values = fitted.values,
                 percentiles   = percentiles)

  if (!show.plot)
    return(answer)



  lp <- length(percentiles)  # Does not include mpv
  tcol.arg <- rep_len(tcol.arg, lp+mpv)
  lcol.arg <- rep_len(lcol.arg, lp+mpv)
  llwd.arg <- rep_len(llwd.arg, lp+mpv)
  llty.arg <- rep_len(llty.arg, lp+mpv)

  X <- model.matrixvlm(object, type = "lm")
  if (is.matrix(X) && length(object@y) && ncol(X)==2 &&
     dimnames(X)[[2]][1] == "(Intercept)") {
    xx <- X[, 2]
    if (!length(xlab))
      xlab <- if (object@misc$nonparametric &&
                 length(object@s.xargument))
                  object@s.xargument else names(object@assign)[2]

    if (!add.arg)
      matplot(x = xx, y = cbind(object@y, fitted.values), main = main,
              xlab = xlab, ylab = ylab, type = "n", ...)

    if (y.arg) {
       matpoints(x = xx, y = object@y, pch = pch, col = pcol.arg)
    }
  } else {
    warning("there is not a single covariate.")
    return(answer)
  }

    for (ii in 1:(lp+mpv)) {
      temp <- cbind(xx, fitted.values[, ii])
      temp <- temp[sort.list(temp[, 1]), ]
      index <- !duplicated(temp[, 1])
      if (spline.fit) {
        lines(spline(temp[index, 1], temp[index, 2]),
              lty = llty.arg[ii], col = lcol.arg[ii], lwd = llwd.arg[ii])
      } else {
        lines(temp[index, 1], temp[index, 2],
              lty = llty.arg[ii], col = lcol.arg[ii], lwd = llwd.arg[ii])
      }
      if (label) {
        mylabel <- (dimnames(answer$fitted)[[2]])[ii]
        text(par()$usr[2], temp[nrow(temp), 2],
             mylabel, adj = tadj, col = tcol.arg[ii], err = -1,
             cex = par()$cex.axis, xpd = par()$xpd)
    }
  }

  invisible(answer)
}





deplot.lms.bcn <- function(object,
                           newdata,
                           y.arg,
                           eta0) {
  if (!any(object@family@vfamily == "lms.bcn"))
    warning("I think you've called the wrong function")

  Zvec <- ((y.arg/eta0[, 2])^(eta0[, 1]) -1) / (eta0[, 1] * eta0[, 3])
  dZ.dy <- ((y.arg/eta0[, 2])^(eta0[, 1]-1)) / (eta0[, 2] * eta0[, 3])
  yvec <- dnorm(Zvec) * abs(dZ.dy)

  list(newdata = newdata, y = y.arg, density = yvec)
}



deplot.lms.bcg <- function(object,
                           newdata,
                           y.arg,
                           eta0) {
  if (!any(object@family@vfamily == "lms.bcg"))
    warning("I think you've called the wrong function")

  Zvec <- (y.arg/eta0[, 2])^(eta0[, 1])  # different from lms.bcn
  dZ.dy <- ((y.arg/eta0[, 2])^(eta0[, 1]-1)) * eta0[, 1] / eta0[, 2]
  lambda <- eta0[, 1]
  sigma <- eta0[, 3]
  shape <- 1 / (lambda * sigma)^2
  yvec <- dgamma(Zvec, shape = shape, rate = shape) * abs(dZ.dy)

  list(newdata = newdata, y = y.arg, density = yvec)
}


deplot.lms.yjn2 <-
deplot.lms.yjn <- function(object,
                           newdata,
                           y.arg,
                           eta0) {

  if (!length(intersect(object@family@vfamily, c("lms.yjn","lms.yjn2"))))
    warning("I think you've called the wrong function")

  lambda <- eta0[, 1]
  Zvec <- (yeo.johnson(y.arg+object@misc$yoffset, lambda = eta0[, 1]) -
               eta0[, 2]) / eta0[, 3]
  dZ.dy <- dyj.dy.yeojohnson(y.arg+object@misc$yoffset,
                             lambda = eta0[, 1]) / eta0[, 3]
  yvec <- dnorm(Zvec) * abs(dZ.dy)

  list(newdata = newdata, y = y.arg, density = yvec)
}


deplot.default <- function(object, ...) {

  warning("no methods function. Returning the object")
  invisible(object)
}




"deplot.vglm" <- function(object, Attach= TRUE, ...) {
  LL <- length(object@family@vfamily)
  newcall <- paste("deplot.", object@family@vfamily[LL],
                   "(object, ...)", sep = "")
  newcall <- parse(text = newcall)[[1]]

  if (Attach) {
    object@post$deplot <- eval(newcall)
    invisible(object)
  } else {
    eval(newcall)
  }
}



"deplot.lmscreg" <- function(object,
                             newdata = NULL,
                             x0,
                             y.arg, show.plot = TRUE, ...) {


  if (!length(newdata)) {
    newdata <- data.frame(x0=x0)
    var1name <- attr(terms(object), "term.labels")[1]
    names(newdata) <- var1name

    ii <- if (object@misc$nonparametric)
          slot(object, "s.xargument") else NULL
    if (length(ii) && any(logic.vec <-
        names(slot(object, "s.xargument")) == var1name))
      names(newdata) <- ii[logic.vec]  # should be the first one
  }

  eta0 <- if (length(newdata)) predict(object, newdata) else
                               predict(object)

  if (!length(double.check.earg <- object@misc$earg))
    double.check.earg <- list(theta = NULL)
  eta0 <- eta2theta(eta0, link = object@misc$link,
                    earg = double.check.earg)  # lambda, mu, sigma

  newcall <- paste("deplot.", object@family@vfamily[1],
                   "(object, newdata, y.arg = y.arg, eta0 = eta0)",
                   sep = "")
  newcall <- parse(text = newcall)[[1]]
  answer <- eval(newcall)

  if (show.plot)
    plotdeplot.lmscreg(answer, y.arg=y.arg, ...)

  invisible(answer)
}



plotdeplot.lmscreg <- function(answer,
                           y.arg,
                           add.arg= FALSE,
                           xlab = "", ylab = "density",
                           xlim = NULL, ylim = NULL,
                           llty.arg = par()$lty, col.arg = par()$col,
                           llwd.arg = par()$lwd, ...) {

  yvec <- answer$density
  xx <- y.arg

  if (!add.arg) {
    if (!is.numeric(xlim))
      xlim <- c(min(xx), max(xx))
    if (!is.numeric(ylim))
      ylim <- c(min(yvec), max(yvec))
    matplot(x = xx, y = yvec,
            xlab = xlab, ylab = ylab, type = "n",
            xlim = xlim, ylim = ylim, ...)
  }

  temp <- cbind(xx, yvec)
  temp <- temp[sort.list(temp[, 1]), ]
  index <- !duplicated(temp[, 1])
  lines(temp[index, 1], temp[index, 2],
        lty = llty.arg, col = col.arg, err = -1, lwd = llwd.arg)

  invisible(answer)
}




if (TRUE) {

  if (!isGeneric("deplot"))
  setGeneric("deplot", function(object, ...) standardGeneric("deplot"))

  setMethod("deplot", signature(object = "vglm"),
            function(object, ...)
            invisible(deplot.vglm(object, ...)))
  setMethod("deplot", signature(object = "vgam"),
            function(object, ...)
            invisible(deplot.vglm(object, ...)))
}




if (TRUE) {

  if (!isGeneric("cdf"))
    setGeneric("cdf", function(object, ...) standardGeneric("cdf"))

  setMethod("cdf", signature(object = "vglm"),
            function(object, ...)
            cdf.vglm(object, ...))

  setMethod("cdf", signature(object = "vgam"),
            function(object, ...)
            cdf.vglm(object, ...))
}


"cdf.vglm" <- function(object, newdata = NULL, Attach = FALSE, ...) {
  LL <- length(object@family@vfamily)
  newcall <- paste("cdf.", object@family@vfamily[LL],
                  "(object, newdata, ...)", sep = "")
  newcall <- parse(text = newcall)[[1]]

  if (Attach) {
    object@post$cdf <- eval(newcall)
    object
  } else {
    eval(newcall)
  }
}



"cdf.lmscreg" <- function(object,
                          newdata = NULL, ...) {



  if (!length(newdata))
    return(object@post$cdf)

  eta0 <- if (length(newdata)) predict(object, newdata) else predict(object)

  if (!length(double.check.earg <- object@misc$earg))
    double.check.earg <- list(theta = NULL)
  eta0 <- eta2theta(eta0, link = object@misc$link,
                    earg = double.check.earg)  # lambda, mu, sigma

  y <- vgety(object, newdata)   # Includes yoffset

  newcall <- paste("cdf.", object@family@vfamily[1],
                   "(y, eta0, ... )", sep = "")
  newcall <- parse(text = newcall)[[1]]
  eval(newcall)
}



cdf.lms.bcn <- function(y, eta0) {
    Zvec <- ((y/eta0[, 2])^(eta0[, 1]) -1) / (eta0[, 1] * eta0[, 3])
    Zvec[abs(eta0[, 3]) < 1e-5] <- log(y/eta0[, 2]) / eta0[, 3]
    ans <- c(pnorm(Zvec))
    names(ans) <- dimnames(eta0)[[1]]
    ans
}


cdf.lms.bcg <- function(y, eta0) {
    shape <- 1 / (eta0[, 1] * eta0[, 3])^2
    Gvec <- shape * (y/eta0[, 2])^(eta0[, 1])
    ans <- c(pgamma(Gvec, shape = shape))
    ans[eta0[, 1] < 0] <- 1-ans
    names(ans) <- dimnames(eta0)[[1]]
    ans
}


cdf.lms.yjn <- function(y, eta0) {


  Zvec <- (yeo.johnson(y, eta0[, 1]) - eta0[, 2])/eta0[, 3]
  ans <- c(pnorm(Zvec))
  names(ans) <- dimnames(eta0)[[1]]
  ans
}


vgety <- function(object, newdata = NULL) {

  y <- if (length(newdata)) {
    yname <- dimnames(attr(terms(object@terms),"factors"))[[1]][1]
    newdata[[yname]]
  } else {
    object@y
  }
  if (length(object@misc$yoffset))
    y <- y + object@misc$yoffset
  y
}






"rlplot.vglm" <- function(object, Attach = TRUE, ...) {

  LL <- length(object@family@vfamily)
  newcall <- paste("rlplot.", object@family@vfamily[LL],
                   "(object, ...)", sep = "")
  newcall <- parse(text = newcall)[[1]]

  if (Attach) {
    object@post$rlplot <- eval(newcall)
    invisible(object)
  } else {
    eval(newcall)
  }
}





"rlplot.vextremes" <- function(object, ...) {



  newcall <- paste("rlplot.", object@family@vfamily[1],
                   "(object = object, ... )", sep = "")
  newcall <- parse(text = newcall)[[1]]
  eval(newcall)
}



rlplot.gevff <-
rlplot.gev <-
  function(object, show.plot = TRUE,
           probability = c((1:9)/100, (1:9)/10, 0.95, 0.99, 0.995, 0.999),
           add.arg = FALSE,
           xlab = if(log.arg) "Return Period (log-scale)" else "Return Period",
           ylab = "Return Level",
           main = "Return Level Plot",
           pch = par()$pch, pcol.arg = par()$col, pcex = par()$cex,
           llty.arg = par()$lty, lcol.arg = par()$col, llwd.arg = par()$lwd,
           slty.arg = par()$lty, scol.arg = par()$col, slwd.arg = par()$lwd,
           ylim = NULL,
           log.arg = TRUE,
           CI = TRUE,
           epsilon = 1.0e-05,
           ...) {

  if (!is.Numeric(epsilon, length.arg = 1) ||
      abs(epsilon) > 0.10)
    stop("bad input for 'epsilon'")
  if (!is.Numeric(probability, positive = TRUE) ||
      max(probability) >= 1 ||
      length(probability) < 5)
    stop("bad input for 'probability'")

  if (!is.logical(log.arg) || length(log.arg) != 1)
    stop("bad input for argument 'log'")
  if (!is.logical(CI) || length(CI) != 1)
    stop("bad input for argument 'CI'")

  if (!object@misc$intercept.only)
    stop("object must be an intercept-only fit, ",
         "i.e., y ~ 1 is the response")

  extra2 <- object@extra
  extra2$percentiles <- 100 * probability  # Overwrite
  zp <- object@family@linkinv(eta = predict(object)[1:2, ],
                              extra = extra2)[1, ]
  yp <- -log(probability)
  ydata <- sort(object@y[, 1])
  n <- object@misc$n
  if (log.arg) {
    if (!add.arg)
      plot(log(1/yp), zp, log = "", type = "n",
           ylim = if (length(ylim)) ylim else
                c(min(c(ydata, zp)), max(c(ydata, zp))),
           xlab = xlab, ylab = ylab, main = main,
           cex.axis = par()$cex.axis,
           cex.main = par()$cex.main,
           cex.lab  = par()$cex.lab,
           ...)
    points(log(-1/log((1:n)/(n+1))), ydata, col = pcol.arg,
           pch = pch, cex = pcex)
    lines(log(1/yp), zp,
          lwd = llwd.arg, col = lcol.arg, lty = llty.arg)
  } else {
    if (!add.arg)
      plot(1/yp, zp, log = "x", type = "n",
           ylim = if (length(ylim)) ylim else
                  c(min(c(ydata, zp)),
                    max(c(ydata, zp))),
           xlab = xlab, ylab = ylab, main = main,
           cex.axis = par()$cex.axis,
           cex.main = par()$cex.main,
           cex.lab  = par()$cex.lab,
           ...)
    points(-1/log((1:n)/(n+1)), ydata, col = pcol.arg,
           pch = pch, cex = pcex)
    lines(1/yp, zp, lwd = llwd.arg, col = lcol.arg, lty = llty.arg)
  }

  if (CI) {
    zpp <- cbind(zp, zp, zp)  # lp x 3
    eta <- predict(object)
    Links <- object@misc$link
    earg <- object@misc$earg
    M <- object@misc$M
    for (ii in 1:M) {
      TTheta <- eta[, ii]
      use.earg <- earg[[ii]]
      newcall <- paste(Links[ii],
                "(theta = TTheta, ",
                "  inverse = TRUE)",
                 sep = "")



      newcall <- parse(text = newcall)[[1]]
      uteta <- eval(newcall)  # Theta, the untransformed parameter
      uteta <- uteta + epsilon  # Perturb it
      newcall <- paste(Links[ii],
                       "(theta = uteta",
                       ")",
                       sep = "")
      newcall <- parse(text = newcall)[[1]]
      teta <- eval(newcall)  # The transformed parameter
      peta <- eta
      peta[, ii] <- teta
      zpp[, ii] <- object@family@linkinv(eta = peta,
                                         extra = extra2)[1, ]
      zpp[, ii] <- (zpp[, ii] - zp) / epsilon  # On the transformed scale
    }




    VCOV <- vcov(object, untransform = TRUE)





    vv <- numeric(nrow(zpp))
    for (ii in 1:nrow(zpp))
      vv[ii] <- t(as.matrix(zpp[ii, ])) %*% VCOV %*% as.matrix(zpp[ii, ])
    if (log.arg) {
      lines(log(1/yp), zp - 1.96 * sqrt(vv),
            lwd = slwd.arg, col = scol.arg, lty = slty.arg)
      lines(log(1/yp), zp + 1.96 * sqrt(vv),
            lwd = slwd.arg, col = scol.arg, lty = slty.arg)
    } else {
      lines(1/yp, zp - 1.96 * sqrt(vv),
            lwd = slwd.arg, col = scol.arg, lty = slty.arg)
      lines(1/yp, zp + 1.96 * sqrt(vv),
            lwd = slwd.arg, col = scol.arg, lty = slty.arg)
    }
  }
  answer <- list(yp = yp,
                zp = zp)
  if (CI) {
    answer$lower <- zp - 1.96 * sqrt(vv)
    answer$upper <- zp + 1.96 * sqrt(vv)
  }
  invisible(answer)
}


if (!isGeneric("rlplot"))
  setGeneric("rlplot",
             function(object, ...) standardGeneric("rlplot"))

setMethod("rlplot",  "vglm", function(object, ...)
          rlplot.vglm(object, ...))







explot.lms.bcn <- function(percentiles = c(25, 50, 75),
                           eta = NULL, yoffset = 0) {

  lp <- length(percentiles)
  answer <- matrix(NA_real_, nrow(eta), lp,
                   dimnames = list(dimnames(eta)[[1]],
                   paste(as.character(percentiles), "%", sep = "")))
  for (ii in 1:lp) {
    answer[, ii] <- eta[, 2] * (1 + eta[, 1] * eta[, 3] *
                    qenorm(percentiles[ii]/100))^(1/eta[, 1])
  }
  answer
}





