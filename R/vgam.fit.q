# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.





vgam.fit <-
  function(x, y, w, mf,
           Xm2 = NULL, Ym2 = NULL,  # Added 20130730
           etastart, mustart, coefstart,
           offset, family, control, criterion = "coefficients",
           constraints = NULL, extra, qr.arg,
           Terms,
           nonparametric, smooth.labels,
           function.name = "vgam", ...) {


  eff.n <- nrow(x)  # + sum(abs(w[1:nrow(x)]))

  specialCM <- NULL
  post <- list()
  check.Rank <- TRUE  # Set this to false for family functions vppr() etc.
  epsilon <- control$epsilon
  maxit <- control$maxit
  save.weight <- control$save.weight
  trace <- control$trace

  bf.maxit <- control$bf.maxit
  bf.epsilon <- control$bf.epsilon
  se.fit <- control$se.fit
  minimize.criterion <- control$min.criterion


  fv <- NULL



  n <- dim(x)[1]




  new.s.call <- expression({
    if (c.list$one.more) {
      fv <- c.list$fit
      new.coeffs <- c.list$coeff

      if (length(family@middle))
        eval(family@middle)

      eta <- fv + offset
      mu <- family@linkinv(eta, extra)

      if (length(family@middle2))
        eval(family@middle2)

      old.crit <- new.crit

      new.crit <- switch(criterion,
                         coefficients = new.coeffs,
                         tfun(mu = mu, y = y, w = w,
                              res = FALSE, eta = eta, extra))
      if (trace) {
        cat("VGAM ", bf, " loop ", iter, ": ", criterion, "= ")

        UUUU <- switch(criterion,
                       coefficients =
                         format(new.crit,
                                dig = round(1 - log10(epsilon))),
                         format(new.crit, 
                                dig = max(4, 
                                          round(-0 - log10(epsilon) +
                                          log10(sqrt(eff.n))))))

        switch(criterion,
               coefficients = {if (length(new.crit) > 2) cat("\n");
               cat(UUUU, fill = TRUE, sep = ", ")},
               cat(UUUU, fill = TRUE, sep = ", "))
      }

      one.more <- eval(control$convergence)

      flush.console()

      if (!is.finite(one.more) || !is.logical(one.more))
        one.more <- FALSE
      if (one.more) {
        iter <- iter + 1
        deriv.mu <- eval(family@deriv)
        wz <- eval(family@weight)
        if (control$checkwz)
          wz <- checkwz(wz, M = M, trace = trace,
                        wzepsilon = control$wzepsilon)

        U <- vchol(wz, M = M, n = n, silent = !trace)
        tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
        z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

        c.list$z <- z
        c.list$wz <- wz
        c.list$U <- U
      }

      c.list$one.more <- one.more
      c.list$coeff <- runif(length(new.coeffs))  # 20030312; twist needed!
      old.coeffs <- new.coeffs
    }
    c.list
  })





  old.coeffs <- coefstart

  intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
  y.names <- predictors.names <- NULL  # May be overwritten in @initialize

  n.save <- n
  if (length(slot(family, "initialize")))
    eval(slot(family, "initialize"))  # Initialize mu & M (& optionally w)

  if (length(etastart)) {
    eta <- etastart
    mu <- if (length(mustart)) mustart else
          if (length(body(slot(family, "linkinv"))))
            slot(family, "linkinv")(eta, extra) else
            warning("argument 'etastart' assigned a value ",
                    "but there is no 'linkinv' slot to use it")
  }

  if (length(mustart)) {
    mu <- mustart
    if (length(body(slot(family, "linkfun")))) {
      eta <- slot(family, "linkfun")(mu, extra)
    } else {
      warning("argument 'mustart' assigned a value ",
              "but there is no 'link' slot to use it")
    }
  }

  M <- if (is.matrix(eta)) ncol(eta) else 1


  if (length(family@constraints))
    eval(family@constraints)
  Blist <- process.constraints(constraints, x, M, specialCM = specialCM)

  ncolBlist <- unlist(lapply(Blist, ncol))
  dimB <- sum(ncolBlist)


    if (nonparametric) {


      smooth.frame <- mf
      assignx <- attr(x, "assign")
      which <- assignx[smooth.labels]

      bf <- "s.vam"
      bf.call <- parse(text = paste(
              "s.vam(x, z, wz, tfit$smomat, which, tfit$smooth.frame,",
              "bf.maxit, bf.epsilon, trace, se = se.fit, X.vlm.save, ",
              "Blist, ncolBlist, M = M, qbig = qbig, Umat = U, ",
              "all.knots = control$all.knots, nk = control$nk)",
              sep = ""))[[1]]

      qbig <- sum(ncolBlist[smooth.labels])  # Number of component funs
      smomat <- matrix(0, n, qbig)
      dy <- if (is.matrix(y)) dimnames(y)[[1]] else names(y)
      d2 <- if (is.null(predictors.names))
              paste("(Additive predictor ",1:M,")", sep = "") else
              predictors.names
      dimnames(smomat) <- list(dy, vlabel(smooth.labels,
                                          ncolBlist[smooth.labels], M))

      tfit <- list(smomat = smomat, smooth.frame = smooth.frame)
    } else {
      bf.call <- expression(vlm.wfit(xmat = X.vlm.save, z,
                                     Blist = NULL, U = U,
                                     matrix.out = FALSE, is.vlmX = TRUE,
                                     qr = qr.arg, xij = NULL))
      bf <- "vlm.wfit"
    }

    X.vlm.save <- lm2vlm.model.matrix(x, Blist, xij = control$xij)


    if (length(coefstart)) {
      eta <- if (ncol(X.vlm.save) > 1) X.vlm.save %*% coefstart +
               offset else X.vlm.save * coefstart + offset
      eta <- if (M > 1) matrix(eta, ncol = M, byrow = TRUE) else c(eta)
      mu <- family@linkinv(eta, extra)
    }


    if (criterion != "coefficients") {
      tfun <- slot(family, criterion)  # Needed 4 R so have to follow suit
    }

    iter <- 1
    new.crit <- switch(criterion,
                       coefficients = 1,
                       tfun(mu = mu, y = y, w = w, res = FALSE,
                            eta = eta, extra))
    old.crit <- ifelse(minimize.criterion,  10 * new.crit + 10,
                                           -10 * new.crit - 10)

    deriv.mu <- eval(family@deriv)
    wz <- eval(family@weight)
    if (control$checkwz)
      wz <- checkwz(wz, M = M, trace = trace,
                    wzepsilon = control$wzepsilon)

    U <- vchol(wz, M = M, n = n, silent = !trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
    z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

    c.list <- list(wz = as.double(wz), z = as.double(z),
                   fit = as.double(t(eta)),
                   one.more = TRUE, U = as.double(U),
                   coeff = as.double(rep(1, ncol(X.vlm.save))))


    dX.vlm <- as.integer(dim(X.vlm.save))
    nrow.X.vlm <- dX.vlm[[1]]
    ncol.X.vlm <- dX.vlm[[2]]
    if (nrow.X.vlm < ncol.X.vlm)
      stop(ncol.X.vlm, " parameters but only ", nrow.X.vlm, " observations")

    while (c.list$one.more) {
      tfit <- eval(bf.call)  # fit$smooth.frame is new

      c.list$coeff <- tfit$coefficients

      tfit$predictors <- tfit$fitted.values + offset

      c.list$fit <- tfit$fitted.values
      c.list <- eval(new.s.call)
      NULL
    }

    if (maxit > 1 && iter >= maxit)
      warning("convergence not obtained in ", maxit, " iterations")


    dnrow.X.vlm <- labels(X.vlm.save)
    xnrow.X.vlm <- dnrow.X.vlm[[2]]
    ynrow.X.vlm <- dnrow.X.vlm[[1]]

    if (length(family@fini))
      eval(family@fini)

    coefs <- tfit$coefficients
    asgn <- attr(X.vlm.save, "assign")    # 20011129 was x 

    names(coefs) <- xnrow.X.vlm
    cnames <- xnrow.X.vlm

    if (!is.null(tfit$rank)) {
      rank <- tfit$rank
      if (rank < ncol(x)) 
        stop("rank < ncol(x) is bad")
    } else rank <- ncol(x)

    R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    attributes(R) <- list(dim = c(ncol.X.vlm, ncol.X.vlm),
                          dimnames = list(cnames, cnames), rank = rank)


    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]



    if (is.matrix(mu)) {
      if (length(dimnames(mu)[[2]])) {
        y.names <- dimnames(mu)[[2]]
      } else
      if (length(dimnames(y)[[2]])) {
        y.names <- dimnames(y)[[2]]
      }
      dimnames(mu) <- list(yn, y.names)
    } else {
      names(mu) <- names(fv)
    }

    tfit$fitted.values <- NULL  # Have to kill it off  3/12/01
    fit <- structure(c(tfit,
           list(assign = asgn,
                constraints = Blist,
                control = control,
                fitted.values = mu,
                formula = as.vector(attr(Terms, "formula")),
                iter = iter,
                offset = offset,
                rank = rank,
                R = R,
                terms = Terms)))

    df.residual <- nrow.X.vlm - rank 

    if (!se.fit) {
      fit$varmat <- NULL
    }

    if (M == 1) {
      wz <- as.vector(wz)  # Convert wz into a vector
    } # else
    fit$weights <- if (save.weight) wz else NULL



    if (M == 1) {
      fit$predictors <- as.vector(fit$predictors)
      fit$residuals <- as.vector(fit$residuals)
      names(fit$residuals) <- names(fit$predictors) <- yn
    } else {
      dimnames(fit$residuals)  <-
      dimnames(fit$predictors) <- list(yn, predictors.names)
    }

    NewBlist <- process.constraints(constraints, x, M,
                                    specialCM = specialCM, by.col = FALSE)

    misc <- list(
        colnames.x = xn,
        colnames.X.vlm = xnrow.X.vlm,
        criterion = criterion,
        function.name = function.name,
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = n,
        new.assign = new.assign(x, NewBlist),
        nonparametric = nonparametric,
        nrow.X.vlm = nrow.X.vlm,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        ncol.X.vlm = ncol.X.vlm,
        ynames = dimnames(y)[[2]])


    if (criterion != "coefficients")
      fit[[criterion]] <- new.crit



    if (se.fit && length(fit$s.xargument)) {
      misc$varassign <- varassign(Blist, names(fit$s.xargument))
    }



    if (nonparametric) {
      misc$smooth.labels <- smooth.labels
    }


    crit.list <- list()
    if (criterion != "coefficients")
      crit.list[[criterion]] <- fit[[criterion]] <- new.crit
    for (ii in names(.min.criterion.VGAM)) {
      if (ii != criterion &&
          any(slotNames(family) == ii) &&
          length(body(slot(family, ii)))) {
        fit[[ii]] <-
        crit.list[[ii]] <-
          (slot(family, ii))(mu = mu, y = y, w = w, res = FALSE,
                             eta = eta, extra)
        }
    }




    if (M == 1) {
      fit$predictors <- as.vector(fit$predictors)
      fit$residuals  <- as.vector(fit$residuals)
      names(fit$residuals)  <-
      names(fit$predictors) <- yn
    } else {
      dimnames(fit$residuals)  <-
      dimnames(fit$predictors) <- list(yn, predictors.names)
    }



 
    if (w[1] != 1 || any(w != w[1]))
      fit$prior.weights <- w

    if (length(family@last))
      eval(family@last)


    if (!is.null(fit$smomat)) {
      fit$nl.chisq <- vgam.nlchisq(fit$qr, fit$resid, wz = wz,
                                   smomat = fit$smomat,
                                   deriv = deriv.mu, U = U,
                                   smooth.labels, attr(x, "assign"),
                                   M = M, n = n, constraints = Blist)
    }


    if (!qr.arg) { 
      fit$qr <- NULL
    }




    fit$misc <- NULL  # 20020608; It is necessary to kill it as it
    structure(c(fit, list(
        contrasts = attr(x, "contrasts"),
        control = control,
        crit.list = crit.list,
        extra = extra,
        family = family,
        iter = iter,
        misc = misc,
        post = post,
        x = x,
        y = y)),
        vclass = family@vfamily)
}





new.assign <- function(X, Blist) {

  M <- nrow(Blist[[1]])
  dn <- labels(X)
  xn <- dn[[2]]

  asgn <- attr(X, "assign")
  nasgn <- names(asgn)
  lasgn <- unlist(lapply(asgn, length))

  ncolBlist <- unlist(lapply(Blist, ncol))
  names(ncolBlist) <- NULL  # This is necessary for below to work 

  temp2 <- vlabel(nasgn, ncolBlist, M)
  L <- length(temp2)
  newasgn <- vector("list", L)

  kk <- 0
  low <- 1
  for (ii in 1:length(asgn)) {
    len  <- low:(low  + ncolBlist[ii] * lasgn[ii] -1)
    temp <- matrix(len, ncolBlist[ii],  lasgn[ii])
    for (mm in 1:ncolBlist[ii])
      newasgn[[kk + mm]] <- temp[mm, ]
    low <- low + ncolBlist[ii] * lasgn[ii]
    kk <- kk + ncolBlist[ii]
  }

  names(newasgn) <- temp2
  newasgn
}




