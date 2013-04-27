# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.






vglm.fit <- function(x, y, w = rep(1, length(x[, 1])),
    X_vlm_arg = NULL,
    Xm2 = NULL, Ym2 = NULL,
    etastart = NULL, mustart = NULL, coefstart = NULL,
    offset = 0, family,
    control = vglm.control(),
    criterion = "coefficients",
    qr.arg = FALSE,
    constraints = NULL,
    extra = NULL,
    Terms = Terms, function.name = "vglm", ...) {

  specialCM <- NULL
  post <- list()
  check.rank <- TRUE # Set this to false for family functions vppr() etc.
  check.rank <- control$Check.rank
  nonparametric <- FALSE
  epsilon <- control$epsilon
  maxit <- control$maxit
  save.weight <- control$save.weight
  trace <- control$trace
  orig.stepsize <- control$stepsize
  minimize.criterion <- control$min.criterion


  fv <- NULL




  n <- dim(x)[1]

  new.s.call <- expression({
    if (c.list$one.more) {
      fv <- c.list$fit
      new.coeffs <- c.list$coeff

      if (length(slot(family, "middle")))
        eval(slot(family, "middle"))

      eta <- fv + offset
      mu <- slot(family, "linkinv")(eta, extra)

      if (length(slot(family, "middle2")))
        eval(slot(family, "middle2"))

      old.crit <- new.crit
      new.crit <- 
          switch(criterion,
              coefficients = new.coeffs,
              tfun(mu = mu, y = y, w = w,
                   res = FALSE, eta = eta, extra))


      if (trace && orig.stepsize == 1) {
        cat("VGLM    linear loop ", iter, ": ", criterion, "= ")
        UUUU <- switch(criterion,
                       coefficients = format(new.crit,
                                      dig = round(2 - log10(epsilon))),
                       format(round(new.crit, 4)))

        switch(criterion,
               coefficients = {if (length(new.crit) > 2) cat("\n");
               cat(UUUU, fill = TRUE, sep = ", ")},
               cat(UUUU, fill = TRUE, sep = ", "))
      }


      {
          take.half.step =
            (control$half.stepsizing &&
             length(old.coeffs)) &&
             ((orig.stepsize != 1) ||
             (criterion != "coefficients" &&
             (if (minimize.criterion) new.crit > old.crit else
             new.crit < old.crit)))
          if (!is.logical(take.half.step))
            take.half.step <- TRUE
          if (take.half.step) {
              stepsize <- 2 * min(orig.stepsize, 2*stepsize)
              new.coeffs.save <- new.coeffs
              if (trace) 
                cat("Taking a modified step")
              repeat {
                  if (trace) {
                    cat(".")
                    flush.console()
                  }
                  stepsize <- stepsize / 2
                  if (too.small <- stepsize < 0.001)
                    break
                  new.coeffs <- (1-stepsize)*old.coeffs +
                                 stepsize*new.coeffs.save

                  if (length(slot(family, "middle")))
                    eval(slot(family, "middle"))

                  fv <- X_vlm_save %*% new.coeffs
                  if (M > 1)
                    fv <- matrix(fv, n, M, byrow = TRUE)

                  eta <- fv + offset
                  mu <- slot(family, "linkinv")(eta, extra)

                  if (length(slot(family, "middle2")))
                    eval(slot(family, "middle2"))


                  new.crit <- 
                      switch(criterion,
                          coefficients = new.coeffs,
                          tfun(mu = mu, y = y, w = w,
                               res = FALSE, eta = eta, extra))

                  if ((criterion == "coefficients") || 
                     ( minimize.criterion && new.crit < old.crit) ||
                     (!minimize.criterion && new.crit > old.crit))
                      break
              } # of repeat

              if (trace) 
                cat("\n")
              if (too.small) {
                warning("iterations terminated because ",
                        "half-step sizes are very small")
                one.more <- FALSE
              } else {
                if (trace) {
                    cat("VGLM    linear loop ",
                        iter, ": ", criterion, "= ")

                    UUUU <- switch(criterion,
                                   coefficients =
                                     format(new.crit,
                                            dig = round(2-log10(epsilon))),
                                   format(round(new.crit, 4)))

                    switch(criterion,
                           coefficients = {
                           if (length(new.crit) > 2) cat("\n");
                           cat(UUUU, fill = TRUE, sep = ", ")},
                           cat(UUUU, fill = TRUE, sep = ", "))
                }

                one.more <- eval(control$convergence)
              }
          } else {
            one.more <- eval(control$convergence)
          }
      }
      flush.console()

      if (!is.logical(one.more))
        one.more <- FALSE
      if (one.more) {
        iter <- iter + 1
        deriv.mu <- eval(slot(family, "deriv"))
        wz <- eval(slot(family, "weight"))
        if (control$checkwz)
          wz <- checkwz(wz, M = M, trace = trace,
                        wzepsilon = control$wzepsilon)

        U <- vchol(wz, M = M, n = n, silent = !trace)
        tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
        z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

        c.list$z <- z
        c.list$U <- U
        if (copy_X_vlm)
          c.list$X_vlm <- X_vlm_save
      }

      c.list$one.more <- one.more
      c.list$coeff <- runif(length(new.coeffs)) # 12/3/03; twist needed!
      old.coeffs <- new.coeffs
    }
    c.list
  })





  copy_X_vlm <- FALSE    # May be overwritten in @initialize
  stepsize <- orig.stepsize
  old.coeffs <- coefstart

  intercept.only <- ncol(x) == 1 &&
                    dimnames(x)[[2]] == "(Intercept)"
  y.names <- predictors.names <- NULL # May be overwritten in @initialize

  n.save <- n 


  if (length(slot(family, "initialize")))
    eval(slot(family, "initialize")) # Initialize mu & M (& optionally w)


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



  if (length(slot(family, "constraints")))
    eval(slot(family, "constraints"))


  Blist <- process.constraints(constraints, x, M,
                               specialCM = specialCM)


  ncolBlist <- unlist(lapply(Blist, ncol))
  dimB <- sum(ncolBlist)





  X_vlm_save <- if (length(X_vlm_arg)) X_vlm_arg else
    lm2vlm.model.matrix(x, Blist, xij = control$xij,
                                  Xm2 = Xm2)





  if (length(coefstart)) {
    eta <- if (ncol(X_vlm_save)>1) X_vlm_save %*% coefstart +
             offset else X_vlm_save * coefstart + offset
    eta <- if (M > 1) matrix(eta, ncol = M, byrow = TRUE) else c(eta) 
    mu <- slot(family, "linkinv")(eta, extra)
  }


  if (criterion != "coefficients") {
    tfun <- slot(family, criterion)   # family[[criterion]]
  }

  iter <- 1
  new.crit <- switch(criterion,
                     coefficients = 1,
                     tfun(mu = mu, y = y, w = w,
                          res = FALSE, eta = eta, extra))
  old.crit <- if (minimize.criterion)
              10*new.crit+10 else
             -10*new.crit-10

  deriv.mu <- eval(slot(family, "deriv"))
  wz <- eval(slot(family, "weight"))
  if (control$checkwz)
    wz <- checkwz(wz, M = M, trace = trace,
                  wzepsilon = control$wzepsilon)

  U <- vchol(wz, M = M, n = n, silent = !trace)
  tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
  z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

  c.list <- list(z = as.double(z), fit = as.double(t(eta)),
                 one.more = TRUE,
                 coeff = as.double(rep(1, ncol(X_vlm_save))),
                 U = as.double(U),
                 copy_X_vlm = copy_X_vlm,
                 X_vlm = if (copy_X_vlm) as.double(X_vlm_save) else
                         double(3))


  dX_vlm <- as.integer(dim(X_vlm_save))
  nrow_X_vlm <- dX_vlm[[1]]
  ncol_X_vlm <- dX_vlm[[2]]

  if (nrow_X_vlm < ncol_X_vlm)
    stop(ncol_X_vlm, "parameters but only ", nrow_X_vlm,
         " observations")



  bf.call <- expression(vlm.wfit(xmat = X_vlm_save, z,
                                 Blist = NULL, U = U,
                                 matrix.out = FALSE,
                                 is.vlmX = TRUE,
                                 qr = qr.arg, xij = NULL))


  while(c.list$one.more) {
      tfit <- eval(bf.call)   # fit$smooth.frame is new
    
      c.list$coeff <- tfit$coefficients 
    
      tfit$predictors <- tfit$fitted.values
    
      c.list$fit <- tfit$fitted.values
      c.list <- eval(new.s.call)
      NULL
  }

  if (maxit > 1 && iter >= maxit && !control$noWarning)
    warning("convergence not obtained in ", maxit, " iterations")




  dnrow_X_vlm <- labels(X_vlm_save)
  xnrow_X_vlm <- dnrow_X_vlm[[2]]
  ynrow_X_vlm <- dnrow_X_vlm[[1]]

  if (length(slot(family, "fini")))
    eval(slot(family, "fini"))

  if (M > 1)
    tfit$predictors <- matrix(tfit$predictors, n, M)

  coefs <- tfit$coefficients
  asgn <- attr(X_vlm_save, "assign")

  names(coefs) <- xnrow_X_vlm

  rank <- tfit$rank
  cnames <- xnrow_X_vlm

  if (check.rank && rank < ncol_X_vlm)
    stop("vglm only handles full-rank models (currently)")

  R <- tfit$qr$qr[1:ncol_X_vlm, 1:ncol_X_vlm, drop = FALSE]
  R[lower.tri(R)] <- 0
  attributes(R) <- list(dim = c(ncol_X_vlm, ncol_X_vlm),
                        dimnames = list(cnames, cnames), rank=rank)

  effects <- tfit$effects
  neff <- rep("", nrow_X_vlm)
  neff[seq(ncol_X_vlm)] <- cnames
  names(effects) <- neff

  dim(tfit$predictors) <- c(n, M)
  dn <- labels(x)
  yn <- dn[[1]]
  xn <- dn[[2]]


  residuals <- z - tfit$predictors
  if (M == 1) {
    tfit$predictors <- as.vector(tfit$predictors)
    residuals <- as.vector(residuals)
    names(residuals) <- names(tfit$predictors) <- yn
  } else {
    dimnames(residuals) <- dimnames(tfit$predictors) <-
                           list(yn, predictors.names)
  }

  if (is.matrix(mu)) {
    if (length(dimnames(y)[[2]])) {
      y.names <- dimnames(y)[[2]]
    }
    if (length(dimnames(mu)[[2]])) {
      y.names <- dimnames(mu)[[2]]
    }
    dimnames(mu) <- list(yn, y.names)
  } else {
    names(mu) <- names(fv)
  }


  df.residual <- nrow_X_vlm - rank
  fit <- list(assign = asgn,
              coefficients = coefs,
              constraints = Blist, 
              df.residual = df.residual,
              df.total = n*M,
              effects = effects, 
              fitted.values = mu,
              offset = offset, 
              rank = rank,
              residuals = residuals,
              R = R,
              terms = Terms) # terms: This used to be done in vglm() 

  if (qr.arg) {
    fit$qr <- tfit$qr
    dimnames(fit$qr$qr) <- dnrow_X_vlm
  }

  if (M == 1) {
    wz <- as.vector(wz)  # Convert wz into a vector
  } # else
  fit$weights <- if (save.weight) wz else NULL


  misc <- list(
      colnames.x = xn,
      colnames.X_vlm = xnrow_X_vlm,
      criterion = criterion,
      function.name = function.name, 
      intercept.only = intercept.only,
      predictors.names = predictors.names,
      M = M,
      n = n,
      nonparametric = nonparametric,
      nrow_X_vlm = nrow_X_vlm,
      orig.assign = attr(x, "assign"),
      p = ncol(x),
      ncol_X_vlm = ncol_X_vlm,
      ynames = dimnames(y)[[2]])


  crit.list <- list()
  if (criterion != "coefficients")
    crit.list[[criterion]] <- fit[[criterion]] <- new.crit

  for(ii in names(.min.criterion.VGAM)) {
    if (ii != criterion &&
        any(slotNames(family) == ii) &&
        length(body(slot(family, ii)))) {
          fit[[ii]] <-
          crit.list[[ii]] <-
          (slot(family, ii))(mu = mu, y = y, w = w,
                             res = FALSE, eta = eta, extra)
    }
  }



  if (w[1] != 1 || any(w != w[1]))
    fit$prior.weights <- w

  if (length(slot(family, "last")))
    eval(slot(family, "last"))

  structure(c(fit,
   list(predictors = tfit$predictors,
        contrasts = attr(x, "contrasts"),
        control = control,
        crit.list = crit.list,
        extra = extra,
        family = family,
        iter = iter,
        misc = misc,
        post = post,
        rss = tfit$rss,
        x = x,
        y = y)),
        vclass = slot(family, "vfamily"))
}

