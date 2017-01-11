# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.






vglm.fit <-
  function(x, y, w = rep_len(1, nrow(x)),
           X.vlm.arg = NULL,
           Xm2 = NULL, Ym2 = NULL,
           etastart = NULL, mustart = NULL, coefstart = NULL,
           offset = 0, family,
           control = vglm.control(),
           qr.arg = FALSE,
           constraints = NULL,
           extra = NULL,
           Terms = Terms, function.name = "vglm", ...) {

  if (is.null(criterion <- control$criterion))
    criterion <- "coefficients"

  eff.n <- nrow(x)  # + sum(abs(w[1:nrow(x)]))

  specialCM <- NULL
  post <- list()
  check.rank <- control$Check.rank
  nonparametric <- FALSE
  epsilon <- control$epsilon
  maxit <- control$maxit
  save.weights <- control$save.weights
  trace <- control$trace

  orig.stepsize <- control$stepsize
  minimize.criterion <- control$min.criterion


  fv <- NULL




  n <- nrow(x)





  stepsize <- orig.stepsize
  old.coeffs <- coefstart  # May be a NULL

  intercept.only <- ncol(x) == 1 && colnames(x) == "(Intercept)"
  y.names <- predictors.names <- NULL  # May be overwritten in @initialize

  n.save <- n


  if (length(slot(family, "initialize")))
    eval(slot(family, "initialize"))  # Initialize mu & M (& optionally w)







  if (length(etastart)) {
    eta <- etastart
    mu <- if (length(mustart)) mustart else
            slot(family, "linkinv")(eta, extra = extra)
  }



  if (length(mustart)) {
    mu <- mustart
    if (length(body(slot(family, "linkfun")))) {
      eta <- slot(family, "linkfun")(mu, extra = extra)
    } else {
      warning("argument 'mustart' assigned a value ",
              "but there is no 'linkfun' slot to use it")
    }
  }


  validparams <- validfitted <- TRUE
  if (length(body(slot(family, "validparams"))))
    validparams <- slot(family, "validparams")(eta, y = y, extra = extra)
  if (length(body(slot(family, "validfitted"))))
    validfitted <- slot(family, "validfitted")(mu, y = y, extra = extra)
  if (!(validparams && validfitted))
    stop("could not obtain valid initial values. ",
         "Try using 'etastart', 'coefstart' or 'mustart', else ",
         "family-specific arguments such as 'imethod'.")




  M <- NCOL(eta)



  if (length(slot(family, "constraints")))
    eval(slot(family, "constraints"))


  Hlist <- process.constraints(constraints, x = x, M = M,
                               specialCM = specialCM,
                               Check.cm.rank = control$Check.cm.rank)


  ncolHlist <- unlist(lapply(Hlist, ncol))





  X.vlm.save <-
    if (length(X.vlm.arg)) {
      X.vlm.arg
    } else {
      lm2vlm.model.matrix(x, Hlist, xij = control$xij, Xm2 = Xm2)
    }



  if (length(coefstart)) {
    eta <- if (ncol(X.vlm.save) > 1) {
             matrix(X.vlm.save %*% coefstart, n, M, byrow = TRUE) + offset
           } else {
             matrix(X.vlm.save  *  coefstart, n, M, byrow = TRUE) + offset
           }
    if (M == 1)
      eta <- c(eta)
    mu <- slot(family, "linkinv")(eta, extra = extra)
  }


  if (criterion != "coefficients") {
    tfun <- slot(family, criterion)  # family[[criterion]]
  }

  iter <- 1
  new.crit <- switch(criterion,
                     coefficients = 1,
                     tfun(mu = mu, y = y, w = w, res = FALSE,
                          eta = eta, extra = extra))

  deriv.mu <- eval(slot(family, "deriv"))
  wz <- eval(slot(family, "weight"))
  if (control$checkwz)
    wz <- checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)

  U <- vchol(wz, M = M, n = n, silent = !trace)
  tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
  z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

  one.more <- TRUE


  nrow.X.vlm <- nrow(X.vlm.save)
  ncol.X.vlm <- ncol(X.vlm.save)
  if (nrow.X.vlm < ncol.X.vlm)
    stop("There are ", ncol.X.vlm, " parameters but only ",
         nrow.X.vlm, " observations")






  while (one.more) {
    tfit <- vlm.wfit(xmat = X.vlm.save, zmat = z,
                     Hlist = NULL, U = U,
                     matrix.out = FALSE, is.vlmX = TRUE,
                     qr = qr.arg, xij = NULL)







    fv <- tfit$fitted.values
    new.coeffs <- tfit$coefficients  # c.list$coeff

    if (length(slot(family, "middle")))
      eval(slot(family, "middle"))

    eta <- fv + offset
    mu <- slot(family, "linkinv")(eta, extra = extra)

    if (length(slot(family, "middle2")))
      eval(slot(family, "middle2"))

    old.crit <- new.crit
    new.crit <-
        switch(criterion,
            coefficients = new.coeffs,
            tfun(mu = mu, y = y, w = w,
                 res = FALSE, eta = eta, extra = extra))


    if (trace && orig.stepsize == 1) {
      cat("VGLM    linear loop ", iter, ": ", criterion, "= ")
      UUUU <- switch(criterion,
                     coefficients =
                       format(new.crit,
                              digits = round(1 - log10(epsilon))),
                       format(new.crit,
                              digits = max(4,
                                           round(-0 - log10(epsilon) +
                                                 log10(sqrt(eff.n))))))
      switch(criterion,
             coefficients = {if (length(new.crit) > 2) cat("\n");
             cat(UUUU, fill = TRUE, sep = ", ")},
             cat(UUUU, fill = TRUE, sep = ", "))
    }  # if (trace && orig.stepsize == 1)


    take.half.step <- (control$half.stepsizing &&
                       length(old.coeffs)) &&
                       ((orig.stepsize != 1) ||
                        (!is.finite(new.crit)) ||  # 20160321
                        (criterion != "coefficients" &&
                        (if (minimize.criterion) new.crit > old.crit else
                                                 new.crit < old.crit)))
    if (!is.logical(take.half.step))
      take.half.step <- TRUE


    if (!take.half.step && length(old.coeffs))  {
      validparams <- validfitted <- TRUE
      if (length(body(slot(family, "validparams"))))
        validparams <- slot(family, "validparams")(eta, y = y, extra = extra)
      if (length(body(slot(family, "validfitted"))))
        validfitted <- slot(family, "validfitted")(mu, y = y, extra = extra)
      take.half.step <- !(validparams && validfitted)


     if (FALSE && take.half.step) {
       stepsize <- orig.stepsize / 4
      }
    }



    if (take.half.step) {
      stepsize <- (1 + (orig.stepsize != 1)) * orig.stepsize
      new.coeffs.save <- new.coeffs
      if (trace)
        cat("Taking a modified step")
      repeat {
        if (trace) {
          cat(".")
          flush.console()
        }
        stepsize <- stepsize / 2
        if (too.small <- stepsize < 1e-6)
          break
        new.coeffs <- (1-stepsize) * old.coeffs +
                         stepsize  * new.coeffs.save

        if (length(slot(family, "middle")))
          eval(slot(family, "middle"))

        fv <- X.vlm.save %*% new.coeffs
        if (M > 1)
          fv <- matrix(fv, n, M, byrow = TRUE)

        eta <- fv + offset
        mu <- slot(family, "linkinv")(eta, extra = extra)

        if (length(slot(family, "middle2")))
          eval(slot(family, "middle2"))

        new.crit <-
          switch(criterion,
                 coefficients = new.coeffs,
                 tfun(mu = mu, y = y, w = w,
                      res = FALSE, eta = eta, extra = extra))


        validparams <- validfitted <- TRUE
        if (length(body(slot(family, "validparams"))))
          validparams <- slot(family, "validparams")(eta, y, extra = extra)
        if (length(body(slot(family, "validfitted"))))
          validfitted <- slot(family, "validfitted")(mu, y, extra = extra)

        if (validparams && validfitted &&
           (is.finite(new.crit)) &&  # 20160321
           (criterion == "coefficients" ||
           (( minimize.criterion && new.crit < old.crit) ||
            (!minimize.criterion && new.crit > old.crit))))
          break
      }  # of repeat

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
                                  digits = round(1 - log10(epsilon))),
                           format(new.crit,
                                  digits = max(4,
                                               round(-0 - log10(epsilon) +
                                                     log10(sqrt(eff.n))))))

          switch(criterion,
                 coefficients = {
                 if (length(new.crit) > 2) cat("\n");
                 cat(UUUU, fill = TRUE, sep = ", ")},
                 cat(UUUU, fill = TRUE, sep = ", "))
        }  # if (trace)

        one.more <- eval(control$convergence)
      }  # Not too.small
    } else {
      one.more <- eval(control$convergence)
    }
    flush.console()

    if (!is.logical(one.more))
      one.more <- FALSE

    if (one.more) {
      iter <- iter + 1
      deriv.mu <- eval(slot(family, "deriv"))
      wz <- eval(slot(family, "weight"))
      if (control$checkwz)
        wz <- checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)

      U <- vchol(wz, M = M, n = n, silent = !trace)
      tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
      z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

    }  # if (one.more)



    if (!one.more && take.half.step && orig.stepsize == 1)
      warning("some quantities such as z, residuals, SEs may ",
              "be inaccurate due to convergence at a half-step")



    old.coeffs <- new.coeffs
  }  # End of while()


  if (maxit > 1 && iter >= maxit && !control$noWarning)
    warning("convergence not obtained in ", maxit, " IRLS iterations")




  dnrow.X.vlm <- labels(X.vlm.save)
  xnrow.X.vlm <- dnrow.X.vlm[[2]]
  ynrow.X.vlm <- dnrow.X.vlm[[1]]

  if (length(slot(family, "fini")))
    eval(slot(family, "fini"))

  if (M > 1)
    fv <- matrix(fv, n, M)

  final.coefs <- new.coeffs  # Was tfit$coefficients prior to 20160317
  asgn <- attr(X.vlm.save, "assign")

  names(final.coefs) <- xnrow.X.vlm

  rank <- tfit$rank
  cnames <- xnrow.X.vlm

  if (check.rank && rank < ncol.X.vlm)
    stop("vglm() only handles full-rank models (currently)")


  R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
  R[lower.tri(R)] <- 0
  attributes(R) <- list(dim = c(ncol.X.vlm, ncol.X.vlm),
                        dimnames = list(cnames, cnames), rank = rank)

  effects <- tfit$effects
  neff <- rep_len("", nrow.X.vlm)
  neff[seq(ncol.X.vlm)] <- cnames
  names(effects) <- neff

  dim(fv) <- c(n, M)
  dn <- labels(x)
  yn <- dn[[1]]
  xn <- dn[[2]]


  wresiduals <- z - fv  # Replaced by fv 20160408
  if (M == 1) {
    fv <- as.vector(fv)
    wresiduals <- as.vector(wresiduals)
    names(wresiduals) <- names(fv) <- yn
  } else {
    dimnames(wresiduals) <-
    dimnames(fv)         <- list(yn, predictors.names)
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


  fit <- list(assign = asgn,
              coefficients = final.coefs,
              constraints = Hlist,
              df.residual = nrow.X.vlm - rank,
              df.total = n * M,
              effects = effects,   # this is good
              fitted.values = mu,   # this is good
              offset = offset,
              rank = rank,   # this is good
              residuals = wresiduals,
              R = R,
              terms = Terms)  # terms: This used to be done in vglm()

  if (qr.arg) {
    fit$qr <- tfit$qr
    dimnames(fit$qr$qr) <- dnrow.X.vlm
  }

  if (M == 1) {
    wz <- as.vector(wz)  # Convert wz into a vector
  } # else
  fit$weights <- if (save.weights) wz else NULL


  misc <- list(
      colnames.x = xn,
      colnames.X.vlm = xnrow.X.vlm,
      criterion = criterion,
      function.name = function.name,
      intercept.only = intercept.only,
      predictors.names = predictors.names,
      M = M,
      n = n,
      nonparametric = nonparametric,
      nrow.X.vlm = nrow.X.vlm,
      orig.assign = attr(x, "assign"),
      p = ncol(x),
      ncol.X.vlm = ncol.X.vlm,
      ynames = colnames(y))


  crit.list <- list()
  if (criterion != "coefficients")
    crit.list[[criterion]] <- fit[[criterion]] <- new.crit

  for (ii in names(.min.criterion.VGAM)) {
    if (ii != criterion &&
        any(slotNames(family) == ii) &&
        length(body(slot(family, ii)))) {
      fit[[ii]] <-
      crit.list[[ii]] <- (slot(family, ii))(mu = mu, y = y, w = w,
                                            res = FALSE, eta = eta,
                                            extra = extra)
    }
  }



  if (w[1] != 1 || any(w != w[1]))
    fit$prior.weights <- w

  if (length(slot(family, "last")))
    eval(slot(family, "last"))

  structure(c(fit,
   list(predictors = fv,  # tfit$predictors,
        contrasts = attr(x, "contrasts"),
        control = control,
        crit.list = crit.list,
        extra = extra,
        family = family,
        iter = iter,
        misc = misc,
        post = post,
        ResSS = tfit$ResSS,
        x = x,
        y = y)),
        vclass = slot(family, "vfamily"))
}  # vglm.fit()


