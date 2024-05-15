# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.









rrvglm.fit <-
  function(x, y, w = rep_len(1, nrow(x)),
           etastart = NULL, mustart = NULL,
           coefstart = NULL,
           offset = 0, family,
           control = rrvglm.control(...),
           criterion = "coefficients",
           qr.arg = FALSE,
           constraints = NULL,
           extra = NULL,
           Terms = Terms,
           function.name = "rrvglm", ...) {

    if (length(slot(family, "start1")))
      eval(slot(family, "start1"))


    eff.n <- nrow(x)  # + sum(abs(w[1:nrow(x)]))

    specialCM <- NULL
    post <- list()
    check.rank <- TRUE  # !control$Quadratic
    check.rank <- control$Check.rank
    is.rrvglm <- control$is.rrvglm  # Default: T
    is.drrvglm <- control$is.drrvglm  # Default: F
    H.A.alt <- control$H.A.alt  # NULL by default
    H.A.thy <- control$H.A.thy  # NULL by default
    H.C <- control$H.C  # list() by default
    nonparametric <- FALSE
    epsilon <- control$epsilon
    maxit <- control$maxit
    save.weights <- control$save.weights
    trace <- control$trace
    orig.stepsize <- control$stepsize
    minimize.criterion <- control$min.criterion
    history <- NULL
    Alt <- list(is.rrvglm = is.rrvglm,
                is.drrvglm = is.drrvglm,
                H.A.thy = list(),
                H.A.alt = list(), H.C = list())


    fv <- one.more <- rrr.expression <-
    modelno <- NULL
    RRR.expression <- paste("rrr",
      control$Algorithm, "expression", sep = ".")

    n <- dim(x)[1]



    copy.X.vlm <- FALSE  # May b overwritten in @initz
    stepsize <- orig.stepsize
    old.coeffs <- coefstart

    intercept.only <- ncol(x) == 1 &&
              dimnames(x)[[2]] == "(Intercept)"
    y.names <-
    predictors.names <- NULL  # Mayb ovrwrttn in @initz


    n.save <- n



    Rank <- control$Rank
    rrcontrol <- control  #

    if (length(slot(family, "initialize")))
      eval(slot(family, "initialize"))  # Initlz mu & M (& opt'ly w)


    eval(rrr.init.expression)


    if (length(etastart)) {
      eta <- etastart
      mu <- if (length(mustart)) mustart else
            if (length(body(slot(family, "linkinv"))))
              slot(family, "linkinv")(eta, extra) else
      warning("arg 'etastart' assigned a value but",
              " there is no @linkinv to use it")
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

    if (is.character(rrcontrol$Dzero)) {
      index <- match(rrcontrol$Dzero,
                     dimnames(as.matrix(y))[[2]])
      if (anyNA(index))
        stop("Dzero arg didn't fully match y-names")
      if (length(index) == M)
        stop("all linear predictors are linear in ",
             "the lat var(s); so set 'Quadratic = F'")
      rrcontrol$Dzero <- control$Dzero <- index
    }





    if (length(family@constraints))
      eval(family@constraints)


    special.matrix <- matrix(-34956.125, M, M)  # Unlikely used matrix

    just.testing <- cm.VGAM(special.matrix, x, rrcontrol$noRRR,
                            constraints)

    findex <- trivial.constraints(just.testing, special.matrix)
    if (is.null(just.testing)) findex <- NULL  # 20100617
    tc1 <- trivial.constraints(constraints)

    if (!is.null(findex) && !control$Quadratic && sum(!tc1)
        && FALSE  # 20240329 Suppress this
       ) {
      for (ii in names(tc1))
        if (!tc1[ii] && !any(ii == names(findex)[findex == 1]))
          warning("'", ii, "' is a non-trivial constraint that ",
                  "will be overwritten by reduced-rank regression")
    }  # if

    if (!is.null(findex) && all(findex == 1))
      stop("use vglm(), not rrvglm()!")
    colx1.index <- names.colx1.index <- NULL
    dx2 <- dimnames(x)[[2]]
    if (sum(findex)) {
      asx <- attr(x, "assign")
      for (ii in names(findex))
        if (findex[ii]) {
          names.colx1.index <- c(names.colx1.index, dx2[asx[[ii]]])
          colx1.index <- c(colx1.index, asx[[ii]])
        }
      names(colx1.index) <- names.colx1.index
    }  # if (sum(findex))
    rrcontrol$colx1.index <-
      control$colx1.index <- colx1.index  # Save it on the object
    colx2.index <- 1:ncol(x)
    names(colx2.index) <- dx2
    if (length(colx1.index))
      colx2.index <- colx2.index[-colx1.index]

    p1 <- length(colx1.index)
    p2 <- length(colx2.index)
    rrcontrol$colx2.index <-
      control$colx2.index <- colx2.index  # Save it on the object
    Index.corner <- control$Index.corner




    skip1vlm <- length(rrcontrol$Ainit) &&
               !length(rrcontrol$Cinit)
      control$skip1vlm <-
    rrcontrol$skip1vlm <- skip1vlm
    Amat <- if (length(rrcontrol$Ainit)) {
      if (length(rrcontrol$str0) &&
          any(is.non0(rrcontrol$Ainit[
              rrcontrol$str0, ]))) {
        warning("Ainit[str0, ] has non-0 values")
        rrcontrol$Ainit[rrcontrol$str0, ] <- 0
      }
      if (rrcontrol$Corner && 
          any(is.non0(c(rrcontrol$Ainit[
              rrcontrol$Index.corner, ]) -
              c(diag(Rank)))))
        stop("Ainit violates corner constraints")
      rrcontrol$Ainit
    } else matrix(
      rnorm(M * Rank, 0, rrcontrol$sd.Cinit),
            M, Rank)
    Cmat <- if (length(rrcontrol$Cinit))
            rrcontrol$Cinit else {
    if (!rrcontrol$Use.Init.Poisson.QO) {
      matrix(rnorm(p2 * Rank, 0,
                   rrcontrol$sd.Cinit), p2, Rank)
    } else {
       .Init.Poisson.QO(ymat = as.matrix(y),
             X1 = if (length(colx1.index))
        x[, colx1.index, drop = FALSE] else NULL,
        X2 = x[, colx2.index, drop = FALSE],
        Rank = rrcontrol$Rank,
        trace = rrcontrol$trace,
        max.ncol.etamat = rrcontrol$Etamat.colmax,
        Crow1positive = rrcontrol$Crow1positive,
        isd.latvar = rrcontrol$isd.latvar)
    }
  }








    if (control$Corner)
      Amat[control$Index.corner,] <- diag(Rank)
    if (length(control$str0))
      Amat[control$str0, ] <- 0

    rrcontrol$Ainit <-
      control$Ainit <- Amat  # Good for valt0()
    rrcontrol$Cinit <-
      control$Cinit <- Cmat  # Good for valt0()

    Hlist <- process.constraints(constraints,
               x, M, specialCM = specialCM)

  if (length(H.C)) {
 setp2 <- names(H.C)  # 1st choice, may be NULL
 if (length(setp2) == p2) {  # keyword: setp2
   match2 <- match(setp2, names(colx2.index))
   if (any(is.na(match2))) {
     warning("names(H.C) & 'colx2.index' do not",
     " match so using H.C[[1]], H.C[[2]], ...",
     " position order for the C matrix")
     setp2 <- 1:p2  # 2nd choice.
   } else {  # Check the order
     if (!all(diff(match2) == 1)) {
       warning("components of H.C will be ",
               "sorted and saved on the object")
       sorted.H.C <- vector("list", p2)
       for (i in seq(p2))
         sorted.H.C[[(match2[i])]] <- H.C[[i]]
       names(sorted.H.C) <- names(colx2.index) 
       rrcontrol$H.C <-  # Overwrite this
       control$H.C <- H.C <- sorted.H.C
     }  # Unsorted
   }
  }
    for (k in setp2)
      Hlist.save <- Hlist <- 
        replaceCMs(Hlist, Amat %*% H.C[[k]],
                   colx2.index[k])  # (Pos456)
  }  # length(H.C)


    nice31 <-  control$Quadratic &&
    (!control$eq.tol || control$I.tolerances) &&
     all(trivial.constraints(Hlist) == 1)

    if (!is.drrvglm)
      Hlist <- Hlist.save <-
        replaceCMs(Hlist, Amat, colx2.index)


    ncolHlist <- unlist(lapply(Hlist, ncol))


    X.vlm.save <- if (control$Quadratic) {
      tmp500 <- lm2qrrvlm.model.matrix(x = x,
                     Hlist = Hlist,
                     C = Cmat, control = control)
      xsmall.qrr <- tmp500$new.latvar.model.matrix
      H.list <- tmp500$constraints
      if (FALSE && modelno == 3) {
  H.list[[1]] <- (H.list[[1]])[, c(TRUE, FALSE),
                       drop = FALSE]  # Amat
  H.list[[2]] <- (H.list[[2]])[, c(TRUE, FALSE),
                       drop = FALSE]  # D
      }

      latvar.mat <- tmp500$latvar.mat
      if (length(tmp500$offset)) {
        offset <- tmp500$offset
      }
      lm2vlm.model.matrix(xsmall.qrr, H.list,
                     label.it = control$label.it,
                     xij = control$xij)
    } else {  # RR-VGLMs:
      latvar.mat <-
        x[, colx2.index, drop = FALSE] %*% Cmat
      lm2vlm.model.matrix(x, Hlist,
                     label.it = control$label.it,
                     xij = control$xij)
    }




    if (length(coefstart)) {
      eta <- if (ncol(X.vlm.save) > 1)
          X.vlm.save %*% coefstart + offset else
          X.vlm.save  *  coefstart + offset
      eta <- if (M > 1)
  matrix(eta, ncol = M, byrow = TRUE) else c(eta)


      mu <- family@linkinv(eta, extra)
    }

    if (criterion != "coefficients") {
      tfun <- slot(family, criterion)  # family[[criterion]]
    }

    iter <- 1
    new.crit <- switch(criterion,
          coefficients = 1,
          tfun(mu = mu, y = y, w = w, res = FALSE,
               eta = eta, extra))
    old.crit <- ifelse(minimize.criterion,
                       10 * new.crit + 10,
                      -10 * new.crit - 10)
    deriv.mu <- eval(family@deriv)

    wz <- eval(family@weight)
    if (control$checkwz)
      wz <- checkwz(wz, M = M, trace = trace,
                    wzepsilon = control$wzepsilon)

    U <- vchol(wz, M = M, n = n, silent = !trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
    z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset

    c.list <-
      list(z = as.double(z),
           fit = as.double(t(eta)),
           one.more = TRUE,
 coeff = as.double(rep_len(1, ncol(X.vlm.save))),
           U = as.double(U),
           copy.X.vlm = copy.X.vlm,
           X.vlm = if (copy.X.vlm)
                     as.double(X.vlm.save) else
                     double(3))



    dX.vlm <- as.integer(dim(X.vlm.save))
    nrow.X.vlm <- dX.vlm[[1]]
    ncol.X.vlm <- dX.vlm[[2]]

    if (nrow.X.vlm < ncol.X.vlm)
      stop(ncol.X.vlm, " parameters but only ",
           nrow.X.vlm, " observations")

    bf.call <- expression(vlm.wfit(X.vlm.save,
        zedd,  # restored
        Hlist = if (control$Quadratic)
                    H.list else Hlist,
        ncolx = ncol(x), U = U,
        Eta.range = control$Eta.range,
        matrix.out = if(control$Quadratic)
                         FALSE else TRUE,
        is.vlmX = TRUE, qr = qr.arg,
        label.it = control$label.it,
        xij = control$xij))

    while (c.list$one.more) {
      if (control$Quadratic) {
        zedd <- as.matrix(z)
        if (control$Corner)
          zedd[, Index.corner] <-
          zedd[, Index.corner] - latvar.mat
      } else {
        zedd <- z
      }

      if (!nice31)
        tfit <- eval(bf.call)  # tfit$fitted.values is n x M

      if (!control$Quadratic) {
        Cmat <- tfit$mat.coef[colx2.index, , drop = FALSE] %*%
                Amat %*% solve(t(Amat) %*% Amat)
        rrcontrol$Ainit <- control$Ainit <- Amat
        rrcontrol$Cinit <- control$Cinit <- Cmat
      }  # Quadratic 

      if (!nice31)
        c.list$coeff <- tfit$coefficients

      if (control$Quadratic) {
        if (control$Corner)
   tfit$fitted.values[, Index.corner] <-
   tfit$fitted.values[, Index.corner] + latvar.mat
      }

      if (!nice31)
tfit$predictors <- tfit$fitted.values  # No offset
      if (!nice31)
        c.list$fit <- tfit$fitted.values

      if (!c.list$one.more) {
        break
      }



      fv <- c.list$fit
      new.coeffs <- c.list$coeff

      if (length(family@middle1))
        eval(family@middle1)

      eta <- fv + offset

      mu <- family@linkinv(eta, extra)

      if (length(family@middle2))
        eval(family@middle2)

      old.crit <- new.crit
      new.crit <-
        switch(criterion,
               coefficients = new.coeffs,
               tfun(mu = mu, y = y, w = w,
                    res = FALSE, eta = eta, extra))


    if (is.null(history))
      history <- matrix(NA_real_, maxit,
      if (criterion == "coefficients")
          length(new.coeffs) else 1)
    history[iter, ] <- new.crit  # Imperfect (e.g., step-halving).



      if (trace && orig.stepsize == 1) {
        cat(if (control$Quadratic) "QRR-VGLM" else "RR-VGLM",
            "   linear loop ", iter, ": ", criterion, "= ")
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

      take.half.step <-
        (control$half.stepsizing &&
         length(old.coeffs)) &&
        !control$Quadratic &&
        ((orig.stepsize != 1) ||
        (criterion != "coefficients" &&
        (if (minimize.criterion)
           new.crit > old.crit else
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
          new.coeffs <- (1 - stepsize) * old.coeffs +
                             stepsize  * new.coeffs.save

          if (length(family@middle1))
            eval(family@middle1)

          fv <- X.vlm.save %*% new.coeffs
          if (M > 1)
            fv <- matrix(fv, n, M, byrow = TRUE)

          eta <- fv + offset

          mu <- family@linkinv(eta, extra)

          if (length(family@middle2))
            eval(family@middle2)

          new.crit <-
            switch(criterion,
                   coefficients = new.coeffs,
                   tfun(mu = mu,y = y,w = w,res = FALSE,
                        eta = eta,extra))

          if ((criterion == "coefficients") ||
              ( minimize.criterion && new.crit < old.crit) ||
              (!minimize.criterion && new.crit > old.crit))
            break
          }

          if (trace)
            cat("\n")
          if (too.small) {
            warning("iterations terminated because ",
                    "half-step sizes are very small")
            one.more <- FALSE
          } else {
            if (trace) {
              cat(if (control$Quadratic) "QRR-VGLM" else "RR-VGLM",
                  "   linear loop ", iter, ": ", criterion, "= ")
              UUUU <-
                switch(criterion,
                       coefficients =
                         format(new.crit,
                                dig = round(1 - log10(epsilon))),
                         format(new.crit,
                                dig = max(4,
                                      round(-0 - log10(epsilon) +
                                      log10(sqrt(eff.n))))))

              switch(criterion,
                     coefficients = {if (length(new.crit) > 2)
                                       cat("\n");
                             cat(UUUU, fill = TRUE, sep = ", ")},
                     cat(UUUU, fill = TRUE, sep = ", "))
            }

            one.more <- eval(control$convergence)
          }
        } else {
          one.more <- eval(control$convergence)
        }
      flush.console()

    if (!is.logical(one.more))
      one.more <- FALSE

      if (one.more) {
        iter <- iter + 1
        deriv.mu <- eval(family@deriv)
        wz <- eval(family@weight)
        if (control$checkwz)
          wz <- checkwz(wz, M = M, trace = trace,
                    wzepsilon = control$wzepsilon)


          wz <- matrix(wz, nrow = n)
        U <- vchol(wz, M, n, silent = !trace)
        tvfor <- vforsub(U, as.matrix(deriv.mu), M, n)
        z <- eta + vbacksub(U, tvfor, M, n) -
            offset  # Contains \bI \bnu

        rrr.expression <- get(RRR.expression)
        eval(rrr.expression)


        c.list$z <- z  # contains \bI_{Rank} \bnu
        c.list$U <- U
        if (copy.X.vlm) c.list$X.vlm <- X.vlm.save
      }

      c.list$one.more <- one.more
      c.list$coeff <- runif(length(
        new.coeffs))  # 20030312; twist needed!
      old.coeffs <- new.coeffs

    }  # End of while()  ....................


  if (maxit > 1 && iter >= maxit && !control$noWarning)
    warning("convergence not obtained in ", maxit, " iterations")





    dnrow.X.vlm <- labels(X.vlm.save)
    xnrow.X.vlm <- dnrow.X.vlm[[2]]
    ynrow.X.vlm <- dnrow.X.vlm[[1]]

    if (length(family@fini1))
      eval(family@fini1)

    if (M > 1 && !nice31)
      tfit$predictors <- matrix(tfit$predictors, n, M)

    asgn <- attr(X.vlm.save, "assign")
    if (nice31) {
      coefs <- rep_len(0, length(xnrow.X.vlm))
        rank <- ncol.X.vlm
    } else {
      coefs <- tfit$coefficients
      names(coefs) <- xnrow.X.vlm
      rank <- tfit$rank
    }

    cnames <- xnrow.X.vlm


  colnames(history) <- if (criterion == "coefficients")
                       xnrow.X.vlm else criterion

    if (check.rank && rank < ncol.X.vlm)
      stop("rrvglm only handles full-rank models")

    if (nice31) {
      R <- matrix(NA_real_, 5, 5)
    } else {
      R <- tfit$qr$qr[1:ncol.X.vlm,
                      1:ncol.X.vlm, drop = FALSE]
      R[lower.tri(R)] <- 0
      attributes(R) <-
          list(dim = c(ncol.X.vlm, ncol.X.vlm),
               dimnames = list(cnames, cnames),
               rank = rank)
    }

    if (nice31) {
      effects <- rep_len(0, 77)
    } else {
      effects <- tfit$effects
      neff <- rep_len("", nrow.X.vlm)
      neff[seq(ncol.X.vlm)] <- cnames
      names(effects) <- neff

      dim(tfit$predictors) <- c(n, M)
    }
    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]

    if (nice31) {
      residuals <- z - fv
      if (M == 1) {
        residuals <- as.vector(residuals)
        names(residuals) <- yn
      } else {
        dimnames(residuals) <- list(yn, predictors.names)
      }
    } else {
        residuals <- z - tfit$predictors
        if (M == 1) {
          tfit$predictors <- as.vector(tfit$predictors)
          residuals <- as.vector(residuals)
          names(residuals) <- names(tfit$predictors) <- yn
        } else {
          dimnames(residuals) <-
          dimnames(tfit$predictors) <- list(yn, predictors.names)
        }
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




  elts.tildeA <- if (is.drrvglm)
     sum(unlist(lapply(H.A.alt, ncol))) else
     (M - Rank - length(control$str0)) * Rank
    no.dpar <- 0
    df.residual <- nrow.X.vlm - rank -
         ifelse(control$Quadratic, Rank*p2, 0) -
         no.dpar - elts.tildeA


  fit <-
    list(assign = asgn,
         coefficients = coefs,
constraints = if (control$Quadratic) H.list else
                                     Hlist,
         df.residual = df.residual,
         df.total = n*M,
         effects = effects,
         fitted.values = mu,
         offset = offset,
         rank = rank,
         residuals = residuals,
         R = R,
         terms = Terms)  # terms was in vglm()

    if (qr.arg && !nice31) {
      fit$qr <- tfit$qr
      dimnames(fit$qr$qr) <- dnrow.X.vlm
    }

    if (M == 1) {
      wz <- as.vector(wz)
    } # else
    fit$weights <- if (save.weights) wz else NULL

    misc <- list(
        colnames.x = xn,
        colnames.X.vlm = xnrow.X.vlm,
        criterion = criterion,
        function.name = function.name,
   history = history[seq(iter), , drop = FALSE],
        intercept.only = intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = n,
        nonparametric = nonparametric,
        nrow.X.vlm = nrow.X.vlm,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        ncol.X.vlm = ncol.X.vlm,
        ynames = dimnames(y)[[2]])

    if (!is.logical(one.more))
      one.more <- FALSE
    if (one.more)
      misc$rrr.expression <- rrr.expression  #


    crit.list <- list()
    if (criterion != "coefficients")
      crit.list[[criterion]] <-
      fit[[criterion]] <- new.crit

    for (ii in names( .min.criterion.VGAM )) {
      if (ii != criterion &&
       any(slotNames(family) == ii) &&
           length(body(slot(family, ii)))) {
            fit[[ii]] <- crit.list[[ii]] <-
            (slot(family, ii))(mu = mu, y = y,
             w = w, res = FALSE, eta = eta, extra)
      }
    }

    if (w[1] != 1 || any(w != w[1]))
      fit$prior.weights <- w

    if (length(family@last))
      eval(family@last)


    structure(c(fit,
                list(predictors = if (nice31)
                matrix(eta, n, M) else
                tfit$predictors,
        contrasts = attr(x, "contrasts"),
  control = rrcontrol,  # Was control < 20231118
        crit.list = crit.list,
        extra = extra,
        family = family,
        iter = iter,
        misc = misc,
        post = post,
        ResSS = if (nice31) 000 else tfit$ResSS,
        is.rrvglm  = is.rrvglm,   # papertrail
        is.drrvglm = is.drrvglm,  # papertrail
        use.H.A = Alt$use.H.A,  # 20231229
        H.A.alt = Alt$H.A.alt,  # 20231114
        H.A.thy = Alt$H.A.thy,  # 20231230
        H.C     = Alt$H.C,      # 20231114
        A.est   = Alt$A.est,    # Transfer these
        C.est   = Alt$C.est,    #
        Amask   = Alt$Amask,    #
        Avec    = Alt$Avec,     #
        B1Cvec  = Alt$B1Cvec,   #
        clist1  = Alt$clist1,   #
        RAvcov  = Alt$RAvcov,
        RCvcov  = Alt$RCvcov,
        valt0.ResSS = Alt$valt0.ResSS,
        x = x,
        y = y)),
        vclass = family@vfamily)
}


