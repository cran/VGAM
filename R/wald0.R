# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.
















 wald.stat.vlm <-
  function(object,
           values0 = 0,
           subset = NULL,  # Useful for Cox model as a poissonff().
           omit1s = TRUE,
           all.out = FALSE,  # If TRUE then lots of output returned
           iterate = TRUE,  # If FALSE then other coeffs at MLE
           trace = FALSE,  # NULL,
           as.summary = FALSE,  # If TRUE then ordinary Wald statistic
           ...) {


  foo1 <- function(score.really = FALSE, ...)
    score.really
  score.really <- foo1(...)

  foo2 <- function(lrt.really = FALSE, ...)
    lrt.really
  lrt.really <- foo2(...)

  if (score.really && lrt.really)
    stop("cannot have both 'score.really' and 'lrt.really'")


  M <- npred(object)  # Some constraints span across responses
  all.Hk <- constraints(object, matrix = TRUE)
  X.lm  <- model.matrix(object, type =  "lm")
  X.vlm.save <- model.matrix(object, type = "vlm")
  eta.mat.orig <- predict(object)
  n.LM <- NROW(eta.mat.orig)
  p.VLM <- ncol(all.Hk)

  vc2 <- vcov(object)
  cobj <- coef(object)
  signed.Lrt.0 <-
  Lrt.0 <-
  Score.0 <-
  SE2.0 <- rep_len(NA_real_, p.VLM)  # More than enough storage


  Pnames <- names(B0 <- coef(object))
  if (any(is.na(B0)))
    stop("currently cannot handle NA-valued regression coefficients")


  if (is.character(subset))
    subset <- match(subset, Pnames)
  if (is.null(subset))
    subset <- 1:p.VLM


  Xm2 <- model.matrix(object, type = "lm2")  # Could be a 0 x 0 matrix
  if (!length(Xm2))
     Xm2 <- NULL  # Make sure. This is safer
  clist <- constraints(object, type = "term")  # type = c("lm", "term")
  H1 <- clist[["(Intercept)"]]
  if (omit1s && length(H1) && any(subset <= ncol(H1))) {
    if (length(clist) == 1)
      return(NULL)  # Regressed against intercept only
    subset <- subset[subset > ncol(H1)]
  }






    if (is.logical(trace))
      object@control$trace <- trace
  mf <- model.frame(object)
  Y <- model.response(mf)
  if (!is.factor(Y))
    Y <- as.matrix(Y)
  OOO.orig <- object@offset
  if (!length(OOO.orig) || all(OOO.orig == 0))
    OOO.orig <- matrix(0, n.LM, M)
  mt <- attr(mf, "terms")
  Wts <- model.weights(mf)
  if (length(Wts) == 0L)
    Wts <- rep(1, n.LM)  # Safest (uses recycling and is a vector)
  summ <- summary(object)
  DispersionParameter <- summ@dispersion
  if (!all(DispersionParameter == 1))
    stop("Currently can only handle dispersion parameters ",
         "that are equal to 1")
  Fam <- object@family
  if (lrt.really) {
    Original.de <- deviance(object)  # Could be NULL
    if (!(use.de <- is.Numeric(Original.de)))
      Original.ll <- logLik(object)

    quasi.type <- if (length(tmp3 <- Fam@infos()$quasi.type))
      tmp3 else FALSE
    if (quasi.type)
      stop("currently this function cannot handle quasi-type",
           " models or models with an estimated dispersion parameter")
  }  # lrt.really








  kvec.use <- subset
  values0.use <- cobj * NA  # A vector of NAs of length == p.VLM
  values0.use[kvec.use] <- values0  # Recycle and put in right place

  if (as.summary) {  # For Wald-type statistics only
    csobj <- coef(summary(object))[kvec.use, , drop = FALSE]
    wald.stat <- csobj[, "z value"]
    SE0 <- csobj[, "Std. Error"]
    cobj <- cobj[kvec.use]
    values0.use <- values0.use[kvec.use]
    names(SE0) <- names(cobj)
    names(values0.use) <- names(cobj)
    names(wald.stat) <- names(cobj)
    if (all.out) return(
      list(wald.stat  = wald.stat,
           SE0        = SE0,
           values0    = values0.use)) else
      return(wald.stat)
  }  # as.summary



  temp1 <- object
  for (kay in kvec.use) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,



    if (iterate) {
      if (NCOL(X.vlm.save) == 1)
        stop("The large model matrix only has one column")
      X.vlm.mk <- X.vlm.save[, -kay, drop = FALSE]
 # This is needed by vglm.fit():
      attr(X.vlm.mk, "assign") <- attr(X.vlm.save, "assign")  # zz wrong!
    
      ooo <- if (values0.use[kay] == 0) OOO.orig else
             OOO.orig + matrix(X.vlm.save[, kay] * values0.use[kay],
                               n.LM, M, byrow = TRUE)
      fm <- vglm.fit(x = X.lm,  # Try this
                     y = Y, w = Wts,
                     X.vlm.arg = X.vlm.mk,
                     Xm2 = Xm2, Terms = mt,
                     constraints = clist, extra = object@extra,
                     etastart = eta.mat.orig,
                     offset = ooo, family = Fam,
                     control = object@control)
    }  # iterate


    if (lrt.really) {  # +++++++++++++++++++++++++++++++++++++

        zee <- if (use.de) {
          fm$crit.list[["deviance"]] - Original.de
        } else {
          2 * (Original.ll - fm$crit.list[["loglikelihood"]])
        }
        if (zee > -1e-3) {
          zee <- max(zee, 0)
        } else {
          warning("omitting 1 column has found a better solution, ",
                  "so the original fit had not converged")
        }
        zedd <- zee  # sgn * sqrt(zee)
      signed.Lrt.0[kay] <- sqrt(zedd) *
                           sign(cobj[kay] - values0.use[kay])
      Lrt.0[kay] <- zedd
    } else {  # +++++++++++++++++++++++++++++++++++++

    eta.mat.use <- if (iterate) {
      if (is.matrix(fm$predictors))
        fm$predictors else as.matrix(fm$predictors)
    } else {
      eta.mat.orig +
      matrix(X.vlm.save[, kay] * (values0.use[kay] - cobj[kay]),
             n.LM, M, byrow = TRUE)
    }


    temp1@predictors <- eta.mat.use
    temp1@fitted.values <- cbind(
      temp1@family@linkinv(eta = temp1@predictors,
                           extra = temp1@extra))  # Make sure a matrix
    wwt.both <- weights(temp1, type = "working", ignore.slot = TRUE,
                        deriv.arg = score.really)
    if (score.really) {
      deriv.new <- wwt.both$deriv
      wwt.new <- wwt.both$weights
    } else {
      wwt.new <- wwt.both
    }


    U <- vchol(wwt.new, M = M, n = n.LM, silent = TRUE)
    w12X.vlm <- mux111(U, X.vlm.save, M = M)


    qrstr <- qr(w12X.vlm)
    if (!all(qrstr$pivot == 1:length(qrstr$pivot)))
      stop("cannot handle pivoting just yet")
    R <- qr.R(qrstr)  # dim(R) == ncol(w12X.vlm); diags may be negative


    covun <- chol2inv(R)  # This is for (t(X.vlm) %*% W %*% X.vlm)^{-1}

    SE2.0[kay] <- diag(covun)[kay]

    if (score.really)
      Score.0[kay] <-
        sum(deriv.new * matrix(X.vlm.save[, kay], n.LM, M, byrow = TRUE))

    }  # !lrt.really +++++++++++++++++++++++++++++++++++++


  }  # for (kay in kvec.use)  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,






  cobj <- cobj[kvec.use]
  SE0 <- sqrt(SE2.0[kvec.use])  # All NAs if 'wald'
  names(SE0) <- names(cobj)
  values0.use <- values0.use[kvec.use]
  names(values0.use) <- names(cobj)


  if (lrt.really) {
    Lrt.0 <- Lrt.0[kvec.use]
    names(Lrt.0) <- names(cobj)
    signed.Lrt.0 <- signed.Lrt.0[kvec.use]
    names(signed.Lrt.0) <- names(cobj)
    if (all.out)
      list(lrt.stat   = signed.Lrt.0,
           Lrt.stat2  = Lrt.0,
           pvalues    = pchisq(Lrt.0, df = 1, lower.tail = FALSE),
           values0    = values0.use) else
      signed.Lrt.0
  } else if (score.really) {
    Score.0 <- Score.0[kvec.use]
    names(Score.0) <- names(cobj)
    score.stat <- Score.0 * SE0
    if (all.out)
      list(score.stat = score.stat,
           SE0        = SE0,  # Same as Wald
           values0    = values0.use) else
      score.stat
  } else {
    wald.stat <- (cobj - values0.use) / SE0
    if (all.out)
      list(wald.stat  = wald.stat,
           SE0        = SE0,
           values0    = values0.use) else
      wald.stat
  }
}  # wald.stat.vlm





 if (!isGeneric("wald.stat"))
   setGeneric("wald.stat", function(object, ...)
     standardGeneric("wald.stat"))


setMethod("wald.stat", "vlm", function(object, ...)
          wald.stat.vlm(object, ...))

















score.stat.vlm <-
  function(object,
           values0 = 0,
           subset = NULL,  # Useful for Cox model as a poissonff().
           omit1s = TRUE,
           all.out = FALSE,  # If TRUE then lots of output returned
           iterate = TRUE,  # If FALSE then other coeffs at MLE
           trace = FALSE,  # NULL,
           ...) {


  wald.stat.vlm(object, values0 = values0,
                subset = subset, omit1s = omit1s, all.out = all.out,
                iterate = iterate,
                trace = trace,
                as.summary = FALSE,  # Does not make sense if TRUE
                score.really = TRUE, ...)
}  # score.stat.vlm





 if (!isGeneric("score.stat"))
   setGeneric("score.stat", function(object, ...)
       standardGeneric("score.stat"))


setMethod("score.stat", "vlm", function(object, ...)
          score.stat.vlm(object, ...))


