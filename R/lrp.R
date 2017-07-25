# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.













lrp.vglm <-
  function(object, which = NULL,  # 1:p.vlm,
           omit1s = TRUE,
           trace = NULL,
           ...) {




  Pnames <- names(B0 <- coef(object))
  nonA <- !is.na(B0)  # Usually == rep(TRUE, p.vlm)
  if (any(is.na(B0)))
    stop("currently cannot handle NA-valued regression coefficients")
  pv0 <- t(as.matrix(B0))  # 1 x p.vlm



  p.vlm <- length(Pnames)
  if (is.character(which))
    which <- match(which, Pnames)
  if (is.null(which))
    which <- 1:p.vlm
  summ <- summary(object)



  M <- npred(object)
  Xm2 <- model.matrix(object, type = "lm2")  # Could be a 0 x 0 matrix
  if (!length(Xm2))
     Xm2 <- NULL  # Make sure. This is safer
  clist <- constraints(object, type = "lm")  # type = c("lm", "term")
  H1 <- clist[["(Intercept)"]]

  if (omit1s && length(H1) && any(which <= ncol(H1))) {
    if (length(clist) == 1)
      return(NULL)  # Regressed against intercept only
    which <- which[which > ncol(H1)]
  }

  mf <- model.frame(object)

  Y <- model.response(mf)
  if (!is.factor(Y))
    Y <- as.matrix(Y)


  n.lm <- nobs(object, type = "lm")
  OOO <- object@offset
  if (!length(OOO) || all(OOO == 0))
    OOO <- matrix(0, n.lm, M)



  mt <- attr(mf, "terms")




  Wts <- model.weights(mf)
  if (length(Wts) == 0L)
    Wts <- rep(1, n.lm)  # Safest (uses recycling and is a vector)
  Original.de <- deviance(object)  # Could be NULL
  if (!(use.de <- is.Numeric(Original.de)))
    Original.ll <- logLik(object)
  DispersionParameter <- summ@dispersion
  if (!all(DispersionParameter == 1))
    stop("Currently can only handle dispersion parameters ",
         "that are equal to 1")
  X.lm  <- model.matrix(object, type =  "lm")
  X.vlm <- model.matrix(object, type = "vlm")
  fam <- object@family



  LPmat <- predict(object)


  quasi.type <- if (length(tmp3 <- fam@infos()$quasi.type))
    tmp3 else FALSE
  if (quasi.type)
    stop("currently this function cannot handle quasi-type models",
         " or models with an estimated dispersion parameter")


  ansvec <- B0[which]  # Replace these with p-values

  iptr <- 0
  for (i in which) {
    aa <- nonA
    aa[i] <- FALSE
    X.vlm.i <- X.vlm[, aa, drop = FALSE]
    X.lm.i  <-  X.lm  # Try this


 # This is needed by vglm.fit():
    attr(X.vlm.i, "assign") <- attr(X.vlm, "assign")  # zz; this is wrong!
    attr( X.lm.i, "assign") <- attr( X.lm, "assign")


    if (is.logical(trace))
      object@control$trace <- trace



        fm <- vglm.fit(x = X.lm.i,  # Possibly use X.lm.i or else X.lm
                       y = Y, w = Wts,
                       X.vlm.arg = X.vlm.i,  # X.vlm,
                       Xm2 = Xm2, Terms = mt,
                       extra = object@extra,
                       etastart = LPmat,
                       offset = OOO,  # ooo,
                       family = fam,
                       control = object@control)



        zee <- if (use.de) {
          fm$crit.list[["deviance"]] - Original.de
        } else {
          2 * (Original.ll - fm$crit.list[["loglikelihood"]])
        }
        if (zee > -1e-3) {
          zee <- max(zee, 0)
        } else {
          stop("omitting 1 column has found a better solution, ",
               "so the original fit had not converged")
        }
        zedd <- zee  # sgn * sqrt(zee)
     iptr <- iptr + 1
     ansvec[iptr] <- pchisq(zedd, df = 1, lower.tail = FALSE)
  }  # for i

  ansvec
}






if (!isGeneric("lrp"))
    setGeneric("lrp",
               function(object, ...)
               standardGeneric("lrp"),
           package = "VGAM")


setMethod("lrp", "vglm",
          function(object, ...)
          lrp.vglm(object = object, ...))









