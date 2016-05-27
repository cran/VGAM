# These functions are
# Copyright (C) 1998-2016 T.W. Yee, University of Auckland.
# All rights reserved.









vgam <- function(formula, 
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = vgam.control(...),
                 offset = NULL, 
                 method = "vgam.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL,
                 constraints = NULL,
                 extra = list(),
                 form2 = NULL,  # Added 20130730
                 qr.arg = FALSE, smart = TRUE, ...) {
  dataname <- as.character(substitute(data))  # "list" if no data= 
  function.name <- "vgam"

  ocall <- match.call()

  if (smart)
    setup.smart("write")

  if (missing(data))
    data <- environment(formula)

  mtsave <- terms(formula, specials = c("s", "ps"), data = data)



  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
      "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method,
         model.frame = return(mf),
         vgam.fit = 1,
         stop("invalid 'method': ", method))
  mt <- attr(mf, "terms")

  xlev <- .getXlevels(mt, mf)
  y <- model.response(mf, "any")  # model.extract(mf, "response")
  x <- if (!is.empty.model(mt))
         model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
  attr(x, "assign") <- attrassigndefault(x, mt)






  if (!is.null(form2)) {
    if (!is.null(subset))
      stop("argument 'subset' cannot be used when ",
            "argument 'form2' is used")
    retlist <- shadowvgam(formula =
                 form2,
                 family = family, data = data,
                 na.action = na.action,
                 control = vgam.control(...),
                 method = method,
                 model = model, x.arg = x.arg, y.arg = y.arg,
                 contrasts = contrasts,
                 constraints = constraints,
                 extra = extra,
                 qr.arg = qr.arg)
    Ym2 <- retlist$Ym2
    Xm2 <- retlist$Xm2

    if (length(Ym2)) {
      if (nrow(as.matrix(Ym2)) != nrow(as.matrix(y)))
        stop("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2)) {
      if (nrow(as.matrix(Xm2)) != nrow(as.matrix(x)))
        stop("number of rows of 'x' and 'Xm2' are unequal")
    }
  } else {
    Xm2 <- Ym2 <- NULL
  }







  offset <- model.offset(mf)
  if (is.null(offset))
    offset <- 0  # yyy ???




  mf2 <- mf
  if (!missing(subset)) {
    mf2$subset <- NULL 
    mf2 <- eval(mf2, parent.frame())  # mf2 is the full data frame. 
    spars2 <-  lapply(mf2, attr, "spar") 
    dfs2   <-  lapply(mf2, attr, "df") 
    sx2 <-  lapply(mf2, attr, "s.xargument") 
    for (ii in seq_along(mf)) {
      if (length(sx2[[ii]])) {
        attr(mf[[ii]], "spar") <- spars2[[ii]]
        attr(mf[[ii]], "dfs2") <- dfs2[[ii]]
        attr(mf[[ii]], "s.xargument") <- sx2[[ii]]
      }
    }
    rm(mf2) 
  }



  w <- model.weights(mf)
  if (!length(w)) {
    w <- rep_len(1, nrow(mf))
  } else if (ncol(as.matrix(w)) == 1 && any(w < 0))
    stop("negative weights not allowed")




  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }

  eval(vcontrol.expression)

  n <- dim(x)[1]


  if (length(slot(family, "first")))
    eval(slot(family, "first"))


  aa <- attributes(mtsave)
  smoothers <- aa$specials

  mgcv.ps <- length(smoothers$ps) > 0
  mgcv.PS <- length(smoothers$PS) > 0
  any.ps.terms <- mgcv.ps || mgcv.PS
  mgcv.s <- length(smoothers$s) > 0
  if (any.ps.terms && mgcv.s)
    stop("cannot include both s() and ps() (or PS()) terms in the formula")




  nonparametric <- length(smoothers$s) > 0
  if (nonparametric) {

      ff <- apply(aa$factors[smoothers[["s"]],,drop = FALSE], 2, any)
      smoothers[["s"]] <-
        if (any(ff)) seq(along = ff)[aa$order == 1 & ff] else NULL

    smooth.labels <- aa$term.labels[unlist(smoothers)]
  } else {
    function.name <- "vglm"  # This is effectively so 
  }
  






  are.ps.terms <- (length(smoothers$ps) + length(smoothers$PS)) > 0
  if (are.ps.terms) {
    control$criterion <- "coefficients"  # Overwrite if necessary

    if (length(smoothers$ps) > 0) {
      ff.ps <- apply(aa$factors[smoothers[["ps"]],,drop = FALSE], 2, any)
      smoothers[["ps"]] <-
        if (any(ff.ps)) seq(along = ff.ps)[aa$order == 1 & ff.ps] else NULL
      smooth.labels <- aa$term.labels[unlist(smoothers)]
    }






    assignx <- attr(x, "assign")
    which.X.ps <- assignx[smooth.labels]
    data <- mf[, names(which.X.ps), drop = FALSE]
    attr(data, "class") <- NULL
    S.arg <- lapply(data, attr, "S.arg")
    lambdalist <- lapply(data, attr, "lambda")
    ridge.adj <- lapply(data, attr, "ridge.adj")
    term.labels <- aa$term.labels


  }



  ps.list <- if (any.ps.terms)
               list(indexterms = ff.ps,
                    intercept = aa$intercept,
                    which.X.ps = which.X.ps,
                    S.arg = S.arg,
                    lambdalist = lambdalist,
                    ridge.adj = ridge.adj,
                    term.labels = term.labels,
                    assignx = assignx) else
               NULL


  fit <- vgam.fit(x = x, y = y, w = w, mf = mf,
                  Xm2 = Xm2, Ym2 = Ym2,  # Added 20130730
      etastart = etastart, mustart = mustart, coefstart = coefstart,
      offset = offset, family = family, control = control,
      criterion = control$criterion,
      constraints = constraints, extra = extra, qr.arg = qr.arg,
      Terms = mtsave,
      nonparametric = nonparametric, smooth.labels = smooth.labels,
      function.name = function.name,
      ps.list = ps.list,
      ...)


  if (is.Numeric(fit$nl.df) && any(fit$nl.df < 0)) {
    fit$nl.df[fit$nl.df < 0] <- 0
  }




  if (!is.null(fit[["smooth.frame"]])) {
    fit <- fit[-1]  # Strip off smooth.frame
  } else {
  }

  fit$smomat <- NULL  # Not needed

  fit$call <- ocall 
  if (model)
    fit$model <- mf 
  if (!x.arg)
    fit$x <- NULL
  if (!y.arg)
    fit$y <- NULL

  if (nonparametric)
    fit$misc$smooth.labels <- smooth.labels


  fit$misc$dataname <- dataname


  if (smart)
    fit$smart.prediction <- get.smart.prediction()






  answer <-
  new(
    if (any.ps.terms) "psvgam" else "vgam",

    "assign"       = attr(x, "assign"),
    "call"         = fit$call,
    "coefficients" = fit$coefficients,
    "constraints"  = fit$constraints,
    "criterion"    = fit$crit.list,
    "df.residual"  = fit$df.residual,
    "dispersion"   = 1,
    "family"       = fit$family,
    "misc"         = fit$misc,
    "model"        = if (model) mf else data.frame(),
    "R"            = fit$R,
    "rank"         = fit$rank,
    "residuals"    = as.matrix(fit$residuals),
    "ResSS"       = fit$ResSS,
    "smart.prediction" = as.list(fit$smart.prediction),
    "terms"        = list(terms = fit$terms))

  if (!smart)
    answer@smart.prediction <- list(smart.arg = FALSE)

  if (qr.arg) {
    class(fit$qr) <- "list"
    slot(answer, "qr") <- fit$qr
  }
  if (length(attr(x, "contrasts")))
    slot(answer, "contrasts") <- attr(x, "contrasts")
  if (length(fit$fitted.values))
    slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action")))
    list(aaa) else list()
  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)
  if (length(fit$weights))
    slot(answer, "weights") <- as.matrix(fit$weights)


  if (x.arg)
    slot(answer, "x") <- x  # The 'small' design matrix



  if (length(fit$misc$Xvlm.aug)) {
    slot(answer, "psslot") <-
      list(Xvlm.aug = fit$misc$Xvlm.aug,
           ps.list  = fit$misc$ps.list,
           magicfit = fit$misc$magicfit)
    fit$misc$Xvlm.aug <- NULL
    fit$misc$ps.list  <- NULL
    fit$misc$magicfit <- NULL
  }



  if (x.arg && length(Xm2))
    slot(answer, "Xm2") <- Xm2  # The second (lm) design matrix
  if (y.arg && length(Ym2))
    slot(answer, "Ym2") <- as.matrix(Ym2)  # The second response
  if (!is.null(form2))
    slot(answer, "callXm2") <- retlist$call
  answer@misc$formula <- formula
  answer@misc$form2   <- form2



  if (length(xlev))
    slot(answer, "xlevels") <- xlev
  if (y.arg)
    slot(answer, "y") <- as.matrix(fit$y)
  answer@misc$formula <- formula





  slot(answer, "control") <- fit$control

  if (length(fit$extra)) {
    slot(answer, "extra") <- fit$extra
  }
  slot(answer, "iter")   <- fit$iter
  slot(answer, "post")   <- fit$post

  fit$predictors <- as.matrix(fit$predictors)  # Must be a matrix
  dimnames(fit$predictors) <- list(dimnames(fit$predictors)[[1]],
                                   fit$misc$predictors.names)
  slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights))
    slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)


  if (nonparametric) {
    slot(answer, "Bspline") <- fit$Bspline
    slot(answer, "nl.chisq") <- fit$nl.chisq
    if (is.Numeric(fit$nl.df))
      slot(answer, "nl.df") <- fit$nl.df
    slot(answer, "spar") <- fit$spar
    slot(answer, "s.xargument") <- fit$s.xargument
    if (length(fit$varmat)) {
      slot(answer, "var") <- fit$varmat
    }




  }
  if (length(fit$effects))
    slot(answer, "effects") <- fit$effects





  if (nonparametric && is.buggy.vlm(answer)) {
    warning("some s() terms have constraint matrices that have columns",
            " which are not orthogonal;",
            " try using ps() instead of s().")
  } else {
  }




  answer
}
attr(vgam, "smart") <- TRUE 





shadowvgam <-
        function(formula,
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = vgam.control(...), 
                 offset = NULL, 
                 method = "vgam.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL, 
                 constraints = NULL,
                 extra = list(), 
                 qr.arg = FALSE, ...) {
    dataname <- as.character(substitute(data))  # "list" if no data=
    function.name <- "shadowvgam"

    ocall <- match.call()

    if (missing(data)) 
        data <- environment(formula)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), vgam.fit = 1,
           stop("invalid 'method': ", method))
    mt <- attr(mf, "terms")

    x <- y <- NULL 

    xlev <- .getXlevels(mt, mf)
    y <- model.response(mf, "any")  # model.extract(mf, "response")
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
    attr(x, "assign") <- attrassigndefault(x, mt)

    list(Xm2 = x, Ym2 = y, call = ocall)
}








is.buggy.vlm <- function(object, each.term = FALSE, ...) {


    
  Hk.list <- constraints(object)
  ncl <- names(Hk.list)
  TFvec <- rep_len(FALSE, length(ncl))
  names(TFvec) <- ncl



  if (!is(object, "vgam")) {
    return(if (each.term) TFvec else any(TFvec))
  }
  if (!length(object@nl.chisq)) {
    return(if (each.term) TFvec else any(TFvec))
  }

  for (kay in seq_along(ncl)) {
    cmat <- Hk.list[[kay]]
    if (ncol(cmat) > 1 && substring(ncl[kay], 1, 2) == "s(") {
      CMat <- crossprod(cmat)  # t(cmat) %*% cmat
      TFvec[kay] <- any(CMat[lower.tri(CMat)] != 0 |
                        CMat[upper.tri(CMat)] != 0)
    }
  }
  if (each.term) TFvec else any(TFvec)
}



if (!isGeneric("is.buggy"))
  setGeneric("is.buggy", function(object, ...)
             standardGeneric("is.buggy"),
             package = "VGAM")



setMethod("is.buggy", signature(object = "vlm"),
          function(object, ...)
          is.buggy.vlm(object, ...))



