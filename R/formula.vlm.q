# These functions are
# Copyright (C) 1998-2023 T.W. Yee, University of Auckland.
# All rights reserved.









assign2assign <- function(vasgn, named = TRUE) {
  asgn <- vasgn
  for (jay in seq_len(length(vasgn)))
    asgn[[jay]] <- asgn[[jay]] * 0 + jay - 1
  if (named) unlist(asgn) else as.vector(unlist(asgn))
}






findterms <- function(Usex.lm, asgn) {

  if (length(Usex.lm) != length(unlist(asgn)))
    stop("two quantities have different number of elements: ",
         length(Usex.lm), " and ",  length(unlist(asgn)))
  nasgn <- names(asgn)
  Col.Usex.lm <- seq_len(length(Usex.lm))[Usex.lm]
  terms.lm <- NULL  # character(0)
  for (jay in seq_len(length(asgn))) {
    if (any(is.element(Col.Usex.lm, asgn[[jay]])))
      terms.lm <- c(terms.lm, nasgn[jay])
  }
  unique(terms.lm)
}






subsetassign <- function(asgn, oTerms) {
  colxptr <- 1
  suba.ptr <- 1
  nasgn <- names(asgn)
  sub.assign <- vector("list", length(oTerms))
  names(sub.assign) <- oTerms
  for (jay in seq_len(length(asgn))) {
    if (is.element(nasgn[jay], oTerms)) {
      lajay <- length(asgn[[jay]])
      sub.assign[[suba.ptr]] <- colxptr:(colxptr + lajay - 1)
      colxptr <- colxptr + lajay  # Next one
      suba.ptr <- suba.ptr + 1
    }
  }
  sub.assign
}  # subsetassign






formula.vlm <- function(x, ...)
  formulavlm(x, ...)



formulavlm <- function(x, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")

  if (!any(slotNames(x) == "misc"))
    stop("cannot find slot 'misc'")

  if (form.number == 1) x@misc$formula else x@misc$form2
}



formulaNA.VGAM <- function(x, ...) {
  stop("a formula does not make sense for object 'x'")
}





setMethod("formula", "vlm",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "vglm",
          function(x, ...)
              formulavlm(x = x, ...))



setMethod("formula", "vgam",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "rrvglm",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "qrrvglm",
          function(x, ...)
              formulavlm(x = x, ...))

setMethod("formula", "grc",
          function(x, ...)
              formulavlm(x = x, ...))










variable.namesvlm <- function(object, full = FALSE, ...) {
  qrslot <- object@qr
  if (!length(qrslot$qr)) {
    use.this <- object@x
    if (!length(use.this))
      stop("argument 'object' has empty 'qr' and 'x' slots.")
  } else {
    use.this <- qrslot$qr
  }
  if (full) dimnames(use.this)[[2]] else
  if (object@rank)
    dimnames(use.this)[[2]][seq_len(object@rank)] else
    character(0)
}




variable.namesrrvglm <- function(object, ...) {

  qrslot <- object@qr
  if (!length(qrslot$qr)) {
    use.this <- object@x
    if (!length(use.this))
      stop("argument 'object' has empty 'qr' and 'x' slots.")
  } else {
    use.this <- qrslot$qr
  }
  dimnames(use.this)[[2]]
}



case.namesvlm <- function(object, full = FALSE, ...) {
  w <- weights(object, type="prior")
  use.this <- residuals(object, type = "working")
  if (!length(use.this))
    use.this <- object@x
  if (!length(use.this))
    use.this <- object@y
  if (!length(use.this))
    stop("argument 'object' has empty 'x' and 'y' slots.")
  dn <- dimnames(use.this)[[1]]
  if (full || is.null(w) || NCOL(w) != 1)
    dn else dn[w != 0]
}


setMethod("variable.names", "vlm",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "vglm",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "vgam",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "rrvglm",
          function(object, ...)
              variable.namesrrvglm(object = object, ...))

setMethod("variable.names", "qrrvglm",
          function(object, ...)
              variable.namesvlm(object = object, ...))

setMethod("variable.names", "grc",
          function(object, ...)
              variable.namesvlm(object = object, ...))






setMethod("case.names", "vlm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "vglm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "vgam",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "rrvglm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "qrrvglm",
          function(object, ...)
              case.namesvlm(object = object, ...))

setMethod("case.names", "grc",
          function(object, ...)
              case.namesvlm(object = object, ...))








has.interceptvlm <- function(object, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")


  if (form.number == 1) {
    if (is.numeric(aa <- attr(terms(object), "intercept")))
      as.logical(aa) else
      FALSE
  } else if (form.number == 2) {
    if (is.numeric(aa <- attr(terms(object, form.number = 2),
                              "intercept")))
      as.logical(aa) else
      FALSE
  }
}



if (!isGeneric("has.intercept"))
    setGeneric("has.intercept", function(object, ...)
               standardGeneric("has.intercept"),
               package = "VGAM")


setMethod("has.intercept",  "vlm", function(object, ...)
           has.interceptvlm(object, ...))







term.namesvlm <- function(model, form.number = 1, ...) {
  if (!is.Numeric(form.number, integer.valued = TRUE,
                  length.arg = 1, positive = TRUE) ||
      form.number > 2)
    stop("argument 'form.number' must be 1 or 2")

    aa <- if (has.intercept(model, form.number = form.number))
          "(Intercept)" else NULL
    bb <- attr(terms(model, form.number = form.number),
               "term.labels")
    c(aa, bb)
}  # term.namesvlm


if (!isGeneric("term.names"))
    setGeneric("term.names", function(model, ...)
               standardGeneric("term.names"),
               package = "VGAM")


setMethod("term.names",  "vlm", function(model, ...)
           term.namesvlm(model, ...))







responseNamevlm <- function(model, form.number = 1, ...) {
  TERMS.MODEL <-terms(model, form.number = form.number)
  if (length(aa <- attr(TERMS.MODEL, "dataClasses")) &&
      length(bb <- attr(TERMS.MODEL, "response"   )) &&
      bb == 1) {
    names(aa)[1]
  } else {
    NULL
  }
}


if (!isGeneric("responseName"))
  setGeneric("responseName", function(model, ...)
             standardGeneric("responseName"),
             package = "VGAM")


setMethod("responseName",  "vlm", function(model, ...)
           responseNamevlm(model, ...))









dftermsvglm <- function(model, term, ...) {

  if (!missing(term) && 1 == length(term)) {
    assign <- attr(model.matrix(model, type =  "lm"), "assign")
    assign <- unlist(lapply(assign, length))
    ind5 <- which(names(assign) == "(Intercept)")
    if (length(ind5) > 0)
      assign <- assign[-ind5]
    which.term <- which(term == labels(terms(model)))
    if (0 == length(which.term))
      stop(paste(term, "is not in the model."))
    Mbyterm <- unlist(lapply(constraints(model, type = "term"),
                             ncol))  # Includes any intercept
    ind5 <- which(names(Mbyterm) == "(Intercept)")
    if (length(ind5) > 0)
      Mbyterm <- Mbyterm[-ind5]
    answer <- assign[which.term] * Mbyterm[which.term]
    answer
  } else {
    terms <- if (missing(term)) labels(terms(model)) else term
    result <- numeric(0)
    for (use.term in terms)
      result <- c(result, Recall(model, term = use.term))
    names(result) <- terms
    result
  }
}  # dftermsvglm



if (!isGeneric("dfterms"))
  setGeneric("dfterms", function(model, term, ...)
             standardGeneric("dfterms"),
             package = "VGAM")


setMethod("dfterms",  "vglm", function(model, term, ...)
           dftermsvglm(model, term, ...))


















drop1.vglm <-
  function(object, scope,
           test = c("none", "LRT"), k = 2, ...) {
  test <- match.arg(test)
  x.lm <- model.matrix(object, type = "lm")
  x.vlm <- model.matrix(object, type = "vlm")
  p.lm <- ncol(x.lm)
  p.vlm <- ncol(x.vlm)
  n.lm <- nobs(object, type = "lm")
  asgn0 <- attr(x.lm, "orig.assign.lm")  # attr(x.lm, "assign")
  if (!length(asgn0))
    stop("could not obtain attribute 'orig.assign.lm' from ",
         "the model matrix; try vglm(..., x = TRUE) and rerun")

  tlab <- attr(terms(object), "term.labels")
  if (missing(scope)) 
    scope <- drop.scope(object) else {
      if (!is.character(scope))
        scope <- attr(terms(update.formula(object, scope)),
                      "term.labels")
      if (!all(match(scope, tlab, 0L) > 0L)) 
        stop("scope is not a subset of term labels")
  }
  ndrop <- match(scope, tlab)
  ns <- length(scope)
  rdf <- df.residual(object)
  chisq <- deviance(object)  # Might be NULL for VGAM
  has.deviance <- !is.null(chisq) && is.finite(chisq)
  dfs <- dev <- numeric(ns)
  dev <- rep(NA_real_, ns)
  llv <- numeric(ns)  # This is new; aka llvec


  M <- npred(object)
  mf <- model.frame(object)
  mt <- attr(mf, "terms")
  OOO <- object@offset
  if (!length(OOO) || all(OOO == 0))
    OOO <- matrix(0, n.lm, M)
  Xm2 <- model.matrix(object, type = "lm2")  # May be 0 x 0
  if (!length(Xm2))
     Xm2 <- NULL  # Make sure. This is safer
  LPmat <- predict(object)
  Fam <- object@family


  Y <- model.response(mf)
  if (!is.factor(Y))
    Y <- as.matrix(Y)
  Wts <- model.weights(mf)
  if (length(Wts) == 0L)
    Wts <- rep(1, n.lm)  # Safest (uses recycling and is a vector)


  big.clist.lm <- constraints(object, type = "lm")
  big.clist.term <- constraints(object, type = "term")
  ncolHlist.lm <- unlist(lapply(big.clist.lm, ncol))
  big.x.lm <- x.lm
  big.x.vlm <- x.vlm
  asgn <- attr(big.x.lm, "assign")  # \pkg{VGAM}
  vasgn <- attr(big.x.vlm, "vassign")  # \pkg{VGAM}
  for (i in seq_len(ns)) {
    ii <- seq_along(asgn0)[asgn0 == ndrop[i]]
    kay.lm <- setdiff(seq(p.lm), ii)
   vecTF <- rep_len(FALSE, p.lm)
   vecTF[kay.lm] <- TRUE



    fit1 <- NULL  # To avoid an warning on CRAN
    oTerms.wint <- voTerms.wint <- ousex.lm <- ousex.vlm <- NULL
    vecTF <- vecTF
    eval(fitmodel.VGAM.expression)


    dfs[i] <- fit1$rank  # Okay
    if (length(tmp5 <- fit1$crit.list$deviance))
      dev[i] <- tmp5
    if (length(tmp5 <- fit1$crit.list$loglikelihood))
      llv[i] <- tmp5  # Almost always okay
  }  # for i


  scope <- c("<none>", scope)
  dfs <- c(object@rank, dfs)
  dev <- c(if (has.deviance) chisq else NA_real_, dev)
  llv <- c(logLik(object), llv)
  dispersion <- 1
  logLIK <- if (has.deviance) dev / dispersion else -2 * llv
  aIC <- logLIK + k * dfs
  dfs <- dfs[1L] - dfs
  dfs[1L] <- NA

  aIC <- aIC + (extractAIC(object, k = k)[2L] - aIC[1L])
   aod <- if (has.deviance)
   data.frame(Df = dfs, Deviance = dev, AIC = aIC,
              row.names = scope, check.names = FALSE) else
   data.frame(Df = dfs, logLik = llv,   AIC = aIC,
                row.names = scope, check.names = FALSE)
                                                            
  if (all(is.na(aIC)))
    aod <- aod[, -3]

  if (test == "LRT") {
    devchange <- pmax(0, logLIK - logLIK[1L])
    devchange[1L] <- NA

    safe_pchisq <- function(q, df, ...) {  # From \pkg{stats}
      df[df <= 0] <- NA
      pchisq(q = q, df = df, ...)
    }

    nas <- !is.na(devchange)
    LRT <- if (dispersion == 1) "LRT" else "scaled dev."
    aod[, LRT] <- devchange
    devchange[nas] <- safe_pchisq(devchange[nas], aod$Df[nas],
                                  lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- devchange
  }

  head <- c("Single term deletions", "\nModel:",
            deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}  # drop1.vglm








extractAIC.vglm <- function(fit, scale = 0, k = 2, ...) {
  n.vlm <- nobs(fit, type = "vlm")
  edf <- n.vlm - df.residual(fit)
  aic <- AIC(fit)
  c(edf, aic + (k - 2) * edf)
}








add1.vglm <-
  function(object, scope,  # scale = 0,
           test = c("none", "LRT"),
           k = 2,
           ...) {
  x <- NULL  # Make life easy for now
  big.x.vlm <- x  # Input could be a very large x.vlm
  rm(x)

  test <- match.arg(test)
  if (!is.character(scope)) 
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  oTerms <- attr(terms(object), "term.labels")
  int <- attr(terms(object), "intercept")
  oTerms.wint <- c(if (int) "(Intercept)" else NULL, oTerms)
  ns <- length(scope)
  llv <- dfs <- dev <- rep(NA_real_, ns + 1)  # numeric(ns + 1)
  names(llv) <- names(dfs) <- names(dev) <- c("<none>", scope)
  add.rhs <- paste(scope, collapse = "+")
  add.rhs <- eval(parse(text = paste("~ . +", add.rhs),
                        keep.source = FALSE))
  new.form <- update.formula(object, add.rhs)
  Terms <- terms(new.form)  # Defines a big model





  n.lm <- nobs(object, type = "lm")
  M <- npred(object)
  mf <- model.frame(object)
  mt <- attr(mf, "terms")
  OOO <- object@offset
  if (!length(OOO) || all(OOO == 0))
    OOO <- matrix(0, n.lm, M)
  Xm2 <- model.matrix(object, type = "lm2")  # May be 0 x 0
  if (!length(Xm2))
    Xm2 <- NULL  # Make sure. This is safer
  LPmat <- predict(object)
  Fam <- object@family

  
  Y <- model.response(mf)
  if (!is.factor(Y))
    Y <- as.matrix(Y)
  Wts <- model.weights(mf)
  if (length(Wts) == 0) 
    Wts <- rep.int(1, n.lm)

  if (is.null(big.x.vlm)) {  # Usually the case

  listcall <- as.list(object@call)
  argslist <- vector("list", length(listcall) - 1)
  for (kay in 2:(length(listcall)))
    argslist[[kay - 1]] <- listcall[[kay]]
  names(argslist) <- c(names(listcall)[-1])
  argslist$formula <- Terms  # A big revised model
    bigfit <- do.call(what = "vglm", args = argslist)

    big.clist.lm <- constraints(bigfit, type = "lm")
    big.clist.term <- constraints(bigfit, type = "term")
    big.x.lm <- model.matrix(bigfit, type = "lm")
    big.x.vlm <- model.matrix(bigfit, type = "vlm")
  ncolHlist.lm <- unlist(lapply(big.clist.lm, ncol))
  ncolHlist.term <- unlist(lapply(big.clist.term, ncol))
  } else {  # Not often the case
  }

  asgn <- attr(big.x.lm, "assign")  # \pkg{VGAM}
  vasgn <- attr(big.x.vlm, "vassign")  # \pkg{VGAM}

  voTerms.wint <-
    vlabel(oTerms.wint, M = M,
           unlist(lapply(big.clist.term[oTerms.wint], ncol)))

  
  tlab <- attr(Terms, "term.labels")  # Terms <-
 

  ousex.lm <- unlist(asgn[oTerms])
  if (int)
    ousex.lm <- c("(Intercept)" = 1, ousex.lm)
  ousex.vlm <- unlist(vasgn[voTerms.wint])  # oTerms zz


  assign0.lm <- subsetassign(asgn, oTerms.wint)
  assign0.vlm <- subsetassign(vasgn, voTerms.wint)
  use.x.lm <- big.x.lm[, ousex.lm, drop = FALSE]
  attr(use.x.lm, "assign") <- assign0.lm
  use.x.vlm <- big.x.vlm[, ousex.vlm, drop = FALSE]
  attr(use.x.vlm, "vassign") <- assign0.vlm


  R.asgn.lm <- assign2assign(asgn)  # length(coef(bigfit)) elts
  orig.assign.lm.lm <- unlist(attr(big.x.lm, "orig.assign.lm"))
  if (length(orig.assign.lm.lm) &&  # It could be NULL, e.g., x = F
      max(abs(orig.assign.lm.lm - as.vector(R.asgn.lm))) > 0) {
    warning("difference found between the original ",  # FYI
            "'assign' attributes. Using the safest choice.")
    R.asgn.lm <- orig.assign.lm.lm  # Safest choice
  }



   if (is.logical(object@control$trace))
      object@control$trace <- FALSE  # Supress 'trace'; keep silent
    prewarn <- options("warn")
    options(warn = -1)  # Supress warnings
    fit0 <- vglm.fit(x = use.x.lm,  # Not really used much
                     y = Y, w = c(Wts),
                     X.vlm.arg = use.x.vlm,
                     Xm2 = Xm2, Terms = mt,
                constraints = big.clist.term[oTerms.wint],  # zz
                     extra = object@extra,
                     etastart = LPmat,
                     offset = OOO, family = Fam,
                     control = object@control)
    options(warn = prewarn[["warn"]])  # Restore warnings

  
  dfs[1L] <- fit0$rank
  if (length(tmpdev <- tmp5 <- fit0$crit.list$deviance))
    dev[1L] <- tmp5
  if (length(tmp5 <- fit0$crit.list$loglikelihood))
    llv[1L] <- tmp5  # Almost always okay
  has.deviance <- !is.null(tmpdev) && is.finite(tmpdev)






  sTerms <- sapply(strsplit(tlab, ":", fixed = TRUE),
                   function(x) paste(sort(x), collapse = ":"))

  for (tt in scope) {
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
    vecTF <- match(R.asgn.lm, match(stt, sTerms), 0L) > 0L





    fit1 <- NULL  # To avoid an warning on CRAN
    vecTF <- vecTF  # Should be called vecTF.lm really
    eval(fitmodel.VGAM.expression)




    dfs[tt] <- fit1$rank
    if (length(tmpdev <- fit1$crit.list$deviance))
      dev[tt] <- tmpdev
    if (length(tmp5 <- fit1$crit.list$loglikelihood))
      llv[tt] <- tmp5  # Almost always okay
  }  # for (tt in scope)


  dispersion <- 1
  loglik <- if (has.deviance) dev / dispersion else -2 * llv

  aic <- loglik + k * dfs
  aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
  dfs <- dfs - dfs[1L]
  dfs[1L] <- NA
  aod <- if (has.deviance)
   data.frame(Df = dfs, Deviance = dev, AIC = aic,
              row.names = names(dfs), check.names = FALSE) else
   data.frame(Df = dfs, logLik   = llv, AIC = aic,
              row.names = names(dfs), check.names = FALSE)
  if (all(is.na(aic)))  # || !has.deviance
    aod <- aod[, -3]

  safe_pchisq <- function(q, df, ...) {  # From \pkg{stats}
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }

  if (test == "LRT") {
    devchange <- pmax(0, loglik[1L] - loglik)
    devchange[1L] <- NA

    LRT <- if (dispersion == 1) "LRT" else "scaled dev."
    aod[, LRT] <- devchange
    nas <- !is.na(devchange)
    devchange[nas] <- safe_pchisq(devchange[nas], aod$Df[nas],
                                  lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- devchange
  }  # test == "LRT"

  head <- c("Single term additions", "\nModel:",
            deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}  # add1.vglm















