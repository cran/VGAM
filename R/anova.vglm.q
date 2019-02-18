# These functions are
# Copyright (C) 1998-2019 T.W. Yee, University of Auckland.
# All rights reserved.














car.relatives <- function(term, names, factors) {
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2])))
  }
  if(length(names) == 1) return(NULL)
  which.term <- which(term == names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
          function(term2) is.relative(term, term2))]
}













anova.vglm <-
  function(object, ...,
           type = c("II", "I", "III", 2, 1, 3),
           test = c("LRT", "none"),  #
           trydev = TRUE,  # Use where possible?
           silent = TRUE
          ) {


  type <- as.character(type)
  type <- match.arg(type, c("II", "I", "III", "2", "1", "3"))
  type[type == "1"] <- "I"
  type[type == "2"] <- "II"
  type[type == "3"] <- "III"
  if (has.intercept(object) && length(constraints(object)) == 1 && 
     names(constraints(object)) == "(Intercept)" && 
     type == "II") {
    type <- "III"
    warning("the model contains only an intercept; ",
            "Type III test substituted")
  }

  if (length( list(...) ) && type != "I")
    stop("argument 'type' must 'I' or 1 for multiple fits")



  dispersion <- 1


  if (mode(test) != "character" && mode(test) != "name")
    test <- as.character(substitute(test))
  test <- match.arg(test, c("LRT", "none"))[1]  # , "Rao"
  test.null <- if (test == "none") NULL else test



  if (!has.intercept(object))
    stop("argument 'object' must have an intercept term")
  object@control$trace <- FALSE
  has.deviance <- !is.null(deviance(object)) && trydev
  if (silent) {
    warn.save <- options()$warn
    options(warn = -1)  # Negative means ignore all warnings
  }


  dotargs <- list(...)
  named <- if (is.null(names(dotargs))) 
    rep_len(FALSE, length(dotargs)) else
    (names(dotargs) != "")
  if (any(named)) 
    warning("the following arguments to 'anova.vglm' are ",
            "invalid and dropped: ", 
        paste(deparse(dotargs[named]), collapse = ", "))



  dotargs <- dotargs[!named]
  is.vglm <- vapply(dotargs, function(x) is(x, "vglm"), NA)
  dotargs <- dotargs[is.vglm]
  if (length(dotargs)) 
    return(anova.vglmlist(c(list(object), dotargs),
                            dispersion = dispersion, 
                            test = test,
                            type = type,
                            .has.deviance = has.deviance,
                            .trydev = trydev
                          ))




  doscore <- !is.null(test.null) && test == "Rao"
  if (doscore) {
    warning("cannot have doscore == TRUE. Setting it FALSE")
    doscore <- FALSE
  }
  varlist <- attr(terms(object), "variables")
  x.lm <- model.matrix(object, type = "lm")
  x.vlm <- model.matrix(object, type = "vlm")
  p.vlm <- ncol(x.vlm)  # Needed for type == "III"
  varseq <- rep(seq(length(attr(x.vlm, "assign"))),
                times = unlist(lapply(attr(x.vlm, "assign"), length))) - 1
  nvars <- max(0, varseq)
  resdev  <- resdf  <- reslogLik  <- NULL
  resdev2 <- resdf2 <- reslogLik2 <- NULL  # For type = "II"



  n.lm <- nobs(object, type = "lm")
  M <- npred(object)
  mf <- model.frame(object)
  mt <- attr(mf, "terms")
  Y <- model.response(mf)
  if (!is.factor(Y))
    Y <- as.matrix(Y)
  Wts <- model.weights(mf)
  if (length(Wts) == 0L)
    Wts <- rep(1, n.lm)  # Safest (uses recycling and is a vector)
  OOO <- object@offset
  if (!length(OOO) || all(OOO == 0))
    OOO <- matrix(0, n.lm, M)
  Xm2 <- model.matrix(object, type = "lm2")  # Could be a 0 x 0 matrix
  if (!length(Xm2))
     Xm2 <- NULL  # Make sure. This is safer
  clist <- constraints(object, type = "term")  # type = c("lm", "term")
  ooo <- OOO  # + matrix(X.vlm[, i] * betai, n.lm, M, byrow = TRUE)
  LPmat <- predict(object)
  Fam <- object@family






  junkval <- NA_real_  # Any junk value will do (makes no difference)
  if (is.element(type, c("I"))) {
    fit0 <- vglm.fit(x = x.lm,  # [, 1, drop = FALSE],
                     y = Y, w = Wts,
                     X.vlm.arg = x.vlm[, varseq <= 0, drop = FALSE],
                     Xm2 = Xm2, Terms = mt,
                     constraints = clist, extra = object@extra,
                     etastart = LPmat,
                     offset = ooo, family = Fam,
                     control = object@control)
    object.df.null <- fit0$df.residual
    if (has.deviance) {
      object.null.deviance <- fit0$crit.list$deviance
      object.null.logLik   <- fit0$crit.list$loglikelihood
    } else {
      object.null.logLik <- fit0$crit.list$loglikelihood
    }
  }  # type == "I"


  if (type == "II") {
    x.lm.vglm <- model.matrix(formula(object),
                              data = get(object@misc$dataname))
    attr.x.lm.vglm <- attr(x.lm.vglm, "assign")
    Names <- term.names(object)[-1]  # Delete the intercept.
  }  # type == "II"



  if (doscore) {
    score <- numeric(nvars)
    method <- "vglm.fit"
    y <- Y <- depvar(object)
    if (FALSE)
    fit <- eval(call(if (is.function(method)) "method" else method, 
                     x = x.vlm[, varseq == 0, drop = FALSE],
                     y = y, weights = object$prior.weights, 
                     start = object$start,
                     offset = object$offset,
                     family = object$family, 
                     control = object$control))

       fit <- vglm.fit(x = x.lm,  # Possibly use X.lm.i or else X.lm
                       y = Y, w = Wts,
                       X.vlm.arg = x.vlm,  # X.vlm,
                       Xm2 = Xm2, Terms = mt,
                       constraints = clist, extra = object@extra,
                       etastart = LPmat,
                       offset = ooo, family = Fam,
                       control = object@control)

    r <- fit$residuals
    w <- fit$weights
  }  # doscore


  if (doscore || nvars > 1) {
    method <- "vglm.fit"  # object$method
    y <- depvar(object)


    if (type == "II") {
    }


    upp.bnd <- switch(type, "I" = nvars - 1L, "II" = , "III" = nvars)


    for (ii in seq_len(upp.bnd)) {
      if (type == "II") {

        asgn <- attr.x.lm.vglm
        which.nms <- function(name)
          which(asgn == which(Names == name))
        index3 <-
          car.relatives(term    = Names[ii],
                        names   = Names,
                        factors = attr(terms(object), "factors"))
        rels <- Names[index3]
        exclude.1 <- as.vector(unlist(sapply(c(Names[ii], rels),
                                             which.nms)))
        exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))

        attr.assign.object <- attr.assign.x.vglm(object)  # Easy!
        vecTF1 <- is.element(attr.assign.object, asgn[exclude.1])
        mod.1 <- colnames(x.vlm)[!vecTF1]
        vecTF2 <- is.element(attr.assign.object, asgn[exclude.2])
        if (length(rels) == 0)
          vecTF2[] <- FALSE  # i.e., so fit2 == object


        mod.2 <- colnames(x.vlm)[!vecTF2]
        fit2 <- vglm.fit(x = x.lm,
                         y = Y, w = Wts,
                         X.vlm.arg = x.vlm[, !vecTF2, drop = FALSE],
                         Xm2 = Xm2, Terms = mt,
                         constraints = clist, extra = object@extra,
                         etastart = LPmat,
                         offset = ooo, family = Fam,
                         control = object@control)
      }  # type == "II"


      vecTF <- switch(type,
                      "I"   = (varseq <= ii),
                      "II"  = !vecTF1,
                      "III" = (varseq != ii))
      fit <- vglm.fit(x = x.lm,  # Possibly use X.lm.i or else X.lm
                      y = Y, w = Wts,
                      X.vlm.arg = x.vlm[, vecTF, drop = FALSE],
                      Xm2 = Xm2, Terms = mt,
                      constraints = clist, extra = object@extra,
                      etastart = LPmat,
                      offset = ooo, family = Fam,
                      control = object@control)


      if (doscore) {
        zedd <- eval(call(if (is.function(method)) "method" else method, 
                   x = x.vlm[, vecTF, drop = FALSE], y = r, 
                   weights = w))
        score[ii] <- zedd$null.deviance - zedd$deviance
        r <- fit$residuals
        w <- fit$weights
      }  # doscore


      if (type == "I") {
        reslogLik <- c(reslogLik, fit$crit.list$loglik)
        resdev <- c(resdev, fit$crit.list$deviance)  # May be NULL
        resdf <- c(resdf, fit$df.residual)
      } else
      if (type == "II") {
        reslogLik2 <- c(reslogLik2, fit2$crit.list$loglik)
        resdev2 <- c(resdev2, fit2$crit.list$deviance)  # May be NULL
        resdf2 <- c(resdf2, fit2$df.residual)
        reslogLik <- c(reslogLik, fit$crit.list$loglik)
        resdev <- c(resdev, fit$crit.list$deviance)  # May be NULL
        resdf <- c(resdf, fit$df.residual)
      } else {
        reslogLik <- c(reslogLik, fit$crit.list$loglik)
        resdev <- c(resdev, fit$crit.list$deviance)  # May be NULL
        resdf <- c(resdf, fit$df.residual)
      } 
    }  # for ii




    if (doscore) {
      zedd <- eval(call(if (is.function(method)) "method" else method, 
                 x = x.vlm, y = r, weights = w))
                 score[nvars] <- zedd$null.deviance - zedd$deviance
    }  # doscore
  }  # if (nvars > 1 || doscore)




  if (is.element(type, c("I"))) {
    resdf <- c(object.df.null, resdf, df.residual(object))
    if (has.deviance) {
      resdev <- c(object.null.deviance, resdev, deviance(object))
    } else {
      reslogLik <- c(object.null.logLik, reslogLik, logLik(object))
    }
  }  # type == "I"

  if (is.element(type, c("II", "III"))) {
    resdf <- c(object.df.null = junkval, resdf)
    if (has.deviance) {
      resdev <- c(object.null.deviance = junkval, resdev)
    } else {
      reslogLik <- c(object.null.logLik = junkval, reslogLik)
    }
  }  # type == "II" or "III"

  if (is.element(type, c("II"))) {
    resdf2 <- c(object.df.null = junkval, resdf2)
    if (has.deviance) {
      resdev2 <- c(object.null.deviance = junkval, resdev2)
    } else {
      reslogLik2 <- c(object.null.logLik = junkval, reslogLik2)
    }
  }  # type == "II"





  table <- if (has.deviance) {
             if (type == "I")
               data.frame(c(NA, -diff(resdf)),
                          c(NA, pmax(0, -diff(resdev))), 
                          resdf, resdev) else
             if (type == "II")
               data.frame(c(NA, resdf[-1] - resdf2[-1]),
                          c(NA, pmax(0, resdev[-1] - resdev2[-1])), 
                          resdf, resdev) else
               data.frame(c(NA, resdf[-1] - df.residual(object)),
                          c(NA, pmax(0, resdev[-1] - deviance(object))), 
                          resdf, resdev)
           } else {
             if (type == "I")
               data.frame(c(NA, -diff(resdf)),
                          c(NA, pmax(0, 2 * diff(reslogLik))), 
                          resdf, reslogLik) else
             if (type == "II")
               data.frame(c(NA, resdf[-1] - resdf2[-1]),
                          c(NA, pmax(0, 2 * (reslogLik2[-1] -
                                             reslogLik[-1]))), 
                          resdf, reslogLik) else
               data.frame(c(NA, resdf[-1] - df.residual(object)),
                          c(NA, pmax(0,2*(logLik(object)-reslogLik[-1]))),
                          resdf, reslogLik)
           }
  tl <- attr(terms(object), "term.labels")




  if (length(tl) == 0L)
    table <- table[1, , drop = FALSE]
  if (1 + length(tl) != nrow(table))
    stop("cannot apply dimnames to 'table'")
  dimnames(table) <-
    list(c("NULL", tl),
         c("Df", "Deviance", "Resid. Df", "Resid. Dev"))




  if (doscore) 
    table <- cbind(table, Rao = c(NA, score))



  suptitle <- if (type == "I") paste("Type I tests: terms added ",
              "sequentially from\nfirst to last", sep = "") else
     if (type == "II")
       "Type II tests" else
       "Type III tests: each term added last"
  title <- paste0("Analysis of Deviance Table (", suptitle, ")",
     "\n\nModel: ",
     paste(paste("'", Fam@vfamily, "'", sep = ""),
           collapse = ", "),
     "\n\nLinks: ",
     paste(paste("'", linkfun(object), "'", sep = ""),
           collapse = ", "),
     "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n\n")
  df.dispersion <- Inf
  if (is.null(dispersion)) {
    stop("cannot handle a NULL 'dispersion'")
  }




  if (!is.null(test.null)) {
    table <- stat.anova(table = table, test = test.null,
                        scale = dispersion, 
                        df.scale = df.dispersion)  #, n = NROW(x)
  }  # (!is.null(test.null))
  if (!has.deviance)
  dimnames(table) <-
    list(c("NULL", tl),
         c("Df", "2 * LogLik Diff.", "Resid. Df", "LogLik",
           if (!is.null(test.null) && doscore) "Rao" else NULL,
           if (is.null(test.null)) NULL else "Pr(>Chi)"))



  if (is.element(type, c("II", "III")) && has.intercept(object)) {
    table <- table[-1, , drop = FALSE]

    index.col <- which(is.element(colnames(table),
      c("Resid. Df", if (has.deviance) "Resid. Dev" else "LogLik")))
    table <- table[, -index.col, drop = FALSE]
  }
  if (silent)
    options(warn = warn.save)  # Restore 'warn'.


  structure(table, heading = title, class = c("anova", "data.frame"))
}  # anova.vglm


    



anova.vglmlist <-
  function(object, ...,
           type = "I",  # c("I", "II","III", 1, 2, 3),
           test = c("LRT", "none"),  #
           .has.deviance = FALSE,
           .trydev = TRUE
          ) {



  type <- as.character(type)
  type <- match.arg(type, c("I", "II","III", "1", "2", "3"))
  type[type == "1"] <- "I"
  type[type == "2"] <- "II"
  type[type == "3"] <- "III"
  if (type != "I")
    stop("argument 'type' must be 'I' since there are several fits")







  if (mode(test) != "character" && mode(test) != "name")
    test <- as.character(substitute(test))
  test <- match.arg(test, c("LRT", "none"))[1]  # , "Rao"
  test.null <- if (test == "none") NULL else test


  doscore <- !is.null(test.null) && test == "Rao"
  if (doscore) {
    warning("cannot have doscore == TRUE. Setting it FALSE")
    doscore <- FALSE
  }


  responses <- as.character(lapply(object, function(x) {
      deparse(formula(x)[[2L]])
  }))
  sameresp <- responses == responses[1L]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(gettextf("models with response %s removed because ",
                     "response differs from model 1", 
        sQuote(deparse(responses[!sameresp]))), domain = NA)
  }



  ns <- as.numeric(lapply(object, function(x) nobs(x)))
  if (any(ns != ns[1L])) 
    stop("models were not all fitted to the same size of dataset")


  nmodels <- length(object)
  if (nmodels == 1) 
    return(anova.vglm(object[[1L]],
                      test = test.null))
  resdf <- as.numeric(lapply(object, function(x) df.residual(x)))
  reslogLik <- as.numeric(lapply(object, function(x) logLik(x)))
  if (.has.deviance && .trydev)
    resdev <- as.numeric(lapply(object, function(x) deviance(x)))


  if (doscore) {
    score <- numeric(nmodels)
    score[1] <- NA
    df <- -diff(resdf)
    for (ii in seq_len(nmodels - 1)) {
      m1 <- if (df[ii] > 0) object[[ii]] else object[[ii + 1]]
      m2 <- if (df[ii] > 0) object[[ii + 1]] else object[[ii]]
      r <- m1$residuals
      w <- m1$weights
      method <- m2$method
      zedd <- eval(call(if (is.function(method)) "method" else method, 
            x = model.matrix(m2), y = r, weights = w))
      score[ii + 1] <- zedd$null.deviance - zedd$deviance
      if (df < 0)
        score[ii + 1] <- -score[ii + 1]
    }
  }  # doscore


  table <- if (.has.deviance)
             data.frame(resdf, resdev, c(NA, -diff(resdf)),
                        c(NA, -diff(resdev))) else
             data.frame(resdf, reslogLik, c(NA, -diff(resdf)),
                        c(NA, 2 * diff(reslogLik)))




  variables <- lapply(object, function(x) paste(deparse(formula(x)), 
      collapse = "\n"))
  dimnames(table) <- list(1L:nmodels,
                          c("Resid. Df", "Resid. Dev",  "Df", "Deviance"))



  if (doscore) 
    table <- cbind(table, Rao = score)
  title <- "Analysis of Deviance Table\n"
  topnote <- paste("Model ", format(1L:nmodels), ": ", variables, 
      sep = "", collapse = "\n")




  if (!is.null(test.null)) {
    bigmodel <- object[[order(resdf)[1L]]]
    dispersion <- 1
    df.dispersion <- if (dispersion == 1) Inf else min(resdf)


    table <- stat.anova(table = table, test = test.null,
                        scale = dispersion,
                        df.scale = df.dispersion)
  }  # !is.null(test.null)



  if (! .has.deviance)
  dimnames(table) <-
    list(1L:nmodels,
         c("Resid. Df", "LogLik", "Df", "2 * LogLik Diff.",
           if (!is.null(test.null) && doscore) "Rao" else NULL,
           if (is.null(test.null)) NULL else "Pr(>Chi)"))



  structure(table, heading = c(title, topnote), class = c("anova", 
            "data.frame"))
}  # anova.vglmlist











if (!isGeneric("anova"))
 setGeneric("anova", function(object, ...)
            standardGeneric("anova"),
            package = "VGAM")



setMethod("anova",
          "vglm", function(object, ...)
          anova.vglm(object, ...))






