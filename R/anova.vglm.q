# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.















car.relatives <- function(term, names, factors) {
  # This function is car:::relatives
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2])))
  }
  if(length(names) == 1) return(NULL)
  which.term <- which(term == names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
          function(term2) is.relative(term, term2))]
}










fitmodel.VGAM.expression <- expression({


    

    Usex.lm <- vecTF  # Should be called vecTF.lm really
    terms.lm <- findterms(Usex.lm, asgn)
    Usex.vlm <- rep(Usex.lm, times = ncolHlist.lm)
    terms.vlm <- findterms(Usex.vlm, vasgn)
    assign1.lm <- subsetassign(asgn, union(oTerms.wint, terms.lm))
    assign1.vlm <- subsetassign(vasgn, union(voTerms.wint, terms.vlm))
    Col.Usex.lm <- seq_len(length(Usex.lm))[Usex.lm]
    Col.Usex.lm <- unique(sort(c(Col.Usex.lm, ousex.lm)))  # dddd
    X.lm <- big.x.lm[, Col.Usex.lm, drop = FALSE]
    attr(X.lm, "assign") <- assign1.lm
    Col.Usex.vlm <- seq_len(length(Usex.vlm))[Usex.vlm]
    Col.Usex.vlm <- unique(sort(c(Col.Usex.vlm, ousex.vlm)))  # dddd
    X.vlm <- big.x.vlm[, Col.Usex.vlm, drop = FALSE]
    attr(X.vlm, "vassign") <- assign1.vlm




   if (is.logical(object@control$trace))
      object@control$trace <- FALSE  # Supress 'trace'; keep silent
    prewarn <- options("warn")
    options(warn = -1)  # Supress warnings
      fit1 <- vglm.fit(x = X.lm, y = Y, w = Wts,
                       X.vlm.arg = X.vlm,
                       Xm2 = Xm2, Terms = mt,
        constraints = big.clist.term[unique(c(oTerms.wint,
                        terms.lm))],  # dddd; Unsorted okay
                       extra = object@extra,
                       etastart = LPmat,
                       offset = OOO, family = Fam,
                       control = object@control)
    options(warn = prewarn[["warn"]])  # Restore warnings
})  # fitmodel.VGAM












anova.vglm <-
  function(object, ...,
           type = c("II", "I", "III", 2, 1, 3),
           test = c("LRT", "none"),  # yettodo: "Rao"
           trydev = TRUE,  # Use where possible?
           silent = TRUE) {
  type <- as.character(type)
  type <- match.arg(type, c("II", "I", "III", "2", "1", "3"))
  type[type == "1"] <- "I"
  type[type == "2"] <- "II"
  type[type == "3"] <- "III"
  if ((int2 <- has.intercept(object)) &&
     length(constraints(object)) == 1 && 
     names(constraints(object)) == "(Intercept)" && type == "II") {
    type <- "III"
    warning("the model contains only an intercept; ",
            "Type III test substituted")
  }

  if (length(list(...)) && type != "I")
    stop("argument 'type' must 'I' or 1 for multiple fits")


  dispersion <- 1
  if ((int <- attr(terms(object), "intercept")) != int2)
    stop("cannot determine whether there is an intercept or not")

  if (mode(test) != "character" && mode(test) != "name")
    test <- as.character(substitute(test))
  test <- match.arg(test, c("LRT", "none"))[1]  # , "Rao"
  test.null <- if (test == "none") NULL else test


  if (!int2)
    stop("argument 'object' must have an intercept term")
  object@control$trace <- FALSE
  has.deviance <- !is.null(dev.object <- deviance(object)) && trydev
  if (silent) {
    warn.save <- unlist(options("warn"))
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
                          test = test, type = type,
                          .has.deviance = has.deviance,
                          .trydev = trydev))


  varlist <- attr(terms(object), "variables")
  x.lm <- model.matrix(object, type = "lm")
  x.vlm <- model.matrix(object, type = "vlm")
  p.lm <- ncol(x.lm)  # Needed for type == "II"
  p.vlm <- ncol(x.vlm)  # Needed for type == "III"
  orig.assign.lm <- varseq <- attr(x.lm, "orig.assign.lm")
  if (!length(varseq))
    stop("could not obtain attribute 'orig.assign.lm' from ",
         "the model matrix; try vglm(..., x = TRUE) and rerun")
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
  Xm2 <- model.matrix(object, type = "lm2")  # Could be 0 x 0
  if (!length(Xm2))
     Xm2 <- NULL  # Make sure. This is safer
  LPmat <- predict(object)
  Fam <- object@family


  big.clist.lm <- constraints(object, type = "lm")
  big.clist.term <- constraints(object, type = "term")
  ncolHlist.lm <- unlist(lapply(big.clist.lm, ncol))
  big.x.lm <- x.lm
  big.x.vlm <- x.vlm
  asgn <- attr(big.x.lm, "assign")  # \pkg{VGAM}
  vasgn <- attr(big.x.vlm, "vassign")  # \pkg{VGAM}


  if (type == "I") {
    if (!int)
      stop("an intercept is needed to fit a null model")

    vecTF <- varseq == 0


    fit1 <- NULL  # To avoid an warning on CRAN
    vecTF <- vecTF
    oTerms.wint <- voTerms.wint <- ousex.lm <- ousex.vlm <- NULL
    eval(fitmodel.VGAM.expression)
    fit0 <- fit1


    object.df.null <- fit0$df.residual
    object.null.deviance <- fit0$crit.list$deviance
    object.null.logLik   <- fit0$crit.list$loglikelihood
  }  # TRUE && is.element(type, c("I", "II", "III"))




    tlab <- attr(terms(object), "term.labels")  # Omits any intercept



  upp.bnd <- switch(type, "I" = nvars - 1L, "II" = , "III" = nvars)

  if (upp.bnd > 0) {  # nvars > 1 for type = "I"
    if (type == "II") {
      which.nms <- function(name)
        which(orig.assign.lm == which(Names == name))
      Fac <- attr(terms(object), "factors")
      Names <- term.names(object)
      if (Names[1] == "(Intercept)")
        Names <- Names[-1]
      if (!all(tlab == Names))
        stop("'tlab' not identical to 'Names'")
    }  # type == "II"



    for (ii in seq_len(upp.bnd)) {
      if (type == "II") {
        index3 <- car.relatives(term    = Names[ii],
                                names   = Names,
                                factors = Fac)
        rels <- Names[index3]
        exclude.1 <- as.vector(unlist(sapply(c(Names[ii], rels),
                                             which.nms)))
        exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))


        vecTF1 <- vecTF2 <- rep(TRUE, p.lm)  # For type == "II"
        vecTF1[exclude.1] <- FALSE
        if (length(rels) > 0)
          vecTF2[exclude.2] <- FALSE


        vecTF <- vecTF2
    oTerms.wint <- voTerms.wint <- ousex.lm <- ousex.vlm <- NULL
        eval(fitmodel.VGAM.expression)
        fit2 <- fit1

      }  # type == "II"




      vecTF <- switch(type,  # Wrt x.lm columns
                      "I"   = (varseq <= ii),
                      "II"  = vecTF1,  # !vecTF1,
                      "III" = (varseq != ii))

      vecTF <- vecTF
      oTerms.wint <- voTerms.wint <- ousex.lm <- ousex.vlm <- NULL
      eval(fitmodel.VGAM.expression)


      reslogLik <- c(reslogLik, fit1$crit.list$loglik)
      resdev <- c(resdev, fit1$crit.list$deviance)  # May be NULL
      resdf <- c(resdf, fit1$df.residual)

      if (type == "II") {
        reslogLik2 <- c(reslogLik2, fit2$crit.list$loglik)
        resdev2 <- c(resdev2, fit2$crit.list$deviance)  # May be NULL
        resdf2 <- c(resdf2, fit2$df.residual)
      }  # "II"
    }  # for ii ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  }  # if (upp.bnd > 0)  # nvars > 1 for type = "I"



  if (type == "I") {
    resdf <- c(object.df.null, resdf, df.residual(object))
    resdev <- c(object.null.deviance, resdev, deviance(object))
    reslogLik <- c(object.null.logLik, reslogLik, logLik(object))
  }  # type == "I"




 
  table <-
    if (has.deviance) {
      if (type == "I")
        data.frame(col1 = c(NA, -diff(resdf)),
                   col2 = c(NA, pmax(0, -diff(resdev))), 
                   resdf, resdev) else
      if (type == "II")
        data.frame(col1 = resdf - resdf2,
                   col2 = pmax(0, resdev - resdev2),
                   resdf, resdev) else
        data.frame(col1 = resdf - df.residual(object),
                   col2 = pmax(0, resdev - deviance(object)), 
                   resdf, resdev)
      } else {
      if (type == "I")
        data.frame(col1 = c(NA, -diff(resdf)),
                   col2 = c(NA, pmax(0, 2 * diff(reslogLik))), 
                   resdf, reslogLik) else
      if (type == "II")
        data.frame(col1 = resdf - resdf2,
                   col2 = pmax(0, 2 * (reslogLik2 - reslogLik)),
                   resdf, reslogLik) else
        data.frame(col1 = resdf - df.residual(object),
                   col2 = pmax(0, 2 * (logLik(object) - reslogLik)),
                   resdf, reslogLik)
      }
  if (length(tlab) == 0L)
    table <- table[1, , drop = FALSE]
  dn1 <- c(if (is.element(type, c("II", "III"))) NULL else "NULL",
           tlab)
  dn2.before <- c("Df", "Deviance",         "Resid. Df", "Resid. Dev")
  dn2.after  <- c("Df", "2 * LogLik Diff.", "Resid. Df", "LogLik")
  dimnames(table) <- list(dn1, dn2.before)  # For stat.anova()

  lfuns <- linkfun(object)
  suptitle <- if (type == "I") paste("Type I tests: terms added ",
              "sequentially from\nfirst to last", sep = "") else
     if (type == "II")
       "Type II tests" else
       "Type III tests: each term added last"
  title <- paste0("Analysis of Deviance Table (", suptitle, ")",
     "\n\nModel: ",
     paste(paste("'", Fam@vfamily, "'", sep = ""), collapse = ", "),
     if (length(lfuns) > 1) "\n\nLinks: " else "\n\nLink: ",
     if (length(unique(lfuns)) == 1)
       paste0("'", lfuns[1], "'") else
       paste(paste("'", lfuns, "'", sep = ""), collapse = ", "),
     "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n")
  df.dispersion <- Inf



  if (!is.null(test.null)) {
    table <- stat.anova(table = table, test = test.null,
                        scale = dispersion, 
                        df.scale = df.dispersion)
  }  # (!is.null(test.null))
  if (!has.deviance)
    dimnames(table) <-
      list(dn1, c(dn2.after,
                  if (is.null(test.null)) NULL else "Pr(>Chi)"))
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
  test <- match.arg(test, c("LRT", "none"))[1]
  test.null <- if (test == "none") NULL else test


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



  ns1 <- as.numeric(lapply(object, function(x) nobs(x, type =  "lm")))
  ns2 <- as.numeric(lapply(object, function(x) nobs(x, type = "vlm")))
  if (any(ns1 != ns1[1L]) || any(ns2 != ns2[1L]))
    stop("models were not all fitted to the same size of dataset")


  nmodels <- length(object)
  if (nmodels == 1) 
    return(anova.vglm(object[[1L]],
                      test = test.null))



  resdf <- as.numeric(lapply(object, function(x) df.residual(x)))
  reslogLik <- as.numeric(lapply(object, function(x) logLik(x)))
  if (.has.deviance && .trydev)
    resdev <- as.numeric(lapply(object, function(x) deviance(x)))


  table <- if (.has.deviance)
             data.frame(resdf, resdev, c(NA, -diff(resdf)),
                        c(NA, -diff(resdev))) else
             data.frame(resdf, reslogLik, c(NA, -diff(resdf)),
                        c(NA, 2 * diff(reslogLik)))




  variables <- lapply(object, function(x) paste(deparse(formula(x)), 
      collapse = "\n"))
  dimnames(table) <-
    list(1L:nmodels,
         c("Resid. Df", "Resid. Dev",  "Df", "Deviance"))


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
           if (is.null(test.null)) NULL else "Pr(>Chi)"))

  structure(table, heading = c(title, topnote), class = c("anova", 
            "data.frame"))
}  # anova.vglmlist













setMethod("anova",
          "vglm", function(object, ...)
          anova.vglm(object, ...))






