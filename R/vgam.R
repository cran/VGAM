# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.







vgam <- function(formula, 
                 family, 
                 data=list(), 
                 weights=NULL,
                 subset=NULL,
                 na.action=na.fail,
                 etastart=NULL, mustart=NULL, coefstart=NULL,
                 control=vgam.control(...),
                 offset=NULL, 
                 method="vgam.fit",
                 model=FALSE, x.arg=TRUE, y.arg=TRUE,
                 contrasts=NULL,
                 constraints=NULL,
                 extra=list(),
                 qr.arg=FALSE, smart=TRUE,
                 ...)
{
    dataname <- as.character(substitute(data))  # "list" if no data= 
    function.name <- "vgam"

    ocall <- match.call()

    if (smart)
        setup.smart("write")

    if (missing(data))
        data <- environment(formula)

    mtsave <- terms(formula, "s", data = data)

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

    xlev = .getXlevels(mt, mf)
    y <- model.response(mf, "any") # model.extract(mf, "response")
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
    attr(x, "assign") = attrassigndefault(x, mt)

    offset <- model.offset(mf)
    if (is.null(offset))
        offset <- 0 # yyy ???




    mf2 = mf
    if (!missing(subset)) {
        mf2$subset <- NULL 
        mf2 <- eval(mf2, parent.frame())   # mf2 is the full data frame. 
        spars2 =  lapply(mf2, attr, "spar") 
        dfs2 =  lapply(mf2, attr, "df") 
        sx2 =  lapply(mf2, attr, "s.xargument") 
        for (ii in 1:length(mf)) {
            if (length(sx2[[ii]])) {
                attr(mf[[ii]], "spar") = spars2[[ii]]
                attr(mf[[ii]], "dfs2") = dfs2[[ii]]
                attr(mf[[ii]], "s.xargument") = sx2[[ii]]
            }
        }
        rm(mf2) 
    }



    w <- model.weights(mf)
    if (!length(w))
        w <- rep(1, nrow(mf))
    else if (ncol(as.matrix(w))==1 && any(w < 0))
        stop("negative weights not allowed")




    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (!inherits(family, "vglmff")) {
        stop("'family=", family, "' is not a VGAM family function")
    }

    eval(vcontrol.expression)

    n <- dim(x)[1]

    if (FALSE && is.R()) {
        family@inverse <- eval(family@inverse)
        family@link <- eval(family@link)

        for (ii in names(.min.criterion.VGAM)) 
            if (length(family[[ii]])) family[[ii]] <- eval(family[[ii]])
    }

    if (length(slot(family, "first")))
        eval(slot(family, "first"))

    if (method != "vgam.fit")
        stop("method must be \"model.frame\" or \"vgam.fit\"")

    # --------------------------------------------------------------

    aa <- attributes(mtsave)
    smoothers <- aa$specials



    nonparametric <- length(smoothers$s) > 0
    if (nonparametric) {

        ff <- apply(aa$factors[smoothers[["s"]],,drop=FALSE], 2, any)
        smoothers[["s"]] <- if (any(ff))
            seq(along=ff)[aa$order==1 & ff] else NULL

        smooth.labels <- aa$term.labels[unlist(smoothers)]
    } else 
        function.name = "vglm"       # This is effectively so 



    fit <- vgam.fit(x=x, y=y, w=w, m=mf,
        etastart=etastart, mustart=mustart, coefstart=coefstart,
        offset=offset, family=family, control=control,
        criterion=control$criterion,
        constraints=constraints, extra=extra, qr.arg=qr.arg,
        Terms=mtsave,
        nonparametric=nonparametric, smooth.labels=smooth.labels,
        function.name=function.name, ...)


    if (is.Numeric(fit$nl.df) && any(fit$nl.df < 0)) {
        fit$nl.df[fit$nl.df < 0] = 0
    }

    # --------------------------------------------------------------

    if (!is.null(fit[["smooth.frame"]])) {
        fit <- fit[-1]       # Strip off smooth.frame
    } else {
    }

    fit$smomat <- NULL          # Not needed

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
    new("vgam", 
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
      "rss"          = fit$rss,
      "smart.prediction" = as.list(fit$smart.prediction),
      "terms"        = list(terms=fit$terms))

    if (!smart) answer@smart.prediction <- list(smart.arg=FALSE)

    if (qr.arg) {
        class(fit$qr) = "list"
        slot(answer, "qr") = fit$qr
    }
    if (length(attr(x, "contrasts")))
        slot(answer, "contrasts") = attr(x, "contrasts")
    if (length(fit$fitted.values))
        slot(answer, "fitted.values") = as.matrix(fit$fitted.values)
    slot(answer, "na.action") = if (length(aaa <- attr(mf, "na.action")))
        list(aaa) else list()
    if (length(offset))
        slot(answer, "offset") = as.matrix(offset)
    if (length(fit$weights))
        slot(answer, "weights") = as.matrix(fit$weights)
    if (x.arg)
        slot(answer, "x") = x # The 'small' design matrix
    if (length(xlev))
        slot(answer, "xlevels") = xlev
    if (y.arg)
        slot(answer, "y") = as.matrix(fit$y)
    answer@misc$formula = formula


    slot(answer, "control") = fit$control

    if (length(fit$extra)) {
        slot(answer, "extra") = fit$extra
    }
    slot(answer, "iter") = fit$iter
    slot(answer, "post") = fit$post

    fit$predictors = as.matrix(fit$predictors)  # Must be a matrix
    dimnames(fit$predictors) = list(dimnames(fit$predictors)[[1]],
                                    fit$misc$predictors.names)
    slot(answer, "predictors") = fit$predictors
    if (length(fit$prior.weights))
        slot(answer, "prior.weights") = as.matrix(fit$prior.weights)


    if (nonparametric) {
        slot(answer, "Bspline") = fit$Bspline
        slot(answer, "nl.chisq") = fit$nl.chisq
        if (is.Numeric(fit$nl.df))
            slot(answer, "nl.df") = fit$nl.df
        slot(answer, "spar") = fit$spar
        slot(answer, "s.xargument") = fit$s.xargument
        if (length(fit$varmat)) {
            slot(answer, "var") = fit$varmat
        }






    }
    if (length(fit$effects))
        slot(answer, "effects") = fit$effects


    answer
}
attr(vgam, "smart") <- TRUE 




care.exp <- function(x, thresh = -log(.Machine$double.eps)) {
    x[x >   thresh]  <-  thresh
    x[x < (-thresh)] <- -thresh
    exp(x)
}





