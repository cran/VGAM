# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.




vlm <- function(formula,
                data=list(), 
                weights=NULL, subset, na.action,
                prior.weights=NULL, 
                control=vlm.control(...), 
                method="qr",
                model=FALSE, x.arg=FALSE, y.arg=TRUE, qr.arg=TRUE,
                contrasts=NULL, 
                constraints=NULL,
                extra=NULL, offset=NULL,  
                smart=TRUE, ...)
{
    dataname <- as.character(substitute(data))  # "list" if no data=
    function.name <- "vlm"

    ocall <- match.call()

    if(smart)
        setup.smart("write")

    mt <- terms(formula, data = data)  # attr(m, "terms")
    if (missing(data))
        data <- sys.frame(sys.parent())

    mf <- match.call(expand=FALSE)
    mf$method <- mf$model <- mf$x.arg <- mf$y.arg <- mf$control <- 
        mf$contrasts <- mf$constraints <- mf$extra <- 
        mf$qr.arg <- mf$smart <- mf$... <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame()) 
    if(method == "model.frame")
        return(mf)
    if(method != "qr")
        stop("only method=\"qr\" is implemented")

    na.act <- attr(mf, "na.action")

    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0)
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }

    y <- model.response(mf, "numeric") # model.extract(mf, "response")
    x <- model.matrix(mt, mf, contrasts)
    attr(x, "assign") <- attrassigndefault(x, mt) # So as to make it like Splus
    offset <- model.offset(mf)
    if(is.null(offset))
        offset <- 0 # yyy ???
    if(length(offset) && any(offset!=0))
        stop("offsets are redundant for (vector) linear models")
    wz <- model.weights(mf)

    y = as.matrix(y)
    M <- ncol(as.matrix(y))
    n <- nrow(x)
    dy <- dimnames(y)
    dy1 <- if(length(dy[[1]])) dy[[1]] else dimnames(mf)[[1]]
    dy2 <- if(length(dy[[2]])) dy[[2]] else paste("Y", 1:M, sep="")
    dimnames(y) <- list(dy1, dy2)
    predictors.names = dy2

    if(!length(prior.weights)) {
        prior.weights = rep(1, len=n)
        names(prior.weights) = dy1
    }
    if(any(prior.weights <= 0))
        stop("only positive weights allowed")
    if(!length(wz)) {
        wz <- matrix(prior.weights, n, M)
        identity.wts <- TRUE
    } else {
        identity.wts <- FALSE
        temp = ncol(as.matrix(wz))
        if(temp < M || temp > M*(M+1)/2)
            stop(paste("input w must have at between", M, "and",
                       M*(M+1)/2, "columns"))
        wz <- prior.weights * wz
    }

    control = control
    Blist <- process.constraints(constraints, x, M)
    intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"

    fit <- vlm.wfit(x=x, z=y, Blist=Blist, wz=wz, U=NULL,
                    matrix.out=FALSE, XBIG=FALSE, rss=TRUE, qr=qr.arg,
                    x.ret=TRUE, offset = offset)

    p.big <- fit$rank
    fit$R <- fit$qr$qr[1:p.big, 1:p.big, drop=FALSE]
    fit$R[lower.tri(fit$R)] <- 0




    fit$constraints <- Blist

    dn.big <- labels(fit$xbig)
    xn.big <- dn.big[[2]]
    dn <- labels(x)
    xn <- dn[[2]]
    dxbig <- as.integer(dim(fit$xbig))
    n.big <- dxbig[[1]]
    p.big <- dxbig[[2]]

    misc <- list(
        colnames.x = xn,
        colnames.xbig = xn.big,
        function.name = function.name,
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = nrow(x),
        n.big = n.big,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        p.big = p.big,
        ynames = dimnames(y)[[2]])
    
    fit$misc <- misc

    fit$misc$dataname <- dataname
    

    
    if(smart) {
        fit$smart.prediction <- get.smart.prediction()
        wrapup.smart()
    }

    answer <-
    new("vlm", 
      "assign"       = attr(x, "assign"),
      "call"         = ocall,
      "coefficients" = fit$coefficients,
      "constraints"  = fit$constraints,
      "control"      = control, 
      "criterion"    = list(deviance=fit$rss),
      "dispersion"   = 1,
      "df.residual"  = fit$df.residual,
      "df.total"     = n*M,
      "effects"      = fit$effects,
      "fitted.values"= as.matrix(fit$fitted.values),
      "misc"         = fit$misc,
      "model"        = if(model) mf else data.frame(),
      "R"            = fit$R,
      "rank"         = fit$rank,
      "residuals"    = as.matrix(fit$residuals),
      "rss"          = fit$rss,
      "smart.prediction" = as.list(fit$smart.prediction),
      "terms"        = list(terms=mt))

    if(!smart) answer@smart.prediction <- list(smart.arg=FALSE)

    slot(answer, "prior.weights") = prior.weights

    if(length(attr(x, "contrasts")))
        slot(answer, "contrasts") = attr(x, "contrasts")
    slot(answer, "na.action") = if(length(na.act)) list(na.act) else list()
    if(length(offset))
        slot(answer, "offset") = as.matrix(offset)
    if(qr.arg) {
        class(fit$qr) = "list"
        slot(answer, "qr") = fit$qr
    }
    if(x.arg)
        slot(answer, "x") = x # The 'small' design matrix
    if(control$save.weight)
        slot(answer, "weights") = wz
    if(length(xlev))
        slot(answer, "xlevels") = xlev
    if(y.arg)
        slot(answer, "y") = as.matrix(y)

    answer
}
attr(vlm, "smart") <- TRUE    




