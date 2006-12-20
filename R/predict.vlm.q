# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.





if(!exists("is.R")) is.R <- function()
    exists("version") && !is.null(version$language) && version$language=="R"


if(FALSE) vmodel.matrix.lm <- function(object, ...) {
    x <- slot(object, "x")
    if(!length(x)) {
        data <- model.framevlm(object, xlev = object@xlevels, ...)


        x <- if(is.R()) {
                 kill.con <- if(length(object@contrasts)) object@contrasts else
                             NULL 
                 vmodel.matrix.default(object, data = data,
                                      contrasts = kill.con)
             } else 
             model.matrix.default(terms(object), data = data, 
                                   contrasts = object@contrasts)
    } 
    x
}


if(!is.R()) {
if(FALSE)    vmodel.frame.lm <- function(formula, data = NULL, 
                                na.action = NULL, ...) {
        m <- formula@model
        if(!is.null(m))
            return(m)
        oc <- formula@call
        oc$method <- "model.frame"
        oc[[1.]] <- as.name("lm")
        if(length(data)) {
            oc$data <- substitute(data)
            eval(oc, sys.parent())
        }
        else eval(oc, list())
    }
}


if(is.R()) {
if(FALSE) vmodel.frame.lm = function (formula, data, na.action, ...) {
        if (is.null(formula@model)) {
            fcall <- formula@call
            fcall$method <- "model.frame"
            fcall[[1]] <- as.name("lm")
            env <- environment(fcall$formula)
            if (is.null(env))
                env <- parent.frame()
            eval(fcall, env)
        } else formula@model
    }
}


predict.vlm <- function(object, newdata=NULL, type=c("response","terms"),
                        se.fit = FALSE, scale = NULL,
                        terms.arg=NULL,
                        raw=FALSE,
                        dispersion = NULL, ...)
{

    if(mode(type) != "character" && mode(type) != "name")
        type <- as.character(substitute(type))
    type <- match.arg(type, c("response","terms"))[1]

    na.act = object@na.action
    object@na.action = list()

    if(raw && type != "terms")
        stop("sorry, raw=TRUE only works when type=\"terms\"")


    if(!length(newdata) && type=="response" && !se.fit &&
       length(object@fitted.values)) {
        if(length(na.act)) {
            return(napredict(na.act[[1]], object@fitted.values))
        } else {
            return(object@fitted.values)
        }
    }


    tt <- terms(object)  # 11/8/03; object@terms$terms
    if(!length(newdata)) {
        offset <- object@offset

        if(FALSE) {

            X <- model.matrixvlm(object, type="lm")



            if(is.R() && !length(object@x)) {
                attrassigndefault <- function(mmat, tt) {
                  if (!inherits(tt, "terms"))
                    stop("need terms object")
                  aa <- attr(mmat, "assign")
                  if (is.null(aa))
                    stop("argument is not really a model matrix")
                  ll <- attr(tt, "term.labels")
                  if (attr(tt, "intercept") > 0)
                    ll <- c("(Intercept)", ll)
                  aaa <- factor(aa, labels = ll)
                  split(order(aa), aaa)
                }
                attrassignlm <- function(object, ...) 
                    attrassigndefault(model.matrix(object, type="lm"),
                                      terms(object))

                attr(X, "assign") <- attrassignlm(X, tt)
            }
        } else {
            X <- model.matrix(object, type="lm")
        }
    } else {

        if(is.smart(object) && length(object@smart.prediction)) {
            setup.smart("read", smart.prediction=object@smart.prediction)

        }



        X = model.matrix(delete.response(tt), newdata,
            contrasts=if(length(object@contrasts)) object@contrasts else NULL,
            xlev = object@xlevels)


        if(is.R() && object@misc$intercept.only && nrow(X)!=nrow(newdata)) {
            as.save = attr(X, "assign")
            X = X[rep(1, nrow(newdata)),,drop=FALSE] # =matrix(1,nrow(newdata),1)
            dimnames(X) = list(dimnames(newdata)[[1]], "(Intercept)")
            attr(X, "assign") = as.save  # Restored 
        }

        offset <- if (!is.null(off.num<-attr(tt,"offset"))) {
            eval(attr(tt,"variables")[[off.num+1]], newdata)
        } else if (!is.null(object@offset))
            eval(object@call$offset, newdata)

        if(is.smart(object) && length(object@smart.prediction)) {
            wrapup.smart() 
        }

        if(is.R())
            attr(X, "assign") <- attrassigndefault(X, tt)

    }


    hasintercept <- attr(tt, "intercept")

    dx1 <- dimnames(X)[[1]]
    M <- object@misc$M
    Blist <- object@constraints
    ncolBlist <- unlist(lapply(Blist, ncol))
    if(hasintercept)
        ncolBlist <- ncolBlist[-1]

    xbar <- NULL
    if(type == "terms" && hasintercept) {
        if(is.null(newdata)) {
            xbar <- apply(X, 2, mean)
            X <- sweep(X, 2, xbar)
        } else {
            xbar <- apply(model.matrixvlm(object, type="lm"), 2, mean)
            xbar <- apply(X, 2, mean)
            X <- sweep(X, 2, xbar)
        }
        nac <- is.na(object@coefficients)
        if(any(nac)) {
            X <- X[, !nac, drop=FALSE]
            xbar <- xbar[!nac]
        }
    }

    if(!is.null(newdata) && !is.data.frame(newdata))
        newdata <- as.data.frame(newdata)

    nn <- if(!is.null(newdata)) nrow(newdata) else object@misc$n
    if(raw) {
        Blist <- canonical.Blist(Blist)
        object@constraints <- Blist
    }


    X <- lm2vlm.model.matrix(X, Blist=Blist, M=M, xij=object@control$xij)
    attr(X, "constant") <- xbar 

    coefs <- coef(object)
    vasgn <- attr(X, "vassign")

 
    if(type == "terms") {
        nv <- names(vasgn)
        if(hasintercept)
            nv <- nv[-(1:ncol(object@constraints[["(Intercept)"]]))]
        terms.arg <- if(is.null(terms.arg)) nv else terms.arg

        index <- if(is.R()) charmatch(terms.arg, nv) else 
                     match.arg(terms.arg, nv)
        if(all(index==0)) {
            warning("no match found; returning all terms")
            index <- 1:length(nv)
        }
        vasgn <- if(is.R()) vasgn[nv[index]] else {
            ans <- vector("list", length(index))
            for(loop in 1:length(index))
                ans[[loop]] <- vasgn[[index[loop]]]
            names(ans) <- index
            ans
        }
    }

    if(any(is.na(object@coefficients)))
        stop("can't handle NA's in object@coefficients")

    dname2 <- object@misc$predictors.names
    if(se.fit) {
        class(object) <- "vlm"   # coerce 
        fit.summary <- summaryvlm(object, dispersion=dispersion)
        sigma <- if(!is.null(fit.summary@sigma)) fit.summary@sigma else 
                 sqrt(deviance(object) / object@df.resid)  # was @rss 
        pred <- Build.terms.vlm(X, coefs,
                                cov=sigma^2 * fit.summary@cov.unscaled,
                                vasgn,
                                collapse=type!="terms", M=M,
                                dimname=list(dx1, dname2),
                                coefmat=coef(object, matrix=TRUE))
        pred$df <- object@df.residual
        pred$sigma <- sigma
    } else
        pred <- Build.terms.vlm(X, coefs, cov=NULL, assign = vasgn,
                                collapse=type!="terms", M=M,
                                dimname=list(dx1, dname2),
                                coefmat=coef(object, matrix=TRUE))
        constant <- attr(pred, "constant")


    if(type != "terms" && length(offset) && any(offset != 0)) {
        if(se.fit) {
            pred$fitted.values <- pred$fitted.values + offset
        } else {
            pred <- pred + offset
        }
    }

    if(type == "terms") {
        Blist <- subconstraints(object@misc$orig.assign, object@constraints)
        ncolBlist <- unlist(lapply(Blist, ncol))
        if(hasintercept)
            ncolBlist <- ncolBlist[-1]

        cs <- cumsum(c(1, ncolBlist))  # Like a pointer
        for(i in 1:(length(cs)-1))
            if(cs[i+1]-cs[i]>1) 
                for(k in (cs[i]+1):(cs[i+1]-1)) 
                    if(se.fit) {
                      pred$fitted.values[,cs[i]]=pred$fitted.values[,cs[i]] +
                                                 pred$fitted.values[,k]
                      pred$se.fit[,cs[i]] = pred$se.fit[,cs[i]] +
                                            pred$se.fit[,k]
                    } else {
                        pred[,cs[i]] <- pred[,cs[i]] + pred[,k]
                    }


        if(se.fit) {
            pred$fitted.values=pred$fitted.values[,cs[-length(cs)],drop=FALSE]
            pred$se.fit <- pred$se.fit[, cs[-length(cs)], drop=FALSE]
        } else {
            pred <- pred[, cs[-length(cs)], drop=FALSE]
        }
      
        pp <- if(se.fit) ncol(pred$fitted.values) else ncol(pred)
        if(se.fit) {
            dimnames(pred$fitted.values) <- dimnames(pred$se.fit) <- NULL
            dim(pred$fitted.values) <- dim(pred$se.fit) <- c(M, nn, pp)
            pred$fitted.values <- aperm(pred$fitted.values, c(2,1,3))
            pred$se.fit <- aperm(pred$se.fit, c(2,1,3))
            dim(pred$fitted.values) <- dim(pred$se.fit) <- c(nn, M*pp)
        } else {
            dimnames(pred) <- NULL   # Saves a warning
            dim(pred) <- c(M, nn, pp)
            pred <- aperm(pred, c(2,1,3))
            dim(pred) <- c(nn, M*pp)
        }

        if(raw) {
            kindex <- NULL
            for(i in 1:pp) 
                kindex <- c(kindex, (i-1)*M + (1:ncolBlist[i]))
            if(se.fit) {
                pred$fitted.values <- pred$fitted.values[,kindex,drop=FALSE]
                pred$se.fit <- pred$se.fit[,kindex,drop=FALSE]
            } else {
                pred <- pred[,kindex,drop=FALSE]
            }
        } 

        temp <- if(raw) ncolBlist else rep(M, length(ncolBlist))
        dd <- vlabel(names(ncolBlist), temp, M)
        if(se.fit) {
            dimnames(pred$fitted.values) <- 
            dimnames(pred$se.fit) <- list(if(length(newdata))
                                     dimnames(newdata)[[1]] else dx1, dd)
        } else {
            dimnames(pred) <- list(if(length(newdata))
                                   dimnames(newdata)[[1]] else dx1, dd)
        }

        if(!length(newdata) && length(na.act)) {
            if(se.fit) {
                pred$fitted.values = napredict(na.act[[1]], pred$fitted.values)
                pred$se.fit = napredict(na.act[[1]], pred$se.fit)
            } else {
                pred = napredict(na.act[[1]], pred)
            }
        }

        if(!raw) {
            cs <- cumsum(c(1, M + 0*ncolBlist))
        }
        fred <- vector("list", length(ncolBlist))
        for(i in 1:length(fred))
            fred[[i]] <- cs[i]:(cs[i+1]-1)
        names(fred) <- names(ncolBlist)
        if(se.fit) {
            attr(pred$fitted.values, "vterm.assign") <- fred
            attr(pred$se.fit, "vterm.assign") <- fred
        } else {
            attr(pred, "vterm.assign") <- fred
        }
    }

    if(!is.null(xbar))
        if(se.fit) {
            attr(pred$fitted.values, "constant") <- constant
        } else {
            attr(pred, "constant") <- constant
        }

    pred
}


subconstraints <- function(assign, constraints) {


    ans <- vector("list", length(assign))
    if(!length(assign) || !length(constraints))
        stop("assign and/or constraints is empty")
    for(i in 1:length(assign))
        ans[[i]] <- constraints[[assign[[i]][1]]]
    names(ans) <- names(assign)
    ans
}


is.linear.term <- function(ch) {
    lch <- length(ch)
    ans <- rep(FALSE, len=lch)
    for(i in 1:lch) {
        nc <- nchar(ch[i])
        x <- substring(ch[i], 1:nc, 1:nc)
        ans[i] <- all(x!="(" & x!="+" & x!="-" & x!="/" & x!="*" &
                      x!="^")
    }
    names(ans) <- ch
    ans 
}




canonical.Blist <- function(Blist) {
    for(i in 1:length(Blist)) {
        temp <- Blist[[i]] * 0
        r <- ncol(temp)
        temp[cbind(1:r,1:r)] <- 1
        Blist[[i]] <- temp
    }
    Blist
}




    setMethod("predict", "vlm", 
              function(object, ...)
              predict.vlm(object, ...))











