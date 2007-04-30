# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.



predict.vgam <- function(object, newdata=NULL,
                         type=c("link", "response", "terms"),
                         se.fit=FALSE, deriv.arg=0, terms.arg=NULL,
                         raw=FALSE,
                         all=TRUE, offset=0, 
                         untransform = FALSE,
                         dispersion=NULL, ...)
{
    if(missing(newdata)) {
        newdata <- NULL
    } else {
        newdata <- as.data.frame(newdata)
    }
    no.newdata = length(newdata)==0

    na.act = object@na.action
    object@na.action = list()

    if(mode(type) != "character" && mode(type) != "name")
        type <- as.character(substitute(type))
    type <- match.arg(type, c("link", "response", "terms"))[1]


    if(untransform && (type!="link" || se.fit || deriv.arg != 0 || offset != 0))
        stop(paste("argument \"untransform\"=TRUE only if type=\"link\",",
                   "se.fit=FALSE, deriv=0"))

    if(raw && type!="terms")
        stop("raw=TRUE only works when type=\"terms\"")

    if(!is.numeric(deriv.arg) || deriv.arg<0 ||
       deriv.arg!=round(deriv.arg) || length(deriv.arg)>1)
        stop("bad input for the deriv argument")

    if(deriv.arg>0 && type!="terms")
        stop("deriv>0 can only be specified if type=\"terms\"")

    if(deriv.arg != 0 && !(type!="response" && !se.fit))
        stop(paste("deriv= only works with type!=\"response\"and se.fit=FALSE"))

    if(se.fit && length(newdata))
        stop("can't specify se.fit=TRUE when there is newdata")


    tt <- terms(object) # 11/8/03; object@terms$terms

    ttf <- attr(tt, "factors")
    tto <- attr(tt, "order")
    intercept <- attr(tt, "intercept")
    if(!intercept)
        stop("an intercept is assumed")

    M <- object@misc$M
    Blist <- object@constraints
    ncolBlist <- unlist(lapply(Blist, ncol))
    if(intercept)
        ncolBlist <- ncolBlist[-1]
    if(raw) {
        Blist <- canonical.Blist(Blist)
        object@constraints <- Blist
    }

    if(!length(newdata)) {
        if(type=="link") {
            if(se.fit) {
                stop("cannot handle this option (se.fit=TRUE) currently")  # zz
            } else {
                if(length(na.act)) {
                    answer = napredict(na.act[[1]], object@predictors)
                } else {
                    answer = object@predictors
                }
                if(untransform) return(untransformVGAM(object, answer)) else
                    return(answer)
            }
        } else 
        if(type=="response") {
            if(se.fit) {
                stop("cannot handle this option (se.fit=TRUE) currently")  # zz
            } else {
                if(length(na.act)) {
                    return(napredict(na.act[[1]], object@fitted.values))
                } else {
                    return(object@fitted.values)
                }
            }
        }

        predictor <- predict.vlm(object,
                         type="terms",
                         se.fit=se.fit,
                         terms.arg=terms.arg,
                         raw=raw,
                         all=all, offset=offset, 
                         dispersion=dispersion, ...) # deriv.arg=deriv.arg,

        newdata <- model.matrixvlm(object, type="lm")


    } else {

        temp.type <- if(type=="link") "response" else type 


        predictor <- predict.vlm(object, newdata,
                            type=temp.type,
                            se.fit=se.fit,
                            terms.arg=terms.arg,
                            raw=raw,
                            all=all, offset=offset, 
                            dispersion=dispersion, ...) # deriv.arg=deriv.arg,
    }


    if(deriv.arg>0)
        if(se.fit) {
            predictor$fitted.values <- predictor$fitted.values * 0
            predictor$se.fit <- predictor$se.fit * NA
        } else {
            predictor <- predictor * 0
        }


    if(length(s.xargument <- object@s.xargument)) {




        dnames2 <- dimnames(newdata)[[2]]
        index1 <- match(s.xargument, dnames2, nomatch=FALSE)
        index2 <- match(names(s.xargument), dnames2, nomatch=FALSE)
        index <- index1 | index2
        if(!length(index) || any(!index))
            stop("required variables not found in newdata")




        if(is.null(tmp6 <- attr(if(se.fit) predictor$fitted.values else 
                                predictor, "vterm.assign"))) {

            Blist <- subconstraints(object@misc$orig.assign,
                                    object@constraints)
            ncolBlist <- unlist(lapply(Blist, ncol))
            if(intercept)
                ncolBlist <- ncolBlist[-1]
    
            cs <- if(raw) cumsum(c(1, ncolBlist)) else
                          cumsum(c(1, M + 0*ncolBlist))
            tmp6 <- vector("list", length(ncolBlist))
            for(i in 1:length(tmp6))
                tmp6[[i]] <- cs[i]:(cs[i+1]-1)
            names(tmp6) <- names(ncolBlist)
        }

        n.s.xargument <- names(s.xargument)   # e.g., c("s(x)", "s(x2)")
        for(i in n.s.xargument) {

            fred <- s.xargument[i]
            if(!any(dimnames(newdata)[[2]] == fred))
                fred <- i

            xx <- newdata[,fred] # [,s.xargument[i]]   # [,nindex[i]]   
            ox <- order(xx)

            raw.mat <- predictvsmooth.spline.fit(
                                 object@Bspline[[i]],
                                 x=xx,
                                 deriv=deriv.arg)$y


            eta.mat <- if(raw) raw.mat else
                     raw.mat %*% t(Blist[[i]])

            if(type=="terms") {
                ii <- tmp6[[i]]
                if(se.fit) {
                    predictor$fitted.values[,ii] = 
                    predictor$fitted.values[,ii] + eta.mat

                        TS <- predictor$sigma^2

                        temp.var <- if(raw) {
                                        iii <- object@misc$varassign
                                        iii <- iii[[i]]
                                        object@var[,iii,drop=FALSE]
                                    } else
                                        stop("can't handle se's with raw=FALSE")

                        predictor$se.fit[,ii] <- (predictor$se.fit[,ii]^2 +
                           TS * temp.var)^0.5
                } else {
                    predictor[,ii] <- predictor[,ii] + eta.mat
                }
            } else {
                if(se.fit) {
                    predictor$fitted.values <- predictor$fitted.values + eta.mat 

                    TS <- 1  # out$residual.scale^2
                    TS <- predictor$sigma^2

                    TT <- ncol(object@var)
                    predictor$se.fit <- sqrt(predictor$se.fit^2 + TS *
                                             object@var %*% rep(1, TT))


                } else {
                    predictor <- predictor + eta.mat 
                }
            }
        }
    }

    if(type=="link") {
        if(no.newdata && length(na.act)) {
            return(napredict(na.act[[1]], predictor))
        } else {
            return(predictor)
        }
    } else
    if(type=="response") {
        fv <- object@family@inverse(if(se.fit) predictor$fitted.values else
                                    predictor, object@extra)
        if(is.matrix(fv) && is.matrix(object@fitted.values))
            dimnames(fv) <- list(dimnames(fv)[[1]],
                                 dimnames(object@fitted.values)[[2]])
        if(is.matrix(fv) && ncol(fv)==1)
            fv <- c(fv)
        if(no.newdata && length(na.act)) {
            if(se.fit) {
                fv = napredict(na.act[[1]], fv)
            } else {
                fv = napredict(na.act[[1]], fv)
            }
        }
        if(se.fit) {
            return(list(fit=fv, se.fit=fv*NA))
        } else {
            return(fv)
        }
    } else {
        if(deriv.arg >= 1)
            if(se.fit) {
                attr(predictor$fitted.values, "constant") <- NULL
            } else {
                attr(predictor, "constant") <- NULL
            }

        if(deriv.arg >= 1) {
            v = attr(if(se.fit) predictor$fitted.values else 
                predictor, "vterm.assign")
            is.lin <- is.linear.term(names(v))
            coefmat <- coef(object, matrix=TRUE)
            ord <- 0
            for(i in names(v)) {
                ord <- ord + 1
                index <- v[[i]]
                lindex <- length(index)
                if(is.lin[i]) {
                    if(tto[ord]>1 || (length(ttf) && ttf[i,i])) {
                        if(se.fit) {
                            predictor$fitted.values[,index] = 
                                if(tto[ord]>1) NA else NA
                        } else {
                            predictor[,index] <- if(tto[ord]>1) NA else NA
                        }
                    } else {
                        ans <- coefmat[i, 1:lindex]
                        if(se.fit) {
                            predictor$fitted.values[,index] = if(deriv.arg==1)
                                matrix(ans, ncol=lindex, byrow=TRUE) else 0
                        } else {
                            predictor[,index] <- if(deriv.arg==1)
                                matrix(ans, ncol=lindex, byrow=TRUE) else 0
                        }
                    }
                } else
                if(length(s.xargument) && any(n.s.xargument == i)) {
                    ans <- coefmat[i, 1:lindex]
                    if(se.fit) {
                        predictor$fitted.values[,index] =
                        predictor$fitted.values[,index] + 
                             (if(deriv.arg==1)
                              matrix(ans, nrow=nrow(predictor$fitted.values),
                               ncol=lindex, byrow=TRUE) else 0)
                    } else {
                        predictor[,index] <- predictor[,index] +
                             (if(deriv.arg==1)
                              matrix(ans, nrow=nrow(predictor), 
                               ncol=lindex, byrow=TRUE) else 0)
                    }
                } else {
                    cat("Derivatives of term ", i, "are unknown\n")
                    if(se.fit) {
                        predictor$fitted.values[,index] <- NA
                    } else {
                        predictor[,index] <- NA
                    }
                }
            }
        }

        if(no.newdata && length(na.act)) {
            if(se.fit) {
                predictor$fitted.values = napredict(na.act[[1]],
                                                    predictor$fitted.values)
                predictor$se.fit = napredict(na.act[[1]], predictor$se.fit)
            } else {
                predictor = napredict(na.act[[1]], predictor)
            }
        }

        if(se.fit) {
            attr(predictor$fitted.values, "derivative") <- deriv.arg
        } else {
            attr(predictor, "derivative") <- deriv.arg
        }

        return(predictor)
    }
}


    setMethod("predict", "vgam",
              function(object, ...)
              predict.vgam(object, ...))



varassign <- function(constraints, n.s.xargument) { 

    if(!length(n.s.xargument))
        stop("length(n.s.xargument) must be > 0")

    ans <- vector("list", length(n.s.xargument))

    ncolBlist <- unlist(lapply(constraints, ncol))

    names(ans) <- n.s.xargument
    ptr <- 1
    for(i in n.s.xargument) {
        temp <- ncolBlist[[i]]
        ans[[i]] <- ptr:(ptr+temp-1)
        ptr <- ptr + temp
    }
    ans 
}




