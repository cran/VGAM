# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.



untransformVGAM = function(object, pred) {
    M = object@misc$M
    Links = object@misc$link
    if(length(Links) != M && length(Links) != 1)
       stop("cannot obtain the link functions to untransform the object")
    upred = pred
    earg = object@misc$earg
    for(ii in 1:M) {
        TTheta = pred[,ii]  # Transformed theta
        newcall = paste(Links[ii], "(theta=TTheta, earg=earg, inverse=TRUE)", sep="")
        newcall = parse(text=newcall)[[1]]
        Theta = eval(newcall) # Theta, the untransformed parameter
        upred[,ii] = Theta
    }
    dmn2 = if(length(names(object@misc$link))) names(object@misc$link) else {
        if(length(object@misc$parameters)) object@misc$parameters else NULL}
    dimnames(upred) = list(dimnames(upred)[[1]], dmn2)
    upred
}




predict.vglm <- function(object,
                         newdata=NULL,
                         type=c("link", "response", "terms"),
                         se.fit=FALSE,
                         deriv=0,
                         dispersion=NULL,
                         untransform=FALSE,
                         extra=object@extra, ...)
{

    na.act = object@na.action
    object@na.action = list()

    mextra <- missing(extra)
    if(mextra) {
    }

    if(deriv!=0)
        stop("deriv must be 0 for predict.vglm")

    if(mode(type) != "character" && mode(type) != "name")
        type <- as.character(substitute(type))
    type <- match.arg(type, c("link", "response", "terms"))[1]

    if(untransform && (type!="link" || se.fit || deriv != 0))
        stop(paste("argument \"untransform\"=TRUE only if type=\"link\",",
                   "se.fit=FALSE, deriv=0"))


    pred = 
    if(se.fit)
        switch(type,
               response = {
                  warning(paste("type=\"response\" and se.fit=TRUE not valid",
                                "together; setting se.fit=FALSE"))
                  se.fit <- FALSE
                    predictor <- predict.vlm(object, newdata=newdata,
                                             type=type, se.fit=se.fit,
                                             deriv=deriv, 
                                             dispersion=dispersion, ...) 
                  fv <- object@family@inverse(predictor, extra)
                  dimnames(fv) <- list(dimnames(fv)[[1]],
                                       dimnames(object@fitted.values)[[2]])
                  fv
               },
               link = {
                       type <- "response"
                       predict.vlm(object, newdata=newdata,
                                   type=type, se.fit=se.fit,
                                   deriv=deriv, dispersion=dispersion, ...) 
               },
                terms={
                    predict.vlm(object, newdata=newdata,
                                type=type, se.fit=se.fit,
                                deriv=deriv, dispersion=dispersion, ...) 
                }) else {
        if(is.null(newdata)) {
            switch(type, 
                   link=object@predictors, 
                   response=object@fitted.values,
                   terms={
                       predict.vlm(object, newdata=newdata,
                                   type=type, se.fit=se.fit,
                                   deriv=deriv, dispersion=dispersion, ...) 
                   })
        } else {
            if(!(length(object@offset)==1 && object@offset==0))
                warning("zero offset used") 
            switch(type, 
                response={
                    predictor <- predict.vlm(object, newdata=newdata,
                                             type=type, se.fit=se.fit,
                                             deriv=deriv, 
                                             dispersion=dispersion, ...) 
                    M <- object@misc$M



                    fv <- object@family@inverse(predictor, extra)
                    if(M > 1 && is.matrix(fv)) {
                        dimnames(fv) <- list(dimnames(fv)[[1]],
                                             dimnames(object@fitted.values)[[2]])
                    } else {
                    }
                    fv
                }, 
                link = {
                    type <- "response"
                    predict.vlm(object, newdata=newdata,
                                type=type, se.fit=se.fit,
                                deriv=deriv, dispersion=dispersion, ...) 
                },
                terms={
                    predict.vlm(object, newdata=newdata,
                                type=type, se.fit=se.fit,
                                deriv=deriv, dispersion=dispersion, ...) 
                }
            )
        }
    }

    if(!length(newdata) && length(na.act)) {
        if(se.fit) {
            pred$fitted.values = napredict(na.act[[1]], pred$fitted.values)
            pred$se.fit = napredict(na.act[[1]], pred$se.fit)
        } else {
            pred = napredict(na.act[[1]], pred)
        }
    }
    
    if(untransform) untransformVGAM(object, pred) else pred
}


setMethod("predict", "vglm", function(object, ...) 
    predict.vglm(object, ...))




predict.rrvglm = function(object, 
                         newdata=NULL, 
                         type=c("link", "response", "terms"),
                         se.fit=FALSE, 
                         deriv=0,
                         dispersion=NULL, 
                         extra=object@extra, ...) 
{

    if(se.fit) {
        stop("11/8/03; predict.rrvglm(..., se.fit=TRUE) not complete yet") 
        pred = 
        switch(type,
               response = {
                  warning(paste("type=\"response\" and se.fit=TRUE not valid",
                                "together; setting se.fit=FALSE"))
                  se.fit <- FALSE
                    predictor <- predict.vlm(object, newdata=newdata,
                                             type=type, se.fit=se.fit,
                                             deriv=deriv, 
                                             dispersion=dispersion, ...) 
                  fv <- object@family@inverse(predictor, extra)
                  dimnames(fv) <- list(dimnames(fv)[[1]],
                                       dimnames(object@fitted.values)[[2]])
                  fv
               },
               link = {
                       type <- "response"
                       predict.vlm(object, newdata=newdata,
                                   type=type, se.fit=se.fit,
                                   deriv=deriv, dispersion=dispersion, ...) 
               },
                terms={
                    predict.vlm(object, newdata=newdata,
                                type=type, se.fit=se.fit,
                                deriv=deriv, dispersion=dispersion, ...) 
                }
              )
    } else {
        return(predict.vglm(object, newdata=newdata,
                            type=type, se.fit=se.fit,
                            deriv=deriv, 
                            dispersion=dispersion, ...))
    }

    na.act = object@na.action

    if(!length(newdata) && length(na.act)) {
        if(se.fit) {
            pred$fitted.values = napredict(na.act[[1]], pred$fitted.values)
            pred$se.fit = napredict(na.act[[1]], pred$se.fit)
        } else {
            pred = napredict(na.act[[1]], pred)
        }
    }

    pred
}


setMethod("predict", "rrvglm", function(object, ...) 
    predict.rrvglm(object, ...))



if(FALSE) {
predict.rrvlm = function(object, 
                         newdata=NULL, 
                         type=c("link", "response", "terms"),
                         se.fit=FALSE, 
                         deriv=0,
                         dispersion=NULL, 
                         extra=object@extra, ...) 
{
    stop("this function hasn't been written yet")
}


setMethod("predict", "rrvlm", function(object, ...) 
    predict.rrvlm(object, ...))
}



