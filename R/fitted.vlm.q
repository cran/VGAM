# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.



fitted.vlm <- function(object, matrix=TRUE, ...)
{

    answer = 
    if(matrix)
        object@fitted.values else
    {
        if(!is.matrix(object@fitted.values) || !length(object@fitted.values))
            stop("object@fitted.values is not a matrix or is empty")
        if(ncol(object@fitted.values) == 1)
            c(object@fitted.values) else {
                warning("ncol(object@fitted.values) is not 1")
                c(object@fitted.values)
            }
    }

    if(length(answer) && length(object@na.action)) {
        napredict(object@na.action[[1]], answer)
    } else {
        answer
    }
}

setMethod("fitted.values",  "vlm",
    function(object, ...)
    fitted.vlm(object, ...))

setMethod("fitted",  "vlm",
    function(object, ...)
    fitted.vlm(object, ...))

setMethod("fitted.values",  "vglm",
    function(object, ...)
    fitted.vlm(object, ...))

setMethod("fitted",  "vglm",
    function(object, ...)
    fitted.vlm(object, ...))


predictors.vglm <- function(object, matrix=TRUE, ...)
{
    answer = 
    if(matrix)
        object@predictors else
    {
        if(!is.matrix(object@predictors) || !length(object@predictors))
            stop("object@predictors is not a matrix or is empty")
        if(ncol(object@predictors) == 1)
            c(object@predictors) else {
                warning("ncol(object@predictors) is not 1")
                c(object@predictors)
            }
    }

    if(length(answer) && length(object@na.action)) {
        napredict(object@na.action[[1]], answer)
    } else {
        answer
    }
}


if(!isGeneric("predictors")) 
    setGeneric("predictors", function(object, ...) standardGeneric("predictors"))

setMethod("predictors",  "vglm",
    function(object, ...)
    predictors.vglm(object, ...))





