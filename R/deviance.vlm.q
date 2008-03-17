# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.



deviance.vlm <- function(object, ...)
    object@criterion$deviance

deviance.vglm <- function(object, ...)
    object@criterion$deviance




if(!isGeneric("deviance"))
    setGeneric("deviance", function(object, ...) standardGeneric("deviance"))


setMethod("deviance", "vlm", function(object, ...)
           deviance.vlm(object, ...))

if(is.R()) {


    setMethod("deviance", "vglm", function(object, ...)
               deviance.vglm(object, ...))
} else {
    setMethod("deviance", "vglm", function(object, ...)
               deviance.vglm(object, ...))
}




df.residual.vlm <- function(object, ...)
    object@df.residual

if(is.R()) {


    setMethod("df.residual", "vlm", function(object, ...)
               df.residual.vlm(object, ...))
} else {
    if(!isGeneric("df.residual"))
    setGeneric("df.residual", function(object, ...)
               standardGeneric("df.residual"))
    setMethod("df.residual", "vlm", function(object, ...)
               df.residual.vlm(object, ...))
}





