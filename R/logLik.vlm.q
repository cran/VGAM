# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.



logLik.vlm <- function(object, ...)
        object@criterion$loglikelihood

if(!isGeneric("logLik"))
    setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

setMethod("logLik",  "vlm", function(object, ...)
    logLik.vlm(object, ...))




if(TRUE) {
constraints.vlm <- function(object, all=TRUE, which, ...) 
    if(all) slot(object, "constraints") else 
    slot(object, "constraints")[[which]]

if(!isGeneric("constraints"))
setGeneric("constraints", function(object, ...) standardGeneric("constraints"))
setMethod("constraints",  "vlm", function(object, ...)
    constraints.vlm(object, ...))

}




