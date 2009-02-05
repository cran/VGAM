# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.








if(TRUE) {


AICvlm = function(object, ..., k=2) {
    estdisp = object@misc$estimated.dispersion
    no.dpar = if(length(estdisp) && is.logical(estdisp) && estdisp)
        length(object@misc$dispersion) else 0
    -2 * logLik.vlm(object, ...) + k * (length(coefvlm(object)) + no.dpar)
}

AICvgam = function(object, ..., k=2) {
    estdisp = object@misc$estimated.dispersion
    no.dpar = if(length(estdisp) && is.logical(estdisp) && estdisp)
        length(object@misc$dispersion) else 0 
    nldf = if(is.Numeric(object@nl.df)) sum(object@nl.df) else 0
    -2 * logLik.vlm(object, ...) + k * (length(coefvlm(object)) + no.dpar + nldf)
}

AICrrvglm = function(object, ..., k=2) {
    estdisp = object@misc$estimated.dispersion
    no.dpar = if(length(estdisp) && is.logical(estdisp) && estdisp)
        length(object@misc$dispersion) else 0 
    elts.tildeA = (object@misc$M - object@control$Rank) * object@control$Rank
    -2 * logLik.vlm(object, ...) + k * (length(coefvlm(object)) +
    no.dpar + elts.tildeA)
}

AICqrrgvlm = function(object, ..., k=2) {
    stop("this function not written yet")
}

setMethod("AIC", "vlm",
         function(object, ..., k=2)
         AICvlm(object, ..., k=k))

setMethod("AIC", "vglm",
         function(object, ..., k=2)
         AICvlm(object, ..., k=k))

setMethod("AIC", "vgam",
         function(object, ..., k=2)
         AICvgam(object, ..., k=k))

setMethod("AIC", "rrvglm",
         function(object, ..., k=2)
         AICrrvglm(object, ..., k=k))

setMethod("AIC", "qrrvglm",
         function(object, ..., k=2)
         AICqrrvglm(object, ..., k=k))
}


if(FALSE) {



AICvglm = function(object, ..., k=2) {
    crit = logLik.vlm(object, ...)
    -2 * crit + k * length(coef(object))
}





AICrrvglm = function(object, ..., k=2) {
    stop("not working yet")
    crit = logLik.vlm(object)
    sign = -2
    if(!length(crit) || !is.numeric(crit)) {
        crit = deviance(object)
        sign = 1
    }
    if(!length(crit) || !is.numeric(crit))
        stop("can't get at the deviance or loglikelihood of the object")

    sign * crit + 2 * (length(coef(object)) +
    object@control$rank * (object@misc$M - object@control$rank))
}



    # setGeneric("AIC", function(object, ..., k = 2) standardGeneric("AIC"))

setMethod("AIC", signature(object="vglm"),
           function(object, ..., k=2)
           AICvglm(object, ..., k=k))



setMethod("AIC", signature(object="rrvglm"),
           function(object, ..., k=2)
           AICrrvglm(object, ..., k=k))


}


