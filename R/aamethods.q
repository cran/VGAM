# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.



if(!exists("is.R")) is.R <- function()
    exists("version") && !is.null(version$language) && version$language=="R"

is.Numeric <- function(x, allowable.length=Inf, integer.valued=FALSE, positive=FALSE)
    if(all(is.numeric(x)) && all(is.finite(x)) &&
    (if(is.finite(allowable.length)) length(x)==allowable.length else TRUE) &&
    (if(integer.valued) all(x==round(x)) else TRUE) &&
    (if(positive) all(x>0) else TRUE)) TRUE else FALSE


if(is.R())
    VGAMenv = new.env()





if(is.R()) {
    .onLoad <- function(lib, pkg) require(methods)  # 25/1/05
 


    if(!any(search()=="package:methods"))
        library(methods)

    if(!any(search()=="package:splines"))
        require(splines)

}








.VGAM.prototype.list = if(is.R())
list(
      "constraints"  = expression({}),
      "fini"         = expression({}),
      "first"        = expression({}),
      "initialize"   = expression({}),
      "last"         = expression({}),
      "middle"       = expression({}),
      "middle2"      = expression({}),
      "deriv"        = expression({}),
      "weight"       = expression({})) else 
list(
      "blurb"        = "",
      "constraints"  = expression({}),
      "deviance"     = function() {},
      "fini"         = expression({}),
      "first"        = expression({}),
      "initialize"   = expression({}),
      "inverse"      = function() {},
      "last"         = expression({}),
      "link"         = function() {},
      "loglikelihood"= function() {},
      "middle"       = expression({}),
      "middle2"      = expression({}),
      "summary.dispersion"  = FALSE,
      "vfamily"      = "",
      "deriv"        = expression({}),
      "weight"       = expression({})) # Splus doesn't use prototypes 


if(is.R())
setClass("vglmff", representation(
      "blurb"        = "character",
      "constraints"  = "expression",
      "deviance"     = "function",
      "fini"         = "expression",
      "first"        = "expression",
      "initialize"   = "expression",
      "inverse"      = "function",
      "last"         = "expression",
      "link"         = "function",
      "loglikelihood"= "function",   # problem: zz function() NULL if unspecified
      "middle"       = "expression",
      "middle2"      = "expression",
      "summary.dispersion"  = "logical",
      "vfamily"      = "character",
      "deriv"        = "expression",
      "weight"       = "expression"),  #  "call"
prototype = .VGAM.prototype.list) else 
setClass("vglmff", representation(
      "blurb"        = "character",
      "constraints"  = "expression",
      "deviance"     = "function",
      "fini"         = "expression",
      "first"        = "expression",
      "initialize"   = "expression",
      "inverse"      = "function",
      "last"         = "expression",
      "link"         = "function",
      "loglikelihood"= "function",   # problem: zz function() NULL if unspecified
      "middle"       = "expression",
      "middle2"      = "expression",
      "summary.dispersion"  = "logical",
      "vfamily"      = "character",
      "deriv"        = "expression",
      "weight"       = "expression"))


valid.vglmff = function(object) {

    compulsory = c("initialize", "weight", "deriv", "inverse")
    for(i in compulsory) {
        if(!length(slot(object, i)))
            stop(paste("slot \"", i, "\" is empty"))
    }

    if(length(as.list(object@inverse)) != 3)
        stop("wrong number of arguments in object@inverse")
}

if(FALSE) 
setValidity("vglmff", valid.vglmff)



print.vglmff <- function(x, ...)
{
    f <- x@vfamily
    if(is.null(f))
        stop("not a VGAM family function")

    nn <- x@blurb
    if(!length(nn))
        invisible(return(x))

    cat("Family: ", f[1], "\n")
    if(length(f)>1) cat("Informal classes:", paste(f, collapse=", "), "\n")
    cat("\n")

    for(i in 1:length(nn))
        cat(nn[i])
    cat("\n")
    invisible(return(x))
}


setMethod("print", "vglmff",
         function(x, ...)
         invisible(print.vglmff(x, ...)))

setMethod("show", "vglmff",
          function(object)
              print.vglmff(x=object))









if(is.R())
setClass("vlm", representation(
      "assign"       = "list",
      "call"         = "call",
      "coefficients" = if(is.R()) "numeric" else "named",
      "constraints"  = "list",
      "contrasts"    = "list",
      "control"      = "list",
      "criterion"    = "list",
      "df.residual"  = "numeric",
      "df.total"     = "numeric",
      "dispersion"   = "numeric",
      "effects"      = "numeric",
      "fitted.values"= "matrix",
      "misc"         = "list",
      "model"        = "data.frame",
      "na.action"    = "list",    # ' if(is.R()) "omit" else ' 
      "offset"       = "matrix",
      "post"         = "list",
      "preplot"      = "list",
      "prior.weights"= if(is.R()) "numeric" else "named",
      "qr"           = "list",
      "R"            = if(is.R()) "matrix" else "upper",
      "rank"         = "integer",
      "residuals"    = "matrix",
      "rss"          = "numeric",
      "smart.prediction" = "list",
      "terms"        = "list",
      "weights"      = "matrix",
      "x"            = if(is.R()) "matrix" else "model.matrix",
      "xlevels"      = "list",
      "y"            = "matrix")
) else 
setClass("vlm", representation(
      "assign"       = "list",
      "call"         = "call",
      "coefficients" = if(is.R()) "numeric" else "named",
      "constraints"  = "list",
      "contrasts"    = "list",
      "control"      = "list",
      "criterion"    = "list",
      "df.residual"  = "numeric",
      "df.total"     = "numeric",
      "dispersion"   = "numeric",
      "effects"      = "numeric",
      "fitted.values"= "matrix",
      "misc"         = "list",
      "model"        = "data.frame",
      "na.action"    = "list",    # ' if(is.R()) "omit" else ' 
      "offset"       = "matrix",
      "post"         = "list",
      "preplot"      = "list",
      "prior.weights"= if(is.R()) "numeric" else "named",
      "qr"           = "qr",
      "R"            = if(is.R()) "matrix" else "upper",
      "rank"         = "integer",
      "residuals"    = "matrix",
      "rss"          = "numeric",
      "smart.prediction" = "list",
      "terms"        = "list",
      "weights"      = "matrix",
      "x"            = if(is.R()) "matrix" else "model.matrix",
      "xlevels"      = "list",
      "y"            = "matrix")
)


setClass("vglm", representation("vlm",
      "extra"            = "list",
      "family"           = "vglmff",
      "iter"             = if(is.R()) "numeric" else "integer",
      "predictors"       = if(is.R()) "matrix" else "matrix"))


setClass("vgam", representation("vglm",
      "Bspline"             = "list", # each [[i]] is a "vsmooth.spline.fit"
      "nl.chisq"            = if(is.R()) "numeric" else "named",
      "nl.df"               = if(is.R()) "numeric" else "named",
      "spar"                = if(is.R()) "numeric" else "named",
      "s.xargument"         = if(is.R()) "character" else "named", 
      "var"                 = "matrix"))





if(is.R()) 
    setClass("summary.vgam",
             representation("vgam",
        anova="data.frame",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"),
prototype(anova=data.frame())) else 
    setClass("summary.vgam",
             representation("vgam",
        anova="data.frame",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"))


    setClass("summary.vglm",
             representation("vglm",
        coef3="matrix",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"))

    setClass("summary.vlm",
             representation("vlm",
        coef3="matrix",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"))



 setClass( "rrvglm", representation("vglm"))

 setClass("qrrvglm", representation(
      "assign"       = "list",
      "call"         = "call",
      "coefficients" = if(is.R()) "numeric" else "named",
      "constraints"  = "list",
      "contrasts"    = "list",
      "control"      = "list",
      "criterion"    = "list",
      "df.residual"  = "numeric",
      "df.total"     = "numeric",
      "dispersion"   = "numeric",
      "extra"        = "list",
      "family"       = "vglmff",
      "fitted.values"= "matrix",
      "iter"         = if(is.R()) "numeric" else "integer",
      "misc"         = "list",
      "model"        = "data.frame",
      "na.action"    = "list",    # ' if(is.R()) "omit" else ' 
      "offset"       = "matrix",
      "post"         = "list",
      "predictors"   = if(is.R()) "matrix" else "matrix",
      "preplot"      = "list",
      "prior.weights"= if(is.R()) "numeric" else "named",
      "residuals"    = "matrix",
      "smart.prediction" = "list",
      "terms"        = "list",
      "weights"      = "matrix",
      "x"            = if(is.R()) "matrix" else "model.matrix",
      "xlevels"      = "list",
      "y"            = "matrix")
)

if(FALSE)
setAs("qrrvglm", "vglm", function(from)
new("vglm", "extra"=from@extra,
 "family"=from@family,
 "iter"=from@iter,
 "predictors"=from@predictors,
 "assign"=from@assign,
 "call"=from@call,
 "coefficients"=from@coefficients,
 "constraints"=from@constraints,
 "contrasts"=from@contrasts,
 "control"=from@control,
 "criterion"=from@criterion,
 "df.residual"=from@df.residual,
 "df.total"=from@df.total,
 "dispersion"=from@dispersion,
 "effects"=from@effects,
 "fitted.values"=from@fitted.values,
 "misc"=from@misc,
 "model"=from@model,
 "na.action"=from@na.action,
 "offset"=from@offset,
 "post"=from@post,
 "preplot"=from@preplot,
 "prior.weights"=from@prior.weights,
 "qr"=from@qr,
 "R"=from@R,
 "rank"=from@rank,
 "residuals"=from@residuals,
 "rss"=from@rss,
 "smart.prediction"=from@smart.prediction,
 "terms"=from@terms,
 "weights"=from@weights,
 "x"=from@x,
 "xlevels"=from@xlevels,
 "y"=from@y))



 setClass("grc", representation("rrvglm", not.needed="numeric"))


setMethod("summary", "grc",
          function(object, ...)
          summary.grc(object, ...))



if(FALSE) {
setClass("vfamily",
  representation("list"))
}




if(!isGeneric("Coef"))
setGeneric("Coef", function(object, ...) standardGeneric("Coef"))
if(!isGeneric("Coefficients"))
setGeneric("Coefficients", function(object, ...)
            standardGeneric("Coefficients"))







if(FALSE) {

if(!isGeneric("AIC"))
    setGeneric("AIC", function(object, ..., k=2) standardGeneric("AIC"))

AIC.vlm = function(object, ..., k=2) {
    ed = object@misc$estimated.dispersion
    no.dpar = if(length(ed) && is.logical(ed) && ed)
        length(object@misc$dispersion) else 0 
    -2 * logLik(object, ...) + k * (length(coef(object)) + no.dpar)
}

AIC.vgam = function(object, ..., k=2) {
    ed = object@misc$estimated.dispersion
    no.dpar = if(length(ed) && is.logical(ed) && ed)
        length(object@misc$dispersion) else 0 
    nldf = if(is.Numeric(object@nl.df)) sum(object@nl.df) else 0
    -2 * logLik(object, ...) + k * (length(coef(object)) + no.dpar + nldf)
}

AIC.rrvglm = function(object, ..., k=2) {
    ed = object@misc$estimated.dispersion
    no.dpar = if(length(ed) && is.logical(ed) && ed)
        length(object@misc$dispersion) else 0 
    elts.tildeA = (object@misc$M - object@control$Rank) * object@control$Rank
    -2 * logLik(object, ...) + k * (length(coef(object)) + no.dpar + elts.tildeA)
}

AIC.qrrgvlm = function(object, ..., k=2) {
    stop("this function not written yet")
}

setMethod("AIC", "vlm",
         function(object, ..., k=2)
         AIC.vlm(object, ..., k=k))

setMethod("AIC", "vglm",
         function(object, ..., k=2)
         AIC.vlm(object, ..., k=k))

setMethod("AIC", "vgam",
         function(object, ..., k=2)
         AIC.vgam(object, ..., k=k))

setMethod("AIC", "rrvglm",
         function(object, ..., k=2)
         AIC.rrvglm(object, ..., k=k))

setMethod("AIC", "qrrvglm",
         function(object, ..., k=2)
         AIC.qrrvglm(object, ..., k=k))
}

if(!isGeneric("logLik"))
    setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

if(!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

if(!isGeneric("vcov"))
    setGeneric("vcov", function(object, ...) standardGeneric("vcov"))








setClass("vlmsmall", representation(
      "call"         = "call",
      "coefficients" = if(is.R()) "numeric" else "named",
      "constraints"  = "list",
      "control"      = "list",
      "criterion"    = "list",
      "fitted.values"= "matrix",
      "misc"         = "list",
      "model"        = "data.frame",
      "na.action"    = "list",    # ' if(is.R()) "omit" else ' 
      "post"         = "list",
      "preplot"      = "list",
      "prior.weights"= if(is.R()) "numeric" else "named",
      "residuals"    = "matrix",
      "weights"      = "matrix",
      "x"            = if(is.R()) "matrix" else "model.matrix",
      "y"            = "matrix"),
)

setClass("uqo", representation("vlmsmall",
      "lv"               = "matrix",
      "extra"            = "list",
      "family"           = "vglmff",
      "iter"             = if(is.R()) "numeric" else "integer",
      "predictors"       = "matrix"))


setClass(Class="cao", repr=representation("vgam", "uqo"))


if(!isGeneric("lvplot"))
setGeneric("lvplot", function(object, ...) standardGeneric("lvplot"))

if(!isGeneric("ccoef"))
    setGeneric("ccoef", function(object, ...) standardGeneric("ccoef")) 





if(!isGeneric("coef"))
    setGeneric("coef", function(object, ...) standardGeneric("coef"))

if(!isGeneric("coefficients"))
    setGeneric("coefficients", function(object, ...)
                               standardGeneric("coefficients"))

if(!isGeneric("df.residual"))
    setGeneric("df.residual", function(object, ...)
                              standardGeneric("df.residual"))

if(!isGeneric("fitted"))
    setGeneric("fitted", function(object, ...) standardGeneric("fitted"))

 if(!isGeneric("fitted.values"))
     setGeneric("fitted.values", function(object, ...)
                                 standardGeneric("fitted.values"))

if(!isGeneric("model.matrix"))
    setGeneric("model.matrix", function(object, ...)
                               standardGeneric("model.matrix"))

if(!isGeneric("model.frame"))
    setGeneric("model.frame", function(formula, ...)
                              standardGeneric("model.frame"))


if(!isGeneric("predict"))
    setGeneric("predict", function(object, ...) standardGeneric("predict"))

if(!isGeneric("resid"))
    setGeneric("resid", function(object, ...) standardGeneric("resid"))

if(!isGeneric("residuals"))
    setGeneric("residuals", function(object, ...) standardGeneric("residuals"))

if(!isGeneric("weights"))
    setGeneric("weights", function(object, ...) standardGeneric("weights"))







