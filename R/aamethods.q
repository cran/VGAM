# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.








is.Numeric <- function(x, allowable.length=Inf, integer.valued=FALSE, positive=FALSE)
    if(all(is.numeric(x)) && all(is.finite(x)) &&
    (if(is.finite(allowable.length)) length(x)==allowable.length else TRUE) &&
    (if(integer.valued) all(x==round(x)) else TRUE) &&
    (if(positive) all(x>0) else TRUE)) TRUE else FALSE


VGAMenv = new.env()


.onLoad <- function(lib, pkg) require(methods)  # 25/1/05
 
 
 
if(!any(search()=="package:methods"))
    library(methods)


if(!any(search()=="package:splines"))
    require(splines)



.VGAM.prototype.list = list(
      "constraints"  = expression({}),
      "fini"         = expression({}),
      "first"        = expression({}),
      "initialize"   = expression({}),
      "last"         = expression({}),
      "middle"       = expression({}),
      "middle2"      = expression({}),
      "deriv"        = expression({}),
      "weight"       = expression({}))


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
      "loglikelihood"= "function",
      "middle"       = "expression",
      "middle2"      = "expression",
      "summary.dispersion"  = "logical",
      "vfamily"      = "character",
      "deriv"        = "expression",
      "weight"       = "expression"),  #  "call"
prototype = .VGAM.prototype.list)


valid.vglmff = function(object) {
    compulsory = c("initialize", "weight", "deriv", "inverse")
    for(ii in compulsory) {
        if(!length(slot(object, ii)))
            stop("slot ", ii, " is empty")
    }

    if(length(as.list(object@inverse)) != 3)
        stop("wrong number of arguments in object@inverse")
}

if(FALSE) 
    setValidity("vglmff", valid.vglmff)





if(!isGeneric("print"))
    setGeneric("print", function(x, ...) standardGeneric("print"),
               package="VGAM")



print.vglmff <- function(x, ...) {
    f <- x@vfamily
    if(is.null(f))
        stop("not a VGAM family function")

    nn <- x@blurb
    if(!length(nn))
        invisible(return(x))

    cat("Family: ", f[1], "\n")
    if(length(f)>1) cat("Informal classes:", paste(f, collapse=", "), "\n")
    cat("\n")

    for(ii in 1:length(nn))
        cat(nn[ii])
    cat("\n")
    invisible(return(x))
}


setMethod("print", "vglmff",
         function(x, ...)
         invisible(print.vglmff(x, ...)))

setMethod("show", "vglmff",
          function(object)
              print.vglmff(x=object))










setClass("vlmsmall", representation(
      "call"         = "call",
      "coefficients" = "numeric",
      "constraints"  = "list",
      "control"      = "list",
      "criterion"    = "list",
      "fitted.values"= "matrix",
      "misc"         = "list",
      "model"        = "data.frame",
      "na.action"    = "list",
      "post"         = "list",
      "preplot"      = "list",
      "prior.weights"= "numeric",
      "residuals"    = "matrix",
      "weights"      = "matrix",
      "x"            = "matrix",
      "y"            = "matrix"),
)


setClass("vlm", representation(
      "assign"       = "list",
      "callXm2"      = "call",
      "contrasts"    = "list",
      "df.residual"  = "numeric",
      "df.total"     = "numeric",
      "dispersion"   = "numeric",
      "effects"      = "numeric",
      "offset"       = "matrix",
      "qr"           = "list",
      "R"            = "matrix",
      "rank"         = "integer",
      "rss"          = "numeric",
      "smart.prediction" = "list",
      "terms"        = "list",
      "Xm2"          = "matrix",
      "Ym2"          = "matrix",
      "xlevels"      = "list"
      ),
    contains = "vlmsmall"
)


setClass("vglm", representation(
      "extra"            = "list",
      "family"           = "vglmff",
      "iter"             = "numeric",
      "predictors"       = "matrix"),
    contains = "vlm")


setClass("vgam", representation(
      "Bspline"             = "list", # each [[i]] is a "vsmooth.spline.fit"
      "nl.chisq"            = "numeric",
      "nl.df"               = "numeric",
      "spar"                = "numeric",
      "s.xargument"         = "character",
      "var"                 = "matrix"),
    contains = "vglm")





setClass("summary.vgam", representation(
        anova="data.frame",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"),
    prototype(anova=data.frame()),
    contains = "vgam")


setClass("summary.vglm", representation(
        coef3="matrix",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"),
    contains = "vglm")


setClass("summary.vlm", representation(
        coef3="matrix",
        cov.unscaled="matrix",
        correlation="matrix",
        df="numeric",
        pearson.resid="matrix",
        sigma="numeric"),
    contains = "vlm")



 setClass(Class="rrvglm",
          contains="vglm")



if(FALSE)
 setClass("qrrvglm", representation(
      "assign"       = "list",
      "call"         = "call",
      "coefficients" = "numeric",
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
      "iter"         = "numeric",
      "misc"         = "list",
      "model"        = "data.frame",
      "na.action"    = "list",
      "offset"       = "matrix",
      "post"         = "list",
      "predictors"   = "matrix",
      "preplot"      = "list",
      "prior.weights"= "numeric",
      "residuals"    = "matrix",
      "smart.prediction" = "list",
      "terms"        = "list",
      "weights"      = "matrix",
      "x"            = "matrix",
      "Xm2"          = "matrix",
      "Ym2"          = "matrix",
      "xlevels"      = "list",
      "y"            = "matrix")
)




 setClass(Class="qrrvglm",
          contains = "rrvglm")






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



 setClass("grc", representation(not.needed="numeric"),
          contains = "rrvglm")


setMethod("summary", "grc",
          function(object, ...)
          summary.grc(object, ...))


if(FALSE) {
    setClass("vfamily",
        representation("list"))
}



if(!isGeneric("Coef"))
setGeneric("Coef", function(object, ...) standardGeneric("Coef"),
           package="VGAM")
if(!isGeneric("Coefficients"))
setGeneric("Coefficients", function(object, ...)
            standardGeneric("Coefficients"),
           package="VGAM")








if(!isGeneric("logLik"))
    setGeneric("logLik", function(object, ...) standardGeneric("logLik"),
           package="VGAM")

if(!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"),
           package="VGAM")

if(!isGeneric("vcov"))
    setGeneric("vcov", function(object, ...) standardGeneric("vcov"),
           package="VGAM")








setClass("uqo", representation(
      "lv"               = "matrix",
      "extra"            = "list",
      "family"           = "vglmff",
      "iter"             = "numeric",
      "predictors"       = "matrix"),
    contains = "vlmsmall")


setClass(Class="cao",
         contains="vgam")


if(!isGeneric("lvplot"))
setGeneric("lvplot", function(object, ...) standardGeneric("lvplot"),
           package="VGAM")

if(!isGeneric("ccoef"))
    setGeneric("ccoef", function(object, ...) standardGeneric("ccoef"),
           package="VGAM")





if(!isGeneric("coef"))
    setGeneric("coef", function(object, ...) standardGeneric("coef"),
           package="VGAM")

if(!isGeneric("coefficients"))
    setGeneric("coefficients", function(object, ...)
                               standardGeneric("coefficients"),
               package="VGAM")

if(!isGeneric("df.residual"))
    setGeneric("df.residual", function(object, ...)
                              standardGeneric("df.residual"),
           package="VGAM")

if(!isGeneric("fitted"))
    setGeneric("fitted", function(object, ...) standardGeneric("fitted"),
           package="VGAM")

 if(!isGeneric("fitted.values"))
     setGeneric("fitted.values", function(object, ...)
                                 standardGeneric("fitted.values"),
           package="VGAM")

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
    setGeneric("residuals", function(object, ...) standardGeneric("residuals"),
           package="VGAM")

if(!isGeneric("weights"))
    setGeneric("weights", function(object, ...) standardGeneric("weights"),
           package="VGAM")





if(!isGeneric("AIC"))
    setGeneric("AIC", function(object, ..., k=2) standardGeneric("AIC"),
           package="VGAM")



  if(!isGeneric("formula"))
      setGeneric("formula", function(x, ...) standardGeneric("formula"),
           package="VGAM")


  if(!isGeneric("case.names"))
      setGeneric("case.names", function(object, ...)
           standardGeneric("case.names"),
           package="VGAM")

  if(!isGeneric("variable.names"))
      setGeneric("variable.names", function(object, ...)
           standardGeneric("variable.names"),
           package="VGAM")













