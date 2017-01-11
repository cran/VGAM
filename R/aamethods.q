# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.








is.Numeric <- function(x, length.arg = Inf,
                       integer.valued = FALSE, positive = FALSE)
    if (all(is.numeric(x)) && all(is.finite(x)) &&
    (if (is.finite(length.arg))
       length(x) == length.arg else TRUE) &&
    (if (integer.valued) all(x == round(x)) else TRUE) &&
    (if (positive) all(x>0) else TRUE)) TRUE else FALSE


VGAMenv <- new.env()













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
      "blurb"         = "character",
      "constraints"   = "expression",
      "deviance"      = "function",
      "fini"          = "expression",
      "first"         = "expression",
      "infos"         = "function",  # Added 20101203
      "initialize"    = "expression",
      "last"          = "expression",
      "linkfun"       = "function",
      "linkinv"       = "function",
      "loglikelihood" = "function",
      "middle"        = "expression",
      "middle2"       = "expression",
      "summary.dispersion"  = "logical",
      "vfamily"       = "character",
      "validparams"   = "function",  # Added 20160305
      "validfitted"   = "function",  # Added 20160305
      "simslot"       = "function",
      "deriv"         = "expression",
      "weight"        = "expression"),  #  "call"
         prototype = .VGAM.prototype.list )


valid.vglmff <- function(object) {
  compulsory <- c("initialize", "weight", "deriv", "linkinv")
  for (ii in compulsory) {
    if (!length(slot(object, ii)))
      stop("slot ", ii, " is empty")
  }

  if (length(as.list(object@linkinv)) != 3)
    stop("wrong number of arguments in object@linkinv")
}


if (FALSE)
  setValidity("vglmff", valid.vglmff)











show.vglmff <- function(object) {
  f <- object@vfamily
  if (is.null(f))
      stop("not a VGAM family function")

  nn <- object@blurb

  cat("Family: ", f[1], "\n")
  if (length(f) > 1)
    cat("Informal classes:", paste(f, collapse = ", "), "\n")
  cat("\n")

  for (ii in seq_along(nn))
    cat(nn[ii])
  cat("\n")



}






setMethod("show", "vglmff",
          function(object)
              show.vglmff(object = object))

















setClass("vlmsmall", representation(
      "call"          = "call",
      "coefficients"  = "numeric",
      "constraints"   = "list",
      "control"       = "list",
      "criterion"     = "list",
      "fitted.values" = "matrix",
      "misc"          = "list",
      "model"         = "data.frame",
      "na.action"     = "list",
      "post"          = "list",
      "preplot"       = "list",
      "prior.weights" = "matrix",
      "residuals"     = "matrix",
      "weights"       = "matrix",
      "x"             = "matrix",
      "y"             = "matrix"),
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
      "ResSS"       = "numeric",
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
      "Bspline"             = "list",  # each [[i]] is a "vsmooth.spline.fit"
      "nl.chisq"            = "numeric",
      "nl.df"               = "numeric",
      "spar"                = "numeric",
      "s.xargument"         = "character",
      "var"                 = "matrix"),
    contains = "vglm")


setClass("pvgam", representation(
         "ospsslot"         = "list"),
         contains = "vglm")





setClass("summary.vgam", representation(
        anova = "data.frame",
        cov.unscaled = "matrix",
        correlation = "matrix",
        df = "numeric",
        pearson.resid = "matrix",
        sigma = "numeric"),
        prototype(anova = data.frame()),
        contains = "vgam")


setClass("summary.vglm", representation(
         coef3 = "matrix",
         cov.unscaled = "matrix",
         correlation = "matrix",
         df = "numeric",
         pearson.resid = "matrix",
         sigma = "numeric"),
         contains = "vglm")


setClass("summary.vlm", representation(
         coef3 = "matrix",
         cov.unscaled = "matrix",
         correlation = "matrix",
         df = "numeric",
         pearson.resid = "matrix",
         sigma = "numeric"),
         contains = "vlm")






setClass("summary.pvgam", representation(
         anova      = "data.frame",
         "ospsslot" = "list"),
         prototype(anova = data.frame()),
         contains = c("summary.vglm", "pvgam")
         )










 setClass(Class = "rrvglm",
          contains = "vglm")



if (FALSE)
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
      "prior.weights"= "matrix",
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




 setClass(Class = "qrrvglm",
          contains = "rrvglm")






if (FALSE)
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
 "ResSS"=from@ResSS,
 "smart.prediction"=from@smart.prediction,
 "terms"=from@terms,
 "weights"=from@weights,
 "x"=from@x,
 "xlevels"=from@xlevels,
 "y"=from@y))



 setClass("rcim0", representation(not.needed = "numeric"),
          contains = "vglm")  # Added 20110506
 setClass("rcim", representation(not.needed = "numeric"),
          contains = "rrvglm")
 setClass("grc",  representation(not.needed = "numeric"),
          contains = "rrvglm")


setMethod("summary", "rcim",
          function(object, ...)
          summary.rcim(object, ...))

setMethod("summary", "grc",
          function(object, ...)
          summary.grc(object, ...))


if (FALSE) {
    setClass("vfamily",
        representation("list"))
}



if (!isGeneric("Coef"))
  setGeneric("Coef", function(object, ...) standardGeneric("Coef"),
             package = "VGAM")
if (!isGeneric("Coefficients"))
  setGeneric("Coefficients", function(object, ...)
              standardGeneric("Coefficients"),
             package = "VGAM")








if (!isGeneric("logLik"))
    setGeneric("logLik", function(object, ...) standardGeneric("logLik"),
           package = "VGAM")

if (!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"),
           package = "VGAM")


if (!isGeneric("vcov"))
    setGeneric("vcov", function(object, ...) standardGeneric("vcov"),
           package = "VGAM")










setClass("uqo", representation(
      "latvar"           = "matrix",
      "extra"            = "list",
      "family"           = "vglmff",
      "iter"             = "numeric",
      "predictors"       = "matrix"),
    contains = "vlmsmall")


setClass(Class = "rrvgam",
         contains = "vgam")


if (!isGeneric("lvplot"))
setGeneric("lvplot", function(object, ...) standardGeneric("lvplot"),
           package = "VGAM")



 if (FALSE) {
 if (!isGeneric("ccoef"))
    setGeneric("ccoef", function(object, ...) {
    .Deprecated("concoef")

    standardGeneric("ccoef")
    })
}



 if (!isGeneric("concoef"))
    setGeneric("concoef", function(object, ...) {
    standardGeneric("concoef")
    })



















if (!isGeneric("model.matrix"))
    setGeneric("model.matrix", function(object, ...)
                               standardGeneric("model.matrix"))


if (!isGeneric("model.frame"))
    setGeneric("model.frame", function(formula, ...)
                              standardGeneric("model.frame"))







if (!isGeneric("predict"))
  setGeneric("predict", function(object, ...) standardGeneric("predict"))






if (!isGeneric("resid"))
  setGeneric("resid", function(object, ...) standardGeneric("resid"))










if (!isGeneric("AIC"))
  setGeneric("AIC", function(object, ..., k=2) standardGeneric("AIC"),
             package = "VGAM")





















if (!isGeneric("summary"))
  setGeneric("summary", function(object, ...)
             standardGeneric("summary"),
             package = "VGAM")




if (!isGeneric("QR.R"))
  setGeneric("QR.R", function(object, ...)
             standardGeneric("QR.R"),
             package = "VGAM")


setMethod("QR.R", "vglm",
          function(object, ...) {
  if (length(object@R)) object@R else {
    warning("empty 'R' slot on object. Returning a NULL")
    NULL
  }
})



if (!isGeneric("QR.Q"))
  setGeneric("QR.Q", function(object, ...)
             standardGeneric("QR.Q"),
             package = "VGAM")


setMethod("QR.Q", "vglm",
          function(object, ...) {
  qr.list <- object@qr
  if (length(qr.list)) {
    class(qr.list) <- "qr"
    qr.Q(qr.list)
  } else {
    warning("empty 'qr' slot on object. Returning a NULL")
    NULL
  }
})






if (!isGeneric("margeffS4VGAM"))
    setGeneric("margeffS4VGAM",
               function(object, subset = NULL,
                        VGAMff,
                        ...)
                 standardGeneric("margeffS4VGAM"),
               package = "VGAM")





if (!isGeneric("summaryvglmS4VGAM"))
    setGeneric("summaryvglmS4VGAM",
               function(object,
                        VGAMff,
                        ...)
                 standardGeneric("summaryvglmS4VGAM"),
               package = "VGAM")


if (!isGeneric("showsummaryvglmS4VGAM"))
    setGeneric("showsummaryvglmS4VGAM",
               function(object,
                        VGAMff,
                        ...)
                 standardGeneric("showsummaryvglmS4VGAM"),
               package = "VGAM")





if (!isGeneric("showvglmS4VGAM"))
    setGeneric("showvglmS4VGAM",
               function(object,
                        VGAMff,
                        ...)
                 standardGeneric("showvglmS4VGAM"),
               package = "VGAM")

if (!isGeneric("showvgamS4VGAM"))
    setGeneric("showvgamS4VGAM",
               function(object,
                        VGAMff,
                        ...)
                 standardGeneric("showvgamS4VGAM"),
               package = "VGAM")



if (!isGeneric("predictvglmS4VGAM"))
    setGeneric("predictvglmS4VGAM",
               function(object,
                        VGAMff,
                        ...)
                 standardGeneric("predictvglmS4VGAM"),
               package = "VGAM")










