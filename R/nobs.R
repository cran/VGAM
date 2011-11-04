# These functions are Copyright (C) 1998-2011 T. W. Yee    All rights reserved.

# nobs.R 

# Notes.
# 1. 20110711 Looked at "NEWS" and found out about nobs().
#    Adding nvar() too while I am at it.


# ======================================================================
# 20110711

nobs.vlm <- function(object, type = c("lm", "vlm"), ...) {

# Notes:
# 1. with type = "vlm" this is n * M.

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("lm", "vlm"))[1]


  if (type == "lm") {
    object@misc$n
  } else {
    object@misc$nrow_X_vlm
  }
}


setMethod("nobs", "vlm",
         function(object, ...)
         nobs.vlm(object, ...))






# ======================================================================
# 20110711
# Here is the 'nvar' methods functions.
# Tricky for "vgam", "rrvglm", "qrrvglm", "cao", "rcam" objects?

nvar.vlm <- function(object, type = c("vlm", "lm"), ...) {

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("vlm", "lm"))[1]


  if (type == "lm") {
    object@misc$p
  } else {
    object@misc$ncol_X_vlm
  }
}



nvar.vgam <- function(object, type = c("vgam", "zz"), ...) {
# 20110711
# Uses the effective dof, or edof, or edf zz??

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("vgam", "zz"))[1]

  stop("function nvar.vgam() has not been written yet")

  if (type == "vgam") {
    object@misc$p
  } else {
    object@misc$ncol_X_vlm
  }
}


nvar.rrvglm <- function(object, type = c("rrvglm", "zz"), ...) {
# 20110711
# Uses the effective dof, or edof, or edf zz??

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("rrvglm", "zz"))[1]

  stop("function nvar.rrvglm() has not been written yet")

  if (type == "vgam") {
    object@misc$p
  } else {
    object@misc$ncol_X_vlm
  }
}



nvar.qrrvglm <- function(object, type = c("qrrvglm", "zz"), ...) {
# 20110711
# Uses the effective dof, or edof, or edf zz??

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("qrrvglm", "zz"))[1]

  stop("function nvar.qrrvglm() has not been written yet")

  if (type == "qrrvglm") {
    object@misc$p
  } else {
    object@misc$ncol_X_vlm
  }
}



nvar.cao <- function(object, type = c("cao", "zz"), ...) {
# 20110711
# Uses the effective dof, or edof, or edf zz??

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("rrvglm", "zz"))[1]

  stop("function nvar.cao() has not been written yet")

  if (type == "cao") {
    object@misc$p
  } else {
    object@misc$ncol_X_vlm
  }
}



nvar.rcam <- function(object, type = c("rcam", "zz"), ...) {
# 20110711
# Uses the effective dof, or edof, or edf zz??

  if(mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("rcam", "zz"))[1]

  stop("function nvar.rcam() has not been written yet")

  if (type == "rcam") {
    object@misc$p
  } else {
    object@misc$ncol_X_vlm
  }
}





if (!isGeneric("nvar"))
setGeneric("nvar", function(object, ...)
            standardGeneric("nvar"),
           package = "VGAM")


setMethod("nvar", "vlm",
         function(object, ...)
         nvar.vlm(object, ...))



setMethod("nvar", "vgam",
         function(object, ...)
         nvar.vgam(object, ...))


setMethod("nvar", "rrvglm",
         function(object, ...)
         nvar.rrvglm(object, ...))



setMethod("nvar", "qrrvglm",
         function(object, ...)
         nvar.qrrvglm(object, ...))



setMethod("nvar", "cao",
         function(object, ...)
         nvar.cao(object, ...))



setMethod("nvar", "rcam",
         function(object, ...)
         nvar.rcam(object, ...))


# ======================================================================



