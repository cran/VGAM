# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.





deviance.vlm <- function(object, ...)
    object@criterion$deviance


deviance.vglm <- function(object, ...)
    object@criterion$deviance



if(!isGeneric("deviance"))
    setGeneric("deviance", function(object, ...)
    standardGeneric("deviance"))


setMethod("deviance", "vlm", function(object, ...)
           deviance.vlm(object, ...))


setMethod("deviance", "vglm", function(object, ...)
           deviance.vglm(object, ...))




df.residual_vlm <- function(object, type = c("vlm", "lm"), ...) {
  type <- type[1]
  switch(type,
         vlm = object@df.residual,
          lm = nobs(object, type = "lm") - nvar(object, type = "lm"),
         stop("argument 'type' unmatched"))
}




setMethod("df.residual", "vlm", function(object, ...)
           df.residual_vlm(object, ...))





