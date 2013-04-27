# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.







logLik.vlm <- function(object, ...)
  object@criterion$loglikelihood


if (!isGeneric("logLik"))
  setGeneric("logLik", function(object, ...)
             standardGeneric("logLik"),
             package = "VGAM")



setMethod("logLik",  "vlm", function(object, ...)
    logLik.vlm(object, ...))


setMethod("logLik",  "vglm", function(object, ...)
    logLik.vlm(object, ...))


setMethod("logLik",  "vgam", function(object, ...)
    logLik.vlm(object, ...))








constraints.vlm <- function(object,
                            type = c("lm", "term"),
                            all = TRUE, which,
                            matrix.out = FALSE,
                            ...) {


  type <- match.arg(type, c("lm", "term"))[1]


  Hlist <- ans <- slot(object, "constraints") # For "lm" (formerly "vlm")

  if (type == "term") {
    oassign.LM <- object@misc$orig.assign

    x.LM <- model.matrix(object)
    att.x.LM  <- attr(x.LM,  "assign")
    names.att.x.LM <- names(att.x.LM)
    ppp <- length(names.att.x.LM)

    ans <- vector("list", ppp)
    for (ii in 1:ppp) {
      col.ptr <- (oassign.LM[[ii]])[1] # 20110114
      ans[[ii]] <- (Hlist[[col.ptr]])
    }
    names(ans) <- names.att.x.LM
  } # End of "term"

  if (matrix.out) {
    if (all) {
      M <- npred(object)
      mat.ans <- matrix(unlist(ans), nrow = M)
      if (length(object@misc$predictors.names) == M)
        rownames(mat.ans) <- object@misc$predictors.names
      if (length(object@misc$colnames.X_vlm) == ncol(mat.ans))
        colnames(mat.ans) <- object@misc$colnames.X_vlm
      mat.ans
    } else {
      ans[[which]]
    }
  } else {
    if (all) ans else ans[[which]]
  }
}



if (!isGeneric("constraints"))
  setGeneric("constraints", function(object, ...)
             standardGeneric("constraints"))


setMethod("constraints",  "vlm", function(object, ...)
  constraints.vlm(object, ...))







