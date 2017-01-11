# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.









logLik.vlm <- function(object,
                       summation = TRUE,
                       ...) {

  if (summation) {
    object@criterion$loglikelihood
  } else {


    Args <- formals(args(object@family@loglikelihood))
    if (length(Args$summation) == 0)
      stop("there is no 'summation' argument for the function in the ",
           "'loglikelihood' slot of the object.")


    object@family@loglikelihood(mu = fitted(object),
                                y = depvar(object),
                          w = as.vector(weights(object, type = "prior")),
                                residuals = FALSE,
                                eta = predict(object),
                                extra = object@extra,
                                summation = summation)
  }
}




logLik.qrrvglm <- function(object,
                           summation = TRUE,
                           ...) {

  ff.code <- object@family
  ll.ff.code <- ff.code@loglikelihood

  prior.weights <- weights(object, type = "prior")
  if (is.matrix(prior.weights) &&
      ncol(prior.weights) == 1)
    prior.weights <- c(prior.weights)

  loglik.try <-
    ll.ff.code(mu = fitted(object),
               y = depvar(object),
               w = prior.weights,
               residuals = FALSE,
               eta = predict(object),
               extra = object@extra,
               summation = summation)
  if (!is.numeric(loglik.try))
    loglik.try <- NULL

  loglik.try
}







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





setMethod("logLik",  "qrrvglm", function(object, ...)
  logLik.qrrvglm(object, ...))


setMethod("logLik",  "rrvgam", function(object, ...)
  logLik.qrrvglm(object, ...))























constraints.vlm <-
  function(object,
           type = c("lm", "term"),
           all = TRUE, which,
           matrix.out = FALSE,
           colnames.arg = TRUE,  # 20130827
           ...) {


  type <- match.arg(type, c("lm", "term"))[1]


  Hlist <- ans <- slot(object, "constraints")  # For "lm" (formerly "vlm")

  if (type == "term") {
    oassign.LM <- object@misc$orig.assign

    x.LM <- model.matrix(object)
    att.x.LM  <- attr(x.LM,  "assign")
    names.att.x.LM <- names(att.x.LM)
    ppp <- length(names.att.x.LM)


    ans <- vector("list", ppp)
    for (ii in 1:ppp) {
      col.ptr <- (oassign.LM[[ii]])[1]  # 20110114
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


      if (colnames.arg)
        dimnames(mat.ans) <-
          list(NULL,
               colnames(model.matrix(object, type = "vlm")))


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







