# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.












check.omit.constant <- function(object) {



  if (is.logical(object@misc$needto.omit.constant) &&
       object@misc$needto.omit.constant &&
      !object@misc$omit.constant)
    warning("Probably 'omit.constant = TRUE' should have been set. ",
            "See the family function '",
            object@family@vfamily[1],
            "' help file.")

}






if (!isGeneric("nparam"))
  setGeneric("nparam", function(object, ...)
             standardGeneric("nparam"),
             package = "VGAM")



nparam.vlm <- function(object, dpar = TRUE, ...) {

  estdisp <- object@misc$estimated.dispersion

  check.omit.constant(object)

  no.dpar <- if (length(estdisp) && is.logical(estdisp) && estdisp)
    length(object@misc$dispersion) else 0

  tot.par <- length(coefvlm(object)) + as.numeric(dpar) * no.dpar
  tot.par
}





nparam.vgam <- function(object, dpar = TRUE,
                        linear.only = FALSE, ...) {

  estdisp <- object@misc$estimated.dispersion

  check.omit.constant(object)

  no.dpar <- if (length(estdisp) && is.logical(estdisp) && estdisp)
             length(object@misc$dispersion) else 0
  nldf <- if (is.Numeric(object@nl.df)) sum(object@nl.df) else 0

  if (linear.only) {
    length(coefvlm(object)) + as.numeric(dpar) * no.dpar
  } else {
    length(coefvlm(object)) + as.numeric(dpar) * no.dpar + nldf
  }
}



nparam.rrvglm <- function(object, dpar = TRUE, ...) {

  check.omit.constant(object)

  estdisp <- object@misc$estimated.dispersion
  no.dpar <- if (length(estdisp) && is.logical(estdisp) && estdisp)
    length(object@misc$dispersion) else 0
  str0 <- object@control$str0
  MMM <- object@misc$M
  Rank <- object@control$Rank
  elts.tildeA <- (MMM - Rank - length(str0)) * Rank



  length(coefvlm(object)) + as.numeric(dpar) * no.dpar + elts.tildeA
}



nparam.qrrvglm <- function(object, dpar = TRUE, ...) {

  check.omit.constant(object)

  estdisp <- object@misc$estimated.dispersion
  no.dpar <- if (length(estdisp) && is.logical(estdisp) && estdisp)
             length(object@misc$dispersion) else 0
  str0 <- object@control$str0
  MMM <- object@misc$M
  Rank <- object@control$Rank
  elts.tildeA <- (MMM - Rank - length(str0)) * Rank




  eq.tolerances <- object@control$eq.tolerances
  I.tolerances <- object@control$I.tolerances
  if (!(length(eq.tolerances) == 1 && is.logical(eq.tolerances)))
    stop("could not determine whether the fitted object used an ",
         "equal-tolerances assumption based on ",
         "argument 'eq.tolerances'")
  if (!(length(I.tolerances) == 1 && is.logical(I.tolerances)))
    stop("could not determine whether the fitted object used an ",
         "equal-tolerances assumption based on argument 'I.tolerances'")


  NOS <- if (length(object@y)) ncol(object@y) else MMM
  MSratio <- MMM / NOS  # First value is g(mean) = quadratic form in l
  if (round(MSratio) != MSratio)
    stop("variable 'MSratio' is not an integer")
  elts.D <- ifelse(I.tolerances || eq.tolerances, 1, NOS) *
            Rank * (Rank + 1) / 2







  elts.B1 <- length(object@extra$B1)
  elts.C  <- length(object@extra$Cmat)
  num.params <- elts.B1 + elts.tildeA  + elts.D + elts.C

  num.params
}



nparam.rrvgam <- function(object, dpar = TRUE, ...) {



  check.omit.constant(object)


  estdisp <- object@misc$estimated.dispersion
  no.dpar <- if (length(estdisp) && is.logical(estdisp) && estdisp)
             length(object@misc$dispersion) else 0
  str0 <- object@control$str0
  MMM <- object@misc$M
  Rank <- object@control$Rank




  NOS <- if (length(object@y)) ncol(object@y) else MMM
  MSratio <- MMM / NOS  # First value is g(mean) = quadratic form in l
  if (round(MSratio) != MSratio)
    stop("variable 'MSratio' is not an integer")





  elts.B1     <- length(object@extra$B1)  # 0 since a NULL
  elts.C      <- length(object@extra$Cmat)
  elts.df1.nl <-    sum(object@extra$df1.nl)

  num.params <- elts.B1 + elts.C + (
                2 * length(object@extra$df1.nl) + elts.df1.nl) -
                (Rank + length(str0)) * Rank


  num.params
}



setMethod("nparam", "vlm",
           function(object, ...)
           nparam.vlm(object, ...))

setMethod("nparam", "vglm",
           function(object, ...)
           nparam.vlm(object, ...))

setMethod("nparam", "vgam",
           function(object, ...)
           nparam.vgam(object, ...))

setMethod("nparam", "rrvglm",
           function(object, ...)
           nparam.rrvglm(object, ...))

setMethod("nparam", "qrrvglm",
           function(object, ...)
           nparam.qrrvglm(object, ...))


setMethod("nparam", "rrvgam",
           function(object, ...)
           nparam.rrvgam(object, ...))









if (!isGeneric("AIC"))
  setGeneric("AIC", function(object, ..., k = 2)
             standardGeneric("AIC"),
             package = "VGAM")






AICvlm <- function(object, ...,
                   corrected = FALSE,
                   k = 2) {
  estdisp <- object@misc$estimated.dispersion


  tot.par <- nparam.vlm(object, dpar = TRUE)
  ans <- (-2) * logLik.vlm(object, ...) + k * tot.par

  if (corrected) {
    ans <- ans + k * tot.par * (tot.par + 1) / (
           nobs(object) - tot.par - 1)
  }
  ans
}




AICvgam <- function(object, ...,
                    k = 2) {


  sum.lco.no.dpar.nldf <- nparam.vgam(object, dpar = TRUE,
                                      linear.only = FALSE)

  -2 * logLik.vlm(object, ...) + k * sum.lco.no.dpar.nldf
}




AICrrvglm <- function(object, ...,
                      k = 2) {



  sum.lco.no.dpar.A <- nparam.rrvglm(object, dpar = TRUE)
  (-2) * logLik.vlm(object, ...) + k * sum.lco.no.dpar.A
}




AICqrrvglm <- function(object, ...,
                       k = 2) {







  loglik.try <- logLik.qrrvglm(object, ...)
  if (!is.numeric(loglik.try))
    warning("cannot compute the log-likelihood of 'object'. ",
            "Returning NULL")

  num.params <- nparam.qrrvglm(object, dpar = TRUE)


  if (is.numeric(loglik.try)) {
    (-2) * loglik.try     + k * num.params
  } else {

    NULL
  }
}






 AICrrvgam <- function(object, ...,
                       k = 2) {







  loglik.try <- logLik(object, ...)
  if (!is.numeric(loglik.try))
    warning("cannot compute the log-likelihood of 'object'. ",
            "Returning NULL")

  num.params <- nparam.rrvgam(object, dpar = TRUE)


  if (is.numeric(loglik.try)) {
    (-2) * loglik.try     + k * num.params
  } else {

    NULL
  }
}





setMethod("AIC", "vlm",
         function(object, ..., k = 2)
           AICvlm(object, ..., k = k))

setMethod("AIC", "vglm",
         function(object, ..., k = 2)
           AICvlm(object, ..., k = k))

setMethod("AIC", "vgam",
         function(object, ..., k = 2)
          AICvgam(object, ..., k = k))

setMethod("AIC", "rrvglm",
           function(object, ..., k = 2)
          AICrrvglm(object, ..., k = k))

setMethod("AIC", "qrrvglm",
            function(object, ..., k = 2)
          AICqrrvglm(object, ..., k = k))


setMethod("AIC", "rrvgam",
          function(object, ..., k = 2)
            AICrrvgam(object, ..., k = k))




if (!isGeneric("AICc"))
  setGeneric("AICc", function(object, ..., k = 2)
             standardGeneric("AICc"),
             package = "VGAM")


setMethod("AICc", "vlm",
         function(object, ..., k = 2)
         AICvlm(object, ..., corrected = TRUE, k = k))

setMethod("AICc", "vglm",
         function(object, ..., k = 2)
         AICvlm(object, ..., corrected = TRUE, k = k))

















if (!isGeneric("BIC"))
  setGeneric("BIC", function(object, ..., k = log(nobs(object)))
             standardGeneric("BIC"),
             package = "VGAM")


BICvlm <- function(object, ..., k = log(nobs(object))) {
  AICvlm(object, ..., k = k)
}



setMethod("BIC", "vlm",
          function(object, ..., k = log(nobs(object)))
            BICvlm(object, ..., k = k))

setMethod("BIC", "vglm",
          function(object, ..., k = log(nobs(object)))
            BICvlm(object, ..., k = k))

setMethod("BIC", "vgam",
          function(object, ..., k = log(nobs(object)))
           AICvgam(object, ..., k = k))

setMethod("BIC", "rrvglm",
           function(object, ..., k = log(nobs(object)))
          AICrrvglm(object, ..., k = k))

setMethod("BIC", "qrrvglm",
            function(object, ..., k = log(nobs(object)))
          AICqrrvglm(object, ..., k = k))


setMethod("BIC", "rrvgam",
          function(object, ..., k = log(nobs(object)))
            AICrrvgam(object, ..., k = k))










