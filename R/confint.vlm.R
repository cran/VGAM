# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





confintvglm <-
  function(object, parm = "(All)", level = 0.95,
           method = c("wald", "profile"),
           trace = NULL,
           ...) {

  method <- match.arg(method, c("wald", "profile"))[1]


  cf <- coef(object)
  pnames <- names(cf)
  if (is.character(parm) && length(parm) == 1 && parm == "(All)")
    parm <- pnames else
  if (is.numeric(parm))
    parm <- pnames[parm]
  format.perc <- function(probs, digits)
    paste(format(100 * probs, trim = TRUE, scientific = FALSE,
          digits = digits), "%")



  if (method == "wald") {
    aa <- (1 - level) / 2
    aa <- c(aa, 1 - aa)
    pct <- format.perc(aa, 3)
    fac <- qnorm(aa)
    ci <- array(NA, dim = c(length(parm), 2L),
                dimnames = list(parm, pct))
    ses <- sqrt(diag(vcov(object)))[parm]
    ci[] <- cf[parm] + ses %o% fac
    return(ci)
  }  # if (method == "wald")




  ppv <- profilevglm(object, which = parm, alpha = (1 - level) / 4,
                     trace = trace, ...)



  MASSconfint.profile.glm(ppv, parm = parm, level = level,
                          trace = trace, ...)
}  # confintvglm



confintrrvglm <- function(object, parm, level = 0.95, ...) {
  stop("currently this function has not been written")

}


confintvgam <- function(object, parm, level = 0.95, ...) {
  stop("currently this function has not been written")
}








if (!isGeneric("confint"))
    setGeneric("confint",
               function(object, parm, level = 0.95, ...)
               standardGeneric("confint"),
           package = "VGAM")


setMethod("confint", "vglm",
          function(object, parm, level = 0.95, ...)
            confintvglm(object = object,
                        parm = if (missing(parm)) "(All)" else parm,
                        level = level, ...))


setMethod("confint", "rrvglm",
          function(object, parm, level = 0.95, ...)
          confintrrvglm(object = object, parm = parm, level = level, ...))

setMethod("confint", "vgam",
          function(object, parm, level = 0.95, ...)
            confintvgam(object = object, parm = parm, level = level, ...))






MASSconfint.profile.glm <-
function (object, parm = seq_along(pnames), level = 0.95, ...) {
  of <- attr(object, "original.fit")
  pnames <- names(coef(of))
  if (is.character(parm))
    parm <- match(parm, pnames, nomatch = 0L)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(round(100 * a, 1), "%")
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(pnames[parm],
      pct))
  cutoff <- qnorm(a)
  for (pm in parm) {
    pro <- object[[pnames[pm]]]
    if (is.null(pro))
      next
    if (length(pnames) > 1L)
      sp <- spline(x = pro[, "par.vals"][, pm], y = pro[, 1])
    else sp <- spline(x = pro[, "par.vals"], y = pro[, 1])
    ci[pnames[pm], ] <- approx(sp$y, sp$x, xout = cutoff)$y
  }
  drop(ci)
}

