# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.





coef.vlm <- function(object, ...) {
  coefvlm(object, ...)
}



coefvlm <- function(object, matrix.out = FALSE, label = TRUE,
                    colon = FALSE) {
  Ans <- object@coefficients

  if (colon) {
    if (matrix.out)
      stop("cannot have 'matrix.out = TRUE' and 'colon = TRUE'")
    if (!label)
      stop("cannot have 'label = FALSE' and 'colon = TRUE'")

    d1 <- object@misc$colnames.x
    Hlist <- object@constraints
    M <- object@misc$M
    ncolHlist <- unlist(lapply(Hlist, ncol))
    new.labs <- vlabel(xn = d1, ncolHlist, M = M, colon = colon)
    names(Ans) <- new.labs
    return(Ans)
  }

  if (!label)
    names(Ans) <- NULL
  if (!matrix.out)
    return(Ans)


  ncolx <- object@misc$p  # = length(object@constraints)
  M <- object@misc$M

  Hlist <- object@constraints
  if (all(trivial.constraints(Hlist) == 1)) {
    Bmat <- matrix(Ans, nrow = ncolx, ncol = M, byrow = TRUE)
  } else {
    Bmat <- matrix(NA_real_, nrow = ncolx, ncol = M)

    if (!matrix.out)
      return(Ans)

    ncolHlist <- unlist(lapply(Hlist, ncol))
    nasgn <- names(Hlist)
    temp <- c(0, cumsum(ncolHlist))
    for (ii in seq_along(nasgn)) {
      index <- (temp[ii] + 1):temp[ii + 1]
      cmat <- Hlist[[nasgn[ii]]]
      Bmat[ii, ] <- cmat %*% Ans[index]
    }
  }

  if (label) {
    d1 <- object@misc$colnames.x
    d2 <- object@misc$predictors.names  # Could be NULL
    dimnames(Bmat) <- list(d1, d2)
  }

  Bmat
}  # coefvlm



setMethod("coefficients", "vlm", function(object, ...)
           coefvlm(object, ...))
setMethod("coef", "vlm", function(object, ...)
           coefvlm(object, ...))
setMethod("coefficients", "vglm", function(object, ...)
           coefvlm(object, ...))
setMethod("coef", "vglm", function(object, ...)
           coefvlm(object, ...))




setMethod("coefficients", "summary.vglm", function(object, ...)
          object@coef3)
setMethod("coef",         "summary.vglm", function(object, ...)
          object@coef3)




Coef.vlm <- function(object, ...) {

  LL <- length(object@family@vfamily)
  funname <- paste("Coef.", object@family@vfamily[LL], sep = "")

  if (exists(funname)) {
    newcall <- paste("Coef.", object@family@vfamily[LL],
                    "(object, ...)", sep = "")
    newcall <- parse(text = newcall)[[1]]
    return(eval(newcall))
  }

  Answer <-
    if (length(tmp2 <- object@misc$link) != 0 &&
        object@misc$intercept.only &&
        all(as.logical(trivial.constraints(object@constraints)))) {



    if (!is.list(use.earg <- object@misc$earg))
      use.earg <- list()

    Answer <- eta2theta(rbind(coefvlm(object)),
                        link = object@misc$link, earg = use.earg)

    Answer <- c(Answer)
    if (length(ntmp2 <- names(tmp2)) == object@misc$M) {
      special.case <- sum(object@misc$link == "multilogitlink") > 0
      try.this <- object@family@infos()$parameters.names
      names(Answer) <- if (special.case &&
                           length(try.this) == length(Answer))
        try.this else ntmp2
    }
    Answer
  } else {
    coefvlm(object, ... )
  }

  if (length(tmp3 <- object@misc$parameter.names) != 0 &&
      object@misc$intercept.only &&
      all(as.logical(trivial.constraints(object@constraints)))) {
    Answer <- c(Answer)
    if (length(tmp3) == object@misc$M && is.character(tmp3))
      names(Answer) <- tmp3
  }

  Answer
}  # Coef.vlm



setMethod("Coefficients", "vlm", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coef", "vlm", function(object, ...)
               Coef.vlm(object, ...))





coefvgam <-
  function(object, type = c("linear", "nonlinear"), ...) {
  type <- match.arg(type, c("linear", "nonlinear"))[1]

  if (type == "linear") {
    coefvlm(object, ...)
  } else {
    object@Bspline
  }
}


setMethod("coefficients", "vgam",
          function(object, ...)
          coefvgam(object, ...))


setMethod("coef", "vgam",
          function(object, ...)
          coefvgam(object, ...))



