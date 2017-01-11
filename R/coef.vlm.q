# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.






coef.vlm <- function(object, ...) {

  coefvlm(object, ...)
}



coefvlm <- function(object, matrix.out = FALSE, label = TRUE,
                    colon = FALSE) {

  ans <- object@coefficients

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
    names(ans) <- new.labs
    return(ans)
  }


  if (!label)
    names(ans) <- NULL
  if (!matrix.out)
    return(ans)



  ncolx <- object@misc$p  # = length(object@constraints)
  M <- object@misc$M

  Hlist <- object@constraints
  if (all(trivial.constraints(Hlist) == 1)) {
    Bmat <- matrix(ans, nrow = ncolx, ncol = M, byrow = TRUE)
  } else {
    Bmat <- matrix(NA_real_, nrow = ncolx, ncol = M)

    if (!matrix.out)
      return(ans)

    ncolHlist <- unlist(lapply(Hlist, ncol))
    nasgn <- names(Hlist)
    temp <- c(0, cumsum(ncolHlist))
    for (ii in seq_along(nasgn)) {
      index <- (temp[ii] + 1):temp[ii + 1]
      cmat <- Hlist[[nasgn[ii]]]
      Bmat[ii,] <- cmat %*% ans[index]
    }
  }

  if (label) {
    d1 <- object@misc$colnames.x
    d2 <- object@misc$predictors.names  # Could be NULL
    dimnames(Bmat) <- list(d1, d2)
  }

  Bmat
}  # end of coefvlm




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


  answer <-
  if (length(tmp2 <- object@misc$link) &&
    object@misc$intercept.only &&
    trivial.constraints(object@constraints)) {




    if (!is.list(use.earg <- object@misc$earg))
      use.earg <- list()

    answer <-
      eta2theta(rbind(coefvlm(object)),
                link = object@misc$link,
                earg = use.earg)


    answer <- c(answer)
    if (length(ntmp2 <- names(tmp2)) == object@misc$M)
      names(answer) <- ntmp2
    answer
  } else {
    coefvlm(object, ... )
  }



  if (length(tmp3 <- object@misc$parameter.names) &&
    object@misc$intercept.only &&
    trivial.constraints(object@constraints)) {
    answer <- c(answer)
    if (length(tmp3) == object@misc$M &&
        is.character(tmp3))
      names(answer) <- tmp3
  }


  answer
}





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







