# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.










residualsvlm  <-
  function(object,
           type = c("response", "deviance", "pearson", "working")) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("response", "deviance", "pearson", "working"))[1]

  na.act <- object@na.action
  object@na.action <- list()

  pooled.weight <- object@misc$pooled.weight
  if (is.null(pooled.weight))
    pooled.weight <- FALSE

  answer <-
  switch(type,
    working = if (pooled.weight) NULL else object@residuals,
    pearson = {
        if (pooled.weight) return(NULL)
        n <- object@misc$n
        M <- object@misc$M
        wz <- weights(object, type = "work")  # $weights
        if (!length(wz))
          wz <- if (M == 1) rep_len(1, n) else matrix(1, n, M)

        if (M == 1) {
          if (any(wz < 0))
            warning("some weights are negative. ",
                    "Their residual will be assigned NA")
          ans <- sqrt(c(wz)) * c(object@residuals)
          names(ans) <- names(object@residuals)
          ans
        } else {
          wz.sqrt <- matrix.power(wz, M = M, power = 0.5, fast = TRUE)
          ans <- mux22(wz.sqrt, object@residuals,
                       M = M, upper = FALSE)
          dim(ans) <- c(M, n)
          ans <- t(ans)
          dimnames(ans) <- dimnames(object@residuals)  # n x M
          ans
        }
    },
    deviance = {
      M <- object@misc$M
      if (M > 1)
        return(NULL)
      residualsvlm(object, type = "pearson")
    },
    response = object@residuals
  )

  if (length(answer) && length(na.act)) {
    napredict(na.act[[1]], answer)
  } else {
    answer
  }
}




residualsvglm  <-
  function(object,
           type = c("working", "pearson", "response", "deviance", "ldot"),
           matrix.arg = TRUE) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
          c("working", "pearson", "response", "deviance", "ldot"))[1]

  na.act <- object@na.action
  object@na.action <- list()

  pooled.weight <- object@misc$pooled.weight
  if (is.null(pooled.weight))
    pooled.weight <- FALSE

  answer <-
  switch(type,
    working = if (pooled.weight) NULL else object@residuals,
    pearson = {
      if (pooled.weight) return(NULL)

      n <- object@misc$n
      M <- object@misc$M
      wz <- weights(object, type = "work")   # $weights

      if (M == 1) {
        if (any(wz < 0))
          warning("some weights are negative. ",
                  "Their residual will be assigned NA")
        ans <- sqrt(c(wz)) * c(object@residuals)
        names(ans) <- names(object@residuals)
        ans
      } else {
        wz.sqrt <- matrix.power(wz, M = M, power = 0.5, fast = TRUE)
        ans <- mux22(wz.sqrt, object@residuals,
                     M = M, upper = FALSE)
        dim(ans) <- c(M, n)
        ans <- t(ans)
        dimnames(ans) <- dimnames(object@residuals)   # n x M
        ans
      }
    },
    deviance = {
      n <- object@misc$n

      y <- as.matrix(object@y)
      mu <- object@fitted.values


      w <- object@prior.weights
      if (!length(w))
        w <- rep_len(1, n)
      eta <- object@predictors

      dev.fn <- object@family@deviance  # May not 'exist' for that model
      if (length(body(dev.fn)) > 0) {
        extra <- object@extra
        ans <- dev.fn(mu = mu,y = y, w = w,
                      residuals = TRUE, eta = eta, extra)
        if (length(ans)) {
          lob <- labels(object@residuals)
          if (is.list(lob)) {
            if (is.matrix(ans))
              dimnames(ans) <- lob else
              names(ans) <- lob[[1]]
          } else {
            names(ans) <- lob
          }
        }
        ans
      } else {
        NULL
      }
    },
    ldot = {
      n <- object@misc$n
      y <- as.matrix(object@y)
      mu <- object@fitted
      w <- object@prior.weights
      if (is.null(w))
          w <- rep_len(1, n)
      eta <- object@predictors
      if (!is.null(ll.fn <- object@family@loglikelihood)) {
        extra <- object@extra
        ans <- ll.fn(mu = mu,y = y,w = w,
                     residuals = TRUE, eta = eta, extra)
        if (!is.null(ans)) {
          ans <- c(ans)  # ldot residuals can only be a vector
          names(ans) <- labels(object@residuals)
        }
        ans
      } else {
        NULL
      }
    },
    response = {
      y <- object@y

      mu <- fitted(object)

      true.mu <- object@misc$true.mu
      if (is.null(true.mu))
        true.mu <- TRUE

      ans <- if (true.mu) y - mu else NULL


      if (!matrix.arg && length(ans)) {
        if (ncol(ans) == 1) {
          names.ans <- dimnames(ans)[[1]]
          ans <- c(ans)
          names(ans) <- names.ans
          ans
        } else {
          warning("ncol(ans) is not 1")
          ans
        }
      } else {
        ans
      }
    })

  if (length(answer) && length(na.act)) {
    napredict(na.act[[1]], answer)
  } else {
    answer
  }
}





residualsqrrvglm  <- function(object,
                              type = c("response"),
                              matrix.arg = TRUE) {


  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("response"))[1]

  na.act <- object@na.action
  object@na.action <- list()

  pooled.weight <- object@misc$pooled.weight
  if (is.null(pooled.weight))
    pooled.weight <- FALSE

  answer <-
  switch(type,
    working = if (pooled.weight) NULL else object@residuals,
    pearson = {
      stop("have not programmed pearson resids yet")
    },
    deviance = {
      stop("have not programmed deviance resids yet")
    },
    ldot = {
      stop("have not programmed ldot resids yet")
    },
    response = {
      y <- object@y
      mu <- fitted(object)

      true.mu <- object@misc$true.mu
      if (is.null(true.mu))
        true.mu <- TRUE

      ans <- if (true.mu) y - mu else NULL


      if (!matrix.arg && length(ans)) {
        if (ncol(ans) == 1) {
          names.ans <- dimnames(ans)[[1]]
          ans <- c(ans)
          names(ans) <- names.ans
          ans
        } else {
          warning("ncol(ans) is not 1")
          ans
        }
      } else {
        ans
      }
  })

  if (length(answer) && length(na.act)) {
    napredict(na.act[[1]], answer)
  } else {
    answer
  }
}




setMethod("residuals",  "vlm",
          function(object, ...)
          residualsvlm(object, ...))
setMethod("residuals",  "vglm",
          function(object, ...)
          residualsvglm(object, ...))
setMethod("residuals",  "vgam",
          function(object, ...)
          residualsvglm(object, ...))
setMethod("residuals",  "qrrvglm",
          function(object, ...)
          residualsqrrvglm(object, ...))

setMethod("resid",  "vlm",
          function(object, ...)
          residualsvlm(object, ...))
setMethod("resid",  "vglm",
          function(object, ...)
          residualsvglm(object, ...))
setMethod("resid",  "vgam",
          function(object, ...)
          residualsvglm(object, ...))
setMethod("resid",  "qrrvglm",
          function(object, ...)
          residualsqrrvglm(object, ...))





