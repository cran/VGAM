# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.







predict.vlm <- function(object,
                        newdata = NULL,
                        type = c("response", "terms"),
                        se.fit = FALSE, scale = NULL,
                        terms.arg = NULL,
                        raw = FALSE,
                        dispersion = NULL, ...) {
  Xm2 <- NULL
  xij.used <- length(form2 <- object@misc$form2) ||
              length(object@control$xij)

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("response", "terms"))[1]

  na.act <- object@na.action
  object@na.action <- list()

  if (raw && type != "terms")
    stop("sorry, 'raw=TRUE' only works when 'type=\"terms\"'")

  if (!length(newdata) && type == "response" && !se.fit &&
    length(object@fitted.values)) {
    if (length(na.act)) {
      return(napredict(na.act[[1]], object@fitted.values))
    } else {
      return(object@fitted.values)
    }
  }

  ttob <- terms(object)  # 20030811; object@terms$terms


  if (!length(newdata)) {
    offset <- object@offset

    if (xij.used) {
      bothList <- model.matrix(object, type = "bothlmlm2")
      X   <- bothList$X
      Xm2 <- bothList$Xm2
    } else {
      X <- model.matrix(object, type = "lm")
    }
  } else {

    if (is.smart(object) && length(object@smart.prediction)) {
      setup.smart("read", smart.prediction = object@smart.prediction)
    }

    X <- model.matrix(delete.response(ttob), newdata,
                      contrasts = if (length(object@contrasts))
                                  object@contrasts else NULL,
                      xlev = object@xlevels)
    if (xij.used) {
      ttXm2 <- terms(form2)
      Xm2 <- model.matrix(delete.response(ttXm2), newdata,
                          contrasts = if (length(object@contrasts))
                                      object@contrasts else NULL,
                          xlev = object@xlevels)
    }

    if (object@misc$intercept.only &&
        nrow(X) != nrow(newdata)) {
      as.save <- attr(X, "assign")
      X <- X[rep_len(1, nrow(newdata)), , drop = FALSE]
      dimnames(X) <- list(dimnames(newdata)[[1]], "(Intercept)")
      attr(X, "assign") <- as.save  # Restored
    }

    offset <- if (!is.null(off.num <- attr(ttob, "offset"))) {
      eval(attr(ttob, "variables")[[off.num + 1]], newdata)
    } else if (!is.null(object@offset))
      eval(object@call$offset, newdata)

    if (is.smart(object) && length(object@smart.prediction)) {
      wrapup.smart()
    }

    attr(X, "assign") <- attrassigndefault(X, ttob)
    if (length(Xm2))
      attr(Xm2, "assign") <- attrassigndefault(Xm2, ttXm2)
  }


  hasintercept <- attr(ttob, "intercept")

  dx1 <- dimnames(X)[[1]]
  M <- object@misc$M
  Hlist <- object@constraints
  ncolHlist <- unlist(lapply(Hlist, ncol))
  if (hasintercept)
    ncolHlist <- ncolHlist[-1]

  xbar <- x2bar <- NULL
  if (type == "terms" && hasintercept) {
    if (length(object@control$xij)) {
      x2bar <- colMeans(Xm2)
      Xm2 <- sweep(Xm2, 2, x2bar)
    }
    xbar <- colMeans(X)
    X <- sweep(X, 2, xbar)
    nac <- is.na(object@coefficients)
    if (any(nac)) {
      if (length(object@control$xij))
        stop("cannot handle 'xij' argument when ",
             "there are NAs in the coefficients")
      X <- X[, !nac, drop = FALSE]
      xbar <- xbar[!nac]
    }
  }

    if (!is.null(newdata) && !is.data.frame(newdata))
        newdata <- as.data.frame(newdata)

    nn <- if (!is.null(newdata)) nrow(newdata) else object@misc$n
    if (raw) {
      Hlist <- canonical.Hlist(Hlist)
      object@constraints <- Hlist
    }



    X_vlm <- lm2vlm.model.matrix(X, Hlist = Hlist, M = M,
                                 xij = object@control$xij, Xm2 = Xm2)


    attr(X_vlm, "constant")  <- xbar
    attr(X_vlm, "constant2") <- x2bar




    coefs <- coefvlm(object)
    vasgn <- attr(X_vlm, "vassign")


    if (type == "terms") {
      nv <- names(vasgn)
      if (hasintercept)
        nv <- nv[-(1:ncol(object@constraints[["(Intercept)"]]))]
      terms.arg <- if (is.null(terms.arg)) nv else terms.arg

      index <- charmatch(terms.arg, nv)
      if (all(index == 0)) {
        warning("no match found; returning all terms")
        index <- seq_along(nv)
      }
      vasgn <- vasgn[nv[index]]
    }

    if (anyNA(object@coefficients))
      stop("cannot handle NAs in 'object@coefficients'")

    dname2 <- object@misc$predictors.names
    if (se.fit) {
      object <- as(object, "vlm")  # Coerce
      fit.summary <- summaryvlm(object, dispersion = dispersion)
      sigma <- if (is.numeric(fit.summary@sigma))
        fit.summary@sigma else
        sqrt(deviance(object) / object@df.residual)  # was @ResSS
      pred <- Build.terms.vlm(x = X_vlm, coefs = coefs,
                              cov = sigma^2 * fit.summary@cov.unscaled,
                              assign = vasgn,
                              collapse = type != "terms", M = M,
                              dimname = list(dx1, dname2),
                              coefmat = coefvlm(object, matrix.out = TRUE))
      pred$df <- object@df.residual
      pred$sigma <- sigma
    } else {
      pred <- Build.terms.vlm(x = X_vlm, coefs = coefs,
                              cov = NULL,
                              assign = vasgn,
                              collapse = type != "terms", M = M,
                              dimname = list(dx1, dname2),
                              coefmat = coefvlm(object, matrix.out = TRUE))
    }

    constant  <- attr(pred, "constant")

  if (type != "terms" && length(offset) && any(offset != 0)) {
    if (se.fit) {
      pred$fitted.values <- pred$fitted.values + offset
    } else {
      pred <- pred + offset
    }
  }



  if (type == "terms") {
    Hlist <- subconstraints(object@misc$orig.assign, object@constraints)
    ncolHlist <- unlist(lapply(Hlist, ncol))
    if (hasintercept)
      ncolHlist <- ncolHlist[-1]

    cs <- cumsum(c(1, ncolHlist))  # Like a pointer
    for (ii in 1:(length(cs)-1))
      if (cs[ii+1] - cs[ii] > 1)
        for (kk in (cs[ii]+1):(cs[ii+1]-1))
          if (se.fit) {
            pred$fitted.values[, cs[ii]] <- pred$fitted.values[, cs[ii]] +
                                            pred$fitted.values[, kk]
            pred$se.fit[, cs[ii]] <- pred$se.fit[, cs[ii]] +
                                     pred$se.fit[, kk]
          } else {
            pred[, cs[ii]] <- pred[, cs[ii]] + pred[, kk]
          }

        if (se.fit) {
          pred$fitted.values <-
          pred$fitted.values[, cs[-length(cs)], drop = FALSE]
          pred$se.fit <- pred$se.fit[, cs[-length(cs)], drop = FALSE]
        } else {
          pred <- pred[, cs[-length(cs)], drop = FALSE]
        }

        pp <- if (se.fit) ncol(pred$fitted.values) else ncol(pred)
        if (se.fit) {
          dimnames(pred$fitted.values) <- dimnames(pred$se.fit) <- NULL
          dim(pred$fitted.values) <- dim(pred$se.fit) <- c(M, nn, pp)
          pred$fitted.values <- aperm(pred$fitted.values, c(2, 1, 3))
          pred$se.fit <- aperm(pred$se.fit, c(2, 1, 3))
          dim(pred$fitted.values) <- dim(pred$se.fit) <- c(nn, M*pp)
        } else {
          dimnames(pred) <- NULL  # Saves a warning
          dim(pred) <- c(M, nn, pp)
          pred <- aperm(pred, c(2, 1, 3))
          dim(pred) <- c(nn, M*pp)
        }

      if (raw) {
        kindex <- NULL
        for (ii in 1:pp)
          kindex <- c(kindex, (ii-1) * M + (1:ncolHlist[ii]))
        if (se.fit) {
          pred$fitted.values <- pred$fitted.values[, kindex, drop = FALSE]
          pred$se.fit <- pred$se.fit[, kindex, drop = FALSE]
        } else {
          pred <- pred[, kindex, drop = FALSE]
        }
      }

      temp <- if (raw) ncolHlist else rep_len(M, length(ncolHlist))
      dd <- vlabel(names(ncolHlist), temp, M)
      if (se.fit) {
        dimnames(pred$fitted.values) <-
        dimnames(pred$se.fit) <-
          list(if (length(newdata)) dimnames(newdata)[[1]] else dx1,
               dd)
      } else {
        dimnames(pred) <-
          list(if (length(newdata)) dimnames(newdata)[[1]] else dx1,
               dd)
      }

      if (!length(newdata) && length(na.act)) {
        if (se.fit) {
          pred$fitted.values <- napredict(na.act[[1]], pred$fitted.values)
          pred$se.fit <- napredict(na.act[[1]], pred$se.fit)
        } else {
          pred <- napredict(na.act[[1]], pred)
        }
      }

    if (!raw)
      cs <- cumsum(c(1, M + 0 * ncolHlist))
    fred <- vector("list", length(ncolHlist))
    for (ii in seq_along(fred))
      fred[[ii]] <- cs[ii]:(cs[ii+1]-1)
    names(fred) <- names(ncolHlist)
    if (se.fit) {
      attr(pred$fitted.values, "vterm.assign") <-
      attr(pred$se.fit,        "vterm.assign") <- fred
    } else {
      attr(pred,               "vterm.assign") <- fred
    }
  }  # End of if (type == "terms")

  if (!is.null(xbar)) {
    if (se.fit) {
      attr(pred$fitted.values, "constant") <- constant
    } else {
      attr(pred,               "constant") <- constant
    }
  }

  pred
}  # predict.vlm()





setMethod("predict", "vlm",
          function(object, ...)
          predict.vlm(object, ...))







predict.vglm.se <- function(fit, ...) {


  H.ss <- hatvalues(fit, type = "centralBlocks")  # diag = FALSE

  M <- npred(fit)
  nn <- nobs(fit, type = "lm")
  U <- vchol(weights(fit, type = "working"), M = M, n = nn)

  Uarray <- array(0, c(M, M, nn))
  ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
  MMp1d2 <- M * (M + 1) / 2
  for (jay in 1:MMp1d2)
    Uarray[ind1$row.index[jay],
           ind1$col.index[jay], ] <- U[jay, ]

  Uinv.array <- apply(Uarray, 3, backsolve, x = diag(M))
  dim(Uinv.array) <- c(M, M, nn)

  Utinv.array <- Uinv.array
  if (M > 1)
    for (jay in 1:(M-1)) {
      for (kay in (jay+1):M) {
        Utinv.array[kay, jay, ] <- Uinv.array[jay, kay, ]
        Utinv.array[jay, kay, ] <- 0
      }
    }

  var.boldeta.i <- mux5(H.ss, Utinv.array, M = M,
                        matrix.arg = TRUE)  # First M cols are SE^2

  sqrt(var.boldeta.i[, 1:M])  # SE(linear.predictor)




  sqrt(var.boldeta.i[, 1:M])
}










subconstraints <- function(assign, constraints) {


  ans <- vector("list", length(assign))
  if (!length(assign) || !length(constraints))
    stop("assign and/or constraints is empty")
  for (ii in seq_along(assign))
    ans[[ii]] <- constraints[[assign[[ii]][1]]]
  names(ans) <- names(assign)
  ans
}



is.linear.term <- function(ch) {
  lchar <- length(ch)
  ans <- rep_len(FALSE, lchar)
  for (ii in 1:lchar) {
    nc <- nchar(ch[ii])
    x <- substring(ch[ii], 1:nc, 1:nc)
    ans[ii] <- all(x != "(" & x != "+" & x != "-" &
                   x != "/" & x != "*" & x != "^")
  }
  names(ans) <- ch
  ans
}



canonical.Hlist <- function(Hlist) {
  for (ii in seq_along(Hlist)) {
    temp <- Hlist[[ii]] * 0
    temp[cbind(1:ncol(temp), 1:ncol(temp))] <- 1
    Hlist[[ii]] <- temp
  }
  Hlist
}



