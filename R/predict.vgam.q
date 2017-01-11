# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.






predict.vgam <-
  function(object, newdata = NULL,
           type = c("link", "response", "terms"),
           se.fit = FALSE, deriv.arg = 0, terms.arg = NULL,
           raw = FALSE,
           all = TRUE, offset = 0,
           untransform = FALSE,
           dispersion = NULL, ...) {
  newdata <- if (missing(newdata)) {
    NULL
  } else {
    as.data.frame(newdata)
  }
  no.newdata <- (length(newdata) == 0)

  na.act <- object@na.action
  object@na.action <- list()

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("link", "response", "terms"))[1]


  if (untransform &&
   (type != "link" || se.fit || deriv.arg != 0 || offset != 0))
    stop("argument 'untransform = TRUE' only if type='link', ",
         "se.fit = FALSE, deriv = 0")

  if (raw && type!="terms")
    stop("'raw = TRUE' only works when 'type = \"terms\"'")

  if (!is.numeric(deriv.arg) || deriv.arg < 0 ||
     deriv.arg != round(deriv.arg) || length(deriv.arg) > 1)
    stop("bad input for the 'deriv' argument")

  if (deriv.arg > 0 && type != "terms")
    stop("'deriv>0' can only be specified if 'type=\"terms\"'")

  if (deriv.arg != 0 && !(type != "response" && !se.fit))
    stop("argument 'deriv' only works with type != 'response' and ",
         "se.fit = FALSE")

  if (se.fit && length(newdata))
    stop("cannot specify 'se.fit = TRUE' when argument 'newdata' ",
         "is assigned")


  tt <- terms(object)  # 20030811; object@terms$terms

  ttf <- attr(tt, "factors")
  tto <- attr(tt, "order")
  intercept <- attr(tt, "intercept")
  if (!intercept)
    stop("an intercept is assumed")

  M <- object@misc$M
  Hlist <- object@constraints
  ncolHlist <- unlist(lapply(Hlist, ncol))
  if (intercept)
    ncolHlist <- ncolHlist[-1]
  if (raw) {
    Hlist <- canonical.Hlist(Hlist)
    object@constraints <- Hlist
  }

  if (!length(newdata)) {
    if (type == "link") {
      if (se.fit) {
        stop("cannot handle this option (se.fit = TRUE) currently")
      } else {
        answer <- if (length(na.act)) {
          napredict(na.act[[1]], object@predictors)
        } else {
          object@predictors
        }
        if (untransform)
          return(untransformVGAM(object, answer)) else
          return(answer)
      }
    } else
    if (type == "response") {
      if (se.fit) {
        stop("cannot handle this option (se.fit = TRUE) currently")
      } else {
        if (length(na.act)) {
          return(napredict(na.act[[1]], object@fitted.values))
        } else {
          return(object@fitted.values)
        }
      }
    }

    predictor <- predict.vlm(object,
                             type = "terms",
                             se.fit = se.fit,
                             terms.arg = terms.arg,
                             raw = raw,
                             all = all, offset = offset,
                             dispersion = dispersion,
                             ...)  # deriv.arg = deriv.arg,

    newdata <- model.matrixvlm(object, type = "lm")


  } else {

    temp.type <- if (type == "link") "response" else type


    predictor <- predict.vlm(object, newdata,
                             type = temp.type,
                             se.fit = se.fit,
                             terms.arg = terms.arg,
                             raw = raw,
                             all = all, offset = offset,
                             dispersion = dispersion,
                             ...)  # deriv.arg = deriv.arg,
  }


  if (deriv.arg > 0)
    if (se.fit) {
      predictor$fitted.values <- predictor$fitted.values * 0
      predictor$se.fit <- predictor$se.fit * NA
    } else {
      predictor <- predictor * 0
    }


  if (length(s.xargument <- object@s.xargument)) {




    dnames2 <- dimnames(newdata)[[2]]
    index1 <- match(s.xargument, dnames2, nomatch = FALSE)
    index2 <- match(names(s.xargument), dnames2, nomatch = FALSE)
    index <- index1 | index2
    if (!length(index) || any(!index))
      stop("required variables not found in newdata")




    if (is.null(tmp6 <- attr(if (se.fit) predictor$fitted.values else
                            predictor, "vterm.assign"))) {

      Hlist <- subconstraints(object@misc$orig.assign,
                              object@constraints)
      ncolHlist <- unlist(lapply(Hlist, ncol))
      if (intercept)
        ncolHlist <- ncolHlist[-1]

      cs <- if (raw) cumsum(c(1, ncolHlist)) else
                     cumsum(c(1, M + 0 * ncolHlist))
      tmp6 <- vector("list", length(ncolHlist))
      for (ii in seq_along(tmp6))
        tmp6[[ii]] <- cs[ii]:(cs[ii+1]-1)
      names(tmp6) <- names(ncolHlist)
    }

    n.s.xargument <- names(s.xargument)  # e.g., c("s(x)", "s(x2)")
      for (ii in n.s.xargument) {

        fred <- s.xargument[ii]
        if (!any(dimnames(newdata)[[2]] == fred))
          fred <- ii

        xx <- newdata[, fred]  # [, s.xargument[ii]]  # [, nindex[ii]]

        rawMat <- predictvsmooth.spline.fit(object@Bspline[[ii]],
                                            x = xx,
                                            deriv = deriv.arg)$y


        eta.mat <- if (raw) rawMat else (rawMat %*% t(Hlist[[ii]]))

        if (type == "terms") {
          hhh <- tmp6[[ii]]
          if (se.fit) {
            predictor$fitted.values[, hhh] <-
            predictor$fitted.values[, hhh] + eta.mat

            TS <- predictor$sigma^2

            temp.var <- if (raw) {
                          tmp7 <- object@misc$varassign
                          tmp7 <- tmp7[[ii]]
                          object@var[, tmp7, drop = FALSE]
                        } else {
                          stop("cannot handle se's with raw = FALSE")
                        }

                        predictor$se.fit[, hhh] <-
                       (predictor$se.fit[, hhh]^2 + TS * temp.var)^0.5
                } else {
                  predictor[, hhh] <- predictor[, hhh] + eta.mat
                }
        } else {
          if (se.fit) {
            predictor$fitted.values <- predictor$fitted.values + eta.mat

            TS <- 1  # out$residual.scale^2
            TS <- predictor$sigma^2

            TT <- ncol(object@var)
            predictor$se.fit <- sqrt(predictor$se.fit^2 +
                                     TS * object@var %*% rep_len(1, TT))
          } else {
            predictor <- predictor + eta.mat
          }
        }
      }
  }

  if (type == "link") {
    if (no.newdata && length(na.act)) {
      return(napredict(na.act[[1]], predictor))
    } else {
      return(predictor)
    }
  } else
  if (type == "response") {
    fv <- object@family@linkinv(if (se.fit) predictor$fitted.values else
                                predictor, object@extra)
    if (is.matrix(fv) && is.matrix(object@fitted.values))
      dimnames(fv) <- list(dimnames(fv)[[1]],
                           dimnames(object@fitted.values)[[2]])
    if (is.matrix(fv) && ncol(fv) == 1)
      fv <- c(fv)
    if (no.newdata && length(na.act)) {
      fv <- if (se.fit) {
        napredict(na.act[[1]], fv)
      } else {
        napredict(na.act[[1]], fv)
      }
    }
    if (se.fit) {
      return(list(fit = fv, se.fit = fv * NA))
    } else {
      return(fv)
    }
  } else {
    if (deriv.arg >= 1) {
      if (se.fit) {
        attr(predictor$fitted.values, "constant") <- NULL
      } else {
        attr(predictor, "constant") <- NULL
      }
    }


    if (deriv.arg >= 1) {
      v <- attr(if (se.fit) predictor$fitted.values else
          predictor, "vterm.assign")
      is.lin <- is.linear.term(names(v))
        coefmat <- coefvlm(object, matrix.out = TRUE)
        ord <- 0
        for (ii in names(v)) {
          ord <- ord + 1
          index <- v[[ii]]
          lindex <- length(index)
          if (is.lin[ii]) {
            if (tto[ord] > 1 || (length(ttf) && ttf[ii, ii])) {
              if (se.fit) {
                predictor$fitted.values[, index] <-
                  if (tto[ord] > 1) NA else NA
              } else {
                predictor[, index] <- if (tto[ord] > 1) NA else NA
              }
            } else {
              ans <- coefmat[ii, 1:lindex]
                if (se.fit) {
                  predictor$fitted.values[, index] <-
                      if (deriv.arg == 1)
                      matrix(ans, ncol = lindex, byrow = TRUE) else 0
                } else {
                  predictor[, index] <- if (deriv.arg == 1)
                      matrix(ans, ncol = lindex, byrow = TRUE) else 0
                }
            }
          } else
            if (length(s.xargument) && any(n.s.xargument == ii)) {
              ans <- coefmat[ii, 1:lindex]
              if (se.fit) {
                predictor$fitted.values[, index] <-
                predictor$fitted.values[, index] +
                     (if (deriv.arg == 1)
                      matrix(ans, nrow = nrow(predictor$fitted.values),
                             ncol = lindex, byrow = TRUE) else 0)
              } else {
                predictor[, index] <- predictor[, index] +
                     (if (deriv.arg == 1)
                      matrix(ans, nrow = nrow(predictor),
                       ncol = lindex, byrow = TRUE) else 0)
              }
          } else {
              cat("Derivatives of term ", ii, "are unknown\n")
              if (se.fit) {
                predictor$fitted.values[, index] <- NA
              } else {
                predictor[, index] <- NA
              }
          }
        }
    }

    if (no.newdata && length(na.act)) {
      if (se.fit) {
        predictor$fitted.values <- napredict(na.act[[1]],
                                             predictor$fitted.values)
        predictor$se.fit <- napredict(na.act[[1]], predictor$se.fit)
      } else {
        predictor <- napredict(na.act[[1]], predictor)
      }
    }

    if (se.fit) {
      attr(predictor$fitted.values, "derivative") <- deriv.arg
    } else {
      attr(predictor, "derivative") <- deriv.arg
    }

    return(predictor)
  }
}


setMethod("predict", "vgam",
          function(object, ...)
          predict.vgam(object, ...))



varassign <- function(constraints, n.s.xargument) {

  if (!length(n.s.xargument))
    stop("length(n.s.xargument) must be > 0")

  ans <- vector("list", length(n.s.xargument))

  ncolHlist <- unlist(lapply(constraints, ncol))

  names(ans) <- n.s.xargument
  ptr <- 1
  for (ii in n.s.xargument) {
    temp <- ncolHlist[[ii]]
    ans[[ii]] <- ptr:(ptr + temp - 1)
    ptr <- ptr + temp
  }
  ans
}




