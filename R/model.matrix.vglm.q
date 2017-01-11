# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.











attrassigndefault <- function(mmat, tt) {
  if (!inherits(tt, "terms"))
    stop("need terms object")
  aa <- attr(mmat, "assign")
  if (is.null(aa))
    stop("argument is not really a model matrix")
  ll <- attr(tt, "term.labels")
  if (attr(tt, "intercept") > 0)
    ll <- c("(Intercept)", ll)
  aaa <- factor(aa, labels = ll)
  split(order(aa), aaa)
}


attrassignlm <- function(object, ...)
  attrassigndefault(model.matrix(object), object@terms)





 vlabel <- function(xn, ncolHlist, M, separator = ":", colon = FALSE) {

  if (length(xn) != length(ncolHlist))
    stop("length of first two arguments not equal")

  n1 <- rep(xn, ncolHlist)
  if (M == 1)
    return(n1)
  n2 <- as.list(ncolHlist)
  n2 <- lapply(n2, seq)
  n2 <- unlist(n2)
  n2 <- as.character(n2)
  n2 <- paste(separator, n2, sep = "")
  n3 <- rep(ncolHlist, ncolHlist)
  if (!colon)
    n2[n3 == 1] <- ""
  n1n2 <- paste(n1, n2, sep = "")
  n1n2
}






 vlm2lm.model.matrix <-
  function(x.vlm, Hlist = NULL,
           which.linpred = 1,
           M = NULL) {







  if (is.numeric(M)) {
    if (M != nrow(Hlist[[1]]))
      stop("argument 'M' does not match argument 'Hlist'")
  } else {
    M <- nrow(Hlist[[1]])
  }


  Hmatrices <- matrix(c(unlist(Hlist)), nrow = M)
  if (ncol(Hmatrices) != ncol(x.vlm))
    stop("ncol(Hmatrices) != ncol(x.vlm)")


  n.lm <- nrow(x.vlm) / M
  if (round(n.lm) != n.lm)
    stop("'n.lm' does not seem to be an integer")
  linpred.index <- which.linpred
  vecTF <- Hmatrices[linpred.index, ] != 0
  X.lm.jay <- x.vlm[(0:(n.lm - 1)) * M + linpred.index, vecTF,
                    drop = FALSE]
  X.lm.jay
}







 lm2vlm.model.matrix <-
  function(x, Hlist = NULL, assign.attributes = TRUE,
           M = NULL, xij = NULL, Xm2 = NULL) {




  if (length(Hlist) != ncol(x))
    stop("length(Hlist) != ncol(x)")

  if (length(xij)) {
    if (inherits(xij, "formula"))
      xij <- list(xij)
    if (!is.list(xij))
      stop("'xij' is not a list of formulae")
  }

  if (!is.numeric(M))
    M <- nrow(Hlist[[1]])

  nrow.X.lm <- nrow(x)
  if (all(trivial.constraints(Hlist) == 1)) {
    X.vlm <- if (M > 1) kronecker(x, diag(M)) else x
    ncolHlist <- rep(M, ncol(x))
  } else {
    allB <- matrix(unlist(Hlist), nrow = M)
    ncolHlist <- unlist(lapply(Hlist, ncol))
    Rsum <- sum(ncolHlist)

    X1 <- rep(c(t(x)), rep(ncolHlist, nrow.X.lm))
    dim(X1) <- c(Rsum, nrow.X.lm)
    X.vlm <- kronecker(t(X1), matrix(1, M, 1)) *
             kronecker(matrix(1, nrow.X.lm, 1), allB)
    rm(X1)
  }

  dn <- labels(x)
  yn <- dn[[1]]
  xn <- dn[[2]]
  dimnames(X.vlm) <- list(vlabel(yn, rep(M, nrow.X.lm), M),
                          vlabel(xn, ncolHlist, M))

  if (assign.attributes) {
      attr(X.vlm, "contrasts")   <- attr(x, "contrasts")
      attr(X.vlm, "factors")     <- attr(x, "factors")
      attr(X.vlm, "formula")     <- attr(x, "formula")
      attr(X.vlm, "class")       <- attr(x, "class")
      attr(X.vlm, "order")       <- attr(x, "order")
      attr(X.vlm, "term.labels") <- attr(x, "term.labels")

      nasgn <- oasgn <- attr(x, "assign")
      lowind <- 0
      for (ii in seq_along(oasgn)) {
          mylen <- length(oasgn[[ii]]) * ncolHlist[oasgn[[ii]][1]]
          nasgn[[ii]] <- (lowind+1):(lowind+mylen)
          lowind <- lowind + mylen
      } # End of ii
      if (lowind != ncol(X.vlm))
        stop("something gone wrong")
      attr(X.vlm, "assign") <- nasgn


      fred <- unlist(lapply(nasgn, length)) / unlist(lapply(oasgn, length))
      vasgn <- vector("list", sum(fred))
      kk <- 0
      for (ii in seq_along(oasgn)) {
        temp <- matrix(nasgn[[ii]], ncol = length(oasgn[[ii]]))
        for (jloc in 1:nrow(temp)) {
          kk <- kk + 1
          vasgn[[kk]] <- temp[jloc, ]
        }
      }
      names(vasgn) <- vlabel(names(oasgn), fred, M)
      attr(X.vlm, "vassign") <- vasgn

      attr(X.vlm, "constraints") <- Hlist
  } # End of if (assign.attributes)




  if (!length(xij))
    return(X.vlm)







  at.x <- attr(x, "assign")
  at.vlmx <- attr(X.vlm, "assign")
  at.Xm2 <- attr(Xm2, "assign")

  for (ii in seq_along(xij)) {
    form.xij <- xij[[ii]]
    if (length(form.xij) != 3)
      stop("xij[[", ii, "]] is not a formula with a response")
    tform.xij <- terms(form.xij)
    aterm.form <- attr(tform.xij, "term.labels")  # Does not include response
    if (length(aterm.form) != M)
      stop("xij[[", ii, "]] does not contain ", M, " terms")

    name.term.y <- as.character(form.xij)[2]
    cols.X.vlm <- at.vlmx[[name.term.y]]  # May be > 1 in length.

    x.name.term.2 <- aterm.form[1]   # Choose the first one
    One.such.term <- at.Xm2[[x.name.term.2]]
    for (bbb in seq_along(One.such.term)) {
      use.cols.Xm2 <- NULL
      for (sss in 1:M) {
        x.name.term.2 <- aterm.form[sss]
        one.such.term <- at.Xm2[[x.name.term.2]]
        use.cols.Xm2 <- c(use.cols.Xm2, one.such.term[bbb])
      } # End of sss

      allXk <- Xm2[, use.cols.Xm2, drop = FALSE]
      cmat.no <- (at.x[[name.term.y]])[1]  # 1st one will do (all the same).
      cmat <- Hlist[[cmat.no]]
      Rsum.k <- ncol(cmat)
      tmp44 <- kronecker(matrix(1, nrow.X.lm, 1), t(cmat)) *
               kronecker(allXk, matrix(1, ncol(cmat), 1))  # n*Rsum.k x M

      tmp44 <- array(t(tmp44), c(M, Rsum.k, nrow.X.lm))
      tmp44 <- aperm(tmp44, c(1, 3, 2))  # c(M, n, Rsum.k)
      rep.index <- cols.X.vlm[((bbb-1)*Rsum.k+1):(bbb*Rsum.k)]
      X.vlm[, rep.index] <- c(tmp44)
    }  # End of bbb
  }  # End of for (ii in seq_along(xij))

  if (assign.attributes) {
    attr(X.vlm, "vassign") <- vasgn
    attr(X.vlm, "assign") <- nasgn
    attr(X.vlm, "xij") <- xij
  }
  X.vlm
}  # lm2vlm.model.matrix









model.matrix.vlm <- function(object, ...)
  model.matrixvlm(object, ...)





 model.matrixvlm <- function(object,
                             type = c("vlm", "lm", "lm2", "bothlmlm2"),
                             linpred.index = NULL,
                            ...) {



  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("vlm", "lm", "lm2", "bothlmlm2"))[1]

  if (length(linpred.index) &&
      type != "lm")
    stop("Must set 'type = \"lm\"' when 'linpred.index' is ",
         "assigned a value")
  if (length(linpred.index) &&
      length(object@control$xij))
    stop("Currently cannot handle 'xij' models when 'linpred.index' is ",
         "assigned a value")


  x   <- slot(object, "x")


  Xm2 <- if (any(slotNames(object) == "Xm2")) slot(object, "Xm2") else
         numeric(0)


  form2 <- if (any(slotNames(object) == "misc")) object@misc$form2 else NULL
  if (type == "lm2" && !length(form2))
    return(Xm2)


  if (!length(x)) {
    data <- model.frame(object, xlev = object@xlevels, ...)

    kill.con <- if (length(object@contrasts)) object@contrasts else NULL

    x <- vmodel.matrix.default(object, data = data,
                               contrasts.arg = kill.con)
    tt <- terms(object)
    attr(x, "assign") <- attrassigndefault(x, tt)
  }




  if ((type == "lm2" || type == "bothlmlm2") &&
      !length(Xm2)) {
    object.copy2 <- object
    data <- model.frame(object.copy2, xlev = object.copy2@xlevels, ...)

    kill.con <- if (length(object.copy2@contrasts))
                object.copy2@contrasts else NULL

    Xm2 <- vmodel.matrix.default(object.copy2, data = data,
                                 contrasts.arg = kill.con)


    if (length(form2)) {
      attr(Xm2, "assign") <- attrassigndefault(Xm2, terms(form2))
    }
  }





  if (type == "lm" && is.null(linpred.index)) {
    return(x)
  } else if (type == "lm2") {
    return(Xm2)
  } else if (type == "bothlmlm2") {
    return(list(X = x, Xm2 = Xm2))
  }


  M <- object@misc$M
  Hlist <- object@constraints  # == constraints(object, type = "lm")
  X.vlm <- lm2vlm.model.matrix(x = x, Hlist = Hlist,
                               xij = object@control$xij, Xm2 = Xm2)

  if (type == "vlm") {
    return(X.vlm)
  } else if (type == "lm" && length(linpred.index)) {
    if (!is.Numeric(linpred.index, integer.valued = TRUE, positive = TRUE,
                    length.arg = 1))
      stop("bad input for argument 'linpred.index'")
    if (!length(intersect(linpred.index, 1:M)))
      stop("argument 'linpred.index' should have ",
           "a single value from the set 1:", M)

    Hlist <- Hlist
    n.lm <- nobs(object)  # Number of rows of the LM matrix
    M <- object@misc$M  # Number of linear/additive predictors
    Hmatrices <- matrix(c(unlist(Hlist)), nrow = M)
    jay <- linpred.index
    index0 <- Hmatrices[jay, ] != 0
    X.lm.jay <- X.vlm[(0:(n.lm - 1)) * M + jay, index0, drop = FALSE]
    X.lm.jay
  } else {
    stop("am confused. Do not know what to return")
  }
}




setMethod("model.matrix",  "vlm", function(object, ...)
           model.matrixvlm(object, ...))




 model.matrixvgam <-
  function(object,
           type = c("lm", "vlm", "lm", "lm2", "bothlmlm2"),
           linpred.index = NULL,
           ...) {
  model.matrixvlm(object = object,
                  type = type[1],
                  linpred.index = linpred.index, ...)
}
setMethod("model.matrix",  "vgam", function(object, ...)
           model.matrixvgam(object, ...))






 model.framevlm <- function(object,
                            setupsmart = TRUE,
                            wrapupsmart = TRUE, ...) {

  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || !length(object@model)) {
    fcall <- object@call
    fcall$method <- "model.frame"
    fcall[[1]] <- as.name("vlm")

    fcall$smart <- FALSE
    if (setupsmart && length(object@smart.prediction)) {
      setup.smart("read", smart.prediction=object@smart.prediction)
    }

    fcall[names(nargs)] <- nargs
    env <- environment(object@terms$terms)  # @terms or @terms$terms ??
    if (is.null(env))
      env <- parent.frame()
    ans <- eval(fcall, env, parent.frame())

    if (wrapupsmart && length(object@smart.prediction)) {
      wrapup.smart()
    }
    ans
  } else object@model
}


if (!isGeneric("model.frame"))
    setGeneric("model.frame", function(formula, ...)
        standardGeneric("model.frame"))

setMethod("model.frame",  "vlm", function(formula, ...)
           model.framevlm(object = formula, ...))





 vmodel.matrix.default <-
  function(object, data = environment(object),
           contrasts.arg = NULL, xlev = NULL, ...) {

  t <- if (missing(data)) terms(object) else terms(object, data = data)
  if (is.null(attr(data, "terms")))
    data <- model.frame(object, data, xlev = xlev) else {
    reorder <- match(sapply(attr(t, "variables"), deparse,
                     width.cutoff = 500)[-1], names(data))
    if (anyNA(reorder))
      stop("model frame and formula mismatch in model.matrix()")
    if (!identical(reorder, seq_len(ncol(data))))
      data <- data[, reorder, drop = FALSE]
  }
  int <- attr(t, "response")
  if (length(data)) {
    contr.funs <- as.character(getOption("contrasts"))
    namD <- names(data)
    for (i in namD) if (is.character(data[[i]])) {
      data[[i]] <- factor(data[[i]])
      warning(gettextf("variable '%s' converted to a factor", i),
              domain = NA)
    }
    isF <- sapply(data, function(x) is.factor(x) || is.logical(x))
    isF[int] <- FALSE
    isOF <- sapply(data, is.ordered)
    for (nn in namD[isF]) if (is.null(attr(data[[nn]], "contrasts")))
      contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
    if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
      if (is.null(namC <- names(contrasts.arg)))
        stop("invalid 'contrasts.arg' argument")
      for (nn in namC) {
        if (is.na(ni <- match(nn, namD)))
          warning(gettextf(
            "variable '%s' is absent, its contrast will be ignored",
            nn), domain = NA) else {
          ca <- contrasts.arg[[nn]]
          if (is.matrix(ca))
            contrasts(data[[ni]], ncol(ca)) <- ca else
            contrasts(data[[ni]]) <- contrasts.arg[[nn]]
        }
      }
    }
  } else {
    isF <- FALSE
    data <- list(x = rep(0, nrow(data)))
  }


  ans  <-          (model.matrix(t, data))




  cons <- if (any(isF))
    lapply(data[isF], function(x) attr(x, "contrasts")) else NULL
  attr(ans, "contrasts") <- cons
  ans
}





depvar.vlm <-
  function(object,
           type = c("lm", "lm2"),
           drop = FALSE,
           ...) {
  type <- match.arg(type, c("lm", "lm2"))[1]
  ans <- if (type == "lm") {
    object@y
  } else {
    object@Ym2
  }
  ans[, , drop = drop]
}



if (!isGeneric("depvar"))
    setGeneric("depvar",
               function(object, ...)
                 standardGeneric("depvar"),
               package = "VGAM")


setMethod("depvar",  "vlm", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "rrvglm", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "qrrvglm", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "rrvgam", function(object, ...)
           depvar.vlm(object, ...))
setMethod("depvar",  "rcim", function(object, ...)
           depvar.vlm(object, ...))





npred.vlm <- function(object,
                      type = c("total", "one.response"),
                      ...) {
  if (!missing(type))
    type <- as.character(substitute(type))
  type.arg <- match.arg(type, c("total", "one.response"))[1]


  MM <-
    if (length(object@misc$M))
      object@misc$M else
    if (ncol(as.matrix(predict(object))) > 0)
      ncol(as.matrix(predict(object))) else
    stop("cannot seem to obtain 'M'")


  if (type.arg == "one.response") {
    M1.infos <- NULL
    infos.fun <- object@family@infos
    Ans.infos <- infos.fun()
    if (is.list(Ans.infos) && length(Ans.infos$M1))
      M1.infos <- Ans.infos$M1

    Q1 <- Ans.infos$Q1
    if (is.numeric(Q1)) {
      S <- ncol(depvar(object)) / Q1  # Number of (multiple) responses
      if (is.numeric(M1.infos) && M1.infos * S != MM)
        warning("contradiction in values after computing it two ways")
    }


    M1 <- if (is.numeric(M1.infos)) M1.infos else
          if (is.numeric(MM      )) MM       else
          stop("failed to compute 'M'")
    M1
  } else {  # One response is assumed, by default
    MM
  }
}


if (!isGeneric("npred"))
    setGeneric("npred", function(object, ...) standardGeneric("npred"),
               package = "VGAM")


setMethod("npred",  "vlm", function(object, ...)
           npred.vlm(object, ...))
setMethod("npred",  "rrvglm", function(object, ...)
           npred.vlm(object, ...))
setMethod("npred",  "qrrvglm", function(object, ...)
           npred.vlm(object, ...))
setMethod("npred",  "rrvgam", function(object, ...)
           npred.vlm(object, ...))
setMethod("npred",  "rcim", function(object, ...)
           npred.vlm(object, ...))






hatvaluesvlm <-
  function(model,
           type = c("diagonal", "matrix", "centralBlocks"), ...) {


  if (!missing(type))
    type <- as.character(substitute(type))
  type.arg <- match.arg(type, c("diagonal", "matrix", "centralBlocks"))[1]


  qrSlot <- model@qr

  if (!is.list(qrSlot) && class(qrSlot) != "qr")
    stop("slot 'qr' should be a list")

  M  <- npred(model)
  nn <- nobs(model, type = "lm")

  if (is.empty.list(qrSlot)) {

    wzedd <- weights(model, type = "working")
    UU <- vchol(wzedd, M = M, n = nn, silent = TRUE)  # Few rows, many cols
    X.vlm <- model.matrix(model, type = "vlm")
    UU.X.vlm <- mux111(cc = UU, xmat = X.vlm, M = M)
    qrSlot <- qr(UU.X.vlm)
  } else {
    X.vlm <- NULL
    class(qrSlot) <- "qr" # S3 class
  }
  Q.S3 <- qr.Q(qrSlot)



  if (type.arg == "diagonal") {
    Diag.Hat <- rowSums(Q.S3^2)
    Diag.Elts <- matrix(Diag.Hat, nn, M, byrow = TRUE)

    if (length(model@misc$predictors.names) == M)
      colnames(Diag.Elts) <- model@misc$predictors.names
    if (length(rownames(model.matrix(model, type = "lm"))))
      rownames(Diag.Elts) <- rownames(model.matrix(model, type = "lm"))

    attr(Diag.Elts, "predictors.names") <- model@misc$predictors.names
    attr(Diag.Elts, "ncol.X.vlm") <- model@misc$ncol.X.vlm

    Diag.Elts
  } else if (type.arg == "matrix") {
    all.mat <- Q.S3 %*% t(Q.S3)
    if (!length(X.vlm))
      X.vlm <- model.matrix(model, type = "vlm")
    dimnames(all.mat) <- list(rownames(X.vlm), rownames(X.vlm))

    attr(all.mat, "M") <- M
    attr(all.mat, "predictors.names") <- model@misc$predictors.names
    attr(all.mat, "ncol.X.vlm") <- model@misc$ncol.X.vlm

    all.mat
  } else {
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    MMp1d2 <- M * (M + 1) / 2
    all.rows.index <- rep((0:(nn-1)) * M, rep(MMp1d2, nn)) + ind1$row.index
    all.cols.index <- rep((0:(nn-1)) * M, rep(MMp1d2, nn)) + ind1$col.index

    H.ss <- rowSums(Q.S3[all.rows.index, ] *
                    Q.S3[all.cols.index, ])

    H.ss <- matrix(H.ss, nn, MMp1d2, byrow = TRUE)
    H.ss
  }
}  # hatvaluesvlm





setMethod("hatvalues",  "vlm", function(model, ...)
           hatvaluesvlm(model, ...))
setMethod("hatvalues",  "vglm", function(model, ...)
           hatvaluesvlm(model, ...))
setMethod("hatvalues",  "rrvglm", function(model, ...)
           hatvaluesvlm(model, ...))
setMethod("hatvalues",  "qrrvglm", function(model, ...)
           hatvaluesvlm(model, ...))
setMethod("hatvalues",  "rrvgam", function(model, ...)
           hatvaluesvlm(model, ...))
setMethod("hatvalues",  "rcim", function(model, ...)
           hatvaluesvlm(model, ...))








hatplot.vlm <-
  function(model, multiplier = c(2, 3),
           lty = "dashed",
           xlab = "Observation",
           ylab = "Hat values",
           ylim = NULL, ...) {

  if (is(model, "vlm")) {
    hatval <- hatvalues(model, diag = TRUE)
  } else {
    hatval <- model
  }

  if (!is.matrix(hatval))
    stop("argument 'model' seems neither a vglm() object or a matrix")

  ncol.X.vlm <- attr(hatval, "ncol.X.vlm")
  M <- attr(hatval, "M")
  predictors.names <- attr(hatval, "predictors.names")
  if (!length(predictors.names)) {
    predictors.names <- paste("Linear/additive predictor", 1:M)
  }

  if (length(M)) {
    N <- nrow(hatval) / M
    hatval <- matrix(hatval, N, M, byrow = TRUE)
  } else {
    M <- ncol(hatval)
    N <- nrow(hatval)
  }

  if (is.null(ylim))
    ylim <- c(0, max(hatval))
  for (jay in 1:M) {
    plot(hatval[, jay], type = "n", main = predictors.names[jay],
         ylim = ylim, xlab = xlab, ylab = ylab,
         ...)
    points(1:N, hatval[, jay], ...)
    abline(h = multiplier * ncol.X.vlm / (N * M), lty = lty, ...)
  }
}




if (!isGeneric("hatplot"))
    setGeneric("hatplot", function(model, ...)
      standardGeneric("hatplot"), package = "VGAM")


setMethod("hatplot",  "matrix", function(model, ...)
           hatplot.vlm(model, ...))

setMethod("hatplot",  "vlm", function(model, ...)
           hatplot.vlm(model, ...))
setMethod("hatplot",  "vglm", function(model, ...)
           hatplot.vlm(model, ...))

setMethod("hatplot",  "rrvglm", function(model, ...)
           hatplot.vlm(model, ...))
setMethod("hatplot",  "qrrvglm", function(model, ...)
           hatplot.vlm(model, ...))
setMethod("hatplot",  "rrvgam", function(model, ...)
           hatplot.vlm(model, ...))
setMethod("hatplot",  "rcim", function(model, ...)
           hatplot.vlm(model, ...))








dfbetavlm <-
  function(model,
           maxit.new = 1,
           trace.new = FALSE,
           smallno = 1.0e-8,
           ...) {

  if (!is(model, "vlm"))
    stop("argument 'model' does not seem to be a vglm() object")

  n.lm <- nobs(model, type = "lm")
  X.lm <- model.matrix(model, type = "lm")
  X.vlm <- model.matrix(model, type = "vlm")
  p.vlm <- ncol(X.vlm)  # nvar(model, type = "vlm")
  M    <- npred(model)
  etastart <- predict(model)
  offset <- matrix(model@offset, n.lm, M)
  new.control <- model@control
  pweights <- weights(model, type = "prior")
  orig.w <- if (is.numeric(model@extra$orig.w))
              model@extra$orig.w else 1
  y.integer <- if (is.logical(model@extra$y.integer))
                 model@extra$y.integer else FALSE
  coef.model <- coef(model)


  new.control$trace <- trace.new
  new.control$maxit <- maxit.new

  dfbeta <- matrix(0, n.lm, p.vlm)

  Terms.zz <- NULL





  for (ii in 1:n.lm) {
    if (trace.new) {
      cat("\n", "Observation ", ii, "\n")
      flush.console()
    }

    w.orig <- if (length(orig.w) != n.lm)
                rep_len(orig.w, n.lm) else
                orig.w
    w.orig[ii] <- w.orig[ii] * smallno  # Relative

    fit <- vglm.fit(x = X.lm,
                    X.vlm.arg = X.vlm,  # Should be more efficient
                    y = if (y.integer)
                      round(depvar(model) * c(pweights) / c(orig.w)) else
                           (depvar(model) * c(pweights) / c(orig.w)),
                    w = w.orig,  # Set to zero so that it is 'deleted'.
                    Xm2 = NULL, Ym2 = NULL,
                    etastart = etastart,  # coefstart = NULL,
                    offset = offset,
                    family = model@family,
                    control = new.control,
                    criterion =  new.control$criterion,  # "coefficients",
                    qr.arg = FALSE,
                    constraints = constraints(model, type = "term"),
                    extra = model@extra,
                    Terms = Terms.zz,
                    function.name = "vglm")

    dfbeta[ii, ] <- coef.model - fit$coeff
  }


  dimnames(dfbeta) <- list(rownames(X.lm), names(coef.model))
  dfbeta
}






setMethod("dfbeta",  "matrix", function(model, ...)
           dfbetavlm(model, ...))


setMethod("dfbeta",  "vlm", function(model, ...)
           dfbetavlm(model, ...))
setMethod("dfbeta",  "vglm", function(model, ...)
           dfbetavlm(model, ...))


setMethod("dfbeta",  "rrvglm", function(model, ...)
           dfbetavlm(model, ...))
setMethod("dfbeta",  "qrrvglm", function(model, ...)
           dfbetavlm(model, ...))
setMethod("dfbeta",  "rrvgam", function(model, ...)
           dfbetavlm(model, ...))
setMethod("dfbeta",  "rcim", function(model, ...)
           dfbetavlm(model, ...))





hatvaluesbasic <- function(X.vlm,
                           diagWm,
                           M = 1) {


  if (M  > 1)
    stop("currently argument 'M' must be 1")

  nn <- nrow(X.vlm)
  ncol.X.vlm <- ncol(X.vlm)

  XtW <- t(c(diagWm) * X.vlm)


  UU <- sqrt(diagWm)  # Only for M == 1
  UU.X.vlm <- c(UU) * X.vlm # c(UU) okay for M==1

  qrSlot <- qr(UU.X.vlm)
  Rmat <- qr.R(qrSlot)

  rinv <- diag(ncol.X.vlm)
  rinv <- backsolve(Rmat, rinv)


  Diag.Hat <- if (FALSE) {
    covun <- rinv %*% t(rinv)
    rhs.mat <- covun %*% XtW
    colSums(t(X.vlm) * rhs.mat)
  } else {
    mymat <- X.vlm %*% rinv
    rowSums(diagWm * mymat^2)
  }
  Diag.Hat
}






 model.matrixpvgam <-
  function(object,
           type = c("vlm", "lm", "lm2", "bothlmlm2",
                    "augmentedvlm", "penalty"),  # This line is new
           linpred.index = NULL,
           ...) {
  type <- match.arg(type, c("vlm", "lm", "lm2", "bothlmlm2",
                            "augmentedvlm", "penalty"))[1]

  if (type == "augmentedvlm" ||
      type == "penalty") {
    rbind(if (type == "penalty") NULL else
          model.matrixvlm(object, type = "vlm",
                          linpred.index = linpred.index, ...),
          get.X.VLM.aug(constraints  = constraints(object, type = "term"),
                        sm.osps.list = object@ospsslot$sm.osps.list))
  } else {
    model.matrixvlm(object, type = type, linpred.index = linpred.index, ...)
  }
}



setMethod("model.matrix",  "pvgam", function(object, ...)
           model.matrixpvgam(object, ...))






