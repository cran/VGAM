# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.






rrvglm.control <-
  function(Rank = 1,
           Algorithm = "alternating",
           Corner = TRUE,
           Uncorrelated.latvar = FALSE,
           Wmat = NULL,
           Svd.arg = FALSE,
           Index.corner = head(setdiff(seq(
           length(str0) + Rank), str0), Rank),
           Ainit = NULL,
           Alpha = 0.5,
           Bestof = 1,
           Cinit = NULL,
           Etamat.colmax = 10,
           sd.Ainit = 0.02,
           sd.Cinit = 0.02,
           str0 = NULL,

           noRRR = ~ 1,
           Norrr = NA,

           noWarning = FALSE,

           trace = FALSE,
           Use.Init.Poisson.QO = FALSE,
           checkwz = TRUE,
           Check.rank = TRUE,
           Check.cm.rank = TRUE,
           wzepsilon = .Machine$double.eps^0.75,
           H.A.alt = list(),  # 20231121
           H.C = list(),  # 20231113
           scaleA = FALSE,
           Crow1positive = TRUE,  # 20231128
           ...) {





  if (length(Norrr) != 1 || !is.na(Norrr)) {
    stop("argument 'Norrr' should not be used")
  }


  if (mode(Algorithm) != "character" &&
    mode(Algorithm) != "name")
  Algorithm <- as.character(substitute(Algorithm))
  Algorithm <- match.arg(Algorithm,
      c("alternating"))[1]

  if (Svd.arg)
    Corner <- FALSE
  if (scaleA && Corner) {
    warning("scaleA=TRUE so set Corner=FALSE")
  }

  if (!is.Numeric(Rank, positive = TRUE,
          length.arg = 1, integer.valued = TRUE))
    stop("bad input for 'Rank'")
  if (!is.Numeric(Alpha, positive = TRUE,
                  length.arg = 1) || Alpha > 1)
    stop("bad input for 'Alpha'")
  if (!is.Numeric(Bestof, positive = TRUE,
          length.arg = 1, integer.valued = TRUE))
    stop("bad input for 'Bestof'")
  if (!is.Numeric(sd.Ainit, positive = TRUE,
                  length.arg = 1))
    stop("bad input for 'sd.Ainit'")
  if (!is.Numeric(sd.Cinit, positive = TRUE,
                  length.arg = 1))
    stop("bad input for 'sd.Cinit'")
  if (!is.Numeric(Etamat.colmax, positive = TRUE,
                  length.arg = 1) ||
      Etamat.colmax < Rank)
    stop("bad input for 'Etamat.colmax'")

  if (length(str0))
    stopifnot(round(str0) == str0, str0 >= 1,
              anyDuplicated(str0) == 0)


  Quadratic <- FALSE
  if (!Quadratic && Algorithm == "derivative" &&
      !Corner) {
    dd <- paste0("derivative algorithm only ",
                 "supports corner constraints")
    if (length(Wmat) || Uncorrelated.latvar || Svd.arg)
      stop(dd)
    warning(dd)
    Corner <- TRUE
  }
  if (Quadratic && Algorithm != "derivative")
      stop("Quadratic model can only be fitted",
           " using the derivative algorithm")

  if (Corner && (Svd.arg || Uncorrelated.latvar
      || length(Wmat)))
      stop("cannot have 'Corner = TRUE' and ",
           "either 'Svd = TRUE' or ",
           "'Uncorrelated.latvar = TRUE' or Wmat")
  if (Corner &&
    length(intersect(str0, Index.corner)))
    stop("cannot have 'str0' & 'Index.corner'",
         " having common values")

  if (Corner) {
    stopifnot(length(Index.corner) == Rank)
  } else {
    if (length(Index.corner) != Rank)
      warning("length(Index.corner) != Rank")
  }

  if (!is.logical(checkwz) ||
      length(checkwz) != 1)
    stop("bad input for 'checkwz'")

  H.A.thy <- H.A.alt
  if (length(H.A.alt)) {
    stopifnot(length(H.A.alt) == Rank)
    stopifnot(is.list(H.A.alt))
    ncol.H.A.alt <- unlist(lapply(H.A.alt, ncol))
    nrow.H.A.alt <- unlist(lapply(H.A.alt, nrow))
    if (!all(diff(nrow.H.A.alt) == 0))
      stop("In 'H.A.alt' some nrow() values ",
           "are unequal")  # M unknown here
    if (length(str0))
      for (rr in 1:length(H.A.alt)) stopifnot(
      all.equal(c((H.A.alt[[rr]])[str0, ]) == 0))
    if (Corner &&
        any(ncol.H.A.alt <= Rank))
      stop("Cannot have Corner = TRUE when ",
           "ncol(H.A.alt[[r]]) <= ", Rank)

    H.A.thy <- H.A.alt
    for (rr in 1:length(H.A.alt)) {
      Mt <- H.A.alt[[rr]]
      Mt[Index.corner, ] <- 0
      cols2rm <- which(colSums(abs(Mt)) == 0)
      if (length(cols2rm) > 0 &&
          length(cols2rm) == NCOL(Mt))
        stop("No columns left from Corner = T")
      H.A.thy[[rr]] <- if (length(cols2rm)) 
        Mt[, -cols2rm, drop = FALSE] else Mt
    }
  }  # length(H.A.alt)


  if (length(H.C)) { 
    stopifnot(is.list(H.C))
    if (!all(unlist(lapply(H.C, nrow)) == Rank))
      stop("In 'H.C' some nrow() values are ",
           "not ", Rank)
  }

  drrvglm <- length(H.A.alt) || length(H.C)  #defn
  if (drrvglm) {  # Some normalizations not okay
    stopifnot(!Svd.arg, !Uncorrelated.latvar)
  }

  if (!is.Numeric(wzepsilon, length.arg = 1,
                  positive = TRUE))
    stop("bad input for 'wzepsilon'")


  if (!all(Crow1positive))
    stop("currently 'Crow1positive' must all ",
         "be TRUE")

  if (!is(noRRR, "formula") && !is.null(noRRR))
    stop("arg 'noRRR' should be a formula or NULL")


  ans <-
  c(vglm.control(
                 trace = trace,
                 checkwz = checkwz,
                 Check.rank = Check.rank,
                 Check.cm.rank = Check.cm.rank,
                 wzepsilon = wzepsilon,
                 noWarning = noWarning,
                 ...),

    switch(Algorithm,
      "alternating" = valt0.control(...),
      "derivative"  = rrvglm.optim.control(...)),
    list(Rank = Rank,
         Ainit = Ainit,
         Algorithm = Algorithm,
         Alpha = Alpha,
         Bestof = Bestof,
         Cinit = Cinit,
         Index.corner = Index.corner,
         noRRR = noRRR,

         Corner = Corner,
         Uncorrelated.latvar = Uncorrelated.latvar,
         Wmat = Wmat,
         OptimizeWrtC = TRUE,  # OptimizeWrtC,
         Quadratic = FALSE,  # A const now, here.
         sd.Ainit = sd.Ainit,
         sd.Cinit = sd.Cinit,
         Etamat.colmax = Etamat.colmax,
         str0 = str0,
         Svd.arg = Svd.arg,
         drrvglm = drrvglm,
         scaleA = scaleA,
         Crow1positive = rep(Crow1positive, Rank),
         H.A.alt = H.A.alt,
         H.A.thy = H.A.thy,  # Differs if Corner
         H.C = H.C,
      Use.Init.Poisson.QO = Use.Init.Poisson.QO),
    if (Quadratic)
      qrrvglm.control(Rank = Rank, ...) else NULL)


  if (Quadratic && ans$I.tolerances) {
      ans$Svd.arg <- FALSE
      ans$Uncorrelated.latvar <- FALSE
      ans$Corner <- FALSE
  }

  ans$half.stepsizing <- FALSE  # Turn it off
  ans
}  # rrvglm.control









setClass("summary.rrvglm",
         representation("rrvglm",
                        coef3 = "matrix",
                        coef4lrt0   = "matrix",
                        coef4score0 = "matrix",
                        coef4wald0  = "matrix",
                        cov.unscaled = "matrix",
                        correlation = "matrix",
                        df = "numeric",
                        pearson.resid = "matrix",
                        sigma = "numeric"))


setMethod("summary", "rrvglm",
         function(object, ...)
         summary.rrvglm(object, ...))




show.summary.rrvglm <-
  function(x, digits = NULL, quote = TRUE,
           prefix = "",
           signif.stars = NULL) {


  show.summary.vglm(x, digits = digits,
                    quote = quote, prefix = prefix)


  invisible(x)
  NULL
}






 setMethod("show", "summary.rrvglm",
           function(object)
             show.summary.rrvglm(x = object))




setMethod("coefficients", "summary.rrvglm",
          function(object, ...)
          object@coef3)
setMethod("coef",         "summary.rrvglm",
          function(object, ...)
          object@coef3)



