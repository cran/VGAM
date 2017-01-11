# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.






rrvglm.control <-
  function(Rank = 1,
           Algorithm = c("alternating", "derivative"),
           Corner = TRUE,
           Uncorrelated.latvar = FALSE,
           Wmat = NULL,
           Svd.arg = FALSE,
           Index.corner = if (length(str0))
           head((1:1000)[-str0], Rank) else 1:Rank,
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
           ...) {





  if (length(Norrr) != 1 || !is.na(Norrr)) {
    warning("argument 'Norrr' has been replaced by 'noRRR'. ",
            "Assigning the latter but using 'Norrr' will become an ",
            "error in the next VGAM version soon.")
    noRRR <- Norrr
  }


  if (mode(Algorithm) != "character" && mode(Algorithm) != "name")
    Algorithm <- as.character(substitute(Algorithm))
  Algorithm <- match.arg(Algorithm, c("alternating", "derivative"))[1]

  if (Svd.arg)
    Corner <- FALSE

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

  if (length(str0) &&
     (any(round(str0) != str0) || any(str0 < 1)))
    stop("bad input for the argument 'str0'")


  Quadratic <- FALSE
  if (!Quadratic && Algorithm == "derivative" && !Corner) {
    dd <- "derivative algorithm only supports corner constraints"
    if (length(Wmat) || Uncorrelated.latvar || Svd.arg)
      stop(dd)
    warning(dd)
    Corner <- TRUE
  }
  if (Quadratic && Algorithm != "derivative")
   stop("Quadratic model can only be fitted using the derivative algorithm")

  if (Corner && (Svd.arg || Uncorrelated.latvar || length(Wmat)))
      stop("cannot have 'Corner = TRUE' and either 'Svd = TRUE' or ",
           "'Uncorrelated.latvar = TRUE' or Wmat")

  if (Corner && length(intersect(str0, Index.corner)))
    stop("cannot have arguments 'str0' and 'Index.corner' having ",
         "common values")

  if (length(Index.corner) != Rank)
    stop("length(Index.corner) != Rank")

  if (!is.logical(checkwz) ||
      length(checkwz) != 1)
    stop("bad input for 'checkwz'")

  if (!is.Numeric(wzepsilon, length.arg = 1,
                  positive = TRUE))
    stop("bad input for 'wzepsilon'")

  if (class(noRRR) != "formula" && !is.null(noRRR))
    stop("argument 'noRRR' should be a formula or a NULL")


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
           "alternating" = valt.control(...),
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
         Quadratic = FALSE,  # A constant now, here.
         sd.Ainit = sd.Ainit,
         sd.Cinit = sd.Cinit,
         Etamat.colmax = Etamat.colmax,
         str0 = str0,
         Svd.arg = Svd.arg,
         Use.Init.Poisson.QO = Use.Init.Poisson.QO),
    if (Quadratic) qrrvglm.control(Rank = Rank, ...) else NULL)


  if (Quadratic && ans$I.tolerances) {
      ans$Svd.arg <- FALSE
      ans$Uncorrelated.latvar <- FALSE
      ans$Corner <- FALSE
  }


  ans$half.stepsizing <- FALSE  # Turn it off
  ans
}






setClass("summary.rrvglm",
         representation("rrvglm",
                        coef3 = "matrix",
                        cov.unscaled = "matrix",
                        correlation = "matrix",
                        df = "numeric",
                        pearson.resid = "matrix",
                        sigma = "numeric"))


setMethod("summary", "rrvglm",
         function(object, ...)
         summary.rrvglm(object, ...))




show.summary.rrvglm <-
  function(x, digits = NULL, quote = TRUE, prefix = "",
           signif.stars = NULL) {


  show.summary.vglm(x, digits = digits, quote = quote, prefix = prefix)


  invisible(x)
  NULL
}






 setMethod("show", "summary.rrvglm",
           function(object)
             show.summary.rrvglm(x = object))




setMethod("coefficients", "summary.rrvglm", function(object, ...)
          object@coef3)
setMethod("coef",         "summary.rrvglm", function(object, ...)
          object@coef3)



