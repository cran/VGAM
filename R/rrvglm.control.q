# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.







ridargs.rrvglm.control <-
  function(Uncorrelated.latvar = NULL,
           Wmat = NULL,
           Svd.arg = NULL,
           Alpha = NULL,
           scaleA = NULL,
           Norrr = NULL, ...) {
  if (length(Uncorrelated.latvar) +
      length(Wmat) +
      length(Svd.arg) +
      length(Alpha) +
      length(scaleA) +
      length(Norrr))
      stop("the following arguments are no ",
           "longer in use: 'Wmat', 'Norrr'",
           "'Svd.arg', 'Uncorrelated.latvar'",
           "'Alpha', 'scaleA'")
}





rrvglm.control <-
  function(Rank = 1,
           Corner = TRUE,
           Index.corner = head(setdiff(seq(
           length(str0) + Rank), str0), Rank),
           noRRR = ~ 1,
           str0 = NULL,
           Crow1positive = NULL,  # 20231128
           trace = FALSE,
           Bestof = 1,
           H.A.thy = list(),  # 20231121
           H.C = list(),  # 20231113
           Ainit = NULL,
           Cinit = NULL,
           sd.Cinit = 0.02,
           Algorithm = "alternating",
           Etamat.colmax = 10,

           noWarning = FALSE,

           Use.Init.Poisson.QO = FALSE,
           checkwz = TRUE,
           Check.rank = TRUE,
           Check.cm.rank = TRUE,
           wzepsilon = .Machine$double.eps^0.75,
           ...) {


  ridargs.rrvglm.control(...)


           label.it = TRUE   # 20240315
           scaleA = FALSE  # 20240327
           Uncorrelated.latvar = FALSE  # 20240327
           Wmat = NULL  # 20240327
           Svd.arg = FALSE  # 20240327
           Norrr = NA  # 20240327
           Alpha = 0.5  # 202403278






  if (length(Norrr) != 1 || !is.na(Norrr)) {
    stop("argument 'Norrr' should not be used")
  }


  if (mode(Algorithm) != "character" &&
    mode(Algorithm) != "name")
  Algorithm <- as.character(substitute(Algorithm))
  Algorithm <- match.arg(Algorithm,
      c("alternating"))[1]  # Had "derivative"

  if (!is.Numeric(Rank, positive = TRUE,
          length.arg = 1, integer.valued = TRUE))
    stop("bad input for 'Rank'")
  if (!is.Numeric(Alpha, positive = TRUE,
                  length.arg = 1) || Alpha > 1)
    stop("bad input for 'Alpha'")
  if (!is.Numeric(Bestof, positive = TRUE,
          length.arg = 1, integer.valued = TRUE))
    stop("bad input for 'Bestof'")
  if (!is.Numeric(sd.Cinit, positive = TRUE,
                  length.arg = 1))
    stop("bad input for 'sd.Cinit'")
  if (!is.Numeric(Etamat.colmax, positive = TRUE,
                  length.arg = 1) ||
      Etamat.colmax < Rank)
    stop("bad input for 'Etamat.colmax'")
  if (!is.Numeric(wzepsilon, length.arg = 1,
                  positive = TRUE))
    stop("bad input for 'wzepsilon'")

  if (length(str0))
    stopifnot(round(str0) == str0, str0 >= 1,
              anyDuplicated(str0) == 0)

  Quadratic <- FALSE

  if (Corner && (Svd.arg || Uncorrelated.latvar
      || scaleA || length(Wmat)))
    stop("cannot have 'Corner = T' and ",
         "either 'Svd = T' or 'scaleA = T' or ",
         "'Uncorrelated.latvar = T' or Wmat")
  if (Corner &&
      length(intersect(str0, Index.corner)))
    stop("cannot have 'str0' & 'Index.corner'",
         " having common values")

  if (Corner) {
    stopifnot(length(Index.corner) == Rank,
        all(round(Index.corner) == Index.corner),
        Index.corner > 0,  # Index.corner <= M,
        anyDuplicated(Index.corner) == 0)
  } else {
    if (length(Index.corner) != Rank)
      warning("length(Index.corner) != Rank")
  }

  if (!is.logical(checkwz) ||
      length(checkwz) != 1)
    stop("bad input for 'checkwz'")

  if (is.matrix(H.A.thy))
    H.A.thy <- list(H.A.thy)

  H.A.alt <- H.A.thy  # H.A.alt is smaller
  Amask <- NULL  # For "rrvglm"
  H.A.thy.trivial <- TRUE  # For "rrvglm"
  H.C.trivial <- TRUE  # For "rrvglm"
  if (length(H.A.thy)) {
    stopifnot(Corner,  # RCC needed.
              length(H.A.thy) == Rank,
              is.list(H.A.thy))
    for (r in seq(Rank))  # <= 1 non0 elt per row
      stopifnot(all(apply(H.A.thy[[r]], 1,
                          is.non0) <= 1))
    foo4 <- function(A) colSums(abs(A))
    if (any(sapply(H.A.thy, foo4) == 0))
      stop("there's a coln of all 0s in H.A.thy")
    ncol.H.A.thy <- sapply(H.A.thy, ncol)
    nrow.H.A.thy <- sapply(H.A.thy, nrow)
    if (!all(diff(nrow.H.A.thy) == 0))
      stop("In 'H.A.thy' some nrow() values ",
           "are unequal")  # M unknown here
    M <- nrow.H.A.thy[1]
    Amask <- matrix(NA_real_, M, Rank)
    if (length(str0))
      stop("cannot use 'str0' when 'H.A' is used",
          ". Instead, build it into H.A")
    H.A.alt <- H.A.thy
    for (rr in 1:length(H.A.thy)) {
      Mt <- H.A.thy[[rr]]
      ind.col <- which(Mt[Index.corner[rr],] != 0)
      if (length(ind.col) != 1)
        stop("row ", Index.corner[rr], " of ",
             "H.A.thy is bad, e.g., all 0s?")
      ind.row <- which(Mt[, ind.col] != 0)
      Mt.ii <- Mt[Index.corner[rr], ind.col]
      if (Mt.ii == 0)
        stop("element Mt.ii is 0")
      Amask[ind.row, rr] <-
         Mt[ind.row, ind.col] / Mt.ii
      Mt[ind.row, ind.col] <- 0
      H.A.alt[[rr]] <- rm0cols(Mt)  # May stop().
    }  # rr

    bigmat <- matrix(unlist(H.A.thy), nrow = M)
    bigmat[bigmat != 0] <- 1
    if (any(rowSums(bigmat) > 1))
      stop("the reduced rank regression is ",
           "not separable (A is not block-",
           "diagonal after reordering of rows)")

    if (length(str0) &&
        any(rowSums(bigmat)[str0] > 0))
      stop("conflict between str0 and H.A.thy")

    H.A.thy.trivial <- all(sapply(H.A.thy,
                                  is.Identity))
  }  # length(H.A.thy)

  if (length(H.C)) { 
    stopifnot(is.list(H.C), Corner)
    if (!all(sapply(H.C, nrow) == Rank))
      stop("In 'H.C' some nrow() values are ",
           "not ", Rank)
    H.C.trivial <- all(sapply(H.C, is.Identity))
  }




  is.rrvglm <- length(c(H.A.thy, H.C)) == 0 ||
               (H.A.thy.trivial && H.C.trivial)
  is.drrvglm <- !is.rrvglm
  if (is.drrvglm) {  # Some normalizations r bad
    stopifnot(!Svd.arg, !Uncorrelated.latvar)
    if (length(str0))
    warning("cant use both 'str0' & 'H.A.thy'")
  }



  if (length(Crow1positive))  # ... compromise4now
    stop("currently 'Crow1positive' must be NULL")

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
      "alternating" = valt0.control(...)),
    list(Rank = Rank,
         Ainit = Ainit,
         Cinit = Cinit,
         Algorithm = Algorithm,
         Alpha = Alpha,
         Bestof = Bestof,
         Index.corner = Index.corner,
         noRRR = noRRR,

         Corner = Corner,
         Uncorrelated.latvar = Uncorrelated.latvar,
         Wmat = Wmat,
         OptimizeWrtC = TRUE,  # OptimizeWrtC,
         Quadratic = FALSE,  # A const now, here.
         sd.Cinit = sd.Cinit,
         Etamat.colmax = Etamat.colmax,
         str0 = str0,  # NULL for "drrvglm"
         Svd.arg = Svd.arg,
         is.drrvglm = is.drrvglm,
         is.rrvglm = is.rrvglm,
         scaleA = scaleA,
         Crow1positive = rep(Crow1positive, Rank),
         label.it = label.it,
         Amask   = Amask,
         H.A.thy = H.A.thy,  # Bigger;
         H.A.alt = H.A.alt,  # Smaller
         H.C = H.C,  # list() if unused
      Use.Init.Poisson.QO = Use.Init.Poisson.QO))

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



