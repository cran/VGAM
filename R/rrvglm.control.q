# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.




rrvglm.control = function(Rank = 1,
                          Algorithm = c("alternating", "derivative"),
                          Corner = TRUE,
                          Uncorrelated.lv = FALSE,
                          Wmat = NULL,
                          Svd.arg = FALSE,
                          Index.corner = if (length(szero)) 
                          head((1:1000)[-szero], Rank) else 1:Rank,
                          Ainit = NULL,
                          Alpha = 0.5, 
                          Bestof = 1,
                          Cinit = NULL,
                          Etamat.colmax = 10,
                          SD.Ainit = 0.02,
                          SD.Cinit = 0.02,
                          szero = NULL,
                          Norrr = ~ 1, 
                          trace = FALSE,
                          Use.Init.Poisson.QO = FALSE,
                          checkwz = TRUE,
                          wzepsilon = .Machine$double.eps^0.75,
                          ...)
{




    if (mode(Algorithm) != "character" && mode(Algorithm) != "name")
        Algorithm <- as.character(substitute(Algorithm))
    Algorithm <- match.arg(Algorithm, c("alternating", "derivative"))[1]

    if (Svd.arg) Corner = FALSE 

    if (!is.Numeric(Rank, positive = TRUE,
                    allowable.length = 1, integer.valued = TRUE))
      stop("bad input for 'Rank'")
    if (!is.Numeric(Alpha, positive = TRUE,
                    allowable.length = 1) || Alpha > 1)
      stop("bad input for 'Alpha'")
    if (!is.Numeric(Bestof, positive = TRUE,
                    allowable.length = 1, integer.valued = TRUE))
      stop("bad input for 'Bestof'")
    if (!is.Numeric(SD.Ainit, positive = TRUE,
                    allowable.length = 1))
      stop("bad input for 'SD.Ainit'")
    if (!is.Numeric(SD.Cinit, positive = TRUE,
                    allowable.length = 1))
      stop("bad input for 'SD.Cinit'")
    if (!is.Numeric(Etamat.colmax, positive = TRUE,
                    allowable.length = 1) ||
        Etamat.colmax < Rank)
      stop("bad input for 'Etamat.colmax'")

    if (length(szero) &&
       (any(round(szero) != szero) ||
       any(szero < 1)))
      stop("bad input for the argument 'szero'")


    Quadratic = FALSE
    if (!Quadratic && Algorithm == "derivative" && !Corner) {
        dd = "derivative algorithm only supports corner constraints"
        if (length(Wmat) || Uncorrelated.lv || Svd.arg)
            stop(dd)
        warning(dd)
        Corner = TRUE
    }
    if (Quadratic && Algorithm != "derivative")
        stop("Quadratic model can only be fitted using the derivative algorithm")

    if (Corner && (Svd.arg || Uncorrelated.lv || length(Wmat)))
        stop("cannot have 'Corner = TRUE' and either 'Svd = TRUE' or ",
             "'Uncorrelated.lv = TRUE' or Wmat")

    if (Corner && length(intersect(szero, Index.corner)))
      stop("cannot have 'szero' and 'Index.corner' having ",
           "common values")

    if (length(Index.corner) != Rank)
      stop("length(Index.corner) != Rank")

    if (!is.logical(checkwz) ||
        length(checkwz) != 1)
      stop("bad input for 'checkwz'")

    if (!is.Numeric(wzepsilon, allowable.length = 1,
                    positive = TRUE))
      stop("bad input for 'wzepsilon'")

    if (class(Norrr) != "formula" && !is.null(Norrr))
      stop("argument 'Norrr' should be a formula or a NULL")

    ans =
    c(vglm.control(trace = trace, ...),
      switch(Algorithm,
             "alternating" = valt.control(...),
             "derivative" = rrvglm.optim.control(...)),
      list(Rank = Rank,
           Ainit = Ainit,
           Algorithm = Algorithm,
           Alpha = Alpha,
           Bestof = Bestof,
           Cinit = Cinit,
           Index.corner = Index.corner,
           Norrr = Norrr,
           Corner = Corner,
           Uncorrelated.lv = Uncorrelated.lv,
           Wmat = Wmat,
           OptimizeWrtC = TRUE, # OptimizeWrtC,
           Quadratic = FALSE,   # A constant now, here.
           SD.Ainit = SD.Ainit,
           SD.Cinit = SD.Cinit,
           Etamat.colmax = Etamat.colmax,
           szero = szero,
           Svd.arg=Svd.arg,
           Use.Init.Poisson.QO = Use.Init.Poisson.QO),
           checkwz = checkwz,
           wzepsilon = wzepsilon,
      if (Quadratic) qrrvglm.control(Rank = Rank, ...) else NULL)

    if (Quadratic && ans$ITolerances) {
        ans$Svd.arg = FALSE
        ans$Uncorrelated.lv = FALSE
        ans$Corner = FALSE
    }

    ans$half.stepsizing = FALSE   # Turn it off 
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




show.summary.rrvglm <- function(x, digits = NULL,
                                quote= TRUE, prefix = "")
{


    show.summary.vglm(x, digits = digits, quote = quote, prefix = prefix)


    invisible(x)
    NULL
}





 setMethod("show", "summary.rrvglm",
           function(object)
             show.summary.rrvglm(x = object))






