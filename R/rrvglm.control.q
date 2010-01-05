# These functions are
# Copyright (C) 1998-2010 T.W. Yee, University of Auckland. All rights reserved.




rrvglm.control = function(Rank=1, 
                          Algorithm=c("alternating", "derivative"),
                          Corner=TRUE,
                          Uncorrelated.lv=FALSE,
                          Wmat=NULL,
                          Svd.arg=FALSE,
                          Index.corner = if (length(Structural.zero)) 
                          head((1:1000)[-Structural.zero], Rank) else 1:Rank,
                          Alpha=0.5, 
                          Bestof = 1,
                          Cinit=NULL,
                          Etamat.colmax = 10,
                          SD.Cinit = 0.02,
                          Structural.zero = NULL,
                          Norrr = ~ 1, 
                          trace = FALSE,
                          Use.Init.Poisson.QO=FALSE,
                          checkwz=TRUE,
                          wzepsilon = .Machine$double.eps^0.75,
                          ...)
{




    if (mode(Algorithm) != "character" && mode(Algorithm) != "name")
        Algorithm <- as.character(substitute(Algorithm))
    Algorithm <- match.arg(Algorithm, c("alternating", "derivative"))[1]

    if (Svd.arg) Corner = FALSE 

    if (!is.Numeric(Rank, posit=TRUE, allow=1, integer=TRUE))
        stop("bad input for 'Rank'")
    if (!is.Numeric(Alpha, posit=TRUE, allow=1) || Alpha > 1)
        stop("bad input for 'Alpha'")
    if (!is.Numeric(Bestof, posit=TRUE, allow=1, integer=TRUE))
        stop("bad input for 'Bestof'")
    if (!is.Numeric(SD.Cinit, posit=TRUE, allow=1))
        stop("bad input for 'SD.Cinit'")
    if (!is.Numeric(Etamat.colmax, posit=TRUE, allow=1) || Etamat.colmax < Rank)
        stop("bad input for 'Etamat.colmax'")

    if (length(Structural.zero) && (any(round(Structural.zero) != Structural.zero)
       || any(Structural.zero<1)))
        stop("bad input for the argument 'Structural.zero'")


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
        stop("cannot have Corner=TRUE and either Svd=TRUE or Uncorrelated.lv=TRUE or Wmat")

    if (Corner && length(intersect(Structural.zero, Index.corner)))
    stop("cannot have Structural.zero and Index.corner having common values")

    if (length(Index.corner) != Rank)
        stop("length(Index.corner) != Rank")

    if (!is.logical(checkwz) || length(checkwz) != 1)
        stop("bad input for 'checkwz'")
    if (!is.Numeric(wzepsilon, allow=1, positive=TRUE))
        stop("bad input for 'wzepsilon'")

    ans =
    c(vglm.control(trace = trace, ...),
      switch(Algorithm,
             "alternating" = valt.control(...),
             "derivative" = if (is.R()) rrvglm.optim.control(...) else
                                nlminbcontrol(...)),
      list(Rank=Rank,
           Algorithm=Algorithm,
           Alpha=Alpha,
           Bestof = Bestof,
           Cinit=Cinit,
           Index.corner=Index.corner,
           Norrr=Norrr,
           Corner=Corner, Uncorrelated.lv=Uncorrelated.lv, Wmat=Wmat,
           OptimizeWrtC = TRUE, # OptimizeWrtC,
           Quadratic = FALSE,   # A constant now, here.
           SD.Cinit = SD.Cinit,
           Etamat.colmax = Etamat.colmax,
           Structural.zero = Structural.zero,
           Svd.arg=Svd.arg,
           Use.Init.Poisson.QO=Use.Init.Poisson.QO),
           checkwz=checkwz,
           wzepsilon = wzepsilon,
      if (Quadratic) qrrvglm.control(Rank=Rank, ...) else NULL)

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
    coef3="matrix",
    cov.unscaled="matrix",
    correlation="matrix",
    df="numeric",
    pearson.resid="matrix",
    sigma="numeric"))

setMethod("summary", "rrvglm",
         function(object, ...)
         summary.rrvglm(object, ...))




printsummary.rrvglm <- function(x, digits=NULL, quote= TRUE, prefix="")
{


    printsummary.vglm(x, digits = digits, quote = quote, prefix = prefix)


    invisible(x)
}


setMethod("print", "summary.rrvglm",
         function(x, ...)
         printsummary.rrvglm(x=x, ...))

    setMethod("show", "summary.rrvglm",
             function(object)
             printsummary.rrvglm(x=object))

