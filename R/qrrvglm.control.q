# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.


qrrvglm.control = function(Rank=1,
          Bestof = if(length(Cinit)) 1 else 10,
          checkwz=TRUE,
          Cinit = NULL,
          Crow1positive=TRUE,
          epsilon = 1.0e-06,
          EqualTolerances = ITolerances,
          Etamat.colmax = 10,
          FastAlgorithm = TRUE,
          GradientFunction=TRUE,
          Hstep = 0.001,
          isdlv = rep(c(2, 1, rep(0.5, len=Rank)), len=Rank),
          iKvector = 0.1,
          iShape = 0.1,
          ITolerances = TRUE,
          maxitl = 40,
          method.init = 1,
          Maxit.optim = 250,
          MUXfactor = rep(7, length=Rank),
          Norrr = ~ 1,
          optim.maxit = 20,
          Parscale = if(ITolerances) 0.001 else 1.0,
          SD.Cinit = 0.02,
          SmallNo = 5.0e-13,
          trace = TRUE,
          Use.Init.Poisson.QO=TRUE,
          wzepsilon = .Machine$double.eps^0.75,
          ...)
{



    if(!is.Numeric(iShape, posit=TRUE)) stop("bad input for \"iShape\"")
    if(!is.Numeric(iKvector, posit=TRUE)) stop("bad input for \"iKvector\"")
    if(!is.Numeric(isdlv, posit=TRUE)) stop("bad input for \"isdlv\"")
    if(any(isdlv < 0.2 | isdlv > 10))
        stop("isdlv values must lie between 0.2 and 10")
    if(length(isdlv) > 1 && any(diff(isdlv) > 0))
        stop("successive isdlv values must not increase")
    if(!is.Numeric(epsilon, posit=TRUE, allow=1)) 
        stop("bad input for \"epsilon\"")
    if(!is.Numeric(Etamat.colmax, posit=TRUE, allow=1) || Etamat.colmax < Rank)
        stop("bad input for \"Etamat.colmax\"")
    if(!is.Numeric(Hstep, posit=TRUE, allow=1)) 
        stop("bad input for \"Hstep\"")
    if(!is.Numeric(maxitl, posit=TRUE, allow=1, integer=TRUE)) 
        stop("bad input for \"maxitl\"")
    if(!is.Numeric(method.init, posit=TRUE, allow=1, integer=TRUE)) 
        stop("bad input for \"method.init\"")
    if(!is.Numeric(Maxit.optim, integ=TRUE, posit=TRUE))
        stop("Bad input for \"Maxit.optim\"")
    if(!is.Numeric(MUXfactor, posit=TRUE)) 
        stop("bad input for \"MUXfactor\"")
    if(any(MUXfactor < 1 | MUXfactor > 10))
        stop("MUXfactor values must lie between 1 and 10")
    if(!is.Numeric(optim.maxit, allow=1, integ=TRUE, posit=TRUE))
        stop("Bad input for \"optim.maxit\"")
    if(!is.Numeric(Rank, posit=TRUE, allow=1, integer=TRUE)) 
        stop("bad input for \"Rank\"")
    if(!is.Numeric(SD.Cinit, posit=TRUE, allow=1)) 
        stop("bad input for \"SD.Cinit\"")
    if(ITolerances && !EqualTolerances)
        stop("EqualTolerances must be TRUE if ITolerances is TRUE")
    if(!is.Numeric(Bestof, posit=TRUE, allow=1, integer=TRUE)) 
        stop("bad input for \"Bestof\"")


    FastAlgorithm = as.logical(FastAlgorithm)[1]
    if(!FastAlgorithm)
        stop("FastAlgorithm=TRUE is now required")

    if((SmallNo < .Machine$double.eps) ||
       (SmallNo > .0001)) stop("SmallNo is out of range") 
    if(any(Parscale <= 0))
       stop("Parscale must contain positive numbers only") 

    if(!is.logical(checkwz) || length(checkwz) != 1)
        stop("bad input for \"checkwz\"")
    if(!is.Numeric(wzepsilon, allow=1, positive=TRUE))
        stop("bad input for \"wzepsilon\"")

    ans = list(
           Bestof = Bestof,
           checkwz=checkwz,
           Cinit = Cinit,
           Crow1positive=as.logical(rep(Crow1positive, len=Rank)),
           ConstrainedQO = TRUE, # A constant, not a control parameter
           Corner = FALSE, # Needed for valt.1iter()
           Dzero = NULL,
           epsilon = epsilon,
           EqualTolerances = EqualTolerances,
           Etamat.colmax = Etamat.colmax,
           FastAlgorithm = FastAlgorithm,
           GradientFunction = GradientFunction,
           Hstep = Hstep,
           isdlv = rep(isdlv, len=Rank),
           iKvector = as.numeric(iKvector),
           iShape = as.numeric(iShape),
           ITolerances = ITolerances,
           maxitl = maxitl,
           method.init = method.init,
           Maxit.optim = Maxit.optim,
           min.criterion = TRUE, # needed for calibrate 
           MUXfactor = rep(MUXfactor, length=Rank),
           Norrr = Norrr,
           optim.maxit = optim.maxit,
           OptimizeWrtC = TRUE,
           Parscale = Parscale,
           Quadratic = TRUE,
           Rank = Rank,
           save.weight = FALSE,
           SD.Cinit = SD.Cinit,
           SmallNo = SmallNo,
           Structural.zero = NULL,
           Svd.arg = TRUE, Alpha=0.5, Uncor = TRUE,
           trace = trace,
           Use.Init.Poisson.QO=as.logical(Use.Init.Poisson.QO)[1],
           wzepsilon = wzepsilon)
    ans
}

