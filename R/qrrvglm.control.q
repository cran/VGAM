# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





qrrvglm.control <- function(Rank = 1,
          Bestof = if (length(Cinit)) 1 else 10,
          checkwz = TRUE,
          Cinit = NULL,
          Crow1positive = TRUE,
          epsilon = 1.0e-06,
          EqualTolerances = NULL,
          eq.tolerances = TRUE,  # 20140520; replaces EqualTolerances
          Etamat.colmax = 10,
          FastAlgorithm = TRUE,
          GradientFunction = TRUE,
          Hstep = 0.001,
          isd.latvar = rep_len(c(2, 1, rep_len(0.5, Rank)), Rank),
          iKvector = 0.1,
          iShape = 0.1,
          ITolerances = NULL,
          I.tolerances = FALSE,  # 20140520; replaces ITolerances
          maxitl = 40,
          imethod = 1,
          Maxit.optim = 250,
          MUXfactor = rep_len(7, Rank),
          noRRR = ~ 1,
          Norrr = NA,
          optim.maxit = 20,
          Parscale = if (I.tolerances) 0.001 else 1.0,
          sd.Cinit = 0.02,
          SmallNo = 5.0e-13,
          trace = TRUE,
          Use.Init.Poisson.QO = TRUE,
          wzepsilon = .Machine$double.eps^0.75,
          ...) {






  if (!is.null(EqualTolerances)) {
    warning("argument 'EqualTolerances' is depreciated. ",
            "Use argument 'eq.tolerances'")
    if (is.logical(EqualTolerances)) {
      if (eq.tolerances != EqualTolerances)
        stop("arguments 'eq.tolerances' and 'EqualTolerances' differ")
    } else {
      stop("argument 'EqualTolerances' is not a logical")
    }
  }




  if (!is.null(ITolerances)) {
    warning("argument 'ITolerances' is depreciated. ",
            "Use argument 'I.tolerances'")
    if (is.logical(ITolerances)) {
      if (I.tolerances != ITolerances)
        stop("arguments 'I.tolerances' and 'ITolerances' differ")
    } else {
      stop("argument 'ITolerances' is not a logical")
    }
  }






  if (length(Norrr) != 1 || !is.na(Norrr)) {
    warning("argument 'Norrr' has been replaced by 'noRRR'. ",
            "Assigning the latter but using 'Norrr' will become an error in ",
            "the next VGAM version soon.")
    noRRR <- Norrr
  }



    if (!is.Numeric(iShape, positive = TRUE))
      stop("bad input for 'iShape'")
    if (!is.Numeric(iKvector, positive = TRUE))
      stop("bad input for 'iKvector'")
    if (!is.Numeric(isd.latvar, positive = TRUE))
      stop("bad input for 'isd.latvar'")
    if (any(isd.latvar < 0.2 |
            isd.latvar > 10))
        stop("isd.latvar values must lie between 0.2 and 10")
    if (length(isd.latvar) > 1 && any(diff(isd.latvar) > 0))
        stop("successive isd.latvar values must not increase")
    if (!is.Numeric(epsilon, positive = TRUE,
                    length.arg = 1))
        stop("bad input for 'epsilon'")
    if (!is.Numeric(Etamat.colmax, positive = TRUE,
                    length.arg = 1) ||
        Etamat.colmax < Rank)
        stop("bad input for 'Etamat.colmax'")
    if (!is.Numeric(Hstep, positive = TRUE,
                   length.arg = 1))
        stop("bad input for 'Hstep'")
    if (!is.Numeric(maxitl, positive = TRUE,
                    length.arg = 1, integer.valued = TRUE))
        stop("bad input for 'maxitl'")
    if (!is.Numeric(imethod, positive = TRUE,
                    length.arg = 1, integer.valued = TRUE))
        stop("bad input for 'imethod'")
    if (!is.Numeric(Maxit.optim, integer.valued = TRUE, positive = TRUE))
        stop("Bad input for 'Maxit.optim'")
    if (!is.Numeric(MUXfactor, positive = TRUE))
        stop("bad input for 'MUXfactor'")
    if (any(MUXfactor < 1 | MUXfactor > 10))
        stop("MUXfactor values must lie between 1 and 10")
    if (!is.Numeric(optim.maxit, length.arg = 1,
                    integer.valued = TRUE, positive = TRUE))
        stop("Bad input for 'optim.maxit'")
    if (!is.Numeric(Rank, positive = TRUE,
                    length.arg = 1, integer.valued = TRUE))
        stop("bad input for 'Rank'")
    if (!is.Numeric(sd.Cinit, positive = TRUE,
                    length.arg = 1))
        stop("bad input for 'sd.Cinit'")
    if (I.tolerances && !eq.tolerances)
        stop("'eq.tolerances' must be TRUE if 'I.tolerances' is TRUE")
    if (!is.Numeric(Bestof, positive = TRUE,
                    length.arg = 1, integer.valued = TRUE))
        stop("bad input for 'Bestof'")


    FastAlgorithm = as.logical(FastAlgorithm)[1]
    if (!FastAlgorithm)
        stop("FastAlgorithm = TRUE is now required")

    if ((SmallNo < .Machine$double.eps) ||
       (SmallNo > .0001))
      stop("SmallNo is out of range")
    if (any(Parscale <= 0))
       stop("Parscale must contain positive numbers only")

    if (!is.logical(checkwz) ||
        length(checkwz) != 1)
        stop("bad input for 'checkwz'")
    if (!is.Numeric(wzepsilon,
                    length.arg = 1, positive = TRUE))
        stop("bad input for 'wzepsilon'")

    ans <- list(
           Bestof = Bestof,
           checkwz = checkwz,
           Cinit = Cinit,
           Crow1positive=as.logical(rep_len(Crow1positive, Rank)),
           ConstrainedQO = TRUE,  # A constant, not a control parameter
           Corner = FALSE,  # Needed for valt.1iter()
           Dzero = NULL,
           epsilon = epsilon,
           eq.tolerances = eq.tolerances,
           Etamat.colmax = Etamat.colmax,
           FastAlgorithm = FastAlgorithm,
           GradientFunction = GradientFunction,
           Hstep = Hstep,
           isd.latvar = rep_len(isd.latvar, Rank),
           iKvector = as.numeric(iKvector),
           iShape = as.numeric(iShape),
           I.tolerances = I.tolerances,
           maxitl = maxitl,
           imethod = imethod,
           Maxit.optim = Maxit.optim,
           min.criterion = TRUE,  # needed for calibrate
           MUXfactor = rep_len(MUXfactor, Rank),
           noRRR = noRRR,
           optim.maxit = optim.maxit,
           OptimizeWrtC = TRUE,
           Parscale = Parscale,
           Quadratic = TRUE,
           Rank = Rank,
           save.weights = FALSE,
           sd.Cinit = sd.Cinit,
           SmallNo = SmallNo,
           str0 = NULL,
           Svd.arg = TRUE, Alpha = 0.5, Uncorrelated.latvar = TRUE,
           trace = trace,
           Use.Init.Poisson.QO = as.logical(Use.Init.Poisson.QO)[1],
           wzepsilon = wzepsilon)
    ans
}



