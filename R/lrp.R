# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.











lrt.stat.vlm <-
  function(object,
           values0 = 0,
           subset = NULL,  # Useful for Cox model as a poissonff().
           omit1s = TRUE,
           all.out = FALSE,  # If TRUE then lots of output returned
           trace = FALSE,  # NULL,
           ...) {


  wald.stat.vlm(object, values0 = values0,
                subset = subset, omit1s = omit1s, all.out = all.out,
                iterate.SE = TRUE,
                trace = trace,
                orig.SE = FALSE,  # Does not make sense if TRUE
                LR.really = TRUE, ...)
}  # lrt.stat.vlm



if (!isGeneric("lrt.stat"))
    setGeneric("lrt.stat", function(object, ...)
               standardGeneric("lrt.stat"))


setMethod("lrt.stat", "vlm", function(object, ...)
          lrt.stat.vlm(object = object, ...))












