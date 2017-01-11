# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





summaryvgam <-
  function(object, dispersion = NULL,
           digits = options()$digits-2,
           presid = TRUE,
           nopredictors = FALSE) {



  if (length(dispersion) && dispersion == 0 &&
      length(object@family@summary.dispersion) &&
      !object@family@summary.dispersion) {
    stop("cannot use the general VGLM formula (based on a residual ",
         "sum of squares) for computing the dispersion parameter")
  }

  newobject <- object
  class(newobject) <- "vglm"
  stuff <- summaryvglm(newobject, dispersion = dispersion)
  rdf <- stuff@df[2] <- object@df.residual  # NA

  M <- object@misc$M
  nrow.X.vlm <- object@misc$nrow.X.vlm
  rank <- if (is.null(object@qr$rank)) length(object@coefficients) else
          object@qr$rank








  useF <- object@misc$useF
  if (is.null(useF))
    useF <- FALSE

  df <- unlist(lapply(object@misc$new.assign, length))
  nldf <- object@nl.df

  if (length(df)) {
    aod <- as.matrix(round(df, 1))
    dimnames(aod) <- list(names(df), "Df")
    if (!is.null(object@nl.chisq)) {
      aod <- cbind(aod, NA, NA, NA)
      nl.chisq <- object@nl.chisq / object@dispersion

      special <- abs(nldf) < 0.1  # This was the quick fix in s.vam()
      nldf[special] <- 1  # Give it a plausible value for pchisq & pf

      snames <- names(nldf)
      aod[snames, 2] <- round(nldf, 1)
      aod[snames, 3] <- if (useF) nl.chisq/nldf  else nl.chisq
      aod[snames, 4] <- if (useF)
          pf(nl.chisq / nldf, nldf, rdf, lower.tail = FALSE) else
          pchisq(nl.chisq, nldf, lower.tail = FALSE)

      if (any(special)) {
        aod[snames[special], 2:4] <- NA
      }

      rnames <- c("Df", "Npar Df", "Npar Chisq", "P(Chi)")
      if (useF)
          rnames[3:4] <- c("Npar F", "Pr(F)")
      dimnames(aod) <- list(names(df), rnames)
      heading <- if (useF)
  "\nDF for Terms and Approximate F-values for Nonparametric Effects\n"
      else
  "\nDF for Terms and Approximate Chi-squares for Nonparametric Effects\n"
    } else {
      heading <- "DF for Terms\n\n"
    }
    aod <- as.vanova(data.frame(aod, check.names = FALSE), heading)

    class(aod) <- "data.frame"
  } else {
    aod <- data.frame()
  }

  answer <-
  new("summary.vgam",
      object,
      call = stuff@call,
      cov.unscaled = stuff@cov.unscaled,
      correlation = stuff@correlation,
      df = stuff@df,
      sigma = stuff@sigma)

  slot(answer, "coefficients") <- stuff@coefficients  # Replace
  if (is.numeric(stuff@dispersion))
    slot(answer, "dispersion") <- stuff@dispersion

  if (presid) {
    Presid <- residuals(object, type = "pearson")
    if (length(Presid))
      answer@pearson.resid <- as.matrix(Presid)
  }

  answer@misc$nopredictors <- nopredictors

  slot(answer, "anova") <- aod

  answer
}




show.summary.vgam <-
  function(x, quote = TRUE, prefix = "",
           digits = options()$digits-2,
           nopredictors = NULL) {

  M <- x@misc$M


  cat("\nCall:\n", paste(deparse(x@call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")


  Presid <- x@pearson.resid
  rdf <- x@df[2]
  if (FALSE &&
     !is.null(Presid) && all(!is.na(Presid))) {
    if (rdf/M > 5) {
      rq <-  apply(as.matrix(Presid), 2, quantile)  # 5 x M
      dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                           x@misc$predictors.names)
      cat("\nPearson residuals:\n")
      print(t(rq), digits = digits)
    } else
    if (rdf > 0) {
      cat("\nPearson residuals:\n")
      print(Presid, digits = digits)
    }
  }



  use.nopredictors <- if (is.logical(nopredictors))
    nopredictors else x@misc$nopredictors  # 20140728
  if (!is.logical(use.nopredictors)) {
    warning("cannot determine 'nopredictors'; choosing FALSE")
    use.nopredictors <- FALSE
  }



  cat("\nNumber of linear predictors:   ", M, "\n")




  if (!is.null(x@misc$predictors.names) && !use.nopredictors) {
    if (M == 1) {
      cat("\nName of linear predictor:",
          paste(x@misc$predictors.names, collapse = ", "), "\n")
    } else
    if (M <= 5) {
      cat("\nNames of linear predictors:",
        paste(x@misc$predictors.names, collapse = ", "), fill = TRUE)
    }
  }


  prose <- ""
  if (length(x@dispersion)) {
    if (is.logical(x@misc$estimated.dispersion) &&
        x@misc$estimated.dispersion) {
      prose <- "(Estimated) "
    } else {
      if (is.numeric(x@misc$default.dispersion) &&
          x@dispersion == x@misc$default.dispersion)
        prose <- "(Default) "

      if (is.numeric(x@misc$default.dispersion) &&
          x@dispersion != x@misc$default.dispersion)
        prose <- "(Pre-specified) "
    }
    cat(paste("\n", prose, "Dispersion Parameter for ",
        x@family@vfamily[1],
        " family:   ",
        format(round(x@dispersion, digits)), "\n", sep = ""))
  }

  if (length(deviance(x)))
    cat("\nResidual deviance: ", format(round(deviance(x), digits)),
        "on", format(round(rdf, 3)), "degrees of freedom\n")

  if (length(logLik.vlm(x)))
    cat("\nLog-likelihood:", format(round(logLik.vlm(x), digits)),
        "on", format(round(rdf, 3)), "degrees of freedom\n")

  if (length(x@criterion)) {
    ncrit <- names(x@criterion)
    for (ii in ncrit)
      if (ii != "loglikelihood" && ii != "deviance")
        cat(paste(ii, ":", sep = ""), format(x@criterion[[ii]]), "\n")
  }


  cat("\nNumber of iterations: ", x@iter, "\n")

  if (length(x@anova)) {
    show.vanova(x@anova, digits = digits)   # ".vanova" for Splus6
  }

  invisible(NULL)
}




setMethod("summary", "vgam",
          function(object, ...)
          summaryvgam(object, ...))



setMethod("show", "summary.vgam",
          function(object)
          show.summary.vgam(object))





show.vanova <- function(x, digits = .Options$digits, ...) {
  rrr <- row.names(x)
  heading <- attr(x, "heading")
  if (!is.null(heading))
    cat(heading, sep = "\n")
  attr(x, "heading") <- NULL
  for (ii in seq_along(x)) {
    xx <- x[[ii]]
    xna <- is.na(xx)
    xx <- format(zapsmall(xx, digits))
    xx[xna] <- ""
    x[[ii]] <- xx
  }
  print.data.frame(as.data.frame(x, row.names = rrr))
  invisible(x)
}




as.vanova <- function(x, heading) {
  if (!is.data.frame(x))
    stop("x must be a data frame")
  rrr <- row.names(x)
  attr(x, "heading") <- heading
  x <- as.data.frame(x, row.names = rrr)
  x
}







