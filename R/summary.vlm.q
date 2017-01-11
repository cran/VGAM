# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.












summaryvlm <-
  function(object, correlation = FALSE, dispersion = NULL,
           Colnames = c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
           presid = TRUE) {



  if (is.logical(object@misc$BFGS) && object@misc$BFGS)
    warning("the estimated variance-covariance matrix is ",
            "usually inaccurate because the working weight matrices are ",
            "obtained by a crude BFGS quasi-Newton approximation")

  M <- object@misc$M
  n <- object@misc$n
  nrow.X.vlm <- object@misc$nrow.X.vlm
  ncol.X.vlm <- object@misc$ncol.X.vlm  # May be NULL for CQO objects

  Coefs <- object@coefficients
  cnames <- names(Coefs)

  Presid <- if (presid) {
    Presid <- residualsvlm(object, type = "pearson")  # NULL if pooled.weight
    Presid
  } else {
    NULL
  }

  if (anyNA(Coefs)) {
    warning("Some NAs in the coefficients---no summary is ",
            " provided; returning 'object'")
    return(object)
  }
  rdf <- object@df.residual

  if (!length(dispersion)) {
    if (is.numeric(object@misc$dispersion)) {
      dispersion <- object@misc$dispersion
      if (all(dispersion == 0))
        stop("dispersion shouldn't be zero here!")
    } else {
      dispersion <- 1
      object@misc$estimated.dispersion <- FALSE
    }
  } else if (dispersion == 0) {
      dispersion <-
        if (!length(object@ResSS)) {
          stop("object@ResSS is empty")
      } else {
        object@ResSS / object@df.residual
      }
      object@misc$estimated.dispersion <- TRUE
  } else {
    if (is.numeric(object@misc$dispersion) &&
        object@misc$dispersion != dispersion)
      warning("overriding the value of object@misc$dispersion")
    object@misc$estimated.dispersion <- FALSE
  }
  sigma <- sqrt(dispersion)  # Can be a vector

  if (is.Numeric(ncol.X.vlm)) {
    R <- object@R

    if (ncol.X.vlm < max(dim(R)))
      stop("R is rank deficient")





    covun <- chol2inv(R)

    dimnames(covun) <- list(cnames, cnames)
  }
  coef3 <- matrix(rep(Coefs, 4), ncol = 4)
  dimnames(coef3) <- list(cnames, Colnames)
  SEs <- sqrt(diag(covun))
  if (length(sigma) == 1 && is.Numeric(ncol.X.vlm)) {
    coef3[, 2] <- SEs %o% sigma  # Fails here when sigma is a vector
    coef3[, 3] <- coef3[, 1] / coef3[, 2]
    pvalue <- 2 * pnorm(-abs(coef3[, 3]))
    coef3[, 4] <- pvalue

    if (is.logical(object@misc$estimated.dispersion) &&
        object@misc$estimated.dispersion)
      coef3 <- coef3[, -4]  # Delete the pvalues column
  } else {
    coef3[, 1] <- coef3[, 2] <- coef3[, 3] <- coef3[, 4] <- NA
    coef3 <- coef3[, -4]  # Delete the pvalues column
  }

  if (correlation) {
    correl <- covun * outer(1 / SEs, 1 / SEs)

    diag(correl) <- 1.0

    dimnames(correl) <- list(cnames, cnames)
  } else {
    correl <- matrix(0, 0, 0)  # was NULL, but now a special matrix
  }




  answer <-
  new("summary.vlm",
      object,
      coef3       = coef3,
      correlation = correl,
      df          = c(ncol.X.vlm, rdf),
      sigma       = sigma)

  if (is.Numeric(ncol.X.vlm))
    answer@cov.unscaled <- covun
  answer@dispersion <- dispersion  # Overwrite this

  if (length(Presid))
    answer@pearson.resid <- as.matrix(Presid)


  answer
}





show.summary.vlm <- function(x, digits = NULL, quote = TRUE,
                             prefix = "") {


  M <- x@misc$M
  coef3 <- x@coef3 # ficients
  correl <- x@correlation

  if (is.null(digits)) {
    digits <- options()$digits
  } else {
    old.digits <- options(digits = digits)
    on.exit(options(old.digits))
  }

  cat("\nCall:\n")
  dput(x@call)

  Presid <- x@pearson.resid
  rdf <- x@df[2]
  if (length(Presid) && all(!is.na(Presid))) {
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

  if (!all(is.na(coef3))) {
    cat("\nCoefficients:\n")
    print(coef3, digits = digits)
  }

  cat("\nNumber of responses: ", M, "\n")


  if (length(x@misc$predictors.names))
  if (M == 1) {
    cat("\nName of response:",
        paste(x@misc$predictors.names, collapse = ", "), "\n")
  } else {
    UUU <- paste(x@misc$predictors.names, collapse = ", ")
    UUU <- x@misc$predictors.names
    cat("\nNames of responses:\n")
    cat(UUU, fill = TRUE, sep = ", ")
  }


  if (!is.null(x@ResSS))
    cat("\nResidual Sum of Squares:", format(round(x@ResSS, digits)),
        "on", round(rdf, digits), "degrees of freedom\n")


  if (length(correl)) {
    ncol.X.vlm <- dim(correl)[2]
    if (ncol.X.vlm > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1, -ncol.X.vlm, drop = FALSE],
            quote = FALSE, digits = digits)
    }
  }

  invisible(NULL)
}


setMethod("summary", "vlm",
          function(object, ...)
          summaryvlm(object, ...))




setMethod("show", "summary.vlm",
          function(object)
          show.summary.vlm(object))



