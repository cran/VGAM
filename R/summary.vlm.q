# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.









summaryvlm <- function(object, correlation = FALSE, dispersion = NULL) {


  if (is.logical(object@misc$BFGS) && object@misc$BFGS)
      warning(paste("the estimated variance-covariance matrix is",
         "usually inaccurate as the working weight matrices are a",
         "crude BFGS quasi-Newton approximation"))

  M <- object@misc$M
  n <- object@misc$n
  nrow_X_vlm <- object@misc$nrow_X_vlm 
  ncol_X_vlm <- object@misc$ncol_X_vlm   # May be NULL for CQO objects

  coef <- object@coefficients
  cnames <- names(coef)
  presid = residualsvlm(object, type = "pearson") # NULL if pooled.weight

  if (any(is.na(coef))) {
    warning(paste("Some NAs in the coefficients---no summary",
                  " provided; returning object\n"))
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
      dispersion <- if (!length(object@rss)) {
          stop("object@rss is empty")
      } else {
          object@rss / object@df.residual
      }
      object@misc$estimated.dispersion <- TRUE
  } else {
      if (is.numeric(object@misc$dispersion) &&
         object@misc$dispersion != dispersion)
          warning("overriding the value of object@misc$dispersion")
      object@misc$estimated.dispersion <- FALSE
  }
  sigma <- dispersion^0.5     # Can be a vector 

  if (is.Numeric(ncol_X_vlm)) {
    R <- object@R

    if (ncol_X_vlm < max(dim(R)))
      stop("R is rank deficient")

    rinv = diag(ncol_X_vlm)
    rinv = backsolve(R, rinv)
    rowlen = drop(((rinv^2) %*% rep(1, ncol_X_vlm))^0.5)
    covun = rinv %*% t(rinv)
    dimnames(covun) <- list(cnames, cnames)
  }
  coef <- matrix(rep(coef, 3), ncol=3)
  dimnames(coef) <- list(cnames, c("Estimate", "Std. Error", "z value"))
  if (length(sigma) == 1 && is.Numeric(ncol_X_vlm)) {
    coef[, 2] <- rowlen %o% sigma      # Fails here when sigma is a vector 
    coef[, 3] <- coef[, 1] / coef[, 2]
  } else {
    coef[,1] = coef[,2] = coef[,3] = NA
  }
  if (correlation) {
    correl <- covun * outer(1 / rowlen, 1 / rowlen)
    dimnames(correl) <- list(cnames, cnames)
  } else {
    correl <- matrix(0, 0, 0)  # was NULL, but now a special matrix
  }




  answer <-
  new("summary.vlm",
      object,
      coef3 = coef, 
      correlation = correl,
      df = c(ncol_X_vlm, rdf),
      sigma = sigma)

  if (is.Numeric(ncol_X_vlm)) answer@cov.unscaled = covun
  answer@dispersion = dispersion        # Overwrite this 

  if (length(presid))
    answer@pearson.resid = as.matrix(presid)


  answer
}




show.summary.vlm <- function(x, digits = NULL, quote = TRUE,
                             prefix = "") {


    M <- x@misc$M 
    coef3 <- x@coef3  # ficients
    correl <- x@correlation

    if (is.null(digits)) {
        digits <- options()$digits
    } else {
        old.digits <- options(digits = digits)
        on.exit(options(old.digits))
    }

    cat("\nCall:\n")
    dput(x@call)

    presid <- x@pearson.resid
    rdf <- x@df[2]
    if (length(presid) && all(!is.na(presid))) {
        cat("\nPearson residuals:\n")
        if (rdf/M > 5) {
            rq <-  apply(as.matrix(presid), 2, quantile) # 5 x M
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                                 x@misc$predictors.names)
            print(t(rq), digits=digits)
        } else
        if (rdf > 0) {
            print(presid, digits=digits)
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
        UUU = paste(x@misc$predictors.names, collapse = ", ")
        UUU = x@misc$predictors.names
        cat("\nNames of responses:\n") 
        cat(UUU, fill = TRUE, sep = ", ")
    }


    if (!is.null(x@rss))
        cat("\nResidual Sum of Squares:", format(round(x@rss, digits)),
            "on", round(rdf, digits), "degrees of freedom\n")


    if (length(correl)) {
        ncol_X_vlm <- dim(correl)[2]
        if (ncol_X_vlm > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1, -ncol_X_vlm, drop = FALSE], quote = FALSE, digits=digits)
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



