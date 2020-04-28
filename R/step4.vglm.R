# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.














step4vglm <-
  function (object, scope,  # scale = 0,
            direction = c("both", "backward", "forward"),
            trace = 1, keep = NULL, steps = 1000, k = 2,
            ...) {
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    res <- if (is.null(dev)) extractAIC(x, k = 0)[2L] else dev
    res
  }
  cut.string <- function(string) {
    if (length(string) > 1L) 
      string[-1L] <- paste0("\n", string[-1L])
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE),
          c(nr, nc), list(namr, namc))
  }
  step.results <- function(models, fit, object) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    AIC <- sapply(models, "[[", "AIC")
 heading <- c("Stepwise Model Path \nAnalysis of Deviance ",
              "Table", 
              "\nInitial Model:", deparse(formula(object)),
              "\nFinal Model:", 
              deparse(formula(fit)), "\n")
 aod <- data.frame(Step = I(change),
                   Df = ddf,
                   Deviance = dd, 
                   `Resid. Df` = rdf,
                   `Resid. Dev` = rd,
                   AIC = AIC, 
                   check.names = FALSE)
    attr(aod, "heading") <- heading
    fit@post$anova <- aod  # fit$anova <- aod
    fit
  }  # step.results

  Terms <- terms(object)
  object@call$formula <- object@misc$formula <- Terms
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- is.element(direction, c("both", "backward"))
  forward  <- is.element(direction, c("both",  "forward"))
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) forward <- FALSE
  } else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object, fdrop)), "factors") else
      numeric()
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object, fadd)), "factors")
    } else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep)) 
    keep.list <- vector("list", steps)
  n.lm <- nobs(object, type = "lm")
  n.vlm <- nobs(object, type = "vlm")
  fit <- object
  bAIC <- extractAIC(fit, k = k, ...)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  if (is.na(bAIC)) 
    stop("AIC is not defined for this model, so 'step4' ",
         "cannot proceed")
  if (bAIC == -Inf)
    stop("AIC is -infinity for this model, so 'step4' ",
         "cannot proceed")
  nm <- 1
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n",
        cut.string(deparse(formula(fit))), 
        "\n\n", sep = "")
    flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit),
                       df.resid = n.vlm - edf,
                       change = "", AIC = bAIC)
  if (!is.null(keep)) 
    keep.list[[nm]] <- keep(fit, bAIC)

  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL

    if (backward && length(scope$drop)) {
      aod <- drop1(fit, scope$drop, k = k, ...)  # trace = trace,
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L]))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }  # if (backward && length(scope$drop))

    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- add1(fit, scope$add, k = k, ...)  # trace = trace,
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L]))
        aod <- if (is.null(aod)) aodf else
                 rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break

      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      oo <- order(aod[, nc])
      if (trace) 
        print(aod[oo, ])
      if (oo[1L] == 1) 
        break
      change <- rownames(aod)[oo[1L]]
    }  # if (is.null(change))

    
    fit <- update.default(fit, paste("~ .", change),
                          evaluate = FALSE)  # update()


    fit <- eval.parent(fit)
    nnew <- nobs(fit, type = "vlm")  # use.fallback = TRUE
    if (all(is.finite(c(n.vlm, nnew))) && nnew != n.vlm) 
      stop("number of rows in use has changed: ",
           "remove missing values?")


    Terms <- terms(fit)
    bAIC <- extractAIC(fit, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      flush.console()
    }
    if (bAIC >= AIC + 1e-07) 
      break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit),
                         df.resid = n.vlm - edf,
                         change = change, AIC = bAIC)
    if (!is.null(keep)) 
      keep.list[[nm]] <- keep(fit, bAIC)
  }  # while (steps > 0)
  
  if (!is.null(keep)) 
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object)
}  # step4vglm





if (!isGeneric("step4"))
  setGeneric("step4", function(object, ...)
             standardGeneric("step4"),
             package = "VGAM")


setMethod("step4", "vglm",
          function(object, ...)
          step4vglm(object, ...))














