# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.










lm <-
function (formula, data, subset, weights, na.action, method = "qr", 
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, smart=TRUE, ...) 
{
    ret.x <- x
    ret.y <- y

    if(smart) setup.smart("write")

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf$smart <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf) else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y)) else if (length(offset) != NROW(y))
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
            3) else numeric(0), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0, df.residual = if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) 
            z$fitted.values <- offset
    } else {
        x <- model.matrix(mt, mf, contrasts)
        z <- if (is.null(w)) 
            lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
                ...) else
            lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
            ...)
    }
    class(z) <- c(if (is.matrix(y)) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- x
    if (ret.y) 
        z$y <- y
    if (!qr) 
        z$qr <- NULL

    if(smart) {
        z$smart.prediction <- get.smart.prediction()
        wrapup.smart()
    } else
        z$smart.prediction <- list(smart.arg=FALSE)

    z
}
attr(lm, "smart") <- TRUE


predict.lm <-
function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    interval = c("none", "confidence", "prediction"), level = 0.95, 
    type = c("response", "terms"), terms = NULL, na.action = na.pass, 
    pred.var = res.var/weights, weights = 1, ...) 
{
    # Smart prediction: handle the prediction flaw
    if(is.smart(object) && length(object$smart.prediction)) {
        setup.smart("read", smart.prediction=object$smart.prediction)
    }

    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    } else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action, 
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        offset <- if (!is.null(off.num <- attr(tt, "offset"))) 
            eval(attr(tt, "variables")[[off.num + 1]], newdata) else 
        if (!is.null(object$offset)) 
            eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- object$qr$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset)) 
        predictor <- predictor + offset
    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata)) 
            warning("Predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- weights.default(object)
            if (!is.null(w)) {
                weights <- w
                warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) && 
            missing(pred.var)) 
            warning("Assuming constant prediction variance even though model fit is weighted\n")
        if (inherits(weights, "formula")) {
            if (length(weights) != 2) 
                stop("'weights' as formula should be one-sided")
            d <- if (missing(newdata) || is.null(newdata)) 
                model.frame(object) else newdata
            weights <- eval(weights[[2]], d, environment(weights))
        }
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
            r <- object$residuals
            w <- object$weights
            rss <- sum(if (is.null(w)) 
                r^2 else r^2 * w)
            df <- n - p
            rss/df
        } else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w)) 
                  qr.Q(object$qr)[, p1, drop = FALSE] else 
                X[, piv] %*% qr.solve(qr.R(object$qr)[p1, 
                  p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            } else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0
        if (hasintercept) 
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrix(object)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (se.fit || interval != "none") {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
            }
            if (hasintercept) 
                X <- sweep(X, 2, avx)
            unpiv <- rep.int(0, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq(1, nterms, length = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0] <- 0
                predictor[, i] <- if (any(iipiv) > 0) 
                  X[, iipiv, drop = FALSE] %*% beta[iipiv] else 0
                if (se.fit || interval != "none") 
                  ip[, i] <- if (any(iipiv) > 0) 
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii, 
                      , drop = FALSE])^2 %*% rep.int(res.var, 
                      p) else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit) 
                  ip <- ip[, terms, drop = FALSE]
            }
        } else {
            predictor <- ip <- matrix(0, n, 0)
        }
        attr(predictor, "constant") <- if (hasintercept) 
            termsconst else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
            prediction = sqrt(ip + pred.var))
        if (type != "terms") {
            predictor <- cbind(predictor, predictor + hwid %o% 
                c(1, -1))
            colnames(predictor) <- c("fit", "lwr", "upr")
        } else {
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none") 
        se <- sqrt(ip)
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit) 
            se <- napredict(na.act, se)
    }

    if(is.smart(object) && length(object$smart.prediction)) {
        wrapup.smart()
    }

    if (type == "terms" && interval != "none") {
        if (missing(newdata) && !is.null(na.act)) {
            lwr <- napredict(na.act, lwr)
            upr <- napredict(na.act, upr)
        }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr, 
            df = df, residual.scale = sqrt(res.var))
    } else if (se.fit) 
        list(fit = predictor, se.fit = se, df = df,
             residual.scale = sqrt(res.var)) else predictor
}
attr(predict.lm, "smart") <- TRUE


predict.glm <- 
function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, ...) 
{
    # Smart prediction: handle the prediction flaw
    if(is.smart(object) && length(object$smart.prediction)) {
        setup.smart("read", smart.prediction=object$smart.prediction)
    }

    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors, 
                response = object$fitted, terms = predict.lm(object, 
                  se.fit = se.fit, scale = 1, type = "terms", 
                  terms = terms))
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        } else {
            pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                type = ifelse(type == "link", "response", type), 
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    } else {
        if (inherits(object, "survreg")) 
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0) 
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predict.lm(object, newdata, se.fit, scale = residual.scale, 
            type = ifelse(type == "link", "response", type), 
            terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    if(is.smart(object) && length(object$smart.prediction)) {
        wrapup.smart()
    }
    pred
}
attr(predict.glm, "smart") <- TRUE


predict.mlm <- 
function (object, newdata, se.fit = FALSE, na.action = na.pass, 
    ...) 
{
    # Smart prediction: handle the prediction flaw
    if(is.smart(object) && length(object$smart.prediction)) {
        setup.smart("read", smart.prediction=object$smart.prediction)
    }

    if (missing(newdata)) 
        return(object$fitted)
    if (se.fit) 
        stop("the 'se.fit' argument is not yet implemented for \"mlm\" objects")
    if (missing(newdata)) {
        X <- model.matrix(object)
        offset <- object$offset
    } else {
        tt <- terms(object)
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action, 
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        offset <- if (!is.null(off.num <- attr(tt, "offset"))) 
            eval(attr(tt, "variables")[[off.num + 1]], newdata) else 
        if (!is.null(object$offset)) 
            eval(object$call$offset, newdata)
    }
    piv <- object$qr$pivot[seq(object$rank)]
    pred <- X[, piv, drop = FALSE] %*% object$coefficients[piv, 
        ]
    if (!is.null(offset)) 
        pred <- pred + offset

    if(is.smart(object) && length(object$smart.prediction)) {
        wrapup.smart()
    }

    if (inherits(object, "mlm")) pred else pred[, 1]
}
attr(predict.mlm, "smart") <- TRUE



glm <- 
function (formula, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    smart = TRUE, ...) 
{
    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) 
        data <- environment(formula)

    if(smart) setup.smart("write")

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$smart <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ", 
        method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(Y)) else 
        if (length(offset) != NROW(Y)) 
            stop(gettextf(
    "number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- glm.fit(x = X, y = Y, weights = weights, start = start, 
        etastart = etastart, mustart = mustart, offset = offset, 
        family = family, control = control, intercept = attr(mt, 
            "intercept") > 0)
    if (length(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, weights = weights, offset = offset, family = family, 
            control = control, intercept = TRUE)$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c("glm", "lm")

    if(smart) {
        fit$smart.prediction <- get.smart.prediction()
        wrapup.smart()
    } else
        fit$smart.prediction <- list(smart.arg=FALSE)

    fit
}
attr(glm, "smart") <- TRUE





smartpredenv = new.env()


smart.mode.is <- function(mode.arg=NULL) {
    if(!length(mode.arg)) {
        if(exists(".smart.prediction", env=smartpredenv)) {
            get(".smart.prediction.mode", env=smartpredenv)
        } else {
            "neutral"
        }
    } else {
        if(mode.arg != "neutral" && mode.arg != "read" && mode.arg != "write")
stop("argument \"mode.arg\" must be one of \"neutral\", \"read\" or \"write\"")
        if(exists(".smart.prediction", env=smartpredenv)) {
            get(".smart.prediction.mode", env=smartpredenv)==mode.arg
        } else {
            mode.arg=="neutral"
        }
    }
}


setup.smart <- function(mode.arg, smart.prediction=NULL, max.smart=30) {
    actual <- if(mode.arg=="write") vector("list", max.smart) else 
              if(mode.arg=="read") smart.prediction else
              stop("value of mode.arg unrecognized")

    wrapup.smart()  # make sure

    if(length(actual)) {
        # Double check that smart.prediction is not trivial (in "read" mode)
        # If it is trivial then ignore it. This saves testing whether 
        # length(object$smart.prediction) > 0 in the predict methods function


        assign(".smart.prediction", actual, envir = smartpredenv)
        assign(".smart.prediction.counter", 0, envir = smartpredenv)
        assign(".smart.prediction.mode", mode.arg, envir = smartpredenv)
        assign(".max.smart", max.smart, envir = smartpredenv)
        assign(".smart.prediction", actual, envir = smartpredenv)
    }
}

wrapup.smart <- function() {
    if(exists(".smart.prediction", envir = smartpredenv))
        rm(".smart.prediction", envir = smartpredenv)
    if(exists(".smart.prediction.counter", envir = smartpredenv))
        rm(".smart.prediction.counter", envir = smartpredenv)
    if(exists(".smart.prediction.mode", envir = smartpredenv))
        rm(".smart.prediction.mode", envir = smartpredenv)
    if(exists(".max.smart", envir = smartpredenv))
        rm(".max.smart", envir = smartpredenv)
}


get.smart.prediction <- function() {

    smart.prediction.counter <- get(".smart.prediction.counter",
        envir = smartpredenv)
    max.smart <- get(".max.smart", envir = smartpredenv)

    if(smart.prediction.counter > 0) {
        # Save this on the object for smart prediction later
        smart.prediction <- get(".smart.prediction", envir = smartpredenv)
        if(max.smart >= (smart.prediction.counter+1))
            for(i in max.smart:(smart.prediction.counter+1))
                smart.prediction[[i]] <- NULL
        smart.prediction
    } else 
        NULL
}


put.smart <- function(smart) {

    # Puts the info, if possible, in frame 1.
    # Does not returns whether it did it or not. 


        # Write the info to frame 0 as well
    max.smart <- get(".max.smart", envir = smartpredenv)
    smart.prediction.counter <- get(".smart.prediction.counter",
        envir = smartpredenv)
    smart.prediction <- get(".smart.prediction", envir = smartpredenv)
        smart.prediction.counter <- smart.prediction.counter + 1

        if(smart.prediction.counter > max.smart) {
            # if list is too small, make it larger
            max.smart <- max.smart + (inc.smart <- 10) # can change inc.smart
            smart.prediction <- c(smart.prediction, vector("list", inc.smart))
            assign(".max.smart", max.smart, envir = smartpredenv)
        }

        smart.prediction[[smart.prediction.counter]] <- smart
        assign(".smart.prediction", smart.prediction, envir = smartpredenv)
        assign(".smart.prediction.counter", smart.prediction.counter,
               envir = smartpredenv)
}


get.smart <- function() {
    # Returns one list component of information
    smart.prediction <- get(".smart.prediction", envir = smartpredenv)
    smart.prediction.counter <- get(".smart.prediction.counter",
        envir = smartpredenv)
    smart.prediction.counter <- smart.prediction.counter + 1
    assign(".smart.prediction.counter", smart.prediction.counter,
           envir = smartpredenv)
    smart <- smart.prediction[[smart.prediction.counter]]
    smart
}

smart.expression <- expression({

        # This expression only works if the first argument of the smart
        # function is "x", e.g., smartfun(x, ...)
        # Nb. .smart.match.call is the name of the smart function.

        smart  <- get.smart()
        assign(".smart.prediction.mode", "neutral", envir = smartpredenv)

        .smart.match.call = as.character(smart$match.call)
        smart$match.call = NULL  # Kill it off for the do.call 

        ans.smart <- do.call(.smart.match.call[1], c(list(x=x), smart))
        assign(".smart.prediction.mode", "read", envir = smartpredenv)

        ans.smart
})



is.smart <- function(object) {
    if(is.function(object)) {
        if(is.logical(a <- attr(object, "smart"))) a else FALSE
    } else {
        if(length(slotNames(object))) {
            if(length(object@smart.prediction) == 1 &&
                is.logical(object@smart.prediction$smart.arg))
            object@smart.prediction$smart.arg else
                any(slotNames(object) == "smart.prediction")
        } else {
            if(length(object$smart.prediction) == 1 &&
                is.logical(object$smart.prediction$smart.arg))
            object$smart.prediction$smart.arg else
            any(names(object) == "smart.prediction")
        }
    }
}




library(splines) 



bs <-
function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
    Boundary.knots = range(x)) 
{
    x <- x  # Evaluate x
    if(smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1]) | (or <- x > 
            Boundary.knots[2])
    } else outside <- FALSE
    ord <- 1 + (degree <- as.integer(degree))
    if (ord <= 1) 
        stop("'degree' must be integer >= 1")
    if (!missing(df) && missing(knots)) {
        nIknots <- df - ord + (1 - intercept)
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used  ", ord - 
                (1 - intercept))
        }
        knots <- if (nIknots > 0) {
            knots <- seq(from = 0, to = 1, length = nIknots + 
                2)[-c(1, nIknots + 2)]
            stats::quantile(x[!outside], knots)
        }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))
    if (any(outside)) {
        warning(
"some 'x' values beyond boundary knots may cause ill-conditioned bases")
        derivs <- 0:degree
        scalef <- gamma(1:ord)
        basis <- array(0, c(length(x), length(Aknots) - degree - 
            1))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1]
            xl <- cbind(1, outer(x[ol] - k.pivot, 1:degree, "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord, 
                derivs)$design
            basis[ol, ] <- xl %*% (tt/scalef)
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2]
            xr <- cbind(1, outer(x[or] - k.pivot, 1:degree, "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord, 
                derivs)$design
            basis[or, ] <- xr %*% (tt/scalef)
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                ord)$design
    } else basis <- spline.des(Aknots, x, ord)$design
    if (!intercept) 
        basis <- basis[, -1, drop = FALSE]
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1:n.col)
    a <- list(degree = degree, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("bs", "basis")

    if(smart.mode.is("write"))
        put.smart(list(df=df,
                       knots=knots,
                       degree=degree,
                       intercept=intercept,
                       Boundary.knots=Boundary.knots,
                       match.call=match.call()))

    basis
}
attr(bs, "smart") <- TRUE

ns <-
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) 
{
    x <- x  # Evaluate x
    if(smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1]) | (or <- x > 
            Boundary.knots[2])
    } else outside <- FALSE
    if (!missing(df) && missing(knots)) {
        nIknots <- df - 1 - intercept
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used ", 1 + intercept)
        }
        knots <- if (nIknots > 0) {
            knots <- seq(0, 1, length = nIknots + 2)[-c(1, nIknots + 
                2)]
            stats::quantile(x[!outside], knots)
        }
    } else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1]
            xl <- cbind(1, x[ol] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2), 4, c(0, 
                1))$design
            basis[ol, ] <- xl %*% tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2]
            xr <- cbind(1, x[or] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2), 4, c(0, 
                1))$design
            basis[or, ] <- xr %*% tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4)$design
    } else basis <- spline.des(Aknots, x, 4)$design
    const <- spline.des(Aknots, Boundary.knots, 4, c(2, 2))$design
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1:2), 
        drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1:n.col)
    a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("ns", "basis")

    if(smart.mode.is("write"))
        put.smart(list(df=df,
                       knots=knots,
                       intercept=intercept,
                       Boundary.knots=Boundary.knots,
                       match.call=match.call()))

    basis
}
attr(ns, "smart") <- TRUE




poly <-
function (x, ..., degree = 1, coefs = NULL, raw = FALSE) 
{
    x <- x  # Evaluate x
    if(!raw && smart.mode.is("read")) {
        smart <- get.smart()
        degree <- smart$degree
        coefs  <- smart$coefs
        raw  <- smart$raw
    }

    dots <- list(...)
    if (nd <- length(dots)) {
        if (nd == 1 && length(dots[[1]]) == 1) 
            degree <- dots[[1]] else 
        return(polym(x, ..., degree = degree, raw = raw))
    }
    if (is.matrix(x)) {
        m <- unclass(as.data.frame(cbind(x, ...)))
        return(do.call("polym", c(m, degree = degree, raw = raw)))
    }
    if (degree < 1) 
        stop("'degree' must be at least 1")

    # At prediction time x may be less than the degree
    if(smart.mode.is("write") || smart.mode.is("neutral"))
    if (degree >= length(x))
        stop("degree must be less than number of points")

    if (any(is.na(x))) 
        stop("missing values are not allowed in 'poly'")
    n <- degree + 1
    if (raw) {
        if (degree >= length(x)) 
            stop("'degree' must be less than number of points")
        Z <- outer(x, 1:degree, "^")
        colnames(Z) <- 1:degree
        attr(Z, "degree") <- 1:degree
        class(Z) <- c("poly", "matrix")
        return(Z)
    }
    if (is.null(coefs)) {
        if (degree >= length(x)) 
            stop("'degree' must be less than number of points")
        xbar <- mean(x)
        x <- x - xbar
        X <- outer(x, seq_len(n) - 1, "^")
        QR <- qr(X)
        z <- QR$qr
        z <- z * (row(z) == col(z))
        raw <- qr.qy(QR, z)
        norm2 <- colSums(raw^2)
        alpha <- (colSums(x * raw^2)/norm2 + xbar)[1:degree]
        Z <- raw/rep(sqrt(norm2), each = length(x))
        colnames(Z) <- 1:n - 1
        Z <- Z[, -1, drop = FALSE]
        attr(Z, "degree") <- 1:degree
        attr(Z, "coefs") <- list(alpha = alpha, norm2 = c(1, 
            norm2))
        class(Z) <- c("poly", "matrix")
    } else {
        alpha <- coefs$alpha
        norm2 <- coefs$norm2
        Z <- matrix(, length(x), n)
        Z[, 1] <- 1
        Z[, 2] <- x - alpha[1]
        if (degree > 1) 
            for (i in 2:degree) Z[, i + 1] <- (x - alpha[i]) * 
                Z[, i] - (norm2[i + 1]/norm2[i]) * Z[, i - 1]
        Z <- Z/rep(sqrt(norm2[-1]), each = length(x))
        colnames(Z) <- 0:degree
        Z <- Z[, -1, drop = FALSE]
        attr(Z, "degree") <- 1:degree
        attr(Z, "coefs") <- list(alpha = alpha, norm2 = norm2)
        class(Z) <- c("poly", "matrix")
    }

    if(smart.mode.is("write"))
        put.smart(list(degree=degree, coefs=attr(Z, "coefs"),
                       raw=FALSE,  # raw is changed above
                       match.call=match.call()))

    Z
}
attr(poly, "smart") <- TRUE


scale.default <-
function (x, center = TRUE, scale = TRUE) 
{
    x <- as.matrix(x)

    if(smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    nc <- ncol(x)
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(x, na.rm = TRUE)
            x <- sweep(x, 2, center)
        }
    } else if (is.numeric(center) && (length(center) == nc)) 
        x <- sweep(x, 2, center) else 
    stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
        if (scale) {
            f <- function(v) {
                v <- v[!is.na(v)]
                sqrt(sum(v^2)/max(1, length(v) - 1))
            }
            scale <- apply(x, 2, f)
            x <- sweep(x, 2, scale, "/")
        }
    } else if (is.numeric(scale) && length(scale) == nc) 
        x <- sweep(x, 2, scale, "/") else 
    stop("length of 'scale' must equal the number of columns of 'x'")
    if (is.numeric(center)) 
        attr(x, "scaled:center") <- center
    if (is.numeric(scale)) 
        attr(x, "scaled:scale") <- scale

    if(smart.mode.is("write")) {
        put.smart(list(center=center, scale=scale,
                       match.call=match.call()))
    }

    x
}
attr(scale.default, "smart") <- TRUE


attr(scale, "smart") <- TRUE





"my1" <- function(x, minx=min(x)) {

    x <- x   # Evaluate x

    if(smart.mode.is("read")) {
        smart  <- get.smart()
        minx <- smart$minx          # Overwrite its value 
    } else 
    if(smart.mode.is("write"))
        put.smart(list(minx=minx))

    (x-minx)^2
}
attr(my1, "smart") <- TRUE




"my2" <- function(x, minx=min(x)) {

    x <- x   # Evaluate x

    if(smart.mode.is("read")) {
        return(eval(smart.expression))
    } else 
    if(smart.mode.is("write"))
        put.smart(list(minx=minx, match.call=match.call()))

    (x-minx)^2
}

attr(my2, "smart") <- TRUE




"stdze1" <- function(x, center=TRUE, scale=TRUE) {

    x <- x  # Evaluate x

    if(!is.vector(x))
        stop("x must be a vector")

    if(smart.mode.is("read")) {
        smart  <- get.smart()
        return((x-smart$center)/smart$scale)
    }

    if(is.logical(center))
        center <- if(center) mean(x) else 0
    if(is.logical(scale))
        scale <- if(scale) sqrt(var(x)) else 1

    if(smart.mode.is("write"))
        put.smart(list(center=center,
                       scale=scale))
    # Normal use
    (x-center)/scale
}
attr(stdze1, "smart") <- TRUE

"stdze2" <- function(x, center=TRUE, scale=TRUE) {

    x <- x  # Evaluate x

    if(!is.vector(x))
        stop("x must be a vector")

    if(smart.mode.is("read")) {
        return(eval(smart.expression))
    }

    if(is.logical(center))
        center <- if(center) mean(x) else 0
    if(is.logical(scale))
        scale <- if(scale) sqrt(var(x)) else 1

    if(smart.mode.is("write"))
        put.smart(list(center=center,
                       scale=scale,
                       match.call=match.call()))

    (x-center)/scale
}
attr(stdze2, "smart") <- TRUE




