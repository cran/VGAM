Build.terms <- function(x, coefs, cov = NULL, assign, collapse = TRUE)
{

    cov.true <- !is.null(cov)
    if(collapse) {
        fit <- drop(x %*% coefs)
        if(cov.true) {
            var <- ((x %*% cov) * x) %*% rep(1, length(coefs))
            list(fitted.values = fit, se.fit = drop(sqrt(var)))
        }
        else fit
    } else {
        constant <- attr(x, "constant")
        if(!is.null(constant))
            constant <- sum(constant * coefs)
        if(missing(assign))
            assign <- attr(x, "assign")
        if(is.null(assign))
            stop("Need an 'assign' list")
        fit <- array(0, c(nrow(x), length(assign)), list(dimnames(x)[[1
            ]], names(assign)))
        if(cov.true)
            se <- fit
        TL <- sapply(assign, length)
        simple <- TL == 1
        complex <- TL > 1
        if(any(simple)) {
            asss <- unlist(assign[simple])
            ones <- rep(1, nrow(x))
            fit[, simple] <- x[, asss] * outer(ones, coefs[asss])
            if(cov.true)
                se[, simple] <- abs(x[, asss]) * outer(ones, 
                  sqrt(diag(cov))[asss])
        }
        if(any(complex)) {
            assign <- assign[complex]
            for(term in names(assign)) {
                TT <- assign[[term]]
                xt <- x[, TT]
                fit[, term] <- xt %*% coefs[TT]
                if(cov.true)
                  se[, term] <- sqrt(drop(((xt %*% cov[TT, TT]) *
                    xt) %*% rep(1, length(TT))))
            }
        }
        attr(fit, "constant") <- constant
        if(is.null(cov)) fit else list(fitted.values = fit, se.fit = se)
    }
}

