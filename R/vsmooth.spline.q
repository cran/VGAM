# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.





if(!exists("is.R")) is.R <- function()
    exists("version") && !is.null(version$language) && version$language=="R"

setClass("vsmooth.spline.fit", representation(
      "Bcoefficients"= "matrix",
      "knots"        = "numeric",
      "xmin"         = "numeric",
      "xmax"         = "numeric"))

setClass("vsmooth.spline", representation(
      "call"         = "call",
      "constraints"  = "list",
      "df"           = "numeric",
      "nlfit"        = "vsmooth.spline.fit",  # is the nonlinear component
      "lev"          = "matrix",
      "lfit"         = "vlm",     # 6/6/02: was "vlm.wfit"; is the linear component
      "spar"         = "numeric",
      "var"          = "matrix",
      "w"            = "matrix",
      "x"            = "numeric",
      "y"            = "matrix",
      "yin"          = "matrix"))


setMethod("coefficients", signature(object="vsmooth.spline"),
          function(object, ...)
          coefvsmooth.spline(object, ...))
setMethod("coef", signature(object="vsmooth.spline"),
          function(object, ...)
          coefvsmooth.spline(object, ...))

setMethod("coefficients", signature(object="vsmooth.spline.fit"),
          function(object, ...)
          coefvsmooth.spline.fit(object, ...))
setMethod("coef", signature(object="vsmooth.spline.fit"),
          function(object, ...)
          coefvsmooth.spline.fit(object, ...))

setMethod("fitted.values", signature(object="vsmooth.spline"),
          function(object, ...)
          fittedvsmooth.spline(object, ...))
setMethod("fitted", signature(object="vsmooth.spline"),
          function(object, ...)
          fittedvsmooth.spline(object, ...))

setMethod("residuals", signature(object="vsmooth.spline"),
          function(object, ...)
          residvsmooth.spline(object, ...))
setMethod("resid", signature(object="vsmooth.spline"),
          function(object, ...)
          residvsmooth.spline(object, ...))

setMethod("predict", signature(object="vsmooth.spline"),
          function(object, ...)
          predictvsmooth.spline(object, ...))
setMethod("print", "vsmooth.spline",
         function(x, ...)
         invisible(printvsmooth.spline(x, ...)))
setMethod("show",  "vsmooth.spline",
          function(object)
          printvsmooth.spline(object))
setMethod("plot", "vsmooth.spline",
          function(x, y, ...) {
          if(!missing(y)) stop("can't process the \"y\" argument")
          invisible(plotvsmooth.spline(x, ...))})
setMethod("predict",  "vsmooth.spline.fit",
          function(object, ...)
          predictvsmooth.spline.fit(object, ...))


vsmooth.spline <- function(x, y, w, df=rep(5,M), spar=NULL, # rep(0,M),
                      all.knots=FALSE, 
                      iconstraint=diag(M),
                      xconstraint=diag(M),
                      constraints=list("(Intercepts)"=diag(M), x=diag(M)),
                      tol.nl=0.01, var.arg=FALSE,
                      scale.w=TRUE,
                      nk=NULL)
{


    if(var.arg) {
        warning("@var will be returned, but no use will be made of it") 
    }


    missing.constraints <- missing(constraints)

    if(!(missing.spar <- missing(spar)) && !missing(df))
        stop("can't specify both spar and df")


    my.call <- match.call()
    if(missing(y)) {
        if(is.list(x)) {
            if(any(is.na(match(c("x", "y"), names(x)))))
                stop("cannot find x and y in list")
            y <- x$y
            x <- x$x
        } else if(is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
        } else if(is.matrix(x)) {
            y <- x[,-1]
            x <- x[,1]
        } else {
            y <- x
            x <- time(x)
        }
    }

    n <- length(x)
    y <- as.matrix(y)
    ny2 <- dimnames(y)[[2]]  # NULL if vector 
    M <- ncol(y)
    if(n != nrow(y))
        stop("lengths of x and y must match")

    if(any(is.na(x)) || any(is.na(y)))
        stop("NAs not allowed in x or y")

    if(missing(w)) {
        w <- matrix(1, n, M)
    } else {
        if(any(is.na(w)))
            stop("NAs not allowed in w")

        w <- as.matrix(w)

        if(nrow(y) != nrow(w) || ncol(w)>M*(M+1)/2)
            stop("w and y don't match")

        if(scale.w)
            w <- w / mean(w[,1:M])    # 'Average' value is 1
    }
    dimw <- ncol(w)

    if(missing.constraints)
        constraints <- list("(Intercepts)"=eval(iconstraint),
                            x=eval(xconstraint))
    constraints <- eval(constraints)
    if(is.matrix(constraints))
       constraints <- list("(Intercepts)"=constraints, x=constraints)
    if(!is.list(constraints) || length(constraints)!=2)
        stop("constraints must equal a list (of length 2) or a matrix")
    for(i in 1:2) 
        if(!is.numeric(constraints[[i]]) || !is.matrix(constraints[[i]]) || 
           nrow(constraints[[i]])!=M || ncol(constraints[[i]])>M)
            stop("something wrong with the constraints")
    names(constraints) <- c("(Intercepts)", "x")


    sx <- unique(sort(as.vector(x)))
    o <- match(x, sx)             # sx[o]==x
    nef <- length(sx)
    if(nef < 7)
        stop("not enough unique x values (need 7 or more)")


    index <- iam(NA, NA, M, both=TRUE, diagonal=TRUE)
    template1 <- template2 <- matrix(0, nef, M)  # Must have M columns 
    ncb <- M
    dimu <- dimw # 10/1/00; was M*(M+1)/2

    collaps <- dotFortran(name="vsuff9",
                as.integer(n), as.integer(nef), as.integer(o),
                as.double(x), as.double(y), as.double(w),
                xbar=double(nef), ybar=as.double(template1), wbar=double(nef*dimu),
                     uwbar=as.double(0), wz=as.double(template2), 
                as.integer(M), dimw=as.integer(dimw), dimu=as.integer(dimu),
                     as.integer(index$row), as.integer(index$col),
                double(M*(M+1)), double(ncb*(ncb+1)),
                as.double(diag(M)), as.integer(M), 
                triv=as.integer(1), wuwbar=as.integer(0), ok=as.integer(0))


    if(collaps$ok != 1)
       stop("some non-positive-definite weight matrices detected in \"vsuff9\"")
    dim(collaps$ybar) <- dim(collaps$wz) <- c(nef, M)


    if(FALSE) {
    } else {
        yin = collaps$ybar   # Includes both linear and nonlinear parts 
        junk.frame = data.frame(x=collaps$xbar, yin = yin)
        x = collaps$xbar  # Warning: From now on "x" is no longer the original x 

        lfit = vlm(yin ~ 1 + x,
                   constraints = constraints,
                   save.weight=FALSE, qr=FALSE, x=FALSE, y=FALSE,
                   smart = FALSE,
                   weight=matrix(collaps$wbar, nrow=nrow(yin), byrow=FALSE))
    }

    ncb <- ncol(constraints[[2]])    # Of x and not of the intercept
    spar <- if(length(spar)) rep(spar, length=ncb) else rep(0, length=ncb)
    df <- rep(df, length=ncb)

    if(!missing.spar) {
        ispar <- 1
        if(any(spar <= 0) || !is.numeric(spar))
            stop("not allowed non-positive or non-numeric smoothing parameters")


        nonlin <- (spar != Inf)
    } else {
        ispar <- 0
        if(!is.numeric(df) || any(df < 2 | df > nef))
            stop(paste("you must supply 2 <= df <=", nef))
        if(tol.nl <= 0) stop("bad value for tol.nl")
        nonlin <- abs(df-2) > tol.nl
    }


    if(all(!nonlin)) {

        junk.fill = new("vsmooth.spline.fit",
                        "Bcoefficients"= matrix(as.numeric(NA), 1, 1),
                        "knots"        = numeric(0),
                        "xmin"         = numeric(0),
                        "xmax"         = numeric(0)) # 8/11/03
        object =
        new("vsmooth.spline",
           "call"         = my.call,
           "constraints"  = constraints,
           "df"           = if(ispar==0) df else rep(2, length(spar)),
           "lfit"         = lfit,
           "nlfit"        = junk.fill,
           "spar"         = if(ispar==1) spar else rep(Inf, length(df)),
           "w"            = as.matrix(collaps$wbar),
           "x"            = sx,
           "y"            = lfit@fitted.values,
           "yin"          = yin)

    
        return(object)
    }
    


    xbar <- (sx - sx[1]) / (sx[nef] - sx[1])
    noround = TRUE   # Improvement 3/8/02
    if(all.knots) {
        if(noround) {
            knot = valid.vknotl2(c(rep(xbar[1], 3), xbar, rep(xbar[nef], 3)))
        } else { 
            knot <- c(rep(xbar[1], 3), xbar, rep(xbar[nef], 3))
        }
        if(length(nk)) warning("overriding nk by all.knots=TRUE")
        nk <- length(knot) - 4     # No longer nef + 2
    } else {
        chosen = length(nk)
        if(chosen && (nk > nef+2 || nk <= 5))
            stop("bad value for nk")
        if(!chosen) nk = 0
        knot.list <- dotFortran(name="vknotl2", as.double(xbar), as.integer(nef),
                              knot=double(nef+6), k=as.integer(nk+4),
                              chosen=as.integer(chosen))
        if(noround) {
            knot = valid.vknotl2(knot.list$knot[1:(knot.list$k)])
            knot.list$k = length(knot)
        } else {
            knot <- knot.list$knot[1:(knot.list$k)]
        }
        nk <- knot.list$k - 4
    }
    if(nk <= 5) stop("not enough distinct knots found")

    conmat <- (constraints[[2]])[,nonlin,drop=FALSE]
    ncb <- sum(nonlin)
    trivc <- trivial.constraints(conmat)
    resmat <- collaps$ybar - lfit@fitted.values     # nef by M
    spar.nl <- spar[nonlin]
    df.nl <- df[nonlin]

    edimu <- if(trivc) dimw else max(ncb*(ncb+1)/2, dimw) # for wbar's size
    dimu <- if(trivc) dimw else ncb*(ncb+1)/2
    o <- 1:nef   # Already sorted

    collaps <- dotFortran(name="vsuff9",
                as.integer(nef), as.integer(nef), as.integer(o),
                as.double(collaps$xbar), as.double(resmat), as.double(collaps$wbar),
                xbar=double(nef), ybar=as.double(template1),
                    wbar=double(nef*edimu), uwbar=as.double(0), wz=as.double(template2),
                M=as.integer(M), dimw=as.integer(dimw), dimu=as.integer(dimu),
                    as.integer(index$row), as.integer(index$col),
                double(M*(M+1)), double(ncb*(ncb+1)),
                as.double(conmat), as.integer(ncb), 
                as.integer(trivc), wuwbar=as.integer(0), ok=as.integer(0))
    if(collaps$ok != 1)
       stop("some non-positive-definite weight matrices detected in \"vsuff9\"")

    dim(collaps$ybar) <- dim(collaps$wz) <- c(nef, M)
    collaps$ybar = collaps$ybar[,1:ncb,drop=FALSE]
    collaps$wz   = collaps$wz[,1:ncb,drop=FALSE]
    dim(collaps$wbar) <- c(nef, edimu)


    ldk = 3 * ncb + 1     # 10/7/02; Previously 4 * ncb
    lev <- if(ncb > 1) matrix(0, nef, ncb) else rep(0, nef)
    varmat <- if(var.arg) {if(ncb > 1) matrix(0, nef, ncb) else
                           rep(0, nef)} else double(1)
    index <- iam(NA, NA, ncb, both=TRUE, diagonal=TRUE)
    dimwbar <- if(trivc) dimw else ncb*(ncb+1)/2

    vsplin <- dotFortran(name="vsplin",
                     xs=as.double(xbar),  wz=as.double(collaps$wz), 
                     w=as.double(collaps$wbar), n=as.integer(nef), 
                     xknot=as.double(knot),
                     nk=as.integer(nk), as.integer(ldk),
                     M=as.integer(ncb), dimw=as.integer(dimwbar),
                     as.integer(index$row), as.integer(index$col),
                     wkmm=double(ncb*ncb*16), spar.nl=as.double(spar.nl), 
                     info=integer(1), fv=double(nef*ncb), Bcoef=double(nk*ncb),
                     hs=double(ldk*nk*ncb), btwy=double(ncb*nk),
                     sgdub=double(nk * max(4,ncb)),
                     var=as.double(varmat), ifvar=as.integer(var.arg),
                     bmb=double(ncb*ncb),
                     lev=as.double(lev),
                     as.double(df.nl), 
                     scrtch=double(min((17+nk)*nk, nk*17+1)),
                     ier=as.integer(0),
                     truen=as.integer(nef))


    if(vsplin$ier != 0) {
        cat("vsplin$ier ==", vsplin$ier, "\n")
        stop("something gone wrong in \"vsplin\"")
    }
    if(vsplin$info != 0)
        stop(paste("leading minor of order", vsplin$info,
                   "is not positive definite"))

    dim(vsplin$lev) <- c(nef, ncb)   # A matrix even when ncb==1
    if(ncb > 1) {
        dim(vsplin$fv) <- c(nef, ncb)
        if(var.arg)
            dim(vsplin$var) <- c(nef, ncb)
    }

    df.nl <- apply(vsplin$lev, 2, sum)  # Actual EDF used 


    fv <- lfit@fitted.values + vsplin$fv %*% t(conmat)
    if(M > 1)
        dimnames(fv) <- list(NULL, ny2)

    df[!nonlin] = 2
    df[ nonlin] = df.nl
    if(ispar==0) {
        spar[!nonlin] = Inf
        spar[ nonlin] = vsplin$spar.nl   # Actually used
    }

    fit.object = new("vsmooth.spline.fit",
                     "Bcoefficients"  = matrix(vsplin$Bcoef, nrow=nk, ncol=ncb),
                     "knots"          = knot,
                     "xmax"           = sx[nef],
                     "xmin"           = sx[1])
 
    object =
    new("vsmooth.spline",
        "call"         = my.call,
        "constraints"  = constraints,
        "df"           = df,
        "nlfit"        = fit.object,
        "lev"          = vsplin$lev,
        "lfit"         = lfit,
        "spar"         = spar,   # if(ispar==1) spar else vsplin$spar,
        "w"            = collaps$wbar,
        "x"            = sx,
        "y"            = fv, 
        "yin"          = yin)

    if(var.arg) 
        object@var = vsplin$var 

    object
}


printvsmooth.spline <- function(x, ...)
{
    if(!is.null(cl <- x@call)) {
        cat("Call:\n")
        dput(cl)
    }

    ncb <- if(length(x@nlfit)) ncol(x@nlfit@Bcoefficients) else NULL
    cat("\nSmoothing Parameter (Spar):", 
        if(length(ncb) && ncb==1) format(x@spar) else
            paste(format(x@spar), collapse=", "), "\n")

    cat("\nEquivalent Degrees of Freedom (Df):", 
        if(length(ncb) && ncb==1) format(x@df) else
            paste(format(x@df), collapse=", "), "\n")

    if(!all(trivial.constraints(x@constraints))) {
        cat("\nConstraint matrices:\n")
        print(x@constraints)
    }

    invisible(x)
}


coefvsmooth.spline = function(object, matrix=FALSE, ...) {
        list(lfit=coef(object@lfit, matrix=matrix), nlfit=coef(object@nlfit))
}


coefvsmooth.spline.fit = function(object, ...) {
    object@Bcoefficients 
}


fittedvsmooth.spline = function(object, ...) {
    object@y
}

residvsmooth.spline = function(object, ...) {
    as.matrix(object@yin - object@y)
}



plotvsmooth.spline <- function(x, xlab="x", ylab="", points=TRUE, 
                                pcol=par()$col, pcex=par()$cex,
                                pch=par()$pch,
                                lcol=par()$col, lwd=par()$lwd, lty=par()$lty, 
                                add=FALSE, ...)
{
    M = ncol(x@y)
    pcol = rep(pcol, length=M)
    pcex = rep(pcex, length=M)
    pch = rep(pch, length=M)
    lcol = rep(lcol, length=M)
    lwd = rep(lwd, length=M)
    lty = rep(lty, length=M)
    if(!add)
        matplot(x@x, x@yin, type="n", xlab=xlab, ylab=ylab, ...)
    for(i in 1:ncol(x@y)) {
        if(points)
            points(x@x, x@yin[,i], col=pcol[i], pch=pch[i], cex=pcex[i])
        lines(x@x, x@y[,i], col=lcol[i], lwd=lwd[i], lty=lty[i])
    }
    invisible(x)
}



predictvsmooth.spline <- function(object, x, deriv=0, se.fit=FALSE)
{
    if(se.fit)
        warning("se.fit=TRUE is not currently implemented. Using se.fit=FALSE")

    lfit <- object@lfit     # Linear part of the vector spline
    nlfit <- object@nlfit   # Nonlinear part of the vector spline

    if(missing(x)) {
        if(deriv==0) {
            return(list(x=object@x, y=object@y))
        } else {
            x <- object@x
            return(Recall(object, x, deriv))
        }

    }

    mat.coef = coef(lfit, matrix=TRUE)
    coeflfit <- t(mat.coef)   # M x p now
    M <- nrow(coeflfit) # if(is.matrix(object@y)) ncol(object@y) else 1 

    pred <- if(deriv==0) predict(lfit, data.frame(x=x)) else 
            if(deriv==1) matrix(coeflfit[,2], length(x), M, byrow=TRUE) else 
                  matrix(0, length(x), M)
    if(!length(nlfit@knots))
        return(list(x=x, y=pred))


    nonlin <- (object@spar != Inf)

    conmat <- if(!length(lfit@constraints)) diag(M) else lfit@constraints[[2]]
    conmat <- conmat[,nonlin,drop=FALSE] # Of nonlinear functions

    list(x=x, y=pred + predict(nlfit, x, deriv)$y %*% t(conmat))
}


predictvsmooth.spline.fit <- function(object, x, deriv=0)
{
    nk = nrow(object@Bcoefficients)
    drangex <- object@xmax - object@xmin
    if(missing(x))
        x <- seq(from=object@xmin, to=object@xmax, length=nk-4)

    xs <- as.double((x - object@xmin)/drangex)

    bad.left <- xs <  0
    bad.right <- xs >  1
    good <- !(bad.left | bad.right)

    ncb <- ncol(object@Bcoefficients)
    y <- matrix(as.numeric(NA), length(xs), ncb)
    if(any(good)) {
        ngood <- sum(good)
        junk <- dotFortran(name="vbvs", as.integer(ngood),
            as.double(object@knots), as.double(object@Bcoefficients),
            as.integer(nk),
            as.double(xs[good]), s=double(ngood*ncb),
            as.integer(deriv), as.integer(ncb))
        y[good,] <- junk$s

        if(TRUE && deriv > 1) {
            edges <- xs <= 0 | xs >= 1   # Zero the edges & beyond explicitly
            y[edges,] <- 0
        }
   }
    if(any(!good)) {
        xrange <- c(object@xmin, object@xmax)
        if(deriv == 0) {
            end.object <- Recall(object, xrange)$y
            end.slopes <- Recall(object, xrange, 1)$y * drangex

            if(any(bad.left))
                y[bad.left,] = rep(end.object[1,], rep(sum(bad.left), ncb)) +
                               rep(end.slopes[1,], rep(sum(bad.left), ncb)) *
                               xs[bad.left]
            if(any(bad.right))
                y[bad.right,] = rep(end.object[2,], rep(sum(bad.right), ncb)) +
                                rep(end.slopes[2,], rep(sum(bad.right), ncb)) *
                                (xs[bad.right] - 1)
        } else if(deriv == 1) {
            end.slopes <- Recall(object, xrange, 1)$y * drangex
            y[bad.left,] <- rep(end.slopes[1,], rep(sum(bad.left), ncb)) 
            y[bad.right,] <- rep(end.slopes[2,], rep(sum(bad.right), ncb)) 
        } else
            y[!good,] <- 0
    }
    if(deriv > 0)
        y <- y / (drangex^deriv)
    list(x=x, y=y)
}


valid.vknotl2 = function(knot, tol=1/1000) {

    junk = dotFortran(name="pknotl2", knot=as.double(knot), as.integer(length(knot)),
                    keep=integer(length(knot)), as.double(tol))
    keep = as.logical(junk$keep)
    knot = junk$knot[keep]
    if(length(knot) <= 11)
        stop("too few (distinct) knots")
    knot
}





