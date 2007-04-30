# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.








if(!exists("is.R")) is.R <- function()
    exists("version") && !is.null(version$language) && version$language=="R"



plotvgam <- function(x, newdata=NULL, y=NULL, residuals=NULL, rugplot=TRUE,
                 se= FALSE, scale=0, 
                 raw= TRUE, offset.arg=0, deriv.arg=0, overlay= FALSE,
                 type.residuals=c("deviance","working","pearson","response"),
                 plot.arg= TRUE, which.term=NULL, which.cf=NULL, 
                 control=plotvgam.control(...), 
                 ...)
{

    missing.control = missing(control)

    na.act = x@na.action
    x@na.action = list() # Don't want NAs returned from predict() or resid()

    missing.type.residuals = missing(type.residuals)
    if(mode(type.residuals) != "character" && mode(type.residuals) != "name")
        type.residuals <- as.character(substitute(type.residuals))
    if(!missing.type.residuals)
        type.residuals <- match.arg(type.residuals,
            c("deviance","working","pearson","response"))[1]


    if(!is.numeric(deriv.arg) || deriv.arg<0 ||
       deriv.arg!=round(deriv.arg) || length(deriv.arg)>1)
        stop("bad input for the deriv argument")

    if(se && deriv.arg>0) {
    warning("standard errors not available with derivatives. Setting se=FALSE")
        se = FALSE
    }

    preplot.object <- x@preplot
    if(!length(preplot.object)) {
        if(is.R()) {
            preplot.object <- preplotvgam(x, newdata=newdata,
                                           raw=raw, deriv=deriv.arg, se=se)
        } else {
            preplot.object <- preplotvgam(x, raw=raw, deriv=deriv.arg, se=se)
        }
    }

    x@preplot = preplot.object


    if(!is.null(residuals) && length(residuals) == 1) {
        if(residuals) {
            if(missing.type.residuals) {
                for(rtype in type.residuals)
                    if(!is.null(residuals <- resid(x, type=rtype))) break
            } else {
             residuals=resid(x,typ=type.residuals) #Get the prespecified type
                if(!length(residuals))
                    warning("residuals are NULL. Ignoring residuals=T")
            }
        } else {
            residuals <- NULL
        }
    }

    if(!missing.control) {
        control = if(is.R()) c(plotvgam.control(.include.dots= FALSE, ...),
                               control,
                               plotvgam.control(...)) else 
                             c(plotvgam.control(.include.dots= FALSE, ...),
                               control,
                               plotvgam.control(...))
    }

    x@post$plotvgam.control = control # Add it to the object 

    if(plot.arg)
        plotpreplotvgam(preplot.object, residuals=residuals, 
                        rugplot=rugplot, scale=scale, se=se,
                        offset.arg=offset.arg, deriv.arg=deriv.arg,
                        overlay=overlay, 
                        which.term=which.term, which.cf=which.cf, 
                        control=control)

    x@na.action = na.act  # Restore it's original value
    invisible(x)
}




ylim.scale <- function(ylim, scale=0) {
    if(length(ylim) != 2 || ylim[2] < ylim[1])
        stop("error in ylim")
    try <- ylim[2] - ylim[1]
    if(try > scale) ylim else
        c(ylim[1]+ylim[2]-scale, ylim[1]+ylim[2]+scale) / 2 
}






preplotvgam = function(object, newdata=NULL,
              terms=if(is.R()) attr((object@terms)$terms, "term.labels") else
                    v.labels.lm(object),
              raw= TRUE, deriv.arg=deriv.arg, se= FALSE)
{
    Terms <- terms(object)  # 11/8/03; object@terms$terms 
    aa <- attributes(Terms)
        Call <- object@call
    all.terms <- labels(Terms)
    xvars <- if(is.R()) parse(text=all.terms) else as.vector(Terms)

 
    if(is.R()) {
        names(xvars) <- all.terms
        terms <- sapply(terms, match.arg, all.terms)
    } else {
        names(xvars) <- all.terms
        terms <- match.arg(terms, all.terms)
    }

    Interactions <- aa$order > 1
    if(any(Interactions)) {
        if(is.R())
            stop("can't handle interactions") 
        all.terms <- all.terms[!Interactions]
        TM <- match(terms, all.terms, 0)
        if(!all(TM)) {
            terms <- terms[TM > 0]
            warning("No terms saved for \"a:b\" style interaction terms")
        }
    }

    if(is.R()) {
        xvars <- xvars[terms]
        xnames <- as.list(terms)
        names(xnames) <- terms
        modes <- sapply(xvars, mode)
        for(term in terms[modes != "name"]) {
            evars <- all.names(xvars[term], functions= FALSE, unique= TRUE)
            if(!length(evars))
                next
            xnames[[term]] <- evars
            evars <- parse(text=evars)
            if(length(evars) == 1)
                evars <- evars[[1]]
            else {
                evars <- c(as.name("list"), evars)
                mode(evars) <- "call"
            }
            xvars[[term]] <- evars
        }
    
    
        xvars <- c(as.name("list"), xvars)
        mode(xvars) <- "call"
        if(length(newdata)) {
            xvars <- eval(xvars, newdata)
        } else {
            if(!is.null(Call$subset) | !is.null(Call$na.action) |
               !is.null(options("na.action")[[1]])) {
                Rownames <- names(fitted(object))
                if(!(Rl <- length(Rownames)))
                    Rownames <- dimnames(fitted(object))[[1]]

                if(length(object@x) && !(Rl <- length(Rownames)))
                    Rownames <- (dimnames(object@x))[[1]]
                if(length(object@y) && !(Rl <- length(Rownames)))
                    Rownames <- (dimnames(object@y))[[1]]

                if(!(Rl <- length(Rownames)))
                    stop(paste("need to have names for fitted.values",
                               "when call has a subset or na.action argument"))

                form <- paste("~", unlist(xnames), collapse="+")
                Mcall <- c(as.name("model.frame"), list(formula =
                           terms(as.formula(form)),
                           subset = Rownames, na.action = function(x) x))
                mode(Mcall) <- "call"
                Mcall$data <- Call$data
                xvars <- eval(xvars, eval(Mcall))
            } else {
                ecall <- substitute(eval(expression(xvars)))
                ecall$local <- Call$data
                xvars <- eval(ecall)
            }
        }
    } else {
        xvars <- xvars[terms]
        xnames <- as.list(terms)
        names(xnames) <- terms
        modes <- sapply(xvars, mode)
        for(term in terms[modes != "name"]) {
            evars <- all.names(xvars[term], functions= FALSE, unique= TRUE)
            if(!length(evars))
                next
            xnames[[term]] <- evars
            evars <- parse(text=evars)
            if(length(evars) == 1)
                evars <- evars[[1]]
            else {
                evars <- c(as.name("list"), evars)
                mode(evars) <- "call"
            }
            xvars[[term]] <- evars
        }
    
        act.vars <- as.character(xvars)
    
        xvars <- c(as.name("list"), xvars)
        mode(xvars) <- "call"
        if(length(newdata)) {
            xvars <- eval(xvars, newdata)
        } else {
            if(!is.null(Call$subset) | !is.null(Call$na.action) |
               !is.null(options("na.action")[[1]])) {
                Rownames <- names(fitted(object))
                if(!(Rl <- length(Rownames)))
                    Rownames <- dimnames(fitted(object))[[1]]
                if(!(Rl <- length(Rownames)))
                    stop(paste("need to have names for fitted.values",
                               "when call has a subset or na.action argument"))
                Mcall <- c(as.name("model.frame"), list(formula=
                    terms.inner(parse(text=unlist(xnames))),
                    subset=Rownames, na.action=function(x) x))
                mode(Mcall) <- "call"
                Mcall$data <- Call$data
                xvars <- eval(xvars, eval(Mcall))
            } else {
                ecall <- substitute(eval(expression(xvars)))
                ecall$local <- Call$data
                xvars <- eval(ecall)
            }
        }
    }

    if(length(newdata)) {
        pred <- predict(object, newdata, type="terms",
                        raw=raw, se.fit=se, deriv.arg=deriv.arg)
    } else {
        pred <- predict(object, type="terms",
                        raw=raw, se.fit=se, deriv.arg=deriv.arg)
    }

    fits <- if(is.atomic(pred)) NULL else pred$fit
    se.fit <- if(is.atomic(pred)) NULL else pred$se.fit
    if(is.null(fits))
        fits <- pred
    fred <- attr(fits, "vterm.assign")   # NULL for M==1

    if(FALSE && is.R()) {
        xnames <- vector("list", length(fred))
        names(xnames) <- names(fred)
    }

    gamplot <- xnames

    if(FALSE && is.R()) {
        s.x = if(any(slotNames(object)=="s.xargument")) object@s.xargument else
              NULL
        n.s.x = names(s.x)
    }

    loop.var = if(is.R()) names(fred) else terms
    for(term in loop.var) {
        if(FALSE && is.R()) {
            useterm <- term
            if(length(n.s.x) && any(n.s.x == useterm))
                useterm <- s.x[useterm]
            innerx <- parse(text=useterm)
            innerx <- all.vars(innerx)
            if(length(innerx) == 0)
                warning(paste("couldn't extract variable from", useterm, "\n"))
            if(length(innerx) > 1) {
                warning(paste("using the first of \"", innerx,
                              "\" terms\n", sep=""))
                innerx <- innerx[1]
            }
        }

        if(FALSE && is.R()) {
            .VGAM.x <- if(length(newdata)) newdata[[innerx]] else {
                if(( is.R() && object@misc$dataname != "list") ||
                   (!is.R() && object@misc$dataname != "sys.parent")) {
                    mytext <- paste(object@misc$dataname,
                                    "[['", innerx, "']]", sep="")
                } else {
                    mytext <- innerx
                }
                getx <- parse(text=mytext)
                .VGAM.ans = if(exists(x=mytext, envir = .GlobalEnv))
                            eval(getx, envir = .GlobalEnv) else eval(getx)
                .VGAM.ans
            }
        } # else {

        .VGAM.x <- xvars[[term]]


        if(FALSE && is.R()) {
           class(.VGAM.x)=unique(c(class(.VGAM.x),data.class(unclass(.VGAM.x))))
        }

        myylab = if(all(substring(term, 1:nchar(term), 1:nchar(term)) != "("))
            paste("partial for", term) else term

        TT <- list(x = .VGAM.x,
                   y = fits[, if(is.null(fred)) term else fred[[term]]],
                   se.y = if(is.null(se.fit)) NULL else
                         se.fit[, if(is.null(fred)) term else fred[[term]]],
                   xlab = xnames[[term]],
                   ylab = myylab)
        class(TT) <- "preplotvgam"
        gamplot[[term]] <- TT
    }
    if(!is.R())
        class(gamplot) <- "preplotvgam"    # Commented out 8/6/02
    invisible(gamplot) 
}


if(!is.R())
v.labels.lm <- function(object, ...)
{
    TL <- terms(object)  # 11/8/03; object@terms$terms 
    if(!is.null(TL)) {
        TL <- attr(TL, "term.labels")
        TA <- object@assign
        if(!is.null(TA)) {
            TA <- names(TA)
            TL <- TL[match(TA, TL, 0.)]
        }
    }
    TL
}


plotvlm <- function(object, residuals=NULL, rugplot= FALSE, ...)
{
    stop("sorry, this function hasn't been written yet")
}


plotvglm <- function(x, residuals=NULL, smooths= FALSE,
                      rugplot= FALSE, id.n= FALSE, ...)
{
    stop("this function hasn't been written yet")  # zz

    M <- x@misc$M
    true.mu <- object@misc$true.mu
    response <- as.matrix(x@y)
    if(is.null(true.mu))
        true.mu <- T

    Residuals <- resid(x, type="deviance")
    if(!is.null(residuals))
    {
        if(length(residuals) == 1 && residuals)
            residuals <- Residuals else
            Residuals <- residuals
    }

    if(ncol(response)==1 && true.mu && !is.null(Residuals))
        invisible(NextMethod("plot")) else
    invisible(x)
}



if(!is.R()) jitter <- function(x, factor=1)
{

        z <- diff(range(x[!is.na(x)]))
        if(all(z==0))
            return(x)
        z <- factor * (z/50)
        x + runif(length(x),  - z, z)
}



plotpreplotvgam <- function(x, y=NULL, residuals=NULL,
                              rugplot= TRUE, se= FALSE, scale=0,
                              offset.arg=0, deriv.arg=0, overlay= FALSE, 
                              which.term=NULL, which.cf=NULL, 
                              control=NULL)
{
    listof <- inherits(x[[1]], "preplotvgam")
    if(listof) {
        TT <- names(x)
        if(is.null(which.term))
            which.term = TT  # Plot them all
        plot.no = 0
        for(i in TT) {
            plot.no = plot.no + 1 
            if((is.character(which.term) && any(which.term==i)) ||
               (is.numeric(which.term) && any(which.term==plot.no)))
                plotpreplotvgam(x[[i]], y=NULL, 
                                residuals, rugplot, se, scale, 
                                offset.arg=offset.arg,
                                deriv.arg=deriv.arg, overlay=overlay, 
                                which.cf=which.cf,
                                control=control)
        }
    } else {
        dummy <- function(residuals=NULL, rugplot= TRUE, se= FALSE, scale=0, 
                          offset.arg=0, deriv.arg=0, overlay= FALSE, 
                          which.cf=NULL, control=plotvgam.control(...))
            c(list(residuals=residuals, rugplot=rugplot, se=se, scale=scale,
                   offset.arg=offset.arg, deriv.arg=deriv.arg,
                   overlay=overlay,
                   which.cf=which.cf),
                   control)

        d <- dummy(residuals=residuals, rugplot=rugplot, se=se, scale=scale,
                   offset.arg=offset.arg, deriv.arg=deriv.arg,
                   overlay=overlay,
                   which.cf=which.cf, 
                   control=control)

        uniq.comps <- unique(c(names(x), names(d)))
        Call <- c(as.name("vplot"), c(d, x)[uniq.comps])
        mode(Call) <- "call"
        invisible(eval(Call))
    }
}


vplot.default <- function(x, y, se.y=NULL, xlab="", ylab="",
                          residuals=NULL, rugplot= FALSE,
                          scale=0, se= FALSE, 
                          offset.arg=0, deriv.arg=0, overlay= FALSE, 
                          which.cf=NULL, ...) {
    switch(data.class(x)[1],
           logical=vplot.factor(factor(x), y, se.y, xlab, ylab, residuals, 
                                rugplot, scale, se,
                                offset.arg=offset.arg, overlay=overlay, ...),
           if(is.numeric(x)) {
               vplot.numeric(as.vector(x), y, se.y, xlab, ylab, 
                             residuals, rugplot, scale, se,
                             offset.arg=offset.arg, overlay=overlay, ...)
           } else {
                   warning(paste("The \"x\" component of \"", ylab,
                                 "\" has class \"",
                           paste(class(x), collapse="\", \""), 
                                 "\"; no vplot() methods available", sep=""))
           }
    )
}



vplot.list <- function(x, y, se.y=NULL, xlab, ylab, 
                       residuals=NULL, rugplot= FALSE, scale=0, se= FALSE, 
                       offset.arg=0, deriv.arg=0, overlay= FALSE, 
                       which.cf=NULL, ...)
{

    if(is.numeric(x[[1]])) {
        vplot.numeric(x[[1]], y, se.y, xlab, ylab, 
                      residuals, rugplot, scale, se, 
                      offset.arg=offset.arg, deriv.arg=deriv.arg,
                      overlay=overlay, ...)
    } else 
        stop("this function hasn't been written yet") 
}


plotvgam.control = function(
                          which.cf=NULL,
                          xlim=NULL, ylim=NULL,
                          llty=par()$lty,
                          slty=if(is.R()) "dashed" else 3,
                          pcex=par()$cex,
                          pch=par()$pch,
                          pcol=par()$col,
                          lcol=par()$col,
                          rcol=par()$col,
                          scol=par()$col,
                          llwd=par()$lwd,
                          slwd=par()$lwd,
                          add.arg= FALSE,
                          one.at.a.time= FALSE, 
                          .include.dots= TRUE,
                          ...) {


    ans = 
    list(which.cf=which.cf,
         xlim=xlim, ylim=ylim,
         llty=llty, slty=slty,
         pcex=pcex, pch=pch,
         pcol=pcol, lcol=lcol, rcol=rcol, scol=scol,
         llwd=llwd, slwd=slwd,
         add.arg=add.arg,
         one.at.a.time=one.at.a.time)

    if(.include.dots) {
        c(list(...), ans)
    } else {
        default.vals = plotvgam.control()
        return.list = list()
        for(i in names(default.vals)) {
            replace.val = 
            if(is.R())
                !((length(ans[[i]]) == length(default.vals[[i]])) &&
                  (length(default.vals[[i]]) > 0) &&
                  (is.logical(all.equal(ans[[i]], default.vals[[i]]))) &&
                              all.equal(ans[[i]], default.vals[[i]])) else
                !((length(ans[[i]]) == length(default.vals[[i]])) &&
                  (length(default.vals[[i]]) > 0) &&
                  all(ans[[i]] == default.vals[[i]]))

            if(FALSE && replace.val) {
            }
            if(replace.val) 
                return.list[[i]] = ans[[i]]
        }
        if(length(return.list)) {
            names(return.list) = names(return.list)
            return.list
        } else NULL 
    }
}


vplot.numeric <- function(x, y, se.y=NULL, xlab, ylab,
                          residuals=NULL, rugplot= FALSE, se= FALSE, scale=0,
                          offset.arg=0, deriv.arg=0, overlay= FALSE,
                          which.cf=NULL,
                          xlim=NULL, ylim=NULL,

                          llty=par()$lty,
                          slty=if(is.R()) "dashed" else 3,
                          pcex=par()$cex,
                          pch=par()$pch,
                          pcol=par()$col,
                          lcol=par()$col,
                          rcol=par()$col,
                          scol=par()$col,
                          llwd=par()$lwd,
                          slwd=par()$lwd,
                          add.arg= FALSE,
                          one.at.a.time= FALSE, 
                          separator = ":",

                          ...)
{




    ylim0 <- ylim

    if(length(y)/length(x)  != round(length(y)/length(x)))
        stop("length of x and y do not seem to match")
    y <- as.matrix(y) 
    if(!length(which.cf))
        which.cf = 1:ncol(y)  # Added 7/8/04

    if(!is.null(se.y))
        se.y <- as.matrix(se.y)
    if(!is.null(se.y) && any(is.na(se.y)))
        se.y <- NULL

    if(!is.null(residuals))  {
        residuals <- as.matrix(residuals)
        if(ncol(residuals) != ncol(y)) {
            warning("ncol(residuals) != ncol(y) so residuals are not plotted")
            residuals <- NULL
        }
    }

    offset.arg <- matrix(offset.arg, nrow(y), ncol(y), byrow= TRUE)
    y <- y + offset.arg

    ylab <- add.hookey(ylab, deriv.arg)


    if(xmeanAdded <- (se && !is.null(se.y) &&
       all(substring(ylab, 1:nchar(ylab), 1:nchar(ylab)) != "("))) {
            x = c(x, mean(x))
            y = rbind(y, 0 * y[1,])
            se.y = rbind(se.y, 0 * se.y[1,])
            if(!is.null(residuals))
                residuals = rbind(residuals, NA*residuals[1,]) # NAs not plotted
    }

    ux <- unique(sort(x))
    o <- match(ux, x)
    uy <- y[o,,drop= FALSE]
    xlim <- range(xlim, ux)
    ylim <- range(ylim, uy[,which.cf], na.rm= TRUE)
    if(rugplot) {
        usex = if(xmeanAdded) x[-length(x)] else x
        jx <- jitter(usex[!is.na(usex)])
        xlim <- range(c(xlim, jx))
    }

    if(se && !is.null(se.y)) {
        se.upper <- uy[,which.cf] + 2 * se.y[o,which.cf,drop= FALSE]
        se.lower <- uy[,which.cf] - 2 * se.y[o,which.cf,drop= FALSE]
        ylim <- range(c(ylim, se.upper, se.lower))
    }

    if(!is.null(residuals)) {
        if(length(residuals) == length(y)) {
            residuals <- as.matrix(y + residuals)
            ylim <- range(c(ylim, residuals[,which.cf]), na.rm= TRUE)
        } else {
            residuals <- NULL
            warning(paste("Residuals do not match x in \"", ylab, 
                "\" preplot object", sep=""))
        }
    }


    all.missingy <- all(is.na(y))

    if(all.missingy)
        return()

    ylim <- ylim.scale(ylim, scale)

    if(overlay) {
        if(!length(which.cf)) which.cf = 1:ncol(uy)  # Added 7/8/04
        if(!add.arg) {
            if(is.R()) {
                matplot(ux, uy[,which.cf], type="n", 
                        xlim=xlim, ylim=ylim, 
                        xlab=xlab, ylab=ylab, ...) 
            } else {
                matplot(ux, uy[,which.cf], type="n", 
                        xlim=xlim, ylim=ylim, 
                        xlab=xlab, ylab=ylab, ...)
            }
        }
        matlines(ux, uy[,which.cf],
                lwd=llwd, col=lcol, lty=llty)
        if(!is.null(residuals))
            if(ncol(y)==1) {
                points(x, residuals, pch=pch, col=pcol, cex=pcex) 
            } else {
                matpoints(x, residuals[,which.cf],
                          pch=pch, col=pcol, cex=pcex) # add.arg=TRUE,
            }
        if(rugplot)
            if(is.R()) rug(jx, col=rcol) else rug(jx)
        if(se && !is.null(se.y)) {
            matlines(ux, se.upper[,which.cf], lty= slty, lwd=slwd, col=scol)
            matlines(ux, se.lower[,which.cf], lty= slty, lwd=slwd, col=scol)
        }
    } else {
        YLAB <- ylab 

        pcex = rep(pcex, len=ncol(uy))
        pch  = rep(pch , len=ncol(uy))
        pcol = rep(pcol, len=ncol(uy))
        lcol = rep(lcol, len=ncol(uy))
        llty = rep(llty,  len=ncol(uy))
        llwd = rep(llwd,  len=ncol(uy))
        slty = rep(slty, len=ncol(uy))
        rcol = rep(rcol, len=ncol(uy))
        scol = rep(scol, len=ncol(uy))
        slwd = rep(slwd, len=ncol(uy))

        for(i in 1:ncol(uy))
        if(!length(which.cf) ||
           (length(which.cf) && any(which.cf==i))) {

            if(is.Numeric(ylim0, allow=2)) {
                ylim = ylim0
            } else {
                ylim <- range(ylim0, uy[,i], na.rm= TRUE)
                if(se && !is.null(se.y))
                    ylim <- range(ylim0, se.lower[,i], se.upper[,i], na.rm= TRUE)
                if(!is.null(residuals))
                    ylim <- range(c(ylim, residuals[,i]), na.rm= TRUE)
                ylim <- ylim.scale(ylim, scale)
            }
            if(ncol(uy)>1 && length(separator))
                YLAB <- paste(ylab, separator, i, sep="")  
            if(!add.arg) {
                if(one.at.a.time) {
                    readline("Hit return for the next plot ")
                }
                if(is.R()) {
                    plot(ux, uy[,i], type="n", 
                         xlim=xlim, ylim=ylim, 
                         xlab=xlab, ylab=YLAB, ...)
                } else {
                    plot(ux, uy[,i], type="n", 
                         xlim=xlim, ylim=ylim, 
                         xlab=xlab, ylab=YLAB, ...)
                }
            }
            lines(ux, uy[,i], 
                 lwd=llwd[i], col=lcol[i], lty=llty[i])
            if(!is.null(residuals))
                points(x, residuals[,i], pch=pch[i], col=pcol[i], cex=pcex[i]) 
            if(rugplot)
                if(is.R()) rug(jx, col=rcol[i]) else rug(jx)

            if(se && !is.null(se.y)) {
                lines(ux, se.upper[,i], lty=slty[i], lwd=slwd[i], col=scol[i])
                lines(ux, se.lower[,i], lty=slty[i], lwd=slwd[i], col=scol[i])
            }
        }
    }
}



vplot.matrix <- function(x, y, se.y=NULL, xlab, ylab,
                         residuals=NULL, rugplot= FALSE, scale=0, se= FALSE, 
                         offset.arg=0, deriv.arg=0, overlay= FALSE, 
                         which.cf=NULL, ...)
{
    stop("You shouldn't ever call this function!") 
}


add.hookey <- function(ch, deriv.arg=0) {

    if(!is.numeric(deriv.arg) || deriv.arg<0 ||
       deriv.arg!=round(deriv.arg) || length(deriv.arg)>1)
        stop("bad input for the deriv argument")

    if(deriv.arg==0)
        return(ch)

    hookey <- switch(deriv.arg, "'", "''", "'''", "''''",
                                "'''''", stop("too high a derivative"))
    nc <- nchar(ch)
    sub <- substring(ch, 1:nc, 1:nc)
    if(nc >= 2 && sub[1]=="s" && sub[2]=="(") {
        paste("s", hookey, substring(ch, 2, nc), sep="", coll="")
    } else {
        paste(ch, hookey, sep="", collapse="")
    }
}



vplot.factor <- function(x, y, se.y=NULL, xlab, ylab, 
                         residuals=NULL, rugplot= FALSE, scale=0, 
                         se= FALSE, xlim=NULL, ylim=NULL, 
                         offset.arg=0, deriv.arg=0, overlay= FALSE, 
                         which.cf=NULL, ...)
{
    if(deriv.arg>0)
        return(NULL)

    if(length(y)/length(x)  != round(length(y)/length(x)))
        stop("length of x and y do not seem to match")
    y <- as.matrix(y) 

    if(!is.null(se.y))
        se.y <- as.matrix(se.y)
    if(!is.null(se.y) && any(is.na(se.y)))
        se.y <- NULL

    if(!is.null(residuals))  {
        residuals <- as.matrix(residuals)
        if(ncol(residuals) != ncol(y)) {
            warning("ncol(residuals) != ncol(y) so residuals are not plotted")
            residuals <- NULL
        }
    }
    if(overlay) {
        vvplot.factor(x, y,
                     se.y=if(is.null(se.y)) NULL else se.y,
                     xlab=xlab, ylab=ylab,
                     residuals=residuals,
                     rugplot=rugplot, scale=scale,
                     se=se, xlim=xlim, ylim=ylim, ...) 
    } else {
        for(i in 1:ncol(y)) {
            ylab <- rep(ylab, len=ncol(y))
            if(ncol(y)>1)
                ylab <- dimnames(y)[[2]]
            vvplot.factor(x, y[,i,drop= FALSE],
                         se.y=if(is.null(se.y)) NULL else se.y[,i,drop= FALSE], 
                         xlab=xlab, ylab=ylab[i],
                         residuals= if(is.null(residuals))
                             NULL else residuals[,i,drop= FALSE],
                         rugplot=rugplot, scale=scale,
                         se=se, xlim=xlim, ylim=ylim, ...) 

        }
    } 
    invisible(NULL)
}





vvplot.factor <- function(x, y, se.y=NULL, xlab, ylab,
                          residuals=NULL, rugplot= FALSE, scale=0,
                          se= FALSE, xlim=NULL, ylim=NULL, 
                          ...)
{

    M <- ncol(y)
    nn <- as.numeric(table(x))
    codex <- as.numeric(x)
    ucodex <- seq(nn)[nn > 0]
    o <- match(ucodex, codex, 0)

    uy <- y[o,,drop= FALSE]
    ylim <- range(ylim, uy)
    xlim <- range(c(0, sum(nn), xlim))
    rightx <- cumsum(nn)
    leftx <- c(0, rightx[ - length(nn)])
    ux <- (leftx + rightx)/2
    delta <- (rightx - leftx)/8

    jx <- runif(length(codex), (ux - delta)[codex], (ux + delta)[codex])
    nnajx <- jx[!is.na(jx)]

    if(rugplot)
        xlim <- range(c(xlim, nnajx))
    if(se && !is.null(se.y)) {
        se.upper <- uy + 2 * se.y[o,,drop= FALSE]
        se.lower <- uy - 2 * se.y[o,,drop= FALSE]
        ylim <- range(c(ylim, se.upper, se.lower))
    }
    if(!is.null(residuals)) {
        if(length(residuals) == length(y)) {
            residuals <- y + residuals
            ylim <- range(c(ylim, residuals))
        } else {
            residuals <- NULL
            warning(paste("Residuals do not match x in \"", ylab, 
                "\" preplot object", sep=""))
        }
    }
    ylim <- ylim.scale(ylim, scale)
    Levels <- levels(x)
    if(!all(nn)) {
        keep <- nn > 0
        nn <- nn[keep]
        ux <- ux[keep]
        delta <- delta[keep]
        leftx <- leftx[keep]
        rightx <- rightx[keep]
        Levels <- Levels[keep]
    }


    about <- function(ux, M, Delta=1/M) {
        if(M==1) return(cbind(ux))
        ans <- matrix(as.numeric(NA), length(ux), M)
        grid <- seq(-Delta, Delta, len=M)
        for(i in 1:M) {
            ans[,i] <- ux + grid[i]
        }
        ans
    }

    uxx <- about(ux, M, Delta=min(delta))
    xlim <- range(c(xlim, uxx))

    if(is.R()) {
        matplot(ux, uy, ylim=ylim, xlim=xlim, xlab="", type="n", 
                ylab=ylab, axes= FALSE, frame.plot=TRUE, ...)
        mtext(xlab, 1, 2, adj=0.5)
        axis(side=2)

        lpos <- par("mar")[3]
        mtext(Levels, side=3, line=lpos/2, at=ux, adj=0.5, srt=45)
    } else {
        matplot(ux, uy, ylim=ylim, xlim=xlim, xlab="", type="n", 
                ylab=ylab, xaxt="c", ...)   # xaxt="c",  xaxt="n", 
        mtext(xlab, 1, 2, adj=0.5)
        axis(side=3, at=ux, labels=Levels, srt=45, ticks= FALSE, adj=0)
    }

    for(i in 1:M)
        segments(uxx[,i] - 1.0 * delta, uy[,i],
                 uxx[,i] + 1.0 * delta, uy[,i])
    if(!is.null(residuals)) {
        for(i in 1:M) {
            jux <- uxx[,i]
            jux <- jux[codex]
            jux <- jux + runif(length(jux), -0.7*min(delta), 0.7*min(delta))
            if(M==1) points(jux, residuals[,i]) else 
                points(jux, residuals[,i], pch=as.character(i))
        }
    }
    if(rugplot)
        rug(nnajx)
    if(se) {
        for(i in 1:M) {
            segments(uxx[,i]+0.5*delta, se.upper[,i],
                     uxx[,i]-0.5*delta, se.upper[,i])
            segments(uxx[,i]+0.5*delta, se.lower[,i],
                     uxx[,i]-0.5*delta, se.lower[,i])
            segments(uxx[,i], se.lower[,i], uxx[,i], se.upper[,i], lty=2)
        }
    }
    invisible(diff(ylim))
}


if(!isGeneric("vplot"))
setGeneric("vplot", function(x, ...) standardGeneric("vplot"))
setMethod("vplot", "factor", function(x, ...)
         vplot.factor(x, ...))
setMethod("vplot", "list", function(x, ...)
         vplot.list(x, ...))
setMethod("vplot", "matrix", function(x, ...)
         vplot.matrix(x, ...))
setMethod("vplot", "numeric", function(x, ...)
         vplot.numeric(x, ...))



setMethod("plot", "vlm",
           function(x, y, ...) {
           if(!missing(y)) stop("can't process the \"y\" argument")
           invisible(plotvlm(x, y, ...))})
setMethod("plot", "vglm",
           function(x, y, ...) {
           if(!missing(y)) stop("can't process the \"y\" argument")
           invisible(plotvglm(x, y, ...))})
setMethod("plot", "vgam",
           function(x, y, ...) {
           if(!missing(y)) stop("can't process the \"y\" argument")
           invisible(plotvgam(x, ...))})




if(FALSE) 
vmerge.list = function(list1, list2) {




    if(!is.list(list1) || !is.list(list1))
        stop("list1 and list2 must be lists")
    
    n1 = names(list1)
    n2 = names(list2)
    un12 = unique(c(n1, n2)) 
    ans = vector("list", length(un12)) 
    names(ans) = un12
    for(i in un12) {
        ans[[i]] = if(length(list1[[i]])) list1[[i]] else list2[[i]]
    }
    ans 
}



plotqrrvglm = function(object,
               rtype = c("pearson", "response", "deviance", "working"), 
               ask = FALSE,
               main = paste(Rtype, "residuals vs latent variable(s)"),
               xlab="Latent Variable",
               ITolerances = object@control$EqualTolerances,
               ...) {
    M = object@misc$M
    n = object@misc$n
    Rank = object@control$Rank
    Coef.object = Coef(object, ITolerances = ITolerances)
    rtype <- match.arg(rtype, c("pearson", "response", "deviance", "working"))[1]
    res = resid(object, type=rtype)

    my.ylab = if(length(object@misc$ynames)) object@misc$ynames else 
              rep(" ", len=M)
    Rtype = switch(rtype, pearson="Pearson", response="Response",
                   deviance="Deviance", working="Working")

    done = 0
    for(r in 1:Rank)
        for(i in 1:M) {
            plot(Coef.object@lv[,r], res[,i],
                 xlab=paste(xlab, if(Rank==1) "" else r, sep=""),
                 ylab=my.ylab[i],
                 main = main, ...)
            done = done + 1
            if(done >= prod(par()$mfrow) && ask && done != Rank*M) {
                done = 0
                readline("Hit return for the next plot: ")
            }
        }
    object
}

setMethod("plot", "qrrvglm", function(x, y, ...)
         invisible(plotqrrvglm(object=x, ...)))





