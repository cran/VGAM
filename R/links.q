# These functions are
# Copyright (C) 1998-2010 T.W. Yee, University of Auckland. All rights reserved.




  ToString = function(x) paste(x, collapse = ",")











TypicalVGAMfamilyFunction <- function(lsigma="loge", esigma=list(),
                                      isigma=NULL, parallel=TRUE,
                                      shrinkage.init = 0.95,
                                      nointercept = NULL, method.init=1,
                                      nsimEIM=100, zero=NULL) {
    NULL
}

TypicalVGAMlinkFunction <- function(theta,
    earg=list(), inverse=FALSE, deriv=0, short=TRUE, tag=FALSE) {
    NULL
}



namesof <- function(theta,
                    link,
                    earg=list(),
                    tag=FALSE,
                    short=TRUE)
{


        string <- paste(link,
            "(theta=theta, earg=earg, short=short, tag=tag)", sep="")
        calls <- parse(text=string)[[1]]
        ans <- eval(calls) 
        return(ans)
}

theta2eta <- function(theta, link, earg=list()) {
    string <- paste(link, "(theta=theta, earg=earg)", sep="")
    calls <- parse(text=string)[[1]]
    eval(calls) 
}




eta2theta <- function(theta, link="identity", earg=list()) {
    if (is.null(link))
        link <- "identity"



    llink <- length(link)
    if (llink == 1) {
        string <- paste(link, "(theta=theta, earg=earg, inverse=TRUE)", sep="")
        calls <- parse(text=string)[[1]]
        return(eval(calls))
    } else 
    if (llink > 1) {
        if (is.matrix(theta) && llink == ncol(theta)) {



            ans <- NULL
            for(iii in 1:llink) {
                use.earg = if (is.list(earg) && length(earg)==llink &&
                              is.list(earg[[iii]])) earg[[iii]] else earg
                string = paste(link[iii],
                           "(theta=theta[,iii], earg=use.earg, inverse=TRUE)",
                           sep="")
                calls <- parse(text=string)[[1]]
                ans <- cbind(ans, eval(calls))
            }
        } else {
            if (length(theta) < llink)
                theta = rep(theta, len=llink)

            if (length(theta) != llink)
                stop("length of theta and link don't match") 

            ans <- NULL
            for(iii in 1:llink) {
                string = paste(link[iii],
                               "(theta=theta[iii], earg=earg, inverse=TRUE)",
                               sep="")
                calls <- parse(text=string)[[1]]
                ans <- c(ans, eval(calls))
            }
        }
        return(ans)
    } else 
        stop("length(link)==0 not allowed") 
}



dtheta.deta <- function(theta, link, earg=list()) {

    string <- paste(link, "(theta=theta, earg=earg, deriv=1)", sep="")
    calls <- parse(text=string)[[1]]
    eval(calls) 
}


d2theta.deta2 <- function(theta, link, earg=list())
{

    string <- paste(link, "(theta=theta, earg=earg, deriv=2)", sep="")
    calls <- parse(text=string)[[1]]
    eval(calls) 
}





.all.links = c("cloglog",
               "fisherz", "fsqrt", "identity", "inverse", 
               "logc", "loge", "logit", "loglog", 
               "logoff", "nreciprocal", "nloge", 
               "powl", "probit", "reciprocal", "rhobit",
               "golf", "polf", "nbolf", "nbolf2")


loglog <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("loglog(",theta,")", sep="") else
            paste("log(log(",theta,"))", sep="")
        if (tag) 
            string <- paste("Log-Log:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval))
        theta[theta <= 1.0] <- earg$bval
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            exp(exp(theta))
        }
    } else {
        switch(deriv+1, {
            log(log(theta))},
            theta * log(theta),
            {  junk <- log(theta)
               -junk^2 / (1 + junk)
            },
            stop("'deriv' unmatched"))
    }
}




cloglog <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                    short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("cloglog(",theta,")", sep="") else
            paste("log(-log(1-",theta,"))", sep="")
        if (tag) 
            string <- paste("Complementary log-log:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval)) {
        theta[theta <= 0.0] <- earg$bval
        theta[theta >= 1.0] <- 1.0 - earg$bval
    }
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            junk <- exp(theta)
            -expm1(-junk)
        }
    } else {
        switch(deriv+1, {
            log(-log1p(-theta))},
            -(1-theta) * log1p(-theta),
            {  junk <- log1p(-theta)
               -(1-theta) * (1 + junk) * junk
            },
            stop("'deriv' unmatched"))
    }
}




probit <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("probit(",theta,")", sep="") else
            paste("qnorm(", theta, ")", sep="")
        if (tag) 
            string <- paste("Probit:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval)) {
        theta[theta <= 0.0] <- earg$bval
        theta[theta >= 1.0] <- 1-earg$bval
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            ans <- pnorm(theta)
            if (is.matrix(theta))
                dim(ans) <- dim(theta)
            ans
        }
    } else {
        switch(deriv+1,{
            ans <- qnorm(theta)
            if (is.matrix(theta))
                dim(ans) <- dim(theta)
            ans
        },
        {
           if (is.matrix(theta)) {
               ans <- dnorm(qnorm(theta))
               dim(ans) <- dim(theta)
               ans
           } else dnorm(qnorm(as.vector(theta)))
        }, 
        {
            junk <- qnorm(theta)
            ans <- -junk * dnorm(junk)
            if (is.vector(theta)) ans else
            if (is.matrix(theta)) {
                dim(ans) <- dim(theta)
                ans
            } else {
                warning("can only handle vectors and matrices; converting to vector")
                ans
            }
        })
    }
}








loge <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                 short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("log(",theta,")", sep="") else
            paste("log(", theta, ")", sep="")
        if (tag) 
            string <- paste("Log:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval))
        theta[theta <= 0.0] <- earg$bval
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            exp(theta)
        }
    } else {
        switch(deriv+1, {
           log(theta)},
           theta,
           theta)
    }
}




identity <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                     short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- theta 
        if (tag) 
            string <- paste("Identity:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            theta
        }
    } else {
        switch(deriv+1,
           theta,
           theta*0 + 1,
           theta*0)
    }
}

nidentity <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                     short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- paste("-", theta, sep="")
        if (tag) 
            string <- paste("Negative-Identity:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            -theta
        }
    } else {
        switch(deriv+1,
           -theta,
           theta*0 - 1,
           theta*0)
    }
}


reciprocal <- function(theta, earg=list(), inverse.arg=FALSE, deriv=0,
                     short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- paste("1/",theta, sep="")
        if (tag) 
            string <- paste("Reciprocal:", string) 
        return(string)
    }
    if (!inverse.arg && is.list(earg) && length(earg$bval))
        theta[theta == 0.0] <- earg$bval
    if (inverse.arg) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse.arg=FALSE, deriv=deriv)
        } else {
            1/theta
        }
    } else {
        switch(deriv+1,{
           1/theta},
           -theta^2,
           2*theta^3)
    }
}


nloge <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                 short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("-log(",theta,")", sep="") else
            paste("-log(", theta, ")", sep="")
        if (tag) 
            string <- paste("Negative log:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval))
        theta[theta <= 0.0] <- earg$bval
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            exp(-theta)
        }
    } else {
        switch(deriv+1, {
           -log(theta)},
           -theta,
           theta)
    }
}



nreciprocal <- function(theta, earg=list(), inverse.arg=FALSE, deriv=0,
                     short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- paste("-1/",theta, sep="")
        if (tag) 
            string <- paste("Negative reciprocal:", string) 
        return(string)
    }
    if (!inverse.arg && is.list(earg) && length(earg$bval))
        theta[theta == 0.0] <- earg$bval
    if (inverse.arg) {
        if (deriv>0) {
            1 / nreciprocal(theta, earg=earg, inverse.arg=FALSE, deriv)
        } else {
            -1/theta
        }
    } else {
        switch(deriv+1, {
           -1/theta},
           theta^2,
           2*theta^3)
    }
}


natural.ig <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                       short=TRUE, tag=FALSE)
{

    if (is.character(theta)) {
        string <- paste("-1/",theta, sep="")
        if (tag) 
            string <- paste("Negative inverse:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1 / nreciprocal(theta, earg=earg, inverse=FALSE, deriv)
        } else {
            1/ sqrt(-2*theta)
        }
    } else {
        switch(deriv+1,
           -1/(2*theta^2),
           theta^3,
           3*theta^5)
    }
}





rhobit <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("rhobit(",theta,")", sep="") else
            paste("log((1+", theta, ")/(1-", theta, "))", sep="")
        if (tag) 
            string <- paste("Rhobit:", string) 
        return(string)
    }

    if (!inverse && is.list(earg) && length(earg)) {
        bminvalue = if (length(earg$bminval)) earg$bminval else NULL
        bmaxvalue = if (length(earg$bmaxval)) earg$bmaxval else NULL
       if (!inverse && length(bminvalue)) theta[theta <= -1.0] <- bminvalue
       if (!inverse && length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
    }

    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            junk <- exp(theta)
            expm1(theta) / (junk+1.0)
        }
    } else {
        switch(deriv+1,{
            log1p(theta) - log1p(-theta)},
            (1 - theta^2) / 2,
            (1 - theta^2)^2 / (4*theta))
    }
}



fisherz <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("fisherz(",theta,")", sep="") else
            paste("(1/2)log((1+", theta, ")/(1-", theta, "))", sep="")
        if (tag) 
            string <- paste("Fisher's Z transformation:", string) 
        return(string)
    }
 
    if (!inverse && is.list(earg) && length(earg)) {
        bminvalue = if (length(earg$bminval)) earg$bminval else NULL
        bmaxvalue = if (length(earg$bmaxval)) earg$bmaxval else NULL
       if (!inverse && length(bminvalue)) theta[theta <= -1.0] <- bminvalue
       if (!inverse && length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
    }

    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            junk <- exp(2*theta)
            expm1(2*theta) / (junk+1.0)
        }
    } else {
        switch(deriv+1,
           0.5 * log1p(theta) - log1p(-theta),
           1.0 - theta^2,
           (1.0 - theta^2)^2 / (2*theta))
    }
}




fsqrt <- function(theta, earg=list(min=0, max=1, mux=sqrt(2)),
                  inverse=FALSE, deriv=0, short=TRUE, tag=FALSE)
{
    min=0; max=1; mux=sqrt(2)
    if (!is.list(earg)) stop("earg must be a list")
    if (is.Numeric(earg$min)) min = earg$min
    if (is.Numeric(earg$max)) max = earg$max
    if (is.Numeric(earg$mux)) mux = earg$mux
    if (!is.Numeric(min,allow=1)) stop("bad input for 'min' component")
    if (!is.Numeric(max,allow=1)) stop("bad input for 'max' component")
    if (!is.Numeric(mux,allow=1,posit=TRUE)) stop("bad input for 'mux' component")
    if (min >= max) stop("'min' >= 'max' is not allowed")

    if (is.character(theta)) {
        string <- if (short) 
            paste("fsqrt(",theta,")", sep="") else {
            if (abs(mux-sqrt(2)) < 1.0e-10)
                paste("sqrt(2*",theta,") - sqrt(2*(1-",theta,"))", sep="") else
            paste(as.character(mux),
            " * (sqrt(",theta,"-",min,") - sqrt(",max,"-",theta,"))", sep="")
        }
        if (tag) 
            string <- paste("Folded Square Root:", string) 
        return(string)
    }

    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            mid = (min + max) / 2
            boundary = mux * sqrt(max - min)
            temp = pmax(0, (theta/mux)^2 * (2*(max-min) - (theta/mux)^2))
            ans = theta
            if (any(ind5 <- theta <  0))
                ans[ind5] = mid - 0.5 * sqrt(temp[ind5])
            if (any(ind5 <- theta >= 0))
                ans[ind5] = mid + 0.5 * sqrt(temp[ind5])
            ans[theta < -boundary] <- NA
            ans[theta >  boundary] <- NA
            ans
        }
    } else {
        switch(deriv+1,
            mux * (sqrt(theta-min) - sqrt(max-theta)),
           (2 / mux) / (1/sqrt(theta-min) + 1/sqrt(max-theta)),
           -(4 / mux) / ((theta-min)^(-3/2) - (max-theta)^(-3/2)))
    }
}



powl <- function(theta, earg=list(power=1), inverse=FALSE, deriv=0,
                 short=TRUE, tag=FALSE)
{

    if (!length(earg) || is.list(earg)) {
        exponent = if (length(earg$power)) earg$power else 1
        if (exponent == 0)
            stop("use the 'loge' link")
    } else {
        stop("'earg' must be a list or NULL")
    }

    if (is.character(theta)) {
        string <- if (short) 
            paste("powl(",theta,", earg=list(power=", as.character(exponent),
                  "))", sep="") else
            paste(theta, "^(", as.character(exponent), ")", sep="")
        if (tag) 
            string <- paste("Power:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            theta^(1/exponent)
        }
    } else {
        switch(deriv+1,
        {
            theta^exponent
        },
        {
            (theta^(1-exponent)) / exponent
        },
        {
            (theta^(2-exponent)) / (exponent * (exponent-1))
        })
    }
}


elogit <- function(theta, earg=list(min=0, max=1), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (!length(earg) || is.list(earg)) {
        A = if (length(earg$min)) earg$min else 0
        B = if (length(earg$max)) earg$max else 1
        bminvalue = if (length(earg$bminval)) earg$bminval else NULL
        bmaxvalue = if (length(earg$bmaxval)) earg$bmaxval else NULL
       if (!inverse && length(bminvalue)) theta[theta <= A] <- bminvalue
       if (!inverse && length(bmaxvalue)) theta[theta >= B] <- bmaxvalue
    } else {
        stop("'earg' must be a list or NULL")
    }
    if (is.character(theta)) {
        string <- if (short) {
            if (A != 0 || B != 1)
            paste("elogit(",theta,", earg=list(min=",A,
                  ", max=",B,"))",sep="") else
            paste("elogit(",theta,")",sep="")
            } else
            paste("log((",theta,"-min)/(max-",theta,"))", sep="")
        if (tag) 
            string <- paste("Extended logit:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            junk <- if (is.R()) care.exp(theta) else care.exp(theta)
            (A + B*junk) / (1.0 + junk)
        }
    } else {
        switch(deriv+1, {
           log((theta-A)/(B-theta))},
           (theta-A) * (B - theta) / (B-A),
           (theta-A) * (B - theta) * (B - 2 * theta + A) / (B-A)^2)
    }
}





 logit <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("logit(",theta,")", sep="") else
            paste("log(",theta,"/(1-",theta,"))", sep="")
        if (tag) 
            string <- paste("Logit:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval)) {
        theta[theta <= 0.0] <- earg$bval;
        theta[theta >= 1.0] <- 1.0 - earg$bval;
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            eta <- care.exp(theta)
            eta / (1.0 + eta)
        }
    } else {
        switch(deriv+1, {
           temp2 = log(theta) - log1p(-theta)
           if (any(near0.5 <- (abs(theta - 0.5) < 0.000125)))
               temp2[near0.5] = log(theta[near0.5] / (1-theta[near0.5]))
           temp2
           },
           exp(log(theta) + log1p(-theta)),
           exp(log(theta) + log1p(-theta)) * (1 - 2 * theta))
    }
}


logc <- function(theta, earg=list(), inverse=FALSE, deriv=0,
                 short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("logc(",theta,")", sep="") else
            paste("log(1-",theta,")", sep="")
        if (tag) 
            string <- paste("Log Complementary:", string) 
        return(string)
    }


    if (!inverse && is.list(earg) && length(earg$bval)) {
        theta[theta >= 1.0] <- earg$bval;
    }
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            -expm1(theta)
        }
    } else {
        switch(deriv+1,{
            log1p(-theta)},
           -(1.0 - theta),
           -(1.0 - theta)^2)
    }
}



logoff <- function(theta, earg=list(offset=0), inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    if (!length(earg) || is.list(earg)) {
        offset = if (length(earg$offset)) earg$offset else 0
    } else {
        stop("'earg' must be a list or NULL")
    }

    if (!is.Numeric(offset))
        stop("bad input for argument 'earg'")

    if (is.character(theta)) {
        string <- if (short) 
            paste("logoff(",theta,
                  ", list(offset=",as.character(offset),"))", sep="") else
            paste("log(", as.character(offset), "+", theta, ")", sep="")
        if (tag) 
            string <- paste("Log with offset:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            exp(theta) - offset
        }
    } else {
        switch(deriv+1,
           log(theta+offset),
           theta + offset,
           theta + offset)
    }
}


if(FALSE)
nlogoff <- function(theta, earg=0, inverse=FALSE, deriv=0,
                   short=TRUE, tag=FALSE)
{
    offset = earg
    if (!is.Numeric(offset))
        stop("bad input for argument earg")
    if (is.character(theta)) {
        string <- if (short) 
            paste("nlogoff(",theta,",",as.character(offset),")", sep="") else
            paste("log(", as.character(offset), "-", theta, ")", sep="")
        if (tag) 
            string <- paste("Negative-log with offset:", string) 
        return(string)
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            offset - exp(theta)
        }
    } else {
        switch(deriv+1,
           log(-theta+offset),
           theta - offset,
           theta - offset)
    }
}



cauchit <- function(theta, earg=list(bvalue= .Machine$double.eps),
                    inverse=FALSE, deriv=0,
                    short=TRUE, tag=FALSE)
{
    if (is.character(theta)) {
        string <- if (short) 
            paste("cauchit(",theta,")", sep="") else
            paste("tan(pi*(",theta,"-0.5))", sep="")
        if (tag) 
            string <- paste("Cauchit:", string) 
        return(string)
    }
    if (!inverse && is.list(earg) && length(earg$bval)) {
        theta[theta <= 0.0] <- earg$bval
        theta[theta >= 1.0] <- 1.0 - earg$bval
    }
    if (inverse) {
        if (deriv>0) {
            1/Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            0.5 + atan(theta)/pi 
        }
    } else {
        switch(deriv+1, {
           tan(pi * (theta-0.5))},
           cos(pi * (theta-0.5))^2 / pi,
           -sin(2 * pi * (theta-0.5)))
    }
}



golf <- function(theta, earg=list(lambda=1), inverse=FALSE, deriv=0,
                 short=TRUE, tag=FALSE)
{


    cutpoint = lambda = NULL
    if (!length(earg)) {
        lambda = 1
        cutpoint = NULL
    } else if (is.list(earg)) {
        lambda = earg$lambda
        cutpoint = earg$cutpoint # Optional; if so then is a NULL
    } else
        stop("'earg' must be a list")
    if (!is.Numeric(lambda, posit=TRUE))
        stop('could not determine lambda or lambda has negative values')
    if (is.Numeric(cutpoint))
        if (any(cutpoint < 0) || !is.Numeric(cutpoint, integer=TRUE))
            warning("'cutpoint' should contain non-negative integer values")

    if (is.character(theta)) {
        string <- if (short) {
            lenl = length(lambda) > 1
            lenc = length(cutpoint) > 1
            paste("golf(",theta,", earg=list(lambda=",
                  if (lenl) "c(" else "",
                  ToString(lambda),
                  if (lenl) ")" else "",
                  if (is.Numeric(cutpoint))
            paste(", cutpoint=",
                  if (lenc) "c(" else "",
            ToString(cutpoint),
                  if (lenc) ")" else "",
            sep="") else "",
                        "))", sep="") } else {
            if (is.Numeric(cutpoint)) {
                paste("-3*log(1-qnorm(",theta,")/(3*sqrt(lambda)))",
                      " + log(cutpoint)", sep="")
            } else {
                paste("-3*log(1-qnorm(",theta,")/(3*sqrt(lambda)))", sep="")
            }
        }
        if (tag) 
            string <- paste("Gamma-ordinal link function:", string) 
        return(string)
    }

    thmat = cbind(theta)
    lambda = rep(lambda, len=ncol(thmat)) # Allow recycling for lambda
    if (is.Numeric(cutpoint)) cutpoint = rep(cutpoint, len=ncol(thmat))
    if (ncol(thmat) > 1) {
        answer = thmat
        for(ii in 1:ncol(thmat))
            answer[,ii] = Recall(theta=thmat[,ii],
                   earg=list(lambda=lambda[ii],
                   cutpoint = if (is.Numeric(cutpoint)) cutpoint[ii] else NULL),
                   inverse=inverse, deriv=deriv)
        return(answer)
    }

    answer =
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            if (is.Numeric(cutpoint)) {
                pnorm((1-care.exp(-(theta-log(cutpoint))/3)) * 3 * sqrt(lambda))
            } else {
                pnorm((1-care.exp(-theta/3)) * 3 * sqrt(lambda))
            }
        }
    } else {
        smallno = 1 * .Machine$double.eps
        Theta = theta
        Theta = pmin(Theta, 1 - smallno)  # Since theta==1 is a possibility
        Theta = pmax(Theta, smallno) # Since theta==0 is a possibility
        Ql = qnorm(Theta)
        switch(deriv+1, {
            temp = Ql / (3*sqrt(lambda))
            temp = pmin(temp, 1.0 - smallno)  # 100 / .Machine$double.eps
            -3*log(1-temp) + if (is.Numeric(cutpoint)) log(cutpoint) else 0},
            (1 - Ql / (3*sqrt(lambda))) * sqrt(lambda) * dnorm(Ql),
            {  stop('cannot handle deriv=2') },
            stop("'deriv' unmatched"))
    }
    if (!is.Numeric(answer)) stop("the answer contains some NAs")
    answer
}


polf <- function(theta, earg=stop("'earg' must be given"), 
                 inverse=FALSE, deriv=0, short=TRUE, tag=FALSE)
{
    cutpoint = NULL
    if (is.Numeric(earg)) cutpoint = earg
    if (is.list(earg)) cutpoint = earg$cutpoint
    if (!is.Numeric(cutpoint))
        stop('could not determine the cutpoint')
    if (any(cutpoint < 0) || !is.Numeric(cutpoint, integer=TRUE))
        warning("'cutpoint' should contain non-negative integer values")


    if (is.character(theta)) {
        string <- if (short) {
            lenc = length(cutpoint) > 1
            paste("polf(",theta,", earg=list(cutpoint=",
                  if (lenc) "c(" else "",
                  ToString(cutpoint),
                  if (lenc) ")" else "",
                  "))", sep="") 
        } else
            paste("2*log(0.5*qnorm(",theta,") + sqrt(cutpoint+7/8))", sep="")
        if (tag) 
            string <- paste("Poisson-ordinal link function:", string) 
        return(string)
    }


    thmat = cbind(theta)
    if (ncol(thmat) > 1) {
        answer = thmat
        cutpoint = rep(cutpoint, len=ncol(thmat)) # Reqd for the for loop
        for(ii in 1:ncol(thmat))
            answer[,ii] = Recall(theta=thmat[,ii], earg=cutpoint[ii],
                                 inverse=inverse, deriv=deriv)
        return(answer)
    }

    answer =
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            if (cutpoint == 0) {
                cloglog(theta=theta, earg=earg, inverse=inverse, deriv=deriv)
            } else {
                pnorm(2 * exp(theta/2) - 2 * sqrt(cutpoint + 7/8))
            }
        }
    } else {
        if (cutpoint == 0) {
            cloglog(theta=theta, earg=earg, inverse=inverse, deriv=deriv)
        } else {
            smallno = 1 * .Machine$double.eps
            SMALLNO = 1 * .Machine$double.xmin
            Theta = theta
            Theta = pmin(Theta, 1 - smallno)  # Since theta==1 is a possibility
            Theta = pmax(Theta, smallno) # Since theta==0 is a possibility
            Ql = qnorm(Theta)
            switch(deriv+1, {
            temp = 0.5 * Ql + sqrt(cutpoint + 7/8)
            temp = pmax(temp, SMALLNO)
            2 * log(temp)},
            (Ql/2 + sqrt(cutpoint + 7/8)) * dnorm(Ql),
            {  stop('cannot handle deriv=2') },
            stop("'deriv' unmatched"))
        }
    }
    if (!is.Numeric(answer)) stop("the answer contains some NAs")
    answer
}


nbolf <- function(theta, earg=stop("'earg' must be given"), 
                  inverse=FALSE, deriv=0, short=TRUE, tag=FALSE)
{

    cutpoint = kay = NULL
    if (is.list(earg)) {
        cutpoint = earg$cutpoint
        kay = earg$k
    }
    if (!is.Numeric(kay, positive=TRUE))
        stop("could not determine 'k' or it is not positive-valued")
    if (!is.Numeric(cutpoint))
        stop("could not determine the cutpoint")
    if (any(cutpoint < 0) || !is.Numeric(cutpoint, integer=TRUE))
        warning("'cutpoint' should contain non-negative integer values")

    if (is.character(theta)) {
        string <- if (short) {
            lenc = length(cutpoint) > 1
            lenk = length(kay) > 1
            paste("nbolf(",theta,", earg=list(cutpoint=",
                  if (lenc) "c(" else "",
                  ToString(cutpoint),
                  if (lenc) ")" else "",
                  ", k=",
                  if (lenk) "c(" else "",
                  ToString(kay),
                  if (lenk) ")" else "",
                  "))", sep="")
        } else
            paste("2*log(sqrt(k) * sinh(qnorm(",theta,")/(2*sqrt(k)) + ",
                  "asinh(sqrt(cutpoint/k))))", sep="")
        if (tag) 
            string <- paste("Negative binomial-ordinal link function:", string)
        return(string)
    }

    thmat = cbind(theta)
    kay = rep(kay, len=ncol(thmat)) # Allow recycling for kay
    cutpoint = rep(cutpoint, len=ncol(thmat)) # Allow recycling for cutpoint
    if (ncol(thmat) > 1) {
        answer = thmat
        for(ii in 1:ncol(thmat))
            answer[,ii] = Recall(theta=thmat[,ii],
                                 earg=list(cutpoint=cutpoint[ii], k=kay[ii]),
                                 inverse=inverse, deriv=deriv)
        return(answer)
    }

    answer =
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            if (cutpoint == 0) {
                1.0 - (kay / (kay + care.exp(theta)))^kay
            } else {
                pnorm((asinh(exp(theta/2)/sqrt(kay)) -
                       asinh(sqrt(cutpoint/kay))) * 2 * sqrt(kay))
            }
        }
    } else {
        smallno = 1 * .Machine$double.eps
        SMALLNO = 1 * .Machine$double.xmin
        Theta = theta
        Theta = pmin(Theta, 1 - smallno)  # Since theta==1 is a possibility
        Theta = pmax(Theta, smallno) # Since theta==0 is a possibility
        if (cutpoint == 0) {
            switch(deriv+1, {
            temp = (1 - Theta)^(-1/kay) - 1
            temp = pmax(temp, SMALLNO)
            log(kay) + log(temp)},
            (kay / (1 - Theta)^(1/kay) - kay) * (1 - Theta)^(kay+1/kay),
            {  stop('cannot handle deriv=2') },
            stop("'deriv' unmatched"))
        } else {
            Ql = qnorm(Theta)
            switch(deriv+1, {
                temp = sqrt(kay) * sinh(Ql/(2*sqrt(kay)) +
                       asinh(sqrt(cutpoint/kay)))
                temp = pmax(temp, SMALLNO)
                2 * log(temp)}, {
                arg1 = (Ql/(2*sqrt(kay)) + asinh(sqrt(cutpoint/kay)))
                sqrt(kay) * tanh(arg1) * dnorm(Ql) },
                {  stop('cannot handle deriv=2') },
                stop("'deriv' unmatched"))
        }
    }
    if (!is.Numeric(answer)) stop("the answer contains some NAs")
    answer
}





nbolf2 <- function(theta, earg=stop("'earg' must be given"), 
                   inverse=FALSE, deriv=0, short=TRUE, tag=FALSE)
{

    cutpoint = kay = NULL
    if (is.list(earg)) {
        cutpoint = earg$cutpoint
        kay = earg$k
    }
    if (!is.Numeric(kay, positive=TRUE))
        stop("could not determine 'k' or it is not positive-valued")
    if (!is.Numeric(cutpoint))
        stop("could not determine the cutpoint")
    if (any(cutpoint < 0) || !is.Numeric(cutpoint, integer=TRUE))
        warning("'cutpoint' should contain non-negative integer values")

    if (is.character(theta)) {
        string <- if (short) {
            lenc = length(cutpoint) > 1
            lenk = length(kay) > 1
            paste("nbolf2(",theta,", earg=list(cutpoint=",
                  if (lenc) "c(" else "",
                  ToString(cutpoint),
                  if (lenc) ")" else "",
                  ", k=",
                  if (lenk) "c(" else "",
                  ToString(kay),
                  if (lenk) ")" else "",
                  "))", sep="")
       } else
            paste("3*log(<a complicated expression>)", sep="")
        if (tag) 
            string = paste("Negative binomial-ordinal link function 2:", string)
        return(string)
    }

    thmat = cbind(theta)
    kay = rep(kay, len=ncol(thmat)) # Allow recycling for kay
    if (ncol(thmat) > 1) {
        answer = thmat
        for(ii in 1:ncol(thmat))
            answer[,ii] = Recall(theta=thmat[,ii],
                                 earg=list(cutpoint=cutpoint[ii], k=kay[ii]),
                                 inverse=inverse, deriv=deriv)
        return(answer)
    }

    answer =
    if (inverse) {
        if (deriv>0) {
            1 / Recall(theta=theta, earg=earg, inverse=FALSE, deriv=deriv)
        } else {
            if (cutpoint == 0) {
                1.0 - (kay / (kay + care.exp(theta)))^kay
            } else {

            a1 = -(9*cutpoint+8) / (cutpoint+1)
            a2 = (9*kay-1) / (kay * (cutpoint+1)^(1/3))
            a3 = 9 / (kay * (cutpoint+1)^(2/3))
            a4 = 9 / (cutpoint+1)
            B = exp(theta/3)
            mymat = rbind(a1^2*a2^2 + 2*a1*a2^3*B + B^2*a2^4, 0,
                       -2*a1*a2*a3*B - 2*a2^2*a3*B^2 - a1^2*a3 - a2^2*a4, 0,
                       B^2 * a3^2 + a3 * a4)
            ans = Re(t(apply(mymat, 2, polyroot)))
            theta2 = invfun = pnorm(-ans)  # pnorm(-x) = 1-pnorm(x)
            for(ii in 1:4) {
                theta2[,ii] = Recall(theta=theta2[,ii],
                                     earg=list(cutpoint=cutpoint, k=kay),
                                     inverse=FALSE, deriv=deriv)
            }
            rankmat = t(apply(abs(theta2 - theta), 1, rank))
            for(ii in 2:4) {
                if (any(index4 <- (rankmat[,ii] == 1))) {
                    invfun[index4,1] = invfun[index4,ii]
                }
            }
            invfun[,1]
            }
        }
    } else {
        smallno = 1 * .Machine$double.eps
        SMALLNO = 1 * .Machine$double.xmin
        Theta = theta
        Theta = pmin(Theta, 1 - smallno)  # Since theta==1 is a possibility
        Theta = pmax(Theta, smallno) # Since theta==0 is a possibility
        if (cutpoint == 0) {
            switch(deriv+1, {
            temp = (1 - Theta)^(-1/kay) - 1
            temp = pmax(temp, SMALLNO)
            log(kay) + log(temp)},
            (kay / (1 - Theta)^(1/kay) - kay) * (1 - Theta)^(kay+1/kay),
            {  stop('cannot handle deriv=2') },
            stop("'deriv' unmatched"))
        } else {
            Ql = qnorm(Theta)
            a1 = -(9*cutpoint+8) / (cutpoint+1)
            a2 = (9*kay-1) / (kay * (cutpoint+1)^(1/3))
            a3 = 9 / (kay * (cutpoint+1)^(2/3))
            a4 = 9 / (cutpoint+1)
            discrim = a1^2 * a3 + a2^2 * a4 - Ql^2 * a3 * a4
            denomin = Ql^2 * a3 - a2^2
            numerat = (a1*a2 - Ql * sqrt(discrim))
            argmax1 = numerat / denomin
            switch(deriv+1, {
                argmax2 = (a1*a2 + Ql * sqrt(discrim)) / denomin
                temp = ifelse(argmax1 > 0, argmax1, argmax2)
                temp = pmax(temp, SMALLNO)
                3 * log(temp)}, {
                 BB = (sqrt(discrim) - Ql^2 * a3 * a4 / sqrt(discrim)) / dnorm(Ql)
                 CC = 2 * Ql * a3 / dnorm(Ql)
                 dA.dtheta = (-denomin * BB - numerat * CC) / denomin^2
                 argmax1 / (3 * dA.dtheta)
                },
                {  stop('cannot handle deriv=2') },
                stop("'deriv' unmatched"))
        }
    }
    if (!is.Numeric(answer)) stop("the answer contains some NAs")
    answer
}



Cut = function(y, breaks=c(-Inf, quantile(c(y), prob = (1:4)/4))) {
    y = as.matrix(y)


    temp = cut(y, breaks=breaks, labels=FALSE)
    temp = c(temp) # integer vector of integers
    if (any(is.na(temp))) stop("there are NAs")
    answer = if (ncol(y) > 1) matrix(temp, nrow(y), ncol(y)) else temp
    if (ncol(y) > 1) {
        ynames = dimnames(y)[[2]]
        if (!length(ynames)) ynames = paste("Y", 1:ncol(y), sep="")
        xnames = dimnames(y)[[1]]
        if (!length(xnames)) xnames = as.character(1:nrow(y))
        dimnames(answer) = list(xnames, ynames)
    }
    attr(answer, "breaks") = breaks
    answer
}


checkCut = function(y) {
    if (!is.Numeric(y, posi=TRUE, integ=TRUE))
        stop("'y' must contain positive integers only")
    uy = unique(y)
    L = max(uy)
    oklevels = 1:L
    if (L == 1) stop("only one unique value")
    for(ii in oklevels) {
        if (all(ii != uy)) stop("there is no ", ii, " value")
    }
    TRUE
}




