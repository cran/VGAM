# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.





s.vam <- function(x, z, wz, s, which, smooth.frame, bf.maxit=10,
                  bf.epsilon=0.001, trace=FALSE, se.fit=TRUE, 
                  xbig.save, Blist, ncolBlist, M, qbig, U,
                  backchat=if(is.R()) FALSE else TRUE,
                  all.knots=FALSE, nk=NULL,
                  sf.only=FALSE)
{
    nwhich <- names(which)


    dxbig <- as.integer(dim(xbig.save))
    pbig <- dxbig[2]


    if(!length(smooth.frame$first)) {
        data <- smooth.frame[, nwhich, drop=FALSE]
        smooth.frame <- vgam.match(data, all.knots=all.knots, nk=nk)
        smooth.frame$first <- FALSE  # No longer first for next time

        dx <- as.integer(dim(x))
        smooth.frame$n <- dx[1]
        smooth.frame$p <- dx[2]
        attr(data, "class") <- NULL

        spar <- lapply(data, attr, "spar")
        df <- lapply(data, attr, "df")
        s.xargument <- lapply(data, attr, "s.xargument")
    
        for(k in 1:length(nwhich)) {
            i <- nwhich[k]

            temp <- spar[[i]]
            if(!is.numeric(temp) || any(temp<0))
                stop("spar cannot be negative or non-numeric")
            if(length(temp) > ncolBlist[i])
                warning(paste("only the first", ncolBlist[i],
                              "values of spar used for variable \"",
                              s.xargument, "\""))
            spar[[i]] <- rep(temp, length=ncolBlist[i])   # recycle
    
            temp <- df[[i]]
            if(!is.numeric(temp) || any(temp<1))
                stop("df is non-numeric or less than 1")
            if(length(temp) > ncolBlist[i])
                warning(paste("only the first", ncolBlist[i],
                              "values of df used for variable \"",
                              s.xargument, "\""))
            df[[i]] <- rep(temp, length=ncolBlist[i])    # recycle
            if(max(temp) > smooth.frame$nef[k]-1)
                stop(paste("df value too high for variable \"",
                           s.xargument, "\""))
    
            if(any(spar[[i]]!=0) && any(df[[i]]!=4))
                stop("can't specify both spar and df")
        }

        spar <- unlist(spar)
        df <- unlist(df)
        smooth.frame$spar <- spar     # original
        smooth.frame$df <- df         # original
    
        if(sum(smooth.frame$df[smooth.frame$spar==0]) + pbig > 
            smooth.frame$n * sum(ncolBlist[nwhich]))
            stop("too many parameters/dof for data on hand")
    
        xn.big <- labels(xbig.save)[[2]]
        asgn <- attr(xbig.save, "assign")
        aa <- NULL
        for(i in nwhich) {
            aa <- c(aa, xn.big[asgn[[i]]])
        }
        smooth.frame$ndfspar <- aa             # Stored here
        smooth.frame$xn.big <- xn.big          # Stored here
        smooth.frame$s.xargument <- s.xargument    # Stored here
    
        smooth.frame$smap=as.vector(cumsum(
            c(1,ncolBlist[nwhich]))[1:length(nwhich)])
    
        smooth.frame$try.spar <- spar
        smooth.frame$prev.dof <- df


        smooth.frame$bindex <- as.integer(cumsum(c(1,
            smooth.frame$nknots*ncolBlist[nwhich])))
        smooth.frame$kindex = as.integer(cumsum(c(1, 4+smooth.frame$nknots)))
    }
    if(sf.only)
        return(smooth.frame)

    ldk <- 4 * max(ncolBlist[nwhich])   # was M;     # Prior to 11/7/02
    ldk <- 3 * max(ncolBlist[nwhich]) + 1   # 11/7/02


    which <- unlist(which)
    p <- smooth.frame$p
    n <- smooth.frame$n
    dimw <- if(is.matrix(wz)) ncol(wz) else 1

    dimu <- if(is.matrix(U)) nrow(U) else 1

    index <- iam(NA, NA, M, both=TRUE)

    nBlist <- names(Blist)
    for(i in length(nBlist):1) {
        if(!any(nBlist[i] == nwhich))
            Blist[[i]] <- NULL
    }
    trivc <- trivial.constraints(Blist)

    ncbvec <- ncolBlist[nwhich]
    ncolb <- max(ncbvec)

    pmax.mwk <- rep(dimw, length(trivc))
    pmax.mwk <- pmax(ncbvec*(ncbvec+1)/2, dimw)

    size.twk <- max((4+4*smooth.frame$nef)*ncbvec + dimu*smooth.frame$nef)

    size.twk <- max(size.twk, M*smooth.frame$n)

    fit <- dotFortran(name="vbfa", 
        as.integer(backchat), n = as.integer(n), M = as.integer(M),
            npetc = as.integer(c(n, p, length(which), se.fit, 0, 
                                 bf.maxit, 0, M, n*M, pbig, 
                                 qbig, dimw, dimu, ier=0, ldk=ldk)),
        as.double(x), 
            y = as.double(z), w = as.double(wz),
            spar = as.double(smooth.frame$try.spar), 
            df = as.double(smooth.frame$df),
      as.integer(smooth.frame$o),as.integer(smooth.frame$nef),as.integer(which),
        etal = double(M*n), smooth = as.double(s), eta = double(M*n),
            s0 = double((2*M)*(2*M)*2),
        beta = double(pbig), var = if(se.fit) as.double(s) else double(1),
            as.double(bf.epsilon),
        qr = as.double(xbig.save), qraux = double(pbig),
        qpivot = as.integer(1:pbig),
        xbig = if(backchat) as.double(xbig.save) else double(1),
            U = as.double(U),
            as.double(unlist(Blist)),
        as.integer(ncbvec), as.integer(smooth.frame$smap),
            rcind = integer(M*(M+1)), trivc = as.integer(trivc),
        work1 = double(3*qbig + (9+2*4+max(smooth.frame$nknots))*
                     max(smooth.frame$nknots)),
            wk2 = double(n*M*3),
            wkmm = double(M*M*16 + M*pbig),
            work3 = double(max(max(2 * smooth.frame$nef * ncbvec^2), 
                           max(smooth.frame$nknots * ncbvec * (4*ncbvec+1)))),
        sgdub = double(max(smooth.frame$nknots) * max(4,ncolb)),
            bmb = double(M*M),
            lev = double(max(smooth.frame$nef * ncbvec)),
        mwk = double(max(smooth.frame$nef * (1 + 2*M + pmax.mwk)) ),
            twk = double(size.twk), 
        bcoefficients = double(sum(smooth.frame$nknots*ncbvec)),
            knots = as.double(unlist(smooth.frame$knots)),
            resss = double(1),
        bindex = as.integer(smooth.frame$bindex),
            nknots = as.integer(smooth.frame$nknots),
            itwk = integer(2*M),
            kindex = as.integer(smooth.frame$kindex))

    dim(fit$qr) = dim(xbig.save)
    dimnames(fit$qr) = dimnames(xbig.save)
    dim(fit$y) = dim(z)
    dimnames(fit$y) = dimnames(z)
    dim(fit$smooth) = dim(s)
    dimnames(fit$smooth) = dimnames(s)   # Needed for vgam.nlchisq
    if(se.fit) {
        dim(fit$var) = dim(s)
        dimnames(fit$var) = dimnames(s)
    }






    if(fit$npetc[14] != 0)
        stop("something went wrong in the Fortran subroutine vbfa()")

    fit$eta <- if(M>1) matrix(fit$eta,n,M,byrow=TRUE) else c(fit$eta)

    nit <- fit$npetc[5]
    qrank <- fit$npetc[7]


    smooth.frame$try.spar <- fit$spar
    change <- abs(smooth.frame$prev.dof-fit$df)/fit$df > 0.05 &
                  smooth.frame$spar==0
    smooth.frame$try.spar[change] <- 0         # For next time
    smooth.frame$prev.dof <- fit$df

    if((nit == bf.maxit) & bf.maxit > 1)
        warning(paste("s.vam convergence not obtained in", bf.maxit, 
            "iterations"))

    R <- fit$qr[1:pbig, 1:pbig]
    R[lower.tri(R)] <- 0



    Bspline <- vector("list", length(nwhich))
    names(Bspline) <- nwhich
    for(i in 1:length(nwhich)) {
        ans = fit$bcoeff[(smooth.frame$bindex[i]):(smooth.frame$bindex[i+1]-1)]
        ans = matrix(ans, ncol=ncolBlist[nwhich[i]])
        Bspline[[i]] = new("vsmooth.spline.fit",
                           "Bcoefficients" = ans,
                           "xmax"          = smooth.frame$xmax[i],
                           "xmin"          = smooth.frame$xmin[i],
                           "knots"         = as.vector(smooth.frame$knots[[i]]))
    }


    rl <- list(
        Bspline = Bspline,
        coefficients = fit$beta,
        df.residual = n*M - qrank - sum(fit$df - 1),
        fitted.values = fit$eta, 
        nl.df = fit$df - 1,
        qr = list(qr=fit$qr, rank=qrank, qraux=fit$qraux, pivot=fit$qpivot),
        R = R, 
        rank = qrank, 
        residuals = fit$y - fit$eta, 
        rss = fit$resss,
        smooth = fit$smooth,
        spar = fit$spar,
        s.xargument = unlist(smooth.frame$s.xargument))


    names(rl$coefficients) <- smooth.frame$xn.big
    names(rl$spar) <- smooth.frame$ndfspar
    names(rl$nl.df) <- smooth.frame$ndfspar

    if(se.fit)
        rl <- c(rl, list(var=fit$var))
    c(list(smooth.frame=smooth.frame), rl)
}


