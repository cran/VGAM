# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.





s.vam <- function(x, zedd, wz, smomat, which, smooth.frame, bf.maxit = 10,
                  bf.epsilon=0.001, trace=FALSE, se.fit = TRUE,
                  X_vlm_save, Blist, ncolBlist, M, qbig, Umat,
                  all.knots=FALSE, nk=NULL,
                  sf.only=FALSE)
{
    nwhich <- names(which)


    dX_vlm <- as.integer(dim(X_vlm_save))
    pbig <- dX_vlm[2]


    if (!length(smooth.frame$first)) {
        data <- smooth.frame[, nwhich, drop=FALSE]
        smooth.frame <- vgam.match(data, all.knots=all.knots, nk=nk)
        smooth.frame$first <- FALSE  # No longer first for next time

        dx <- as.integer(dim(x))
        smooth.frame$n_lm <- dx[1]
        smooth.frame$p_lm <- dx[2]
        attr(data, "class") <- NULL

        sparv <- lapply(data, attr, "spar")
        dfvec <- lapply(data, attr, "df")
        s.xargument <- lapply(data, attr, "s.xargument")
    
        for(kk in 1:length(nwhich)) {
            ii <- nwhich[kk]

            temp <- sparv[[ii]]
            if (!is.numeric(temp) || any(temp < 0)) {
                stop("spar cannot be negative or non-numeric")
            }
            if (length(temp) > ncolBlist[ii]) {
                warning("only the first ", ncolBlist[ii], " values of ",
                        "'spar' are used for variable '", s.xargument, "'")
            }
            sparv[[ii]] <- rep(temp, length=ncolBlist[ii])   # recycle
    
            temp <- dfvec[[ii]]
            if (!is.numeric(temp) || any(temp < 1)) {
                stop("df is non-numeric or less than 1")
            }
            if (length(temp) > ncolBlist[ii]) {
                warning("only the first", ncolBlist[ii], "values of 'df' ",
                        "are used for variable '", s.xargument, "'")
            }
            dfvec[[ii]] <- rep(temp, length=ncolBlist[ii])    # recycle
            if (max(temp) > smooth.frame$nef[kk]-1) {
                stop("'df' value too high for variable '", s.xargument, "'")
            }
    
            if (any(sparv[[ii]] != 0) && any(dfvec[[ii]] != 4)) {
                stop("cannot specify both 'spar' and 'df'")
            }
        } # End of kk loop

        sparv <- unlist(sparv)
        dfvec <- unlist(dfvec)
        smooth.frame$sparv <- sparv     # original
        smooth.frame$dfvec <- dfvec         # original
    
        if (sum(smooth.frame$dfvec[smooth.frame$sparv == 0]) + pbig >
            smooth.frame$n_lm * sum(ncolBlist[nwhich])) {
            stop("too many parameters/dof for data on hand")
        }
    
        xnrow_X_vlm <- labels(X_vlm_save)[[2]]
        asgn <- attr(X_vlm_save, "assign")
        aa <- NULL
        for(ii in nwhich) {
            aa <- c(aa, xnrow_X_vlm[asgn[[ii]]])
        }
        smooth.frame$ndfsparv <- aa                # Stored here
        smooth.frame$xnrow_X_vlm <- xnrow_X_vlm    # Stored here
        smooth.frame$s.xargument <- s.xargument    # Stored here
    
        smooth.frame$smap=as.vector(cumsum(
            c(1, ncolBlist[nwhich]))[1:length(nwhich)])
    
        smooth.frame$try.sparv <- sparv
        smooth.frame$prev.dof <- dfvec


        smooth.frame$bindex <- as.integer(cumsum(c(1,
            smooth.frame$nknots*ncolBlist[nwhich])))
        smooth.frame$kindex = as.integer(
            cumsum(c(1, 4 + smooth.frame$nknots)))
    } # End of first
    if (sf.only) {
        return(smooth.frame)
    }

    ldk <- 3 * max(ncolBlist[nwhich]) + 1   # 11/7/02


    which <- unlist(which)
    p_lm <- smooth.frame$p_lm
    n_lm <- smooth.frame$n_lm
    dim2wz <- if (is.matrix(wz)) ncol(wz) else 1

    dim1U <- if (is.matrix(Umat)) nrow(Umat) else 1

    nBlist <- names(Blist)
    for(ii in length(nBlist):1) {
        if (!any(nBlist[ii] == nwhich)) {
            Blist[[ii]] <- NULL
        }
    }
    trivc <- trivial.constraints(Blist)

    ncbvec <- ncolBlist[nwhich]
    ncolbmax <- max(ncbvec)




    contr.sp <- list(low = -1.5,## low = 0.      was default till R 1.3.x
                     high = 1.5,
                     tol = 1e-4,## tol = 0.001   was default till R 1.3.x
                     eps = 2e-8,## eps = 0.00244 was default till R 1.3.x
                     maxit = 500 )

  if (FALSE)
    contr.sp <- list(low = -1.5,## low = 0.      was default till R 1.3.x
                     high = 1.5,
                     tol = 0.001,     # was default till R 1.3.x
                     eps = 0.00244,   # was default till R 1.3.x
                     maxit = 500 )

    fit <- dotC(name="Yee_vbfa",  # ---------------------------------
         npetc = as.integer(c(n_lm, p_lm, length(which), se.fit, 0,
               bf.maxit, qrank = 0, M, nbig = n_lm * M, pbig,
               qbig, dim2wz, dim1U, ier=0, ldk=ldk, # ldk may be unused
               contr.sp$maxit, iinfo = 0
               )),
         doubvec = as.double(c(bf.epsilon, resSS=0, unlist(contr.sp[1:4]))),
     as.double(x),
         y = as.double(zedd), wz = as.double(wz),
         dfvec  = as.double(smooth.frame$dfvec),
         lamvec = as.double(smooth.frame$try.sparv),
         sparv  = as.double(smooth.frame$try.sparv),
   as.integer(smooth.frame$o), as.integer(smooth.frame$nef),
         as.integer(which),
   smomat = as.double(smomat), etamat = double(M * n_lm),
   beta = double(pbig),
       varmat = if (se.fit) as.double(smomat) else double(1),
     qr = as.double(X_vlm_save), qraux = double(pbig),
     qpivot = as.integer(1:pbig),
         as.double(Umat),
         as.double(unlist(Blist)),
     as.integer(ncbvec), as.integer(smooth.frame$smap),
      trivc = as.integer(trivc),

         levmat = if (se.fit) as.double(smomat) else double(1), # 20100227

     bcoefficients = double(sum(smooth.frame$nknots * ncbvec)),
         knots = as.double(unlist(smooth.frame$knots)),
     bindex = as.integer(smooth.frame$bindex),
         nknots = as.integer(smooth.frame$nknots),
         kindex = as.integer(smooth.frame$kindex)) # End of dotC

    dim(fit$qr) = dim(X_vlm_save)
    dimnames(fit$qr) = dimnames(X_vlm_save)
    dim(fit$y) = dim(zedd)
    dimnames(fit$y) = dimnames(zedd)
    dim(fit$smomat) = dim(smomat)
    dimnames(fit$smomat) = dimnames(smomat)   # Needed for vgam.nlchisq
    if (se.fit) {
        dim(fit$varmat) = dim(smomat)
        dimnames(fit$varmat) = dimnames(smomat)
        dim(fit$levmat) = dim(smomat)
        dimnames(fit$levmat) = dimnames(smomat)

    }






    if (fit$npetc[14] != 0 || fit$npetc[17] != 0) {
        stop("something went wrong in the C function 'vbfa'")
    }

    fit$etamat = if (M > 1) matrix(fit$etamat, n_lm, M, byrow=TRUE) else
                 c(fit$etamat)  # May no longer be a matrix
    nits <- fit$npetc[5]
    qrank <- fit$npetc[7]


    smooth.frame$try.sparv <- fit$sparv

    change <- abs(smooth.frame$prev.dof - fit$dfvec)/(1+fit$dfvec) > 0.00 &
                  smooth.frame$sparv == 0


    smooth.frame$try.sparv[change] <- 0         # For next time
    smooth.frame$prev.dof <- fit$dfvec

    if ((nits == bf.maxit) & bf.maxit > 1) {
        warning("'s.vam' convergence not obtained in ", bf.maxit,
                " iterations")
    }

    R <- fit$qr[1:pbig, 1:pbig]
    R[lower.tri(R)] <- 0



    Bspline <- vector("list", length(nwhich))
    names(Bspline) <- nwhich
    for(ii in 1:length(nwhich)) {
        ans = fit$bcoeff[(smooth.frame$bindex[ii]):
                         (smooth.frame$bindex[ii+1]-1)]
        ans = matrix(ans, ncol=ncolBlist[nwhich[ii]])
        Bspline[[ii]] =
            new("vsmooth.spline.fit",
                "Bcoefficients" = ans,
                "xmax"          = smooth.frame$xmax[ii],
                "xmin"          = smooth.frame$xmin[ii],
                "knots"         = as.vector(smooth.frame$knots[[ii]]))
    }


    rl <- list(
      Bspline = Bspline,
      coefficients = fit$beta,
      df.residual = n_lm * M - qrank - sum(fit$dfvec - 1),
      fitted.values = fit$etamat,
      nl.df = fit$dfvec - 1,
      qr = list(qr=fit$qr, rank=qrank, qraux=fit$qraux, pivot=fit$qpivot),
      R = R, 
      rank = qrank, 
      residuals = fit$y - fit$etamat,
      rss = fit$doubvec[2],
      smomat = fit$smomat,
      sparv = fit$sparv,
      s.xargument = unlist(smooth.frame$s.xargument))


    names(rl$coefficients) <- smooth.frame$xnrow_X_vlm
    names(rl$sparv) <- smooth.frame$ndfspar
    names(rl$nl.df) <- smooth.frame$ndfspar

    if (se.fit) {
        rl <- c(rl, list(varmat = fit$varmat))
    }
    c(list(smooth.frame = smooth.frame), rl)
}


