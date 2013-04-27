# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.









s.vam <- function(x, zedd, wz, smomat, which, smooth.frame, bf.maxit = 10,
                  bf.epsilon = 0.001, trace = FALSE, se.fit = TRUE,
                  X_vlm_save, Blist, ncolBlist, M, qbig, Umat,
                  all.knots = FALSE, nk = NULL,
                  sf.only = FALSE) {
  nwhich <- names(which)


  dX_vlm <- as.integer(dim(X_vlm_save))
  pbig <- dX_vlm[2]


  if (!length(smooth.frame$first)) {
    data <- smooth.frame[, nwhich, drop = FALSE]
    smooth.frame <- vgam.match(data, all.knots = all.knots, nk = nk)
    smooth.frame$first <- TRUE  # Only executed at the first time

    dx <- as.integer(dim(x))
    smooth.frame$n_lm <- dx[1]
    smooth.frame$p_lm <- dx[2]
    attr(data, "class") <- NULL

    osparv <- lapply(data, attr, "spar")  # "o" for original
    odfvec <- lapply(data, attr, "df")
    s.xargument <- lapply(data, attr, "s.xargument")
    
    for (kk in 1:length(nwhich)) {
      ii <- nwhich[kk]

      temp <- osparv[[ii]]
      if (!is.numeric(temp) || any(temp < 0)) {
        stop("spar cannot be negative or non-numeric")
      }
      if (length(temp) > ncolBlist[ii]) {
        warning("only the first ", ncolBlist[ii], " values of ",
                "'spar' are used for variable '", s.xargument, "'")
      }
      osparv[[ii]] <- rep(temp, length = ncolBlist[ii])   # recycle
    
      temp <- odfvec[[ii]]
      if (!is.numeric(temp) || any(temp < 1)) {
        stop("df is non-numeric or less than 1")
      }
      if (length(temp) > ncolBlist[ii]) {
        warning("only the first ", ncolBlist[ii], " value(s) of 'df' ",
                "are used for variable '", s.xargument, "'")
      }
      odfvec[[ii]] <- rep(temp, length = ncolBlist[ii]) # recycle
      if (max(temp) > smooth.frame$neffec[kk]-1) {
        stop("'df' value too high for variable '", s.xargument, "'")
      }

      if (any(osparv[[ii]] != 0) &&
          any(odfvec[[ii]] != 4)) {
        stop("cannot specify both 'spar' and 'df'")
      }
    }  # End of kk loop


    osparv <- unlist(osparv)
    odfvec <- unlist(odfvec)
    smooth.frame$osparv <- osparv  # Original
    smooth.frame$odfvec <- odfvec  # Original
    
    if (sum(smooth.frame$dfvec[smooth.frame$osparv == 0]) + pbig >
      smooth.frame$n_lm * sum(ncolBlist[nwhich])) {
      stop("too many parameters/dof for data on hand")
    }
    
    xnrow_X_vlm <- labels(X_vlm_save)[[2]]
    asgn <- attr(X_vlm_save, "assign")
    aa <- NULL
    for (ii in nwhich) {
      aa <- c(aa, xnrow_X_vlm[asgn[[ii]]])
    }
    smooth.frame$ndfsparv <- aa                # Stored here
    smooth.frame$xnrow_X_vlm <- xnrow_X_vlm    # Stored here
    smooth.frame$s.xargument <- s.xargument    # Stored here

    smooth.frame$smap <- as.vector(cumsum(
        c(1, ncolBlist[nwhich]))[1:length(nwhich)])

    smooth.frame$try.sparv <- osparv
    smooth.frame$lamvector <- double(length(odfvec))


    smooth.frame$bindex <-
      as.integer(cumsum(c(1, smooth.frame$nknots * ncolBlist[nwhich])))


    smooth.frame$lindex <-
      as.integer(cumsum(c(1, smooth.frame$neffec * ncolBlist[nwhich])))



    smooth.frame$kindex <-
      as.integer(cumsum(c(1, 4 + smooth.frame$nknots)))
  } else {
    smooth.frame$first <- FALSE
  }


  if (sf.only) {
    return(smooth.frame)
  }


  ldk <- 3 * max(ncolBlist[nwhich]) + 1  # 20020711


  which <- unlist(which)
  p_lm <- smooth.frame$p_lm
  n_lm <- smooth.frame$n_lm
  dim2wz <- if (is.matrix(wz)) ncol(wz) else 1

  dim1U <- if (is.matrix(Umat)) nrow(Umat) else 1

  nBlist <- names(Blist)
  for (ii in length(nBlist):1) {
    if (!any(nBlist[ii] == nwhich)) {
      Blist[[ii]] <- NULL
    }
  }
  trivc <- trivial.constraints(Blist)

  ncbvec <- ncolBlist[nwhich]
  ncolbmax <- max(ncbvec)





  contr.sp <- list(low   = -1.5,  ## low = 0.      was default till R 1.3.x
                   high  =  1.5,
                   tol   = 1e-4,  ## tol = 0.001   was default till R 1.3.x
                   eps   = 2e-8,  ## eps = 0.00244 was default till R 1.3.x
                   maxit =  500)




  fit <-
    dotC(name = "Yee_vbfa",  # ---------------------------------
         npetc = as.integer(c(n_lm, p_lm, length(which), se.fit, 0,
                              bf.maxit, qrank = 0, M, nbig = n_lm * M, pbig,
                              qbig, dim2wz, dim1U, ier = 0,
                              ldk = ldk, # ldk may be unused
                              contr.sp$maxit, iinfo = 0
                             )),
       doubvec = as.double(c(bf.epsilon, resSS = 0, unlist(contr.sp[1:4]))),
     as.double(x),
         y = as.double(zedd), wz = as.double(wz),
         dfvec  = as.double(smooth.frame$odfvec + 1),  # 20130427; + 1 added
         lamvec = as.double(smooth.frame$lamvector),
         sparv  = as.double(smooth.frame$try.sparv),
   as.integer(smooth.frame$matcho), as.integer(smooth.frame$neffec),
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






         levmat = double(sum(smooth.frame$neffec * ncbvec)),  # 20130427;


     bcoefficients = double(sum(smooth.frame$nknots * ncbvec)),
         knots = as.double(unlist(smooth.frame$knots)),
     bindex = as.integer(smooth.frame$bindex),
     lindex = as.integer(smooth.frame$lindex),
         nknots = as.integer(smooth.frame$nknots),
         kindex = as.integer(smooth.frame$kindex)) # End of dotC

  if (exists("flush.console")) flush.console()
 

  if (smooth.frame$first) {
  }


  dim(fit$qr) <- dim(X_vlm_save)
  dimnames(fit$qr) <- dimnames(X_vlm_save)
  dim(fit$y) <- dim(zedd)
  dimnames(fit$y) <- dimnames(zedd)
  dim(fit$smomat) <- dim(smomat)
  dimnames(fit$smomat) <- dimnames(smomat)   # Needed for vgam.nlchisq
  if (se.fit) {
    dim(fit$varmat) <- dim(smomat)
    dimnames(fit$varmat) <- dimnames(smomat)
  }






  if (fit$npetc[14] != 0 ||
      fit$npetc[17] != 0) {
    stop("something went wrong in the C function 'vbfa'")
  }

  fit$etamat <- if (M > 1) matrix(fit$etamat, n_lm, M, byrow = TRUE) else
                           c(fit$etamat)  # May no longer be a matrix
  nits <- fit$npetc[5]
  qrank <- fit$npetc[7]













  if (smooth.frame$first) {
    smooth.frame$try.sparv <- fit$sparv
  }



  if ((nits == bf.maxit) & bf.maxit > 1) {
    warning("'s.vam' convergence not obtained in ", bf.maxit,
            " iterations")
  }

  R <- fit$qr[1:pbig, 1:pbig]
  R[lower.tri(R)] <- 0



  Bspline <- vector("list", length(nwhich))
  names(Bspline) <- nwhich
  for (ii in 1:length(nwhich)) {
    b_coefs <- fit$bcoeff[(smooth.frame$bindex[ii]):
                          (smooth.frame$bindex[ii+1]-1)]
    b_coefs <- matrix(b_coefs, ncol = ncolBlist[nwhich[ii]])
    Bspline[[ii]] <-
        new("vsmooth.spline.fit",
            "Bcoefficients" = b_coefs,
            "xmax"          = smooth.frame$xmax[ii],
            "xmin"          = smooth.frame$xmin[ii],
            "knots"         = as.vector(smooth.frame$knots[[ii]]))
  }




  Leverages <- vector("list", length(nwhich))
  names(Leverages) <- nwhich
  for (ii in 1:length(nwhich)) {
    levvec <- fit$levmat[(smooth.frame$lindex[ii]):
                         (smooth.frame$lindex[ii+1]-1)]
    levmat <- matrix(levvec,
                     nrow = smooth.frame$neffec[ii],
                     ncol = ncolBlist[nwhich[ii]])
    Leverages[[ii]] <- levmat
  }



  nl.df <- fit$dfvec - 1  # Used to be -1; Decrement/increment ?


  retlist <- list(
    Bspline = Bspline,
    coefficients = fit$beta,
    df.residual = n_lm * M - qrank - sum(nl.df),  # Decrement/increment ?
    fitted.values = fit$etamat,
    Leverages = Leverages,
    nl.df = nl.df,
    qr = list(qr = fit$qr, rank = qrank,
              qraux = fit$qraux, pivot = fit$qpivot),
    R = R, 
    rank = qrank, 
    residuals = fit$y - fit$etamat,
    rss = fit$doubvec[2],
    smomat = fit$smomat,
    sparv = fit$sparv,
    s.xargument = unlist(smooth.frame$s.xargument))


  names(retlist$coefficients) <- smooth.frame$xnrow_X_vlm
  names(retlist$sparv) <-
  names(retlist$nl.df) <- smooth.frame$ndfspar

  if (se.fit) {
    retlist <- c(retlist, list(varmat = fit$varmat))
  }

  c(list(smooth.frame = smooth.frame), retlist)
}


