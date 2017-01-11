# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.









s.vam <- function(x, zedd, wz, smomat, which, smooth.frame, bf.maxit = 10,
                  bf.epsilon = 0.001, trace = FALSE, se.fit = TRUE,
                  X.vlm.save, Hlist, ncolHlist, M, qbig, Umat,
                  all.knots = FALSE, nk = NULL,
                  sf.only = FALSE) {
  nwhich <- names(which)


  dX.vlm <- as.integer(dim(X.vlm.save))
  pbig <- dX.vlm[2]


  if (!length(smooth.frame$first)) {
    data <- smooth.frame[, nwhich, drop = FALSE]
    smooth.frame <- vgam.match(data, all.knots = all.knots, nk = nk)
    smooth.frame$first <- TRUE  # Only executed at the first time

    dx <- as.integer(dim(x))
    smooth.frame$n.lm <- dx[1]
    smooth.frame$p.lm <- dx[2]
    attr(data, "class") <- NULL

    osparv <- lapply(data, attr, "spar")  # "o" for original
    odfvec <- lapply(data, attr, "df")
    s.xargument <- lapply(data, attr, "s.xargument")

    for (kk in seq_along(nwhich)) {
      ii <- nwhich[kk]

      temp <- osparv[[ii]]
      if (!is.numeric(temp) || any(temp < 0)) {
        stop("spar cannot be negative or non-numeric")
      }
      if (length(temp) > ncolHlist[ii]) {
        warning("only the first ", ncolHlist[ii], " values of ",
                "'spar' are used for variable '", s.xargument, "'")
      }
      osparv[[ii]] <- rep_len(temp, ncolHlist[ii])  # Recycle

      temp <- odfvec[[ii]]
      if (!is.numeric(temp) || any(temp < 1)) {
        stop("argument 'df' is non-numeric or less than 1")
      }
      if (length(temp) > ncolHlist[ii]) {
        warning("only the first ", ncolHlist[ii], " value(s) of 'df' ",
                "are used for variable '", s.xargument, "'")
      }
      odfvec[[ii]] <- rep_len(temp, ncolHlist[ii])  # Recycle
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
      smooth.frame$n.lm * sum(ncolHlist[nwhich])) {
      stop("too many parameters/dof for data on hand")
    }

    xnrow.X.vlm <- labels(X.vlm.save)[[2]]
    asgn <- attr(X.vlm.save, "assign")
    aa <- NULL
    for (ii in nwhich) {
      aa <- c(aa, xnrow.X.vlm[asgn[[ii]]])
    }
    smooth.frame$ndfsparv <- aa                # Stored here
    smooth.frame$xnrow.X.vlm <- xnrow.X.vlm    # Stored here
    smooth.frame$s.xargument <- s.xargument    # Stored here

    smooth.frame$smap <-
      as.vector(cumsum(c(1, ncolHlist[nwhich]))[seq_along(nwhich)])

    smooth.frame$try.sparv <- osparv


    smooth.frame$bindex <-
      as.integer(cumsum(c(1, smooth.frame$nknots * ncolHlist[nwhich])))


    smooth.frame$lindex <-
      as.integer(cumsum(c(1, smooth.frame$neffec * ncolHlist[nwhich])))


    smooth.frame$kindex <-
      as.integer(cumsum(c(1, 4 + smooth.frame$nknots)))
  } else {
    smooth.frame$first <- FALSE
  }




  if (sf.only) {
    return(smooth.frame)
  }



  ldk <- 3 * max(ncolHlist[nwhich]) + 1  # 20020711

  which <- unlist(which)
  p.lm <- smooth.frame$p.lm
  n.lm <- smooth.frame$n.lm
  dim2wz <- if (is.matrix(wz)) ncol(wz) else 1

  dim1U <- if (is.matrix(Umat)) nrow(Umat) else 1

  nHlist <- names(Hlist)
  for (ii in length(nHlist):1) {
    if (!any(nHlist[ii] == nwhich)) {
      Hlist[[ii]] <- NULL
    }
  }
  trivc <- trivial.constraints(Hlist)

  ncbvec <- ncolHlist[nwhich]
  ncolbmax <- max(ncbvec)


  contr.sp <- list(low   = -1.5,  ## low = 0.      was default till R 1.3.x
                   high  =  1.5,
                   tol   = 1e-4,  ## tol = 0.001   was default till R 1.3.x
                   eps   = 2e-8,  ## eps = 0.00244 was default till R 1.3.x
                   maxit =  500)


  fit <-
    .C("Yee_vbfa",  # ---------------------------------
         npetc = as.integer(c(n.lm, p.lm, length(which), se.fit, 0,
                              bf.maxit, qrank = 0, M, nbig = n.lm * M, pbig,
                              qbig, dim2wz, dim1U, ier = 0,
                              ldk = ldk,  # ldk may be unused
                              contr.sp$maxit, iinfo = 0
                             )),
       doubvec = as.double(c(bf.epsilon, resSS = 0, unlist(contr.sp[1:4]))),
     as.double(x),
         y = as.double(zedd), wz = as.double(wz),
         dfvec  = as.double(smooth.frame$odfvec + 1),  # 20130427; + 1 added
         lamvec = double(length(smooth.frame$odfvec)),
         sparv  = as.double(smooth.frame$try.sparv),
   as.integer(smooth.frame$matcho), as.integer(smooth.frame$neffec),
         as.integer(which),
   smomat = as.double(smomat), etamat = double(M * n.lm),
   beta = double(pbig),
       varmat = if (se.fit) as.double(smomat) else double(1),
     qr = as.double(X.vlm.save), qraux = double(pbig),
     qpivot = as.integer(1:pbig),
         as.double(Umat),
         as.double(unlist(Hlist)),
     as.integer(ncbvec), as.integer(smooth.frame$smap),
      trivc = as.integer(trivc),

         levmat = double(sum(smooth.frame$neffec * ncbvec)),  # 20130427;


     bcoefficients = double(sum(smooth.frame$nknots * ncbvec)),
         knots = as.double(unlist(smooth.frame$knots)),
     bindex = as.integer(smooth.frame$bindex),
     lindex = as.integer(smooth.frame$lindex),
         nknots = as.integer(smooth.frame$nknots),
         kindex = as.integer(smooth.frame$kindex))  # End of dotC


  if (exists("flush.console"))
    flush.console()




  dim(fit$qr) <- dim(X.vlm.save)
  dimnames(fit$qr) <- dimnames(X.vlm.save)
  dim(fit$y) <- dim(zedd)
  dimnames(fit$y) <- dimnames(zedd)
  dim(fit$smomat) <- dim(smomat)
  dimnames(fit$smomat) <- dimnames(smomat)  # Needed for vgam.nlchisq
  if (se.fit) {
    dim(fit$varmat) <- dim(smomat)
    dimnames(fit$varmat) <- dimnames(smomat)
  }

  if (fit$npetc[14] != 0 ||
      fit$npetc[17] != 0) {
    stop("something went wrong in the C function 'vbfa'")
  }

  fit$etamat <- if (M > 1)
                matrix(fit$etamat, n.lm, M, byrow = TRUE) else
                c(fit$etamat)  # May no longer be a matrix
  nits <- fit$npetc[5]
  qrank <- fit$npetc[7]




  if (smooth.frame$first) {
    smooth.frame$try.sparv <- fit$sparv
  }

  if ((nits == bf.maxit) && bf.maxit > 1) {
    warning("'s.vam()' convergence not obtained in ", bf.maxit,
            " iterations")
  }

  R <- fit$qr[1:pbig, 1:pbig]
  R[lower.tri(R)] <- 0



  Bspline <- vector("list", length(nwhich))
  names(Bspline) <- nwhich
  for (ii in seq_along(nwhich)) {
    b.coefs <- fit$bcoeff[(smooth.frame$bindex[ii]):
                          (smooth.frame$bindex[ii + 1] - 1)]
    b.coefs <- matrix(b.coefs, ncol = ncolHlist[nwhich[ii]])
    Bspline[[ii]] <-
        new("vsmooth.spline.fit",
            "Bcoefficients" = b.coefs,
            "xmax"          = smooth.frame$xmax[ii],
            "xmin"          = smooth.frame$xmin[ii],
            "knots"         = as.vector(smooth.frame$knots[[ii]]))
  }




  Leverages <- vector("list", length(nwhich))
  names(Leverages) <- nwhich
  for (ii in seq_along(nwhich)) {
    levvec <- fit$levmat[(smooth.frame$lindex[ii]):
                         (smooth.frame$lindex[ii+1]-1)]
    levmat <- matrix(levvec,
                     nrow = smooth.frame$neffec[ii],
                     ncol = ncolHlist[nwhich[ii]])
    Leverages[[ii]] <- levmat
  }



  nl.df <- fit$dfvec - 1  # Decrement/increment ?

  retlist <- list(
    Bspline = Bspline,
    coefficients = fit$beta,
    df.residual = n.lm * M - qrank - sum(nl.df),  # Decrement/increment ?
    fitted.values = fit$etamat,
    Leverages = Leverages,
    nl.df = nl.df,
    qr = list(qr = fit$qr, rank = qrank,
              qraux = fit$qraux, pivot = fit$qpivot),
    R = R,
    rank = qrank,
    residuals = fit$y - fit$etamat,
    ResSS = fit$doubvec[2],
    smomat = fit$smomat,
    sparv = fit$sparv,
    s.xargument = unlist(smooth.frame$s.xargument))


  names(retlist$coefficients) <- smooth.frame$xnrow.X.vlm
  names(retlist$sparv) <-
  names(retlist$nl.df) <- smooth.frame$ndfspar

  if (se.fit) {
    retlist <- c(retlist, list(varmat = fit$varmat))
  }

  c(list(smooth.frame = smooth.frame), retlist)
}


