# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.





























Trunc <- function(Range, mux = 2, location = 0, omits = TRUE) {
  if (!is.finite(mux) || length(mux) != 1 || round(mux) != mux ||
      mux < 1)
    stop("argument 'mux' must be a positive integer")
  if (any(!is.finite(Range)) || length(Range) < 2 ||
      !all(round(Range) == Range))
    stop("bad input in argument 'Range'")
  if (length(Range) > 2)
    Range <- range(Range, na.rm = TRUE)  # For vglm() and na.omit()
  Min <- Range[1]
  Max <- Range[2]

  allx <- location + (mux * Min):(mux * Max)
  multiples <- location + mux * (Min:Max)
  if (omits) {
    ans <- setdiff(allx, multiples)
    if (length(ans)) ans else NULL
  } else {
    multiples
  }
}  # Trunc












 gaitzeta <-
  function(alt.mix = NULL, inf.mix = NULL,  # Unstructured probs are
           alt.mlm = NULL, inf.mlm = NULL,  # contiguous
           truncate = NULL, max.support = Inf,
           zero = c("pobs", "pstr"),  # Pruned later, handles all four
           eq.ap = FALSE,  # TRUE applies to the intercept, g_a(theta_a)
           eq.ip = FALSE,  # TRUE applies to the intercept, g_i(theta_i)
           parallel.ap = FALSE,  # TRUE applies to the intercept
           parallel.ip = FALSE,  # TRUE applies to the intercept
           lshape.p = "loglink",
           lshape.a = "loglink",
           lshape.i = "loglink",
           type.fitted = c("mean", "shapes",
                           "pobs.mlm", "pstr.mlm",
                           "pobs.mix", "pstr.mix",
                           "Pobs.mix", "Pstr.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "ptrunc.p", "cdf.max.s"),
           gshape.p = 1 + exp(-seq(7)),
           gpstr.mix = ppoints(9) / 2,
           gpstr.mlm = ppoints(9) / (2 + length(inf.mlm)),
           imethod = 1,
           imux = 0.5,  # General downward multiplier for init values
           ishape.p = NULL, ishape.a = ishape.p,
           ishape.i = ishape.p,
           ipobs.mix = NULL, ipstr.mix = NULL,  # 0.25, 
           ipobs.mlm = NULL, ipstr.mlm = NULL,  # 0.25, 
           byrow.ai = FALSE,
           ishrinkage = 0.95,
           probs.y = 0.35) {
  lowsup <- 1
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support, min.support = lowsup)
  lalt.mix <- length(alt.mix)
  linf.mix <- length(inf.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mlm <- length(inf.mlm)
  ltruncat <- length(truncate)
  ltrunc.use <- ltruncat > 0 || !is.infinite(max.support) 


  if (is.finite(max.support))
    stop("argument 'max.support' must be 'Inf'.")

  lshape.p <- as.list(substitute(lshape.p))
  eshape.p <- link2list(lshape.p)
  lshape.p <- attr(eshape.p, "function.name")
  lshape.p.save <- lshape.p
  gshape.p.save <- gshape.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lshape.a <- as.list(substitute(lshape.a))
  eshape.a <- link2list(lshape.a)
  lshape.a <- attr(eshape.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lshape.i <- as.list(substitute(lshape.i))
  eshape.i <- link2list(lshape.i)
  lshape.i <- attr(eshape.i, "function.name")


  if (is.vector(zero) && is.character(zero) && length(zero) == 2) {
    if (linf.mix + linf.mlm == 0)
      zero <- setdiff(zero, "pstr")
    if (lalt.mix + lalt.mlm == 0)
      zero <- setdiff(zero, "pobs")
  }


  lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
  if (lall.len + ltruncat == 0 && is.infinite(max.support))
    return(eval(substitute(
           zetaff(lshape = .lshape.p.save ,
                  gshape = .gshape.p.save ,
                  zero = NULL),
            list( .lshape.p.save = lshape.p.save,
                  .gshape.p.save = gshape.p.save  ))))

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")
  if (!is.logical(parallel.ap) || length(parallel.ap) != 1)
    stop("argument 'parallel.ap' must be a single logical")
  if (!is.logical(parallel.ip) || length(parallel.ip) != 1)
    stop("argument 'parallel.ip' must be a single logical")


  if (lalt.mix == 1 && eq.ap)
    warning("Only one unstructured altered value (no 'shape.a')",
            ", so setting 'eq.ap = TRUE' is meaningless")
  if (linf.mix == 1 && eq.ip)
    warning("Only one unstructured inflated value (no 'shape.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (lalt.mlm == 1 && parallel.ap)  # Only \omega_1
    warning("Only one altered mixture probability, 'pobs", alt.mlm,
            "', so setting 'parallel.ap = TRUE' is meaningless")
  if (linf.mlm == 1 && parallel.ip)  # Only \phi_1
    warning("Only one inflated mixture probability, 'pstr", inf.mlm,
            "', so setting 'parallel.ip = TRUE' is meaningless")


  type.fitted.choices <-
            c("mean", "shapes",
              "pobs.mlm", "pstr.mlm",
              "pobs.mix", "pstr.mix",
              "Pobs.mix", "Pstr.mix",
              "nonspecial", "Numer", "Denom.p",
              "sum.mlm.i", "ptrunc.p", "cdf.max.s")
  type.fitted <- match.arg(type.fitted[1], type.fitted.choices)[1]

  tmp7a <- if (lalt.mlm) paste0("pobs.mlm", alt.mlm) else NULL
  tmp7b <- if (linf.mlm) paste0("pstr.mlm", inf.mlm) else NULL
  tmp3 <- c(shape.p = lshape.p,
            pobs.mix = if (lalt.mix) "multilogitlink" else NULL,
            shape.a = if (lalt.mix > 1) lshape.a else NULL,
            pstr.mix = if (linf.mix) "multilogitlink" else NULL,
            shape.i = if (linf.mix > 1) lshape.i else NULL,
            if (lalt.mlm) rep("multilogitlink", lalt.mlm) else NULL,
            if (linf.mlm) rep("multilogitlink", linf.mlm) else NULL)
  Ltmp3 <- length(tmp3) 
  if (lalt.mlm + linf.mlm)
    names(tmp3)[(Ltmp3 - lalt.mlm - linf.mlm + 1):Ltmp3] <-
      c(tmp7a, tmp7b)
  par1or2 <- 1  # 2
  tmp3.TF <- c(TRUE, lalt.mix > 0, lalt.mix > 1,
                     linf.mix > 0, linf.mix > 1,
                     lalt.mlm > 0, linf.mlm > 0)
  indeta.finish <- cumsum(c(par1or2, 1, par1or2,
                                     1, par1or2,
                            lalt.mlm, linf.mlm,
                            linf.mlm + 1) * c(tmp3.TF, 1)) 
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  rownames(indeta) <- c("shape.p", "pobs.mix", "shape.a",
                        "pstr.mix", "shape.i", "pobs.mlm", "pstr.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.
      

  blurb1 <- "Z"
  if (lalt.mlm + lalt.mix) blurb1 <- "Generally-altered Z"
  if (linf.mlm + linf.mix) blurb1 <- "Generally-inflated Z"
  if (ltrunc.use) blurb1 <- "Generally-truncated Z"
  if ( (lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered and -inflated Z"
  if ( (lalt.mlm + lalt.mix) && !(linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered and -truncated Z"
  if (!(lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated and -truncated Z"
  if ( (lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated and -truncated Z"

      
  new("vglmff",
  blurb = c(blurb1, "eta regression\n",
            "(GAIT-Zeta(shape.p)-Zeta(shape.a)-MLM-",
                  "Zeta(shape.i)-MLM generally)\n\n",
            "Links: ",
            namesof("shape.p", lshape.p, earg = eshape.p, tag = FALSE),
            if (lalt.mix > 0) c(", ", "multilogit(pobs.mix)"),
            if (lalt.mix > 1) c(", ",
            namesof("shape.a",  lshape.a, eshape.a, tag = FALSE)),
            if (lalt.mix && linf.mix) ", \n       ",
            if (linf.mix > 0) c(  if (lalt.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (linf.mix > 1) c(", ",
            namesof("shape.i",  lshape.i, eshape.i, tag = FALSE)),
            if (lalt.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7a, collapse = "),\n"),
            ")") else NULL,
            if (linf.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7b, collapse = "),\n"),
            ")") else NULL),
  constraints = eval(substitute(expression({
    M1 <- max(extra$indeta, na.rm = TRUE)
    lalt.mix <- ( .lalt.mix )
    linf.mix <- ( .linf.mix )
    lalt.mlm <- ( .lalt.mlm )
    linf.mlm <- ( .linf.mlm )

    use.mat.mlm.a <- if (lalt.mlm) {
      if ( .parallel.ap ) matrix(1, lalt.mlm, 1) else diag(lalt.mlm)
    } else {
      NULL
    }

    use.mat.mlm.i <- if (linf.mlm) {
       if ( .parallel.ip ) matrix(1, linf.mlm, 1) else diag(linf.mlm)
    } else {
      NULL
    }

    if (lalt.mlm + linf.mlm == 0) {
      use.mat <- use.mat.mlm <- cbind(M)  # shape.p only
    }
    if (lalt.mlm + linf.mlm) {
      nc1 <- if (length(use.mat.mlm.a)) ncol(use.mat.mlm.a) else 0
      nc2 <- if (length(use.mat.mlm.i)) ncol(use.mat.mlm.i) else 0
      use.mat.mlm <- cbind(1, matrix(0, 1, nc1 + nc2))
      if (lalt.mlm)
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, lalt.mlm, 1),
                                   use.mat.mlm.a,
                                   if (length(use.mat.mlm.i) == 0) NULL
                                   else matrix(0, lalt.mlm, nc2)))
      if (linf.mlm )
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, linf.mlm, 1 + nc1),
                                   use.mat.mlm.i))
      use.mat <- use.mat.mlm
    }  # lalt.mlm + linf.mlm






    use.mat.mix <- if ( ( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,1,0,1,  0,1,0,0,0,  0,0,0,1,0), 5, 3) else
    if (!( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,0,0,1,  0,1,0,0,0,  0,0,1,0,0,
               0,0,0,1,0), 5, 4) else
    if ( ( .eq.ap ) && !( .eq.ip ))
      matrix(c(1,0,1,0,0,  0,1,0,0,0,  0,0,0,1,0,
               0,0,0,0,1), 5, 4) else
      diag(5)

    tmp3.TF <- ( .tmp3.TF )
    tmp3.TF <- tmp3.TF[-(6:7)]  # zz for now. temporary only.
    use.mat.mix <- use.mat.mix[tmp3.TF, , drop = FALSE]

    mincl <- apply(use.mat.mix, 2, min)
    maxcl <- apply(use.mat.mix, 2, max)
    use.mat.mix <- use.mat.mix[, mincl != 0 | maxcl != 0, drop = FALSE]

    if (lalt.mix + linf.mix > 0)
      use.mat <- use.mat.mix





    if (lalt.mlm + linf.mlm > 0 &&
        lalt.mix + linf.mix > 0) {
      use.mat <- rbind(use.mat.mix,
                 matrix(0, nrow(use.mat.mlm) - 1, ncol(use.mat.mix)))
      use.mat <- cbind(use.mat, matrix(0, nrow(use.mat),
                                          ncol(use.mat.mlm) - 1))
      use.mat[row(use.mat) > nrow(use.mat.mix) &
              col(use.mat) > ncol(use.mat.mix)] <- use.mat.mlm[-1, -1]
    }  # lalt.mlm + linf.mlm > 0 && lalt.mix + linf.mix > 0



    if (is.null(constraints)) {
      constraints <-
        cm.VGAM(use.mat, x = x, apply.int = TRUE,  # FALSE
                bool = .eq.ap || .eq.ip || .parallel.ap || .parallel.ip ,
                constraints = constraints)  # FALSE
    }



    if (lalt.mix + linf.mix + lalt.mlm + linf.mlm)
      constraints <-
        cm.zero.VGAM(constraints, x = x, .zero , M = M, M1 = M1,
                     predictors.names = paste0(predictors.names,
                                        names(predictors.names)))
  }), list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF,
            .eq.ap = eq.ap, .eq.ip = eq.ip,
            .parallel.ap = parallel.ap, .parallel.ip = parallel.ip,
            .lalt.mlm = lalt.mlm, .linf.mlm = linf.mlm,
            .lalt.mix = lalt.mix, .linf.mix = linf.mix
          ))),
  infos = eval(substitute(function(...) {
    list(M1 = .M1 ,
         Q1 = 1,
         link = .predictors.names ,  # .tmp3, as.vector strips names off
         link1parameter = as.logical( .lall.len <= 2),  # <= 1 safer
         mixture.links  = any(c( .lalt.mlm , .linf.mlm , .lalt.mix ,
                                 .linf.mix ) > 1),  # FALSE if NULL
         alt.mix = as.vector( .alt.mix ),  # Handles NULL
         alt.mlm = as.vector( .alt.mlm ),
         inf.mix = as.vector( .inf.mix ),
         inf.mlm = as.vector( .inf.mlm ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ),
         Support  = c( .lowsup , Inf, 1),  # a(b)c format as a,c,b.
         expected = TRUE,
         multipleResponses = FALSE,  # zetaff might be called if TRUE
         parameters.names = names( .predictors.names ),
         parent.name = c("zetaff", "zeta"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
         baseparams.argnames  = "shape",
         MM1 = 1,  # One parameter for 1 response
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .lshape.p = lshape.p, .eshape.p = eshape.p,
           .lshape.a = lshape.a, .eshape.a = eshape.a,
           .lshape.i = lshape.i, .eshape.i = eshape.i,
           .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
           .alt.mix = alt.mix, .inf.mix = inf.mix,
           .lalt.mlm = lalt.mlm, .linf.mlm = linf.mlm,
           .lalt.mix = lalt.mix, .linf.mix = linf.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len
         ))),
  initialize = eval(substitute(expression({
    extra$indeta <- ( .indeta )  # Avoids recomputing it several times
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
    truncate <- as.vector( .truncate )
    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(y)
    M <- NOS * M1

    tmp3.TF <- ( .tmp3.TF )

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Since max.support = 9 is possible
              ncol.y.max = 1,
              out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gaitcombo.check(y, truncate = truncate,
                               alt.mlm = alt.mlm, alt.mix = alt.mix,
                               inf.mlm = inf.mlm, inf.mix = inf.mix,
                               max.support = .max.support ,
                               min.support = .min.support )
    extra$skip.mix.a <- glist$skip.mix.a
    extra$skip.mix.i <- glist$skip.mix.i
    extra$skip.mlm.a <- glist$skip.mlm.a
    extra$skip.mlm.i <- glist$skip.mlm.i


    
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- as.vector( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1
    extra$index.M <- iam(NA, NA, M, both = TRUE)  # Used in @weight


    predictors.names <- ( .predictors.names )  # Got it, named


    if (!length(etastart)) {
      shape.p.init <- if (length( .ishape.p )) .ishape.p else {
        zetaff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
          sum(c(w) * dzeta(x = y, shape = shapeval, log = TRUE))
        }
        shape.p.grid <- ( .gshape.p )
        grid.search(shape.p.grid, objfun = zetaff.Loglikfun,
                    y = y, w = w)
      }
      shape.p.init <- rep(shape.p.init, length = n)

      shape.a.init <-  shape.i.init <- shape.p.init  # Needed
      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(shape.p.init, .lshape.p , earg = .eshape.p ))


      pobs.mix.init <- numeric(n)
      if (tmp3.TF[2]) {  # lalt.mix > 0
        pobs.mix.init <- if (length( .ipobs.mix )) {
        rep_len( .ipobs.mix , n)
      } else {
        is.alt.mix <- rowSums(glist$y0.mix.a) > 0
        rep(sum(w[is.alt.mix]) / sum(w), n)
      }  
      }  # lalt.mix > 0


      if (tmp3.TF[3]) {  # Assign coln 3; lalt.mix > 1
        shape.a.init <- if (length( .ishape.a ))
          rep_len( .ishape.a , n) else shape.p.init  # A vector
        etastart[, 3] <-
          theta2eta(shape.a.init, .lshape.a , earg = .eshape.a )
      }

      pstr.mix.init <- numeric(n)
      try.gridsearch.pstr.mix <- FALSE
      if (tmp3.TF[4]) {  # linf.mix > 0
        pstr.mix.init <- if (length( .ipstr.mix )) {
          rep_len( .ipstr.mix , n)
        } else {
          try.gridsearch.pstr.mix <- TRUE
          numeric(n)  # Overwritten by gridsearch
        }
      }  # linf.mix > 0


      if (tmp3.TF[5]) {  # linf.mix > 1
        shape.i.init <- if (length( .ishape.i ))
          rep_len( .ishape.i , n) else shape.p.init  # A vector
        etastart[, (extra$indeta[5, 'launch'])] <-
          theta2eta(shape.i.init, .lshape.i , earg = .eshape.i )
      }  # linf.mix > 1


      if (tmp3.TF[6]) {  #  lalt.mlm
        if (length( .ipobs.mlm )) {
          pobs.mlm.init <- matrix( .ipobs.mlm , n, lalt.mlm,
                                  byrow = .byrow.ai )
        } else {
          pobs.mlm.init <- colSums(c(w) * extra$skip.mlm.a) / colSums(w)
          pobs.mlm.init <- pobs.mlm.init * as.vector( .imux )
          pobs.mlm.init <- matrix(pobs.mlm.init, n, lalt.mlm, byrow=TRUE)
        }
      } else {
        pobs.mlm.init <- matrix(0, n, 1)
      }


      try.gridsearch.pstr.mlm <- FALSE
      if (tmp3.TF[7]) {  #  linf.mlm
        try.gridsearch.pstr.mlm <- !(length( .ipstr.mlm ))
        pstr.mlm.init <- 0  # Might be overwritten by gridsearch

        if (length( .ipstr.mlm ))
          pstr.mlm.init <- as.vector( .ipstr.mlm )

        pstr.mlm.init <- matrix(pstr.mlm.init, n, linf.mlm,
                                byrow = .byrow.ai )
      } else {
        pstr.mlm.init <- matrix(0, n, 1)
      }






      gaitzeta.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitzeta(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  shape.p    = extraargs$shape.p,
                  shape.a    = extraargs$shape.a,
                  shape.i    = extraargs$shape.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }

 gaitzeta.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitzeta(y, pstr.mlm = pstr.mlm.val,
                  pstr.mix    = extraargs$pstr.mix,  # Differs here
                  shape.p    = extraargs$shape.p,
                  shape.a    = extraargs$shape.a,
                  shape.i    = extraargs$shape.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }

 gaitzeta.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitzeta(y, pstr.mix = pstr.mix.val, pstr.mlm = pstr.mlm.val,
                  shape.p    = extraargs$shape.p,
                  shape.a    = extraargs$shape.a,
                  shape.i    = extraargs$shape.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }



      if (linf.mix + linf.mlm) {
        extraargs <- list(
              shape.p    = shape.p.init,
              shape.a    = shape.a.init,
              shape.i    = shape.i.init,
              alt.mix     = alt.mix,
              alt.mlm     = alt.mlm,
              inf.mix     = inf.mix,
              inf.mlm     = inf.mlm,
              truncate    = truncate,
              max.support = as.vector( .max.support ),
              pobs.mix    = pobs.mix.init ,
              pobs.mlm    = pobs.mlm.init )
        pre.warn <- options()$warn
        options(warn = -1)  # Ignore warnings during gridsearch 

        try.this <-
          if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
            grid.search2( .gpstr.mix ,  .gpstr.mlm ,
                         objfun = gaitzeta.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- pstr.mlm.init
            grid.search ( .gpstr.mix ,
                         objfun = gaitzeta.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- pstr.mix.init
            grid.search ( .gpstr.mlm ,
                         objfun = gaitzeta.Loglikfun1.mlm,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          }

        options(warn = pre.warn)  # Restore warnings 
        if (any(is.na(try.this)))
          warning("gridsearch returned NAs. It's going to crash.",
                  immediate. = TRUE)
        if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
          pstr.mix.init <- rep_len(try.this["Value1"], n)
          pstr.mlm.init <- matrix(try.this["Value2"], n, linf.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mix=seq(5)/100'",
                 " and/or 'gpstr.mlm = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mix) {
          pstr.mix.init <- rep_len(try.this["Value"], n)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mix=seq(5)/100'.")
        } else if (try.gridsearch.pstr.mlm) {
          pstr.mlm.init <- matrix(try.this["Value"], n, linf.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mlm=seq(5)/100'.")
        }
      }  # lalt.mix + lnf.mix



        


      while (any((vecTF <- pobs.mix.init + pstr.mix.init +
                           rowSums(pobs.mlm.init) +
                           rowSums(pstr.mlm.init) > 0.96875))) {
        pobs.mix.init[vecTF]   <- 0.875 * pobs.mix.init[vecTF]
        pstr.mix.init[vecTF]   <- 0.875 * pstr.mix.init[vecTF]
        pobs.mlm.init[vecTF, ] <- 0.875 * pobs.mlm.init[vecTF, ]
        pstr.mlm.init[vecTF, ] <- 0.875 * pstr.mlm.init[vecTF, ]
      }

      Numer <- 1 - rowSums(pobs.mlm.init) - rowSums(pstr.mlm.init) -
               pobs.mix.init - pstr.mix.init
        
      etastart.z <- if (lall.len == 0) NULL else
        multilogitlink(cbind(if (tmp3.TF[2]) pobs.mix.init else NULL,
                             if (tmp3.TF[4]) pstr.mix.init else NULL,
                             if (tmp3.TF[6]) pobs.mlm.init else NULL,
                             if (tmp3.TF[7]) pstr.mlm.init else NULL,
                             Numer))
      if (!is.matrix(etastart.z)) etastart.z <- cbind(etastart.z)

      nextone <- 1  # Might not be used actually
      if (tmp3.TF[2]) {
        etastart[, 2] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[4]) {  # Coln 2 or 4
        etastart[, (extra$indeta[4, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[6]) {
        ind6 <- (extra$indeta[6, 'launch']):(extra$indeta[6, 'finish'])
        etastart[, ind6] <- etastart.z[, nextone:(nextone+lalt.mlm-1)]
        nextone <- nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind7 <- (extra$indeta[7, 'launch']):(extra$indeta[7, 'finish'])
        etastart[, ind7] <- etastart.z[, nextone:ncol(etastart.z)]
      }
    }
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .ishape.p = ishape.p,
    .ishape.a = ishape.a,
    .ishape.i = ishape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .ipstr.mix = ipstr.mix, .ipobs.mix = ipobs.mix,
    .ipstr.mlm = ipstr.mlm, .ipobs.mlm = ipobs.mlm,
    .byrow.ai = byrow.ai,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support,
    .min.support = lowsup,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .imux = imux,
    .gpstr.mix = gpstr.mix, .gshape.p = gshape.p,
    .gpstr.mlm = gpstr.mlm,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .indeta = indeta,
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"
      }
    type.fitted <-
      match.arg(type.fitted[1],
                c("mean", "shapes",
                  "pobs.mlm", "pstr.mlm",
                  "pobs.mix", "pstr.mix",
                  "Pobs.mix", "Pstr.mix",
                  "nonspecial", "Numer", "Denom.p", "sum.mlm.i",
                  "ptrunc.p", "cdf.max.s"))[1]

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <- shape.p  # Needed; and answer not corrupted

    if (any(tmp3.TF[c(3, 5)])) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      iptr <- 1
      shape.a <- if (!tmp3.TF[3]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[5]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0
    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")

      Nextone <- 0  # Might not be used actually
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitcombo.zeta(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              shape.a = shape.a, shape.i = shape.i,
              truncate = truncate, max.support = max.support)

    n.eta <- nrow(eta)
    Denom.p <- c(Bits[["cdf.max.s"]] - Bits[["SumT0.p"]] -
                 Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    Numer <- c(1 - pobs.mix - pstr.mix -
               (if (lalt.mlm) rowSums(pobs.mlm) else 0) -
               (if (linf.mlm) rowSums(pstr.mlm) else 0))

    if (!lalt.mlm && type.fitted %in% c("pobs.mlm")) {
      warning("No altered MLM values; returning an NA")
      return(NA)
    }
    if (!linf.mlm && type.fitted %in% c("sum.mlm.i", "pstr.mlm")) {
      warning("No inflated MLM values; returning an NA")
      return(NA)
    }
    if (!lalt.mix && type.fitted %in% c("Pobs.mix")) {
      warning("No altered mixture values; returning an NA")
      return(NA)
    }
    if (!linf.mix && type.fitted %in% c("sum.mix.i", "Pstr.mix")) {
      warning("No inflated mixture values; returning an NA")
      return(NA)
    }

    if (lalt.mix) {
      tmp13 <-  # dzeta() does not retain the matrix format
        dzeta(matrix(alt.mix, NROW(eta), lalt.mix, byrow = TRUE),
              matrix(shape.a, NROW(eta), lalt.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(NROW(eta), lalt.mix)
      dimnames(tmp13) <- list(rownames(eta), as.character(alt.mix))
      propn.mat.a <- tmp13
    }

    if (linf.mix) {
      tmp55 <-  # dzeta() does not retain the matrix format
        dzeta(matrix(inf.mix, NROW(eta), linf.mix, byrow = TRUE),
              matrix(shape.i, NROW(eta), linf.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(NROW(eta), linf.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(inf.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "shapes"     = cbind(shape.p,
                           if (tmp3.TF[3]) shape.a else NULL,
                           if (tmp3.TF[5]) shape.i else NULL),
      "pobs.mlm"   = pobs.mlm,  # aka omegamat, n x lalt.mlm
      "pstr.mlm"   = pstr.mlm,  # aka phimat, n x linf.mlm
      "pobs.mix"   = pobs.mix,  # n-vector
      "pstr.mix"   = pstr.mix,  # n-vector
      "Pobs.mix"   = c(pobs.mix) * propn.mat.a,  # matrix
      "Pstr.mix"   = c(pstr.mix) * propn.mat.i,
      "nonspecial" = Numer * (1 -
         (Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) / Denom.p),
      "Numer"      = Numer,
      "Denom.p"    = Denom.p,
      "sum.mlm.i"  = pstr.mlm + Numer *
           dzeta(matrix(inf.mlm, NROW(eta), linf.mlm, byrow = TRUE),
                 matrix(shape.p, NROW(eta), linf.mlm)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)

   ynames.pobs.mlm <- as.character(alt.mlm)  # Works with NULLs
   ynames.pstr.mlm <- as.character(inf.mlm)  # Works with NULLs
   if (length(ans))
     label.cols.y(ans, NOS = NOS, colnames.y =
     switch(type.fitted,
            "shapes"    = c("shape.p", "shape.a",  # Some colns NA
                            "shape.i")[(tmp3.TF[c(1, 3, 5)])],
            "Pobs.mix"  = as.character(alt.mix),
            "Pstr.mix"  = as.character(inf.mix),
            "pobs.mlm"  = ynames.pobs.mlm,
            "sum.mlm.i" = ,  #
            "pstr.mlm"  = ynames.pstr.mlm,
            extra$colnames.y)) else ans
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  last = eval(substitute(expression({
    pred.names <- c( .predictors.names )  # Save it
    link.names <- as.vector( .predictors.names )
    parameter.names <- names(pred.names)
    predictors.names <- NULL
    for (jay in seq(M))
      predictors.names <- c(predictors.names,
        namesof(parameter.names[jay], link.names[jay], tag = FALSE,
                earg = list()))  # This line isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .eshape.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[2])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[3])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.a )
    if (tmp3.TF[4])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[5])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.i )
    if (tmp3.TF[6]) {  # lalt.mlm
      for (ii in seq(lalt.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # lalt.mlm
    if (tmp3.TF[7]) {  # linf.mlm
      for (ii in seq(linf.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # linf.mlm
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .predictors.names = predictors.names,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <- shape.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5)])) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      iptr <- 1
      shape.a <- if (!tmp3.TF[3]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[5]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0
    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitzeta(y, shape.p, log = TRUE,  # byrow.ai = F,
                  alt.mix = alt.mix, inf.mix = inf.mix,
                  alt.mlm = alt.mlm, inf.mlm = inf.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  shape.a = shape.a, shape.i = shape.i, 
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    small. <- 1e-14
    pobs.mix <- pstr.mix <- small.
    pobs.mlm <- pstr.mlm <- matrix(small., NROW(eta), 1) # 4 rowSums()
    shape.a <- shape.i <- 0.5  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    shape.p <-
      cbind(eta2theta(eta[, 1], .lshape.p , earg = .eshape.p ))
    iptr <- 1
    ind.shape.z <- 1  # Points to shape.p only.
    if (lalt.mix + linf.mix > 1) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (lalt.mix <= 1) shape.p else {
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a ,
                  earg = .eshape.a )
      }
      shape.i <- if (linf.mix <= 1) shape.p else {
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i ,
                  earg = .eshape.i )
      }
    }  # lalt.mix + linf.mix > 0

    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm)
    okay.mix <-
      all(is.finite(shape.p)) && all(0 < shape.p) &&
      all(is.finite(shape.a)) && all(0 < shape.a) &&
      all(is.finite(shape.i)) && all(0 < shape.i) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(pobs.mix + pstr.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1) 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    shape.a <- shape.i <- shape.p  # Needed
    ind.shape.z <- 1  # Points to shape.p only.

    if (lalt.mix + linf.mix > 1) {  # At least one shape.[ai]
      ind.shape.z <- object@extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      iptr <- 1
      shape.a <- if (lalt.mix <= 1) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a , .eshape.a )
      shape.i <- if (linf.mix <= 1) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0

    tmp3.TF <- ( .tmp3.TF )
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    rgaitzeta(nsim * length(shape.p), shape.p,
             pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
             pobs.mix = pobs.mix, pstr.mix = pstr.mix,
             shape.a = shape.a, shape.i = shape.i, 
             alt.mlm = alt.mlm, inf.mlm = inf.mlm,
             alt.mix = alt.mix, inf.mix = inf.mix,
             truncate = .truncate , max.support = .max.support )
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .tmp3.TF = tmp3.TF,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({

    tmp3.TF <- ( .tmp3.TF )
    calA.p  <- tmp3.TF[2]
    calI.p  <- tmp3.TF[4]
    calA.np <- tmp3.TF[6]
    calI.np <- tmp3.TF[7]

    Denom1.a <- Denom1.i <- Denom2.i <- 0  # Denom2.a is not needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <- shape.p  # Needed; doesnt corrupt answer

    if (any(tmp3.TF[c(3, 5)])) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      shape.a <- if (!tmp3.TF[3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0

    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
      }
    }  # lall.len


    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")



    is.alt.mixed <- if (tmp3.TF[2])
      rowSums(extra$skip.mix.a) > 0 else rep(FALSE, n)
    is.inf.mixed <- if (tmp3.TF[4])
      rowSums(extra$skip.mix.i) > 0 else rep(FALSE, n)
    is.alt.mlmed <- if (tmp3.TF[6])
      rowSums(extra$skip.mlm.a) > 0 else rep(FALSE, n)
    is.inf.mlmed <- if (tmp3.TF[7])
      rowSums(extra$skip.mlm.i) > 0 else rep(FALSE, n)

    is.ns <- !is.alt.mlmed & !is.inf.mlmed  &
             !is.alt.mixed & !is.inf.mixed  # & !is.truncd


    prob.mlm.a <- if (lalt.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (linf.mlm) rowSums(pstr.mlm) else 0  # scalar okay



   pmf.deriv1 <- function(y, shape) {
      -dzeta(y, shape) * (log(y) +
        zeta(shape + 1, deriv = 1) / zeta(shape + 1))
    }
    pmf.deriv2 <- function(y, shape) {
      tmp2 <- zeta(shape + 1, deriv = 1) / zeta(shape + 1)
      dzeta(y, shape) * ((log(y) + tmp2)^2 -
        zeta(shape + 1, deriv = 2) / zeta(shape + 1) + tmp2^2)
    }
    

    
    sumD.mix.1a.p <- sumD.mix.2a.p <- matrix(0, n, NOS)
    if (lalt.mix > 0) {  # \calA_p
      Da.mix.0mat.a <-  # Matches naming convention further below
      Da.mix.1mat.a <- Da.mix.2mat.a <- matrix(0, n, lalt.mix)
      for (jay in seq(lalt.mix)) {
        aval <- alt.mix[jay]
        sumD.mix.1a.p <- sumD.mix.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mix.2a.p <- sumD.mix.2a.p + pmf.deriv2(aval, shape.p)
        pmf.a <- dzeta(aval, shape.a)
        Da.mix.0mat.a[, jay] <- pmf.a
        Da.mix.1mat.a[, jay] <- pmf.deriv1(aval, shape.a)
      }
      Denom1.a <- rowSums(Da.mix.1mat.a)  # aka sumD.mix.1a.a
      Denom2.a <- rowSums(Da.mix.2mat.a)  # Needed; sumD.mix.2a.a <-
    }  # lalt.mix > 0




    if (linf.mix) {
      Di.mix.0mat.i <-  # wrt inflated distribution
      Di.mix.1mat.i <- Di.mix.2mat.i <- matrix(0, n, linf.mix)
      Dp.mix.0mat.i <-  # wrt parent distribution
      Dp.mix.1mat.i <- Dp.mix.2mat.i <- matrix(0, n, linf.mix)
        for (jay in seq(linf.mix)) {
          ival <- inf.mix[jay]
          pmf.i <- dzeta(ival, shape.i)
          Di.mix.0mat.i[, jay] <- pmf.i
          Di.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.i)
          Di.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.i)
          pmf.p <- dzeta(ival, shape.p)
          Dp.mix.0mat.i[, jay] <- pmf.p
          Dp.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.p)
          Dp.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.p)
        }  # jay
      Denom1.i <- rowSums(Di.mix.1mat.i)
      Denom2.i <- rowSums(Di.mix.2mat.i)
    }  # linf.mix



    bits <- moments.gaitcombo.zeta(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              shape.a = shape.a, shape.i = shape.i,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p <- sumD.mlm.2a.p <- matrix(0, n, NOS)
    if (lalt.mlm)
      for (aval in alt.mlm) {
        sumD.mlm.1a.p <- sumD.mlm.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mlm.2a.p <- sumD.mlm.2a.p + pmf.deriv2(aval, shape.p)
      }


    Denom0.p <- c(bits[["cdf.max.s"]] - bits[["SumT0.p"]] -
                  bits[["SumA0.mix.p"]] - bits[["SumA0.mlm.p"]])
    Numer <- 1 - pobs.mix - pstr.mix - prob.mlm.a - prob.mlm.i
    Denom0.a <- c(bits[["SumA0.mix.a"]])  # Not .p
    Denom0.i <- c(bits[["SumI0.mix.i"]])


 

    Dp.mlm.0Mat.i <-  # wrt parent distribution
    Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, NOS)
    if (linf.mlm > 0) {
      Dp.mlm.0Mat.i <-  # wrt parent distribution
      Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, linf.mlm)
      for (jay in seq(linf.mlm)) {
        ival <- inf.mlm[jay]
        pmf.p <- dzeta(ival, shape.p)
        Dp.mlm.0Mat.i[, jay] <- pmf.p
        Dp.mlm.1Mat.i[, jay] <- pmf.deriv1(ival, shape.p)
        Dp.mlm.2Mat.i[, jay] <- pmf.deriv2(ival, shape.p)
      }  # jay
    }  # linf.mlm


    

    sumD.1t.a <- sumD.2t.a <-
    sumD.1t.i <- sumD.2t.i <-
    sumD.1t.p <- sumD.2t.p <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p <- sumD.1t.p + pmf.deriv1(tval, shape.p)
        sumD.2t.p <- sumD.2t.p + pmf.deriv2(tval, shape.p)
        sumD.1t.a <- sumD.1t.a + pmf.deriv1(tval, shape.a)
        sumD.2t.a <- sumD.2t.a + pmf.deriv2(tval, shape.a)
        sumD.1t.i <- sumD.1t.i + pmf.deriv1(tval, shape.i)
        sumD.2t.i <- sumD.2t.i + pmf.deriv2(tval, shape.i)
      }



      
    if (is.finite(max.support)) {
      stop("argument 'max.support' must be 'Inf'.")
    }  # is.finite(max.support)





    Denom1.p <- c(-sumD.1t.p - sumD.mlm.1a.p - sumD.mix.1a.p)
    Denom2.p <- c(-sumD.2t.p - sumD.mlm.2a.p - sumD.mix.2a.p)


    d0B.pi.mlm <- Dp.mlm.0Mat.i / Denom0.p
    d1B.pi.mlm <- Dp.mlm.1Mat.i / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.i * Denom1.p / Denom0.p^2
    d2B.pi.mlm <-     Dp.mlm.2Mat.i / Denom0.p -
                  2 * Dp.mlm.1Mat.i * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.i * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.i * (Denom1.p^2) / Denom0.p^3


    Delta.mlm <- if (linf.mlm > 0) {
      pstr.mlm + Numer * d0B.pi.mlm  # n x linf.mlm.
    } else {
      matrix(0, n, 1)  # If linf.mlm == 0, for rowSums().
    }


    if (linf.mix > 0) {
      d0A.i <- Di.mix.0mat.i / Denom0.i
      d0B.pi.mix <- Dp.mix.0mat.i / Denom0.p
      Delta.mix <- pstr.mix * d0A.i + Numer * d0B.pi.mix


      d1A.i <- (Di.mix.1mat.i - Di.mix.0mat.i *
                Denom1.i / Denom0.i) / Denom0.i
      d2A.i <- (Di.mix.2mat.i - (2 * Di.mix.1mat.i * Denom1.i +
                Di.mix.0mat.i * Denom2.i) / Denom0.i +
            2 * Di.mix.0mat.i * (Denom1.i / Denom0.i)^2) / Denom0.i

      d1B.pi.mix <- Dp.mix.1mat.i / Denom0.p -
                    Dp.mix.0mat.i * Denom1.p / Denom0.p^2
      d2B.pi.mix <-     Dp.mix.2mat.i / Denom0.p -
                    2 * Dp.mix.1mat.i * Denom1.p / Denom0.p^2 -
                        Dp.mix.0mat.i * Denom2.p / Denom0.p^2 +
                    2 * Dp.mix.0mat.i * (Denom1.p^2) / Denom0.p^3
    }  # linf.mix > 0


    if (lalt.mix) {
      d0A.a <- Da.mix.0mat.a / Denom0.a
      d1A.a <- Da.mix.1mat.a / Denom0.a -
               Da.mix.0mat.a * Denom1.a / Denom0.a^2
      d2A.a <- (Da.mix.2mat.a - (2 * Da.mix.1mat.a * Denom1.a +
                Da.mix.0mat.a * Denom2.a) / Denom0.a +
            2 * Da.mix.0mat.a * (Denom1.a / Denom0.a)^2) / Denom0.a
    }  # lalt.mix




    fred0.p <- zeta(shape.p + 1)
    fred1.p <- zeta(shape.p + 1, deriv = 1)
    dl.dshape.p <- -log(y) - fred1.p / fred0.p  # Usual formula

    dl.dshape.p[!is.ns] <- 0  # For is.alt.mixed & is.alt.mlmed
    dl.dshape.a <-
    dl.dshape.i <- numeric(n)  # Replace some elts below
    dl.dpstr.mix <- (-1) / Numer  # \notin A, I, T
    dl.dpobs.mix <- numeric(n)  # Replace some elts below
    dl.dpobs.mix[is.ns] <- (-1) / Numer[is.ns]
    dl.dpstr.mix[is.alt.mixed] <- 0
    dl.dpstr.mix[is.alt.mlmed] <- 0
    dl.dpobs.mlm <-
    dl.dpstr.mlm <- matrix(0, n, 1)  # Might not be needed



    if (tmp3.TF[6] && lalt.mlm) {  # aka \calA_{np}
      dl.dpobs.mlm <- matrix(-1 / Numer, n, lalt.mlm)  # \notin calS
      dl.dpobs.mlm[!is.ns, ] <- 0  # For alt.mix only really
      for (jay in seq(lalt.mlm)) {
        aval <- alt.mlm[jay]
        is.alt.j.mlm <- extra$skip.mlm.a[, jay]  # Logical vector
        tmp7a <- 1 / pobs.mlm[is.alt.j.mlm, jay]
        dl.dpobs.mlm[is.alt.j.mlm, jay] <- tmp7a
      }  # jay
    }  # lalt.mlm



    dl.dshape.p[is.ns] <- dl.dshape.p[is.ns] -
                           (Denom1.p / Denom0.p)[is.ns]


    
    if (tmp3.TF[7] && linf.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, linf.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For alt.mlm and alt.mix

      for (jay in seq(linf.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7i <- Numer * d1B.pi.mlm[, jay] / Delta.mlm[, jay]
        dl.dshape.p[is.inf.j.mlm] <- tmp7i[is.inf.j.mlm]


        if (tmp3.TF[6] && lalt.mlm) {  # 20200919
          tmp9i <- (-d0B.pi.mlm[, jay] / Delta.mlm[, jay])
          dl.dpobs.mlm[is.inf.j.mlm, ] <- tmp9i[is.inf.j.mlm]
        }


        if (tmp3.TF[2] && lalt.mix) {  # 20200919 linf.mix wrong
          tmp2 <- (-d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
          dl.dpobs.mix[is.inf.j.mlm] <- tmp2[is.inf.j.mlm]
        }


        if (tmp3.TF[4] && linf.mix) {  # 20200916
          tmp2 <- (-d0B.pi.mlm[, jay] / Delta.mlm[, jay])
          dl.dpstr.mix[is.inf.j.mlm] <- tmp2[is.inf.j.mlm]
        }

          
        tmp9 <- (0 - d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
        tmp8 <- (1 - d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
        dl.dpstr.mlm[is.inf.j.mlm, ] <- tmp9[is.inf.j.mlm]
        dl.dpstr.mlm[is.inf.j.mlm, jay] <- tmp8[is.inf.j.mlm]
      }  # jay
    }  # linf.mlm > 0


    



    if (tmp3.TF[2] && lalt.mix) {  # aka \calA_{p}
      dl.dpobs.mix[is.alt.mixed] <- 1 / pobs.mix[is.alt.mixed]

      if (tmp3.TF[3] && lalt.mix > 1)
        for (jay in seq(lalt.mix)) {
          is.alt.j.mix <- extra$skip.mix.a[, jay]  # Logical vector
          tmp2 <- d1A.a[, jay] / d0A.a[, jay]
          dl.dshape.a[is.alt.j.mix] <- tmp2[is.alt.j.mix]  # ccc.
        }  # jay
    }  # lalt.mix


    if (tmp3.TF[4] && linf.mix > 0) {  # aka \calI_{p}
      for (jay in seq(linf.mix)) {
        ival <- inf.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7b <- Numer * d1B.pi.mix[, jay] / Delta.mix[, jay]
        dl.dshape.p[is.inf.j.mix] <- tmp7b[is.inf.j.mix]
        if (tmp3.TF[2] && lalt.mix) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpobs.mix[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }
        tmp8 <- (d0A.i[, jay] - d0B.pi.mix[, jay]) / Delta.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]

        if (linf.mix > 1) {
          tmp2 <- pstr.mix * d1A.i[, jay] / Delta.mix[, jay]
          dl.dshape.i[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }

        if (tmp3.TF[6] && lalt.mlm) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpobs.mlm[is.inf.j.mix, ] <- tmp2[is.inf.j.mix]
        }
        if (tmp3.TF[7] && linf.mlm) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpstr.mlm[is.inf.j.mix, ] <- tmp2[is.inf.j.mix]
        }
      }  # jay
    }  # linf.mix > 0



    tmp3.TF <- !is.na(rowSums(extra$indeta))



    if (lall.len) {  # MLM fitted
      all4.dldp <- cbind(
        if (tmp3.TF[2]) dl.dpobs.mix else NULL,
        if (tmp3.TF[4]) dl.dpstr.mix else NULL,
        if (tmp3.TF[6]) dl.dpobs.mlm else NULL,
        if (tmp3.TF[7]) dl.dpstr.mlm else NULL)

      dl.deta.mlm4 <- matrix(0, n, ncol(all4.dldp))
      ind.quad.mlm <- seq(ncol(all4.dldp))  # Contiguous & starts at 1

      for (jay in ind.quad.mlm) {
        for (sss in ind.quad.mlm) {  # Dont need all4probs, allprobs ok
          dl.deta.mlm4[, jay] <- dl.deta.mlm4[, jay] +
            allprobs[, sss] * ((sss == jay) - allprobs[, jay]) *
            all4.dldp[, sss]
        }  # sss
      }  # jay
    }  # lall.len


    

    dshape.p.deta <- dtheta.deta(shape.p, .lshape.p , .eshape.p )
    if (tmp3.TF[3])
      dshape.a.deta <- dtheta.deta(shape.a, .lshape.a , .eshape.a )
    if (tmp3.TF[5])
      dshape.i.deta <- dtheta.deta(shape.i, .lshape.i , .eshape.i )


    iptr <- 0
    ansd <- cbind(shape.p = c(dl.dshape.p * dshape.p.deta))
    ansd <- cbind(ansd,
      pobs.mix = if (tmp3.TF[2]) dl.deta.mlm4[, (iptr <- iptr+1)],
      shape.a = if (tmp3.TF[3]) dl.dshape.a * dshape.a.deta)
    ansd <- cbind(ansd,
      pstr.mix = if (tmp3.TF[4]) dl.deta.mlm4[, (iptr <- iptr+1)],
      shape.i = if (tmp3.TF[5]) dl.dshape.i * dshape.i.deta)
    if (any(tmp3.TF[6:7])) {  # The remainder
      ansd <- cbind(ansd, dl.deta.mlm4[, (iptr+1):ncol(dl.deta.mlm4)])
      ind3 <- (ncol(ansd) - lalt.mlm - linf.mlm + 1):ncol(ansd)
      colnames(ansd)[ind3] <- as.character(c(alt.mlm, inf.mlm))
    }

    c(w) * ansd
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({

    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    probns <- Numer * (1 - c(bits[["SumI0.mix.p"]] +
                             bits[["SumI0.mlm.p"]]) / Denom0.p)






    zero0n <- numeric(n)
    ned2l.dpobs.mix.shape.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpobs.mix.shape.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.p <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.shape.i <- zero0n  # mB overwritten below


      posn.pobs.mix <- as.vector(extra$indeta[2, 'launch'])
      posn.shape.a <- as.vector(extra$indeta[3, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[4, 'launch'])
      posn.shape.i <- as.vector(extra$indeta[5, 'launch'])
      posn.pobs.mlm <- as.vector(extra$indeta[6, 'launch'])
      posn.pstr.mlm <- as.vector(extra$indeta[7, 'launch'])


    ned2l.dpstr.mix2         <-  # Elt (4, 4)
    ned2l.dpobs.mlm.pstr.mix <-  # Elts (4, >=6)
    ned2l.dpobs.mix.pstr.mix <- probns / Numer^2  # ccc Elt (2, 4)
    if (all(c(lalt.mix, linf.mlm) > 0))
    ned2l.dpobs.mix.pstr.mlm <- matrix(probns / Numer^2, n, linf.mlm)
    if (all(c(linf.mix, linf.mlm) > 0))
    ned2l.dpstr.mix.pstr.mlm <- matrix(probns / Numer^2, n, linf.mlm)
    



    ned2l.dshape.p2 <- probns * (
      zeta(shape.p + 1, deriv = 2) / fred0.p -
      (fred1.p / fred0.p)^2 +  # ccc
      Denom2.p / Denom0.p - (Denom1.p / Denom0.p)^2) + 
      (if (tmp3.TF[4] && linf.mix) Numer *
      rowSums(Numer * (d1B.pi.mix^2) / Delta.mix - d2B.pi.mix) else 0) +
      (if (tmp3.TF[7] && linf.mlm) Numer *
      rowSums(Numer * (d1B.pi.mlm^2) / Delta.mlm - d2B.pi.mlm) else 0)

    wz[, iam(1, 1, M)] <- ned2l.dshape.p2 * dshape.p.deta^2


    ned2l.dpobs.mix2 <- 1 / pobs.mix + probns / Numer^2
    if (tmp3.TF[4] && linf.mix > 0) {
      ned2l.dpobs.mix2 <-  # More just below, ccc
      ned2l.dpobs.mix2 + rowSums(d0B.pi.mix^2 / Delta.mix)
    }
    if (tmp3.TF[7] && linf.mlm > 0) {
      ned2l.dpobs.mix2 <-  # ccc.
      ned2l.dpobs.mix2 + rowSums(d0B.pi.mlm^2 / Delta.mlm)
    }
    if (tmp3.TF[2] && lalt.mix > 0)
      wz[, iam(2, 2, M)] <- ned2l.dpobs.mix2  # Link done later


    if (tmp3.TF[3] && lalt.mix > 1) {
      ned2l.dshape.a2 <- pobs.mix * (
        rowSums((Da.mix.1mat.a^2) / Da.mix.0mat.a) / Denom0.a -
        (Denom1.a / Denom0.a)^2)  # ccc.
      wz[, iam(3, 3, M)] <- ned2l.dshape.a2 * dshape.a.deta^2
    }






    if (tmp3.TF[4] && linf.mix > 0) {

      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.pi.mix)^2 / Delta.mix)


      if (tmp3.TF[2] && lalt.mix > 0)
        ned2l.dpobs.mix.shape.p <-
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ned2l.dpstr.mix.shape.p <-
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.pi.mix * (1 + Numer * (d0A.i - d0B.pi.mix) / Delta.mix))


    if (all(tmp3.TF[c(2, 4)]))
      ned2l.dpobs.mix.pstr.mix <-  # ccc
      ned2l.dpobs.mix.pstr.mix +
        rowSums(-d0B.pi.mix * (d0A.i - d0B.pi.mix) / Delta.mix)

    }  # (tmp3.TF[4] && linf.mix > 0)



    
    if (all(tmp3.TF[c(2, 4, 7)])) {  # was  lalt.mix > 0 & Delta.mix
      ned2l.dpobs.mix.pstr.mix <-  # ccc.
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.pi.mlm^2 / Delta.mlm)
    }

    if (!is.na(posn.pobs.mix) && !is.na(posn.pstr.mix))
      wz[, iam(posn.pobs.mix, posn.pstr.mix, M)] <-
        ned2l.dpobs.mix.pstr.mix  # Link done later




    if (tmp3.TF[5] && linf.mix > 1) {  # \calI_{p}, includes \theta_iota

      ned2l.dshape.p.shape.i <- pstr.mix * Numer *
        rowSums(d1A.i * d1B.pi.mix / Delta.mix)  # ccc.
      wz[, iam(1, posn.shape.i, M)] <- ned2l.dshape.p.shape.i *
        dshape.p.deta * dshape.i.deta  # All links done here

      ned2l.dshape.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i^2) / Delta.mix - d2A.i)  # ccc.
      wz[, iam(posn.shape.i, posn.shape.i, M)] <-
        ned2l.dshape.i2 * dshape.i.deta^2

      if (tmp3.TF[2]) {  # tmp3.TF[4] is TRUE, given tmp3.TF[5]
        ned2l.dpobs.mix.shape.i <-
          rowSums(-pstr.mix * d1A.i * d0B.pi.mix / Delta.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.shape.i, M)] <-  # link done later
          ned2l.dpobs.mix.shape.i  # * dshape.i.deta
      }

      if (tmp3.TF[4]) {
        ned2l.dpstr.mix.shape.i <- rowSums(  # ccc.
          d1A.i * (pstr.mix * (d0A.i - d0B.pi.mix) / Delta.mix - 1))
        wz[, iam(posn.pstr.mix, posn.shape.i, M)] <-  # link done later
          ned2l.dpstr.mix.shape.i  # * dshape.i.deta
      }

      if (tmp3.TF[6]) {
        ned2l.dpobs.mlm.shape.i <- rowSums(
          -pstr.mix * d0B.pi.mix * d1A.i / Delta.mix)  # ccc.
        for (uuu in seq(lalt.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.shape.i, M)] <-
            ned2l.dpobs.mlm.shape.i  # * dshape.i.deta done later
      }
    }  # (tmp3.TF[5] && linf.mix > 1)



        
    if (tmp3.TF[7] && linf.mlm > 0) {  # \calI_{np}, includes \phi_s

      if (lalt.mix && tmp3.TF[2])
        ned2l.dpobs.mix.shape.p <-  # ccc.
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm))

      ned2l.dpstr.mix.shape.p <-  # ccc.
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.pi.mlm^2 / Delta.mlm)
      }
    }  # tmp3.TF[7] && linf.mlm > 0


    if (!is.na(posn.pobs.mix))  # Optional (1, 2) element:
      wz[, iam(1, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.shape.p  # One link done later
    if (!is.na(posn.pstr.mix))  # Optional (1, 4) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.shape.p  # One link done later

    if (!is.na(posn.pstr.mix))  # Optional (4, 4) element
      wz[, iam(posn.pstr.mix,
               posn.pstr.mix, M)] <- ned2l.dpstr.mix2  # Link done later






    if (tmp3.TF[6] && lalt.mlm) {  # \calA_{np}, includes \omega_s
      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (uuu in seq(lalt.mlm)) {  # Diagonal elts only
        wz[, iam(ofset + uuu,
                 ofset + uuu, M)] <- 1 / pobs.mlm[, uuu]
      }  # uuu

      tmp8a <- probns / Numer^2
      if (tmp3.TF[4] && linf.mix)
        tmp8a <- tmp8a + rowSums((d0B.pi.mix^2) / Delta.mix)
      if (tmp3.TF[7] && linf.mlm)
        tmp8a <- tmp8a + rowSums((d0B.pi.mlm^2) / Delta.mlm)
      for (uuu in seq(lalt.mlm))  # All elts
        for (vvv in uuu:lalt.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
          wz[, iam(ofset + uuu, ofset + vvv, M)] + tmp8a  # All elts
    }  # lalt.mlm


 

    if (tmp3.TF[6] && lalt.mlm) {

      init0.val <- if (tmp3.TF[7] && linf.mlm) rowSums(
        d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm)) else zero0n
      ned2l.dpobs.mlm.shape.p <- init0.val  # Vector, not matrix

      if (tmp3.TF[4] && linf.mix)
        ned2l.dpobs.mlm.shape.p <-
        ned2l.dpobs.mlm.shape.p + rowSums(
          d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (vvv in seq(lalt.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.shape.p
    }  # lalt.mlm > 0





    
    if (tmp3.TF[7] && linf.mlm > 0) {  # \calI_{np}, includes \phi_s
      init0.val <- probns / Numer^2
      if (linf.mix)
        init0.val <- init0.val + rowSums((d0B.pi.mix^2) / Delta.mix)
      ned2l.dpstr.mlm2 <-
        matrix(init0.val, n, linf.mlm * (linf.mlm + 1) / 2)

      for (uuu in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mlm2[, iam(uuu, uuu, linf.mlm)] <-
          ned2l.dpstr.mlm2[, iam(uuu, uuu, linf.mlm)] +
            ((sss == uuu) - d0B.pi.mlm[, sss])^2 / Delta.mlm[, sss]
      if (linf.mlm > 1) {
        for (uuu in 1:(linf.mlm-1))
          for (vvv in (uuu+1):linf.mlm)
            for (sss in seq(linf.mlm))
              ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] <-
              ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] +
              ((sss == uuu) - d0B.pi.mlm[, sss]) *
              ((sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss]
      }  # if (linf.mlm > 1)

      ofset <- posn.pstr.mlm - 1
      for (uuu in seq(linf.mlm))
        for (vvv in uuu:linf.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] 
    }  # linf.mlm > 0



    if (tmp3.TF[7] && linf.mlm > 0) {
      ned2l.dpstr.mlm.theta.p <- matrix(0, n, linf.mlm)
      for (vvv in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mlm.theta.p[, vvv] <-
          ned2l.dpstr.mlm.theta.p[, vvv] +
          d1B.pi.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss])

      if (linf.mix && tmp3.TF[4])
        ned2l.dpstr.mlm.theta.p <-
        ned2l.dpstr.mlm.theta.p + rowSums(
          d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(linf.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta.p[, vvv]
    }  # linf.mlm > 0




    if (linf.mlm && linf.mix > 1) {

      ned2l.dpstr.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.pi.mix * d1A.i / Delta.mix)

      for (vvv in seq(linf.mlm))
        wz[, iam(posn.shape.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.i  # ccc.
    }  # linf.mlm && linf.mix > 1



    if (all(c(lalt.mlm, linf.mlm) > 0)) {
      ned2l.dpobs.mlm.pstr.mlm2 <-
        array(probns / Numer^2, c(n, lalt.mlm, linf.mlm))
      for (uuu in seq(lalt.mlm))
        for (vvv in seq(linf.mlm))
          for (sss in seq(linf.mlm))
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] - d0B.pi.mlm[, sss] *
              ((sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss]

      if (tmp3.TF[4] && linf.mix)
        ned2l.dpobs.mlm.pstr.mlm2 <-
        ned2l.dpobs.mlm.pstr.mlm2 + rowSums(d0B.pi.mix^2 / Delta.mix)

      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pstr <- posn.pstr.mlm - 1
      for (uuu in seq(lalt.mlm))
        for (vvv in seq(linf.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pstr + vvv, M)] <-
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] 
    }  # all(c(lalt.mlm, linf.mlm) > 0)






    if (all(c(lalt.mix, lalt.mlm) > 0)) {
      ned2l.dpobs.mix.pobs.mlm <- probns / Numer^2  # Initialize
      if (linf.mix)  # tmp3.TF[4]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.pi.mix^2 / Delta.mix)
      if (linf.mlm)  # tmp3.TF[7]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.pi.mlm^2 / Delta.mlm)

      for (uuu in seq(lalt.mlm))  # ccc.
        wz[, iam(posn.pobs.mix, posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pobs.mlm  # Link done later
    }



    if (all(c(lalt.mix, linf.mlm) > 0)) {  # all(tmp3.TF[c(2, 7)])
      if (linf.mix)  # tmp3.TF[4]
        ned2l.dpobs.mix.pstr.mlm <-
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.pi.mix^2 / Delta.mix)

      for (uuu in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpobs.mix.pstr.mlm[, uuu] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.pi.mlm[, sss]) *
                            d0B.pi.mlm[, sss] / Delta.mlm[, sss]

      for (uuu in seq(linf.mlm))  # ccc.
        wz[, iam(posn.pobs.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu]  # Link done later
    }



    if (all(c(linf.mix, lalt.mlm) > 0)) {  # all(tmp3.TF[c(4, 6)])
      if (linf.mlm)  # tmp3.TF[7]
        ned2l.dpobs.mlm.pstr.mix <-
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.pi.mlm^2 / Delta.mlm)
        ned2l.dpobs.mlm.pstr.mix <-  # tmp3.TF[4] && linf.mix
        ned2l.dpobs.mlm.pstr.mix -
        rowSums((d0A.i - d0B.pi.mix) * d0B.pi.mix / Delta.mix)

      for (uuu in seq(lalt.mlm))  # ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pstr.mix  # Link done later
    }



    if (all(c(linf.mix, linf.mlm) > 0)) {  # all(tmp3.TF[c(4, 7)])

      for (uuu in seq(linf.mlm))  # tmp3.TF[4]
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mix.pstr.mlm[, uuu] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.pi.mlm[, sss]) *
              d0B.pi.mlm[, sss] / Delta.mlm[, sss]
      ned2l.dpstr.mix.pstr.mlm <-
      ned2l.dpstr.mix.pstr.mlm -
      rowSums((d0A.i - d0B.pi.mix) * d0B.pi.mix / Delta.mix)

      for (uuu in seq(linf.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu]  # Link done later
    }

 



 
    if (lall.len) {
      wz.4 <- matrix(0, n, M * (M + 1) / 2)  # Or == 0 * wz
      ind.rc <- setdiff(1:M, ind.shape.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the MLM

 

      for (uuu in ind.shape.z)
        for (sss in seq(M))
          wz.4[, iam(uuu, sss, M)] <- wz[, iam(uuu, sss, M)]





      speed.up <- intercept.only && (
                  length(offset) == 1 || all(offset[1] == offset))

      ind.mlm <- iam(NA, NA, lind.rc, both = TRUE, diag = TRUE)
      n.use <- 2
      bread <- if (speed.up)
               -allprobs[1:n.use, ind.mlm$row, drop = FALSE] *
                allprobs[1:n.use, ind.mlm$col, drop = FALSE] else
               -allprobs[, ind.mlm$row, drop = FALSE] *
                allprobs[, ind.mlm$col, drop = FALSE]
      if (speed.up) {
        bread[, 1:lind.rc] <-
        bread[, 1:lind.rc] + allprobs[1:n.use, 1:lind.rc]
      } else {
        bread[, 1:lind.rc] <-
        bread[, 1:lind.rc] + allprobs[, 1:lind.rc]  # -use.refLevel
      }
      bread <- m2a(bread, M = lind.rc)  # Half wasteful really


      if (!length(extra$ind.wz.match)) {
        imat <- matrix(NA, lind.rc, lind.rc)
        for (jay in seq(lind.rc)) {
          iptr <- jay
          for (kay in (ind.rc[jay]):M) {
            if (!any(kay %in% ind.shape.z)) {
              imat[jay, iptr] <-
                which(extra$index.M$row == ind.rc[jay] &
                      extra$index.M$col == kay)
              iptr <- iptr + 1
            }  # if
          }  # kay
        }  # jay
        ind.wz.match <- imat[cbind(ind.mlm$row.ind, ind.mlm$col.ind)]
        extra$ind.wz.match <- ind.wz.match  # Assign it once
      }  # !length(extra$ind.wz.match)
      filling <- if (speed.up)
        wz[1:n.use, extra$ind.wz.match, drop = FALSE] else
        wz[, extra$ind.wz.match, drop = FALSE]

      wz.5 <- mux5(filling, bread, M = lind.rc, matrix.arg = TRUE)

      wz.4[, extra$ind.wz.match] <- if (speed.up)
        matrix(wz.5[1, ], n, ncol(wz.5), byrow = TRUE) else c(wz.5)



 


      dstar.deta <- cbind(dshape.p.deta,
                          if (tmp3.TF[3]) dshape.a.deta else NULL,
                          if (tmp3.TF[5]) dshape.i.deta else NULL)
      iptr <- 0
      if (length(ind.shape.z))
      for (uuu in ind.shape.z) {  # Could delete 3 for shape.a (orthog)
        iptr <- iptr + 1
        for (ttt in seq(lind.rc)) {
          wz.4[, iam(uuu, ind.rc[ttt], M)] <- 0  # Initialize
          for (sss in seq(lind.rc)) {
            wz.4[, iam(uuu, ind.rc[ttt], M)] <-
            wz.4[, iam(uuu, ind.rc[ttt], M)] +
              allprobs[, sss] * (max(0, sss == ttt) - allprobs[, ttt]) *
              wz[, iam(uuu, ind.rc[sss], M)] * dstar.deta[, iptr]
          }  # sss
        }  # ttt
      }  # uuu


      wz <- wz.4  # Completed
    }  # lall.len



    mytiny <- (allprobs <       sqrt(.Machine$double.eps)) |
              (allprobs > 1.0 - sqrt(.Machine$double.eps))
    atiny <- rowSums(mytiny) > 0
    if (any(atiny)) {
      ind.diags <- setdiff(1:M, ind.shape.z)  # Exclude thetas
      wz[atiny, ind.diags] <- .Machine$double.eps +
      wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
    }



    c(w) * wz
  }), list( .truncate = truncate ))))
}  # gaitzeta












 gaitlog <-
  function(alt.mix = NULL, inf.mix = NULL,  # Unstructured probs are
           alt.mlm = NULL, inf.mlm = NULL,  # contiguous
           truncate = NULL, max.support = Inf,
           zero = c("pobs", "pstr"),  # Pruned later, handles all four
           eq.ap = FALSE,  # TRUE applies to the intercept, g_a(theta_a)
           eq.ip = FALSE,  # TRUE applies to the intercept, g_i(theta_i)
           parallel.ap = FALSE,  # TRUE applies to the intercept
           parallel.ip = FALSE,  # TRUE applies to the intercept
           lshape.p = "logitlink",
           lshape.a = "logitlink",
           lshape.i = "logitlink",
           type.fitted = c("mean", "shapes",
                           "pobs.mlm", "pstr.mlm",
                           "pobs.mix", "pstr.mix",
                           "Pobs.mix", "Pstr.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "ptrunc.p", "cdf.max.s"),
           gshape.p = -expm1(-7 * ppoints(12)),
           gpstr.mix = ppoints(9) / 2,
           gpstr.mlm = ppoints(9) / (2 + length(inf.mlm)),
           imethod = 1,
           imux = 0.5,  # General downward multiplier for init values
           ishape.p = NULL, ishape.a = ishape.p,
           ishape.i = ishape.p,
           ipobs.mix = NULL, ipstr.mix = NULL,  # 0.25, 
           ipobs.mlm = NULL, ipstr.mlm = NULL,  # 0.25, 
           byrow.ai = FALSE,
           ishrinkage = 0.95,
           probs.y = 0.35) {
  lowsup <- 1
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support, min.support = lowsup)
  lalt.mix <- length(alt.mix)
  linf.mix <- length(inf.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mlm <- length(inf.mlm)
  ltruncat <- length(truncate)
  ltrunc.use <- ltruncat > 0 || !is.infinite(max.support) 

  lshape.p <- as.list(substitute(lshape.p))
  eshape.p <- link2list(lshape.p)
  lshape.p <- attr(eshape.p, "function.name")
  lshape.p.save <- lshape.p
  gshape.p.save <- gshape.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lshape.a <- as.list(substitute(lshape.a))
  eshape.a <- link2list(lshape.a)
  lshape.a <- attr(eshape.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lshape.i <- as.list(substitute(lshape.i))
  eshape.i <- link2list(lshape.i)
  lshape.i <- attr(eshape.i, "function.name")


  if (is.vector(zero) && is.character(zero) && length(zero) == 2) {
    if (linf.mix + linf.mlm == 0)
      zero <- setdiff(zero, "pstr")
    if (lalt.mix + lalt.mlm == 0)
      zero <- setdiff(zero, "pobs")
  }


  lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
  if (lall.len + ltruncat == 0 && is.infinite(max.support))
    return(eval(substitute(
           logff(lshape = .lshape.p.save ,
                 gshape = .gshape.p.save ,
                 zero = NULL),
            list( .lshape.p.save = lshape.p.save,
                  .gshape.p.save = gshape.p.save  ))))

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")
  if (!is.logical(parallel.ap) || length(parallel.ap) != 1)
    stop("argument 'parallel.ap' must be a single logical")
  if (!is.logical(parallel.ip) || length(parallel.ip) != 1)
    stop("argument 'parallel.ip' must be a single logical")


  if (lalt.mix == 1 && eq.ap)
    warning("Only one unstructured altered value (no 'shape.a')",
            ", so setting 'eq.ap = TRUE' is meaningless")
  if (linf.mix == 1 && eq.ip)
    warning("Only one unstructured inflated value (no 'shape.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (lalt.mlm == 1 && parallel.ap)  # Only \omega_1
    warning("Only one altered mixture probability, 'pobs", alt.mlm,
            "', so setting 'parallel.ap = TRUE' is meaningless")
  if (linf.mlm == 1 && parallel.ip)  # Only \phi_1
    warning("Only one inflated mixture probability, 'pstr", inf.mlm,
            "', so setting 'parallel.ip = TRUE' is meaningless")


  type.fitted.choices <-
            c("mean", "shapes",
              "pobs.mlm", "pstr.mlm",
              "pobs.mix", "pstr.mix",
              "Pobs.mix", "Pstr.mix",
              "nonspecial", "Numer", "Denom.p",
              "sum.mlm.i", "ptrunc.p", "cdf.max.s")
  type.fitted <- match.arg(type.fitted[1], type.fitted.choices)[1]

  tmp7a <- if (lalt.mlm) paste0("pobs.mlm", alt.mlm) else NULL
  tmp7b <- if (linf.mlm) paste0("pstr.mlm", inf.mlm) else NULL
  tmp3 <- c(shape.p = lshape.p,
            pobs.mix = if (lalt.mix) "multilogitlink" else NULL,
            shape.a = if (lalt.mix > 1) lshape.a else NULL,
            pstr.mix = if (linf.mix) "multilogitlink" else NULL,
            shape.i = if (linf.mix > 1) lshape.i else NULL,
            if (lalt.mlm) rep("multilogitlink", lalt.mlm) else NULL,
            if (linf.mlm) rep("multilogitlink", linf.mlm) else NULL)
  Ltmp3 <- length(tmp3) 
  if (lalt.mlm + linf.mlm)
    names(tmp3)[(Ltmp3 - lalt.mlm - linf.mlm + 1):Ltmp3] <-
      c(tmp7a, tmp7b)
  par1or2 <- 1  # 2
  tmp3.TF <- c(TRUE, lalt.mix > 0, lalt.mix > 1,
                     linf.mix > 0, linf.mix > 1,
                     lalt.mlm > 0, linf.mlm > 0)
  indeta.finish <- cumsum(c(par1or2, 1, par1or2,
                                     1, par1or2,
                            lalt.mlm, linf.mlm,
                            linf.mlm + 1) * c(tmp3.TF, 1)) 
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  rownames(indeta) <- c("shape.p", "pobs.mix", "shape.a",
                        "pstr.mix", "shape.i", "pobs.mlm", "pstr.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.
      

  blurb1 <- "L"
  if (lalt.mlm + lalt.mix) blurb1 <- "Generally-altered L"
  if (linf.mlm + linf.mix) blurb1 <- "Generally-inflated L"
  if (ltrunc.use) blurb1 <- "Generally-truncated L"
  if ( (lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered and -inflated L"
  if ( (lalt.mlm + lalt.mix) && !(linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered and -truncated L"
  if (!(lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated and -truncated L"
  if ( (lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated and -truncated L"

      
  new("vglmff",
  blurb = c(blurb1, "ogarithmic regression\n",
            "(GAIT-Log(shape.p)-Log(shape.a)-MLM-",
                  "Log(shape.i)-MLM generally)\n\n",
            "Links: ",
            namesof("shape.p", lshape.p, earg = eshape.p, tag = FALSE),
            if (lalt.mix > 0) c(", ", "multilogit(pobs.mix)"),
            if (lalt.mix > 1) c(", ",
            namesof("shape.a",  lshape.a, eshape.a, tag = FALSE)),
            if (lalt.mix && linf.mix) ", \n       ",
            if (linf.mix > 0) c(  if (lalt.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (linf.mix > 1) c(", ",
            namesof("shape.i",  lshape.i, eshape.i, tag = FALSE)),
            if (lalt.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7a, collapse = "),\n"),
            ")") else NULL,
            if (linf.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7b, collapse = "),\n"),
            ")") else NULL),
  constraints = eval(substitute(expression({
    M1 <- max(extra$indeta, na.rm = TRUE)
    lalt.mix <- ( .lalt.mix )
    linf.mix <- ( .linf.mix )
    lalt.mlm <- ( .lalt.mlm )
    linf.mlm <- ( .linf.mlm )

    use.mat.mlm.a <- if (lalt.mlm) {
      if ( .parallel.ap ) matrix(1, lalt.mlm, 1) else diag(lalt.mlm)
    } else {
      NULL
    }

    use.mat.mlm.i <- if (linf.mlm) {
       if ( .parallel.ip ) matrix(1, linf.mlm, 1) else diag(linf.mlm)
    } else {
      NULL
    }

    if (lalt.mlm + linf.mlm == 0) {
      use.mat <- use.mat.mlm <- cbind(M)  # shape.p only
    }
    if (lalt.mlm + linf.mlm) {
      nc1 <- if (length(use.mat.mlm.a)) ncol(use.mat.mlm.a) else 0
      nc2 <- if (length(use.mat.mlm.i)) ncol(use.mat.mlm.i) else 0
      use.mat.mlm <- cbind(1, matrix(0, 1, nc1 + nc2))
      if (lalt.mlm)
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, lalt.mlm, 1),
                                   use.mat.mlm.a,
                                   if (length(use.mat.mlm.i) == 0) NULL
                                   else matrix(0, lalt.mlm, nc2)))
      if (linf.mlm )
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, linf.mlm, 1 + nc1),
                                   use.mat.mlm.i))
      use.mat <- use.mat.mlm
    }  # lalt.mlm + linf.mlm






    use.mat.mix <- if ( ( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,1,0,1,  0,1,0,0,0,  0,0,0,1,0), 5, 3) else
    if (!( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,0,0,1,  0,1,0,0,0,  0,0,1,0,0,
               0,0,0,1,0), 5, 4) else
    if ( ( .eq.ap ) && !( .eq.ip ))
      matrix(c(1,0,1,0,0,  0,1,0,0,0,  0,0,0,1,0,
               0,0,0,0,1), 5, 4) else
      diag(5)

    tmp3.TF <- ( .tmp3.TF )
    tmp3.TF <- tmp3.TF[-(6:7)]  # zz for now. temporary only.
    use.mat.mix <- use.mat.mix[tmp3.TF, , drop = FALSE]

    mincl <- apply(use.mat.mix, 2, min)
    maxcl <- apply(use.mat.mix, 2, max)
    use.mat.mix <- use.mat.mix[, mincl != 0 | maxcl != 0, drop = FALSE]

    if (lalt.mix + linf.mix > 0)
      use.mat <- use.mat.mix





    if (lalt.mlm + linf.mlm > 0 &&
        lalt.mix + linf.mix > 0) {
      use.mat <- rbind(use.mat.mix,
                 matrix(0, nrow(use.mat.mlm) - 1, ncol(use.mat.mix)))
      use.mat <- cbind(use.mat, matrix(0, nrow(use.mat),
                                          ncol(use.mat.mlm) - 1))
      use.mat[row(use.mat) > nrow(use.mat.mix) &
              col(use.mat) > ncol(use.mat.mix)] <- use.mat.mlm[-1, -1]
    }  # lalt.mlm + linf.mlm > 0 && lalt.mix + linf.mix > 0



    if (is.null(constraints)) {
      constraints <-
        cm.VGAM(use.mat, x = x, apply.int = TRUE,  # FALSE
                bool = .eq.ap || .eq.ip || .parallel.ap || .parallel.ip ,
                constraints = constraints)  # FALSE
    }



    if (lalt.mix + linf.mix + lalt.mlm + linf.mlm)
      constraints <-
        cm.zero.VGAM(constraints, x = x, .zero , M = M, M1 = M1,
                     predictors.names = paste0(predictors.names,
                                        names(predictors.names)))
  }), list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF,
            .eq.ap = eq.ap, .eq.ip = eq.ip,
            .parallel.ap = parallel.ap, .parallel.ip = parallel.ip,
            .lalt.mlm = lalt.mlm, .linf.mlm = linf.mlm,
            .lalt.mix = lalt.mix, .linf.mix = linf.mix
          ))),
  infos = eval(substitute(function(...) {
    list(M1 = .M1 ,
         Q1 = 1,
         link = .predictors.names ,  # .tmp3, as.vector strips names off
         link1parameter = as.logical( .lall.len <= 2),  # <= 1 safer
         mixture.links  = any(c( .lalt.mlm , .linf.mlm , .lalt.mix ,
                                 .linf.mix ) > 1),  # FALSE if NULL
         alt.mix = as.vector( .alt.mix ),  # Handles NULL
         alt.mlm = as.vector( .alt.mlm ),
         inf.mix = as.vector( .inf.mix ),
         inf.mlm = as.vector( .inf.mlm ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ),
         Support  = c( .lowsup , Inf, 1),  # a(b)c format as a,c,b.
         expected = TRUE,
         multipleResponses = FALSE,  # logff might be called if TRUE
         parameters.names = names( .predictors.names ),
         parent.name = c("logff", "log"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
         baseparams.argnames  = "shape",
         MM1 = 1,  # One parameter for 1 response
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .lshape.p = lshape.p, .eshape.p = eshape.p,
           .lshape.a = lshape.a, .eshape.a = eshape.a,
           .lshape.i = lshape.i, .eshape.i = eshape.i,
           .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
           .alt.mix = alt.mix, .inf.mix = inf.mix,
           .lalt.mlm = lalt.mlm, .linf.mlm = linf.mlm,
           .lalt.mix = lalt.mix, .linf.mix = linf.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len
         ))),
  initialize = eval(substitute(expression({
    extra$indeta <- ( .indeta )  # Avoids recomputing it several times
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
    truncate <- as.vector( .truncate )
    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(y)
    M <- NOS * M1

    tmp3.TF <- ( .tmp3.TF )

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Since max.support = 9 is possible
              ncol.y.max = 1,
              out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gaitcombo.check(y, truncate = truncate,
                               alt.mlm = alt.mlm, alt.mix = alt.mix,
                               inf.mlm = inf.mlm, inf.mix = inf.mix,
                               max.support = .max.support ,
                               min.support = .min.support )
    extra$skip.mix.a <- glist$skip.mix.a
    extra$skip.mix.i <- glist$skip.mix.i
    extra$skip.mlm.a <- glist$skip.mlm.a
    extra$skip.mlm.i <- glist$skip.mlm.i


    
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- as.vector( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1
    extra$index.M <- iam(NA, NA, M, both = TRUE)  # Used in @weight


    predictors.names <- ( .predictors.names )  # Got it, named


    if (!length(etastart)) {
      shape.p.init <- if (length( .ishape.p )) .ishape.p else {
        logff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
          sum(c(w) * dlog(x = y, shape = shapeval, log = TRUE))
        }
        shape.p.grid <- ( .gshape.p )
        grid.search(shape.p.grid, objfun = logff.Loglikfun, y = y, w = w)
      }
      shape.p.init <- rep(shape.p.init, length = n)

      shape.a.init <-  shape.i.init <- shape.p.init  # Needed
      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(shape.p.init, .lshape.p , earg = .eshape.p ))


      pobs.mix.init <- numeric(n)
      if (tmp3.TF[2]) {  # lalt.mix > 0
        pobs.mix.init <- if (length( .ipobs.mix )) {
        rep_len( .ipobs.mix , n)
      } else {
        is.alt.mix <- rowSums(glist$y0.mix.a) > 0
        rep(sum(w[is.alt.mix]) / sum(w), n)
      }  
      }  # lalt.mix > 0


      if (tmp3.TF[3]) {  # Assign coln 3; lalt.mix > 1
        shape.a.init <- if (length( .ishape.a ))
          rep_len( .ishape.a , n) else shape.p.init  # A vector
        etastart[, 3] <-
          theta2eta(shape.a.init, .lshape.a , earg = .eshape.a )
      }

      pstr.mix.init <- numeric(n)
      try.gridsearch.pstr.mix <- FALSE
      if (tmp3.TF[4]) {  # linf.mix > 0
        pstr.mix.init <- if (length( .ipstr.mix )) {
          rep_len( .ipstr.mix , n)
        } else {
          try.gridsearch.pstr.mix <- TRUE
          numeric(n)  # Overwritten by gridsearch
        }
      }  # linf.mix > 0


      if (tmp3.TF[5]) {  # linf.mix > 1
        shape.i.init <- if (length( .ishape.i ))
          rep_len( .ishape.i , n) else shape.p.init  # A vector
        etastart[, (extra$indeta[5, 'launch'])] <-
          theta2eta(shape.i.init, .lshape.i , earg = .eshape.i )
      }  # linf.mix > 1


      if (tmp3.TF[6]) {  #  lalt.mlm
        if (length( .ipobs.mlm )) {
          pobs.mlm.init <- matrix( .ipobs.mlm , n, lalt.mlm,
                                  byrow = .byrow.ai )
        } else {
          pobs.mlm.init <- colSums(c(w) * extra$skip.mlm.a) / colSums(w)
          pobs.mlm.init <- pobs.mlm.init * as.vector( .imux )
          pobs.mlm.init <- matrix(pobs.mlm.init, n, lalt.mlm, byrow=TRUE)
        }
      } else {
        pobs.mlm.init <- matrix(0, n, 1)
      }


      try.gridsearch.pstr.mlm <- FALSE
      if (tmp3.TF[7]) {  #  linf.mlm
        try.gridsearch.pstr.mlm <- !(length( .ipstr.mlm ))
        pstr.mlm.init <- 0  # Might be overwritten by gridsearch

        if (length( .ipstr.mlm ))
          pstr.mlm.init <- as.vector( .ipstr.mlm )

        pstr.mlm.init <- matrix(pstr.mlm.init, n, linf.mlm,
                                byrow = .byrow.ai )
      } else {
        pstr.mlm.init <- matrix(0, n, 1)
      }






      gaitlog.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitlog(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  shape.p    = extraargs$shape.p,
                  shape.a    = extraargs$shape.a,
                  shape.i    = extraargs$shape.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }

 gaitlog.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitlog(y, pstr.mlm = pstr.mlm.val,
                  pstr.mix    = extraargs$pstr.mix,  # Differs here
                  shape.p    = extraargs$shape.p,
                  shape.a    = extraargs$shape.a,
                  shape.i    = extraargs$shape.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }

 gaitlog.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitlog(y, pstr.mix = pstr.mix.val, pstr.mlm = pstr.mlm.val,
                  shape.p    = extraargs$shape.p,
                  shape.a    = extraargs$shape.a,
                  shape.i    = extraargs$shape.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }



      if (linf.mix + linf.mlm) {
        extraargs <- list(
              shape.p    = shape.p.init,
              shape.a    = shape.a.init,
              shape.i    = shape.i.init,
              alt.mix     = alt.mix,
              alt.mlm     = alt.mlm,
              inf.mix     = inf.mix,
              inf.mlm     = inf.mlm,
              truncate    = truncate,
              max.support = as.vector( .max.support ),
              pobs.mix    = pobs.mix.init ,
              pobs.mlm    = pobs.mlm.init )
        pre.warn <- options()$warn
        options(warn = -1)  # Ignore warnings during gridsearch 

        try.this <-
          if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
            grid.search2( .gpstr.mix ,  .gpstr.mlm ,
                         objfun = gaitlog.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- pstr.mlm.init
            grid.search ( .gpstr.mix ,
                         objfun = gaitlog.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- pstr.mix.init
            grid.search ( .gpstr.mlm ,
                         objfun = gaitlog.Loglikfun1.mlm,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          }

        options(warn = pre.warn)  # Restore warnings 
        if (any(is.na(try.this)))
          warning("gridsearch returned NAs. It's going to crash.",
                  immediate. = TRUE)
        if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
          pstr.mix.init <- rep_len(try.this["Value1"], n)
          pstr.mlm.init <- matrix(try.this["Value2"], n, linf.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mix=seq(5)/100'",
                 " and/or 'gpstr.mlm = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mix) {
          pstr.mix.init <- rep_len(try.this["Value"], n)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mix=seq(5)/100'.")
        } else if (try.gridsearch.pstr.mlm) {
          pstr.mlm.init <- matrix(try.this["Value"], n, linf.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mlm=seq(5)/100'.")
        }
      }  # lalt.mix + lnf.mix



        


      while (any((vecTF <- pobs.mix.init + pstr.mix.init +
                           rowSums(pobs.mlm.init) +
                           rowSums(pstr.mlm.init) > 0.96875))) {
        pobs.mix.init[vecTF]   <- 0.875 * pobs.mix.init[vecTF]
        pstr.mix.init[vecTF]   <- 0.875 * pstr.mix.init[vecTF]
        pobs.mlm.init[vecTF, ] <- 0.875 * pobs.mlm.init[vecTF, ]
        pstr.mlm.init[vecTF, ] <- 0.875 * pstr.mlm.init[vecTF, ]
      }

      Numer <- 1 - rowSums(pobs.mlm.init) - rowSums(pstr.mlm.init) -
               pobs.mix.init - pstr.mix.init
        
      etastart.z <- if (lall.len == 0) NULL else
        multilogitlink(cbind(if (tmp3.TF[2]) pobs.mix.init else NULL,
                             if (tmp3.TF[4]) pstr.mix.init else NULL,
                             if (tmp3.TF[6]) pobs.mlm.init else NULL,
                             if (tmp3.TF[7]) pstr.mlm.init else NULL,
                             Numer))
      if (!is.matrix(etastart.z)) etastart.z <- cbind(etastart.z)

      nextone <- 1  # Might not be used actually
      if (tmp3.TF[2]) {
        etastart[, 2] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[4]) {  # Coln 2 or 4
        etastart[, (extra$indeta[4, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[6]) {
        ind6 <- (extra$indeta[6, 'launch']):(extra$indeta[6, 'finish'])
        etastart[, ind6] <- etastart.z[, nextone:(nextone+lalt.mlm-1)]
        nextone <- nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind7 <- (extra$indeta[7, 'launch']):(extra$indeta[7, 'finish'])
        etastart[, ind7] <- etastart.z[, nextone:ncol(etastart.z)]
      }
    }
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .ishape.p = ishape.p,
    .ishape.a = ishape.a,
    .ishape.i = ishape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .ipstr.mix = ipstr.mix, .ipobs.mix = ipobs.mix,
    .ipstr.mlm = ipstr.mlm, .ipobs.mlm = ipobs.mlm,
    .byrow.ai = byrow.ai,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support,
    .min.support = lowsup,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .imux = imux,
    .gpstr.mix = gpstr.mix, .gshape.p = gshape.p,
    .gpstr.mlm = gpstr.mlm,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .indeta = indeta,
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"
      }
    type.fitted <-
      match.arg(type.fitted[1],
                c("mean", "shapes",
                  "pobs.mlm", "pstr.mlm",
                  "pobs.mix", "pstr.mix",
                          "Pobs.mix", "Pstr.mix",
                  "nonspecial", "Numer", "Denom.p", "sum.mlm.i",
                  "ptrunc.p", "cdf.max.s"))[1]

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <- shape.p  # Needed; and answer not corrupted

    if (any(tmp3.TF[c(3, 5)])) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      iptr <- 1
      shape.a <- if (!tmp3.TF[3]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[5]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0
    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")

      Nextone <- 0  # Might not be used actually
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitcombo.log(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              shape.a = shape.a, shape.i = shape.i,
              truncate = truncate, max.support = max.support)

    n.eta <- nrow(eta)
    Denom.p <- c(Bits[["cdf.max.s"]] - Bits[["SumT0.p"]] -
                 Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    Numer <- c(1 - pobs.mix - pstr.mix -
               (if (lalt.mlm) rowSums(pobs.mlm) else 0) -
               (if (linf.mlm) rowSums(pstr.mlm) else 0))

    if (!lalt.mlm && type.fitted %in% c("pobs.mlm")) {
      warning("No altered MLM values; returning an NA")
      return(NA)
    }
    if (!linf.mlm && type.fitted %in% c("sum.mlm.i", "pstr.mlm")) {
      warning("No inflated MLM values; returning an NA")
      return(NA)
    }
    if (!lalt.mix && type.fitted %in% c("Pobs.mix")) {
      warning("No altered mixture values; returning an NA")
      return(NA)
    }
    if (!linf.mix && type.fitted %in% c("sum.mix.i", "Pstr.mix")) {
      warning("No inflated mixture values; returning an NA")
      return(NA)
    }

    if (lalt.mix) {
      tmp13 <-  # dlog() does not retain the matrix format
        dlog(matrix(alt.mix, NROW(eta), lalt.mix, byrow = TRUE),
             matrix(shape.a, NROW(eta), lalt.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(NROW(eta), lalt.mix)
      dimnames(tmp13) <- list(rownames(eta), as.character(alt.mix))
      propn.mat.a <- tmp13
    }

    if (linf.mix) {
      tmp55 <-  # dlog() does not retain the matrix format
        dlog(matrix(inf.mix, NROW(eta), linf.mix, byrow = TRUE),
             matrix(shape.i, NROW(eta), linf.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(NROW(eta), linf.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(inf.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "shapes"     = cbind(shape.p,
                           if (tmp3.TF[3]) shape.a else NULL,
                           if (tmp3.TF[5]) shape.i else NULL),
      "pobs.mlm"   = pobs.mlm,  # aka omegamat, n x lalt.mlm
      "pstr.mlm"   = pstr.mlm,  # aka phimat, n x linf.mlm
      "pobs.mix"   = pobs.mix,  # n-vector
      "pstr.mix"   = pstr.mix,  # n-vector
      "Pobs.mix"   = c(pobs.mix) * propn.mat.a,  # matrix
      "Pstr.mix"   = c(pstr.mix) * propn.mat.i,
      "nonspecial" = Numer * (1 -
         (Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) / Denom.p),
      "Numer"      = Numer,
      "Denom.p"    = Denom.p,
      "sum.mlm.i"  = pstr.mlm + Numer *
             dlog(matrix(inf.mlm,  NROW(eta), linf.mlm, byrow = TRUE),
                  matrix(shape.p, NROW(eta), linf.mlm)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)

   ynames.pobs.mlm <- as.character(alt.mlm)  # Works with NULLs
   ynames.pstr.mlm <- as.character(inf.mlm)  # Works with NULLs
   if (length(ans))
     label.cols.y(ans, NOS = NOS, colnames.y =
     switch(type.fitted,
            "shapes"    = c("shape.p", "shape.a",  # Some colns NA
                            "shape.i")[(tmp3.TF[c(1, 3, 5)])],
            "Pobs.mix"  = as.character(alt.mix),
            "Pstr.mix"  = as.character(inf.mix),
            "pobs.mlm"  = ynames.pobs.mlm,
            "sum.mlm.i" = ,  #
            "pstr.mlm"  = ynames.pstr.mlm,
            extra$colnames.y)) else ans
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  last = eval(substitute(expression({
    pred.names <- c( .predictors.names )  # Save it
    link.names <- as.vector( .predictors.names )
    parameter.names <- names(pred.names)
    predictors.names <- NULL
    for (jay in seq(M))
      predictors.names <- c(predictors.names,
        namesof(parameter.names[jay], link.names[jay], tag = FALSE,
                earg = list()))  # This line isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .eshape.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[2])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[3])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.a )
    if (tmp3.TF[4])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[5])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.i )
    if (tmp3.TF[6]) {  # lalt.mlm
      for (ii in seq(lalt.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # lalt.mlm
    if (tmp3.TF[7]) {  # linf.mlm
      for (ii in seq(linf.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # linf.mlm
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .predictors.names = predictors.names,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <- shape.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5)])) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      iptr <- 1
      shape.a <- if (!tmp3.TF[3]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[5]) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0
    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitlog(y, shape.p, log = TRUE,  # byrow.ai = F,
                 alt.mix = alt.mix, inf.mix = inf.mix,
                 alt.mlm = alt.mlm, inf.mlm = inf.mlm,
                 truncate = truncate,
                 max.support = as.vector( .max.support ),
                 shape.a = shape.a, shape.i = shape.i, 
                 pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                 pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitlog"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    small. <- 1e-14
    pobs.mix <- pstr.mix <- small.
    pobs.mlm <- pstr.mlm <- matrix(small., NROW(eta), 1) # 4 rowSums()
    shape.a <- shape.i <- 0.5  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    shape.p <-
      cbind(eta2theta(eta[, 1], .lshape.p , earg = .eshape.p ))
    iptr <- 1
    ind.shape.z <- 1  # Points to shape.p only.
    if (lalt.mix + linf.mix > 1) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (lalt.mix <= 1) shape.p else {
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a ,
                  earg = .eshape.a )
      }
      shape.i <- if (linf.mix <= 1) shape.p else {
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i ,
                  earg = .eshape.i )
      }
    }  # lalt.mix + linf.mix > 0

    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm)
    okay.mix <-
      all(is.finite(shape.p)) && all(0 < shape.p) && all(shape.p < 1) &&
      all(is.finite(shape.a)) && all(0 < shape.a) && all(shape.a < 1) &&
      all(is.finite(shape.i)) && all(0 < shape.i) && all(shape.i < 1) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(pobs.mix + pstr.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1) 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    shape.a <- shape.i <- shape.p  # Needed
    ind.shape.z <- 1  # Points to shape.p only.

    if (lalt.mix + linf.mix > 1) {  # At least one shape.[ai]
      ind.shape.z <- object@extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      iptr <- 1
      shape.a <- if (lalt.mix <= 1) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.a , .eshape.a )
      shape.i <- if (linf.mix <= 1) shape.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0

    tmp3.TF <- ( .tmp3.TF )
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    rgaitlog(nsim * length(shape.p), shape.p,
             pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
             pobs.mix = pobs.mix, pstr.mix = pstr.mix,
             shape.a = shape.a, shape.i = shape.i, 
             alt.mlm = alt.mlm, inf.mlm = inf.mlm,
             alt.mix = alt.mix, inf.mix = inf.mix,
             truncate = .truncate , max.support = .max.support )
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .tmp3.TF = tmp3.TF,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({

    tmp3.TF <- ( .tmp3.TF )
    calA.p  <- tmp3.TF[2]
    calI.p  <- tmp3.TF[4]
    calA.np <- tmp3.TF[6]
    calI.np <- tmp3.TF[7]

    Denom1.a <- Denom1.i <- Denom2.i <- 0  # Denom2.a is not needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <- shape.p  # Needed; doesnt corrupt answer

    if (any(tmp3.TF[c(3, 5)])) {  # At least one shape.[ai]
      ind.shape.z <- extra$indeta[c(1, 3, 5), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      shape.a <- if (!tmp3.TF[3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
    }  # lalt.mix + linf.mix > 0

    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
      }
    }  # lall.len


    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")



    is.alt.mixed <- if (tmp3.TF[2])
      rowSums(extra$skip.mix.a) > 0 else rep(FALSE, n)
    is.inf.mixed <- if (tmp3.TF[4])
      rowSums(extra$skip.mix.i) > 0 else rep(FALSE, n)
    is.alt.mlmed <- if (tmp3.TF[6])
      rowSums(extra$skip.mlm.a) > 0 else rep(FALSE, n)
    is.inf.mlmed <- if (tmp3.TF[7])
      rowSums(extra$skip.mlm.i) > 0 else rep(FALSE, n)

    is.ns <- !is.alt.mlmed & !is.inf.mlmed  &
             !is.alt.mixed & !is.inf.mixed  # & !is.truncd

    A8.p <- -1 / log1p(-shape.p)
    A8.a <- -1 / log1p(-shape.a)
    A8.i <- -1 / log1p(-shape.i)
    prob.mlm.a <- if (lalt.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (linf.mlm) rowSums(pstr.mlm) else 0  # scalar okay
    pmf.deriv1 <- function(y, shape) {
      A8 <- -1 / log1p(-shape)
      deriv0 <- A8 * (shape^y) / y
      A8 * (shape^(y-1) - deriv0 / (1 - shape))
    }
    pmf.deriv2 <- function(y, shape) {
      A8 <- -1 / log1p(-shape)
      A8prime <- -(A8^2) / (1 - shape)
      deriv0 <- A8 * (shape^y) / y
      deriv1 <- A8 * (shape^(y-1) - deriv0 / (1 - shape))
      A8prime * (shape^(y-1) - deriv0 / (1 - shape)) +
      A8 * ((y - 1) * shape^(y - 2) - deriv0 / (1 - shape)^2 -
            deriv1 / (1 - shape)) 
    }


    sumD.mix.1a.p <- sumD.mix.2a.p <- matrix(0, n, NOS)
    if (lalt.mix > 0) {  # \calA_p
      Da.mix.0mat.a <-  # Matches naming convention further below
      Da.mix.1mat.a <- Da.mix.2mat.a <- matrix(0, n, lalt.mix)
      for (jay in seq(lalt.mix)) {
        aval <- alt.mix[jay]
        sumD.mix.1a.p <- sumD.mix.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mix.2a.p <- sumD.mix.2a.p + pmf.deriv2(aval, shape.p)
        pmf.a <- dlog(aval, shape.a)
        Da.mix.0mat.a[, jay] <- pmf.a
        Da.mix.1mat.a[, jay] <- pmf.deriv1(aval, shape.a)
      }
      Denom1.a <- rowSums(Da.mix.1mat.a)  # aka sumD.mix.1a.a
      Denom2.a <- rowSums(Da.mix.2mat.a)  # Needed; sumD.mix.2a.a <-
    }  # lalt.mix > 0




    if (linf.mix) {
      Di.mix.0mat.i <-  # wrt inflated distribution
      Di.mix.1mat.i <- Di.mix.2mat.i <- matrix(0, n, linf.mix)
      Dp.mix.0mat.i <-  # wrt parent distribution
      Dp.mix.1mat.i <- Dp.mix.2mat.i <- matrix(0, n, linf.mix)
        for (jay in seq(linf.mix)) {
          ival <- inf.mix[jay]
          pmf.i <- dlog(ival, shape.i)
          Di.mix.0mat.i[, jay] <- pmf.i
          Di.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.i)
          Di.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.i)
          pmf.p <- dlog(ival, shape.p)
          Dp.mix.0mat.i[, jay] <- pmf.p
          Dp.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.p)
          Dp.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.p)
        }  # jay
      Denom1.i <- rowSums(Di.mix.1mat.i)
      Denom2.i <- rowSums(Di.mix.2mat.i)
    }  # linf.mix



    bits <- moments.gaitcombo.log(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              shape.a = shape.a, shape.i = shape.i,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p <- sumD.mlm.2a.p <- matrix(0, n, NOS)
    if (lalt.mlm)
      for (aval in alt.mlm) {
        sumD.mlm.1a.p <- sumD.mlm.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mlm.2a.p <- sumD.mlm.2a.p + pmf.deriv2(aval, shape.p)
      }


    Denom0.p <- c(bits[["cdf.max.s"]] - bits[["SumT0.p"]] -
                  bits[["SumA0.mix.p"]] - bits[["SumA0.mlm.p"]])
    Numer <- 1 - pobs.mix - pstr.mix - prob.mlm.a - prob.mlm.i
    Denom0.a <- c(bits[["SumA0.mix.a"]])  # Not .p
    Denom0.i <- c(bits[["SumI0.mix.i"]])


 

    Dp.mlm.0Mat.i <-  # wrt parent distribution
    Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, NOS)
    if (linf.mlm > 0) {
      Dp.mlm.0Mat.i <-  # wrt parent distribution
      Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, linf.mlm)
      for (jay in seq(linf.mlm)) {
        ival <- inf.mlm[jay]
        pmf.p <- dlog(ival, shape.p)
        Dp.mlm.0Mat.i[, jay] <- pmf.p
        Dp.mlm.1Mat.i[, jay] <- pmf.deriv1(ival, shape.p)
        Dp.mlm.2Mat.i[, jay] <- pmf.deriv2(ival, shape.p)
      }  # jay
    }  # linf.mlm


    

    sumD.1t.a <- sumD.2t.a <-
    sumD.1t.i <- sumD.2t.i <-
    sumD.1t.p <- sumD.2t.p <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p <- sumD.1t.p + pmf.deriv1(tval, shape.p)
        sumD.2t.p <- sumD.2t.p + pmf.deriv2(tval, shape.p)
        sumD.1t.a <- sumD.1t.a + pmf.deriv1(tval, shape.a)
        sumD.2t.a <- sumD.2t.a + pmf.deriv2(tval, shape.a)
        sumD.1t.i <- sumD.1t.i + pmf.deriv1(tval, shape.i)
        sumD.2t.i <- sumD.2t.i + pmf.deriv2(tval, shape.i)
      }



      
    if (is.finite(max.support)) {
    tmp1.p <- A8.p * (shape.p^max.support -
      (1 - plog(max.support, shape.p))) / (1 - shape.p)
    sumD.1t.p <- sumD.1t.p + tmp1.p
    sumD.2t.p <- sumD.2t.p + (A8.p / (1 - shape.p)) * (
      (shape.p^max.support) / (1 - shape.p) +
      max.support * shape.p^(max.support - 1) -
      (1 - plog(max.support, shape.p)) / (1 - shape.p) - 2 * tmp1.p)


    tmp1.a <- A8.a * (shape.a^max.support -
      (1 - plog(max.support, shape.a))) / (1 - shape.a)
    sumD.1t.a <- sumD.1t.a + tmp1.a
    sumD.2t.a <- sumD.2t.a + (A8.a / (1 - shape.a)) * (
      (shape.a^max.support) / (1 - shape.a) +
      max.support * shape.a^(max.support - 1) -
      (1 - plog(max.support, shape.a)) / (1 - shape.a) - 2 * tmp1.a)
        
    tmp1.i <- A8.i * (shape.i^max.support -
      (1 - plog(max.support, shape.i))) / (1 - shape.i)
    sumD.1t.i <- sumD.1t.i + tmp1.i
    sumD.2t.i <- sumD.2t.i + (A8.i / (1 - shape.i)) * (
      (shape.i^max.support) / (1 - shape.i) +
      max.support * shape.i^(max.support - 1) -
      (1 - plog(max.support, shape.i)) / (1 - shape.i) - 2 * tmp1.i)
    }  # is.finite(max.support)






    Denom1.p <- c(-sumD.1t.p - sumD.mlm.1a.p - sumD.mix.1a.p)
    Denom2.p <- c(-sumD.2t.p - sumD.mlm.2a.p - sumD.mix.2a.p)


    d0B.pi.mlm <- Dp.mlm.0Mat.i / Denom0.p
    d1B.pi.mlm <- Dp.mlm.1Mat.i / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.i * Denom1.p / Denom0.p^2
    d2B.pi.mlm <-     Dp.mlm.2Mat.i / Denom0.p -
                  2 * Dp.mlm.1Mat.i * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.i * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.i * (Denom1.p^2) / Denom0.p^3


    Delta.mlm <- if (linf.mlm > 0) {
      pstr.mlm + Numer * d0B.pi.mlm  # n x linf.mlm.
    } else {
      matrix(0, n, 1)  # If linf.mlm == 0, for rowSums().
    }


    if (linf.mix > 0) {
      d0A.i <- Di.mix.0mat.i / Denom0.i
      d0B.pi.mix <- Dp.mix.0mat.i / Denom0.p
      Delta.mix <- pstr.mix * d0A.i + Numer * d0B.pi.mix


      d1A.i <- (Di.mix.1mat.i - Di.mix.0mat.i *
                Denom1.i / Denom0.i) / Denom0.i
      d2A.i <- (Di.mix.2mat.i - (2 * Di.mix.1mat.i * Denom1.i +
                Di.mix.0mat.i * Denom2.i) / Denom0.i +
            2 * Di.mix.0mat.i * (Denom1.i / Denom0.i)^2) / Denom0.i

      d1B.pi.mix <- Dp.mix.1mat.i / Denom0.p -
                    Dp.mix.0mat.i * Denom1.p / Denom0.p^2
      d2B.pi.mix <-     Dp.mix.2mat.i / Denom0.p -
                    2 * Dp.mix.1mat.i * Denom1.p / Denom0.p^2 -
                        Dp.mix.0mat.i * Denom2.p / Denom0.p^2 +
                    2 * Dp.mix.0mat.i * (Denom1.p^2) / Denom0.p^3
    }  # linf.mix > 0


    if (lalt.mix) {
      d0A.a <- Da.mix.0mat.a / Denom0.a
      d1A.a <- Da.mix.1mat.a / Denom0.a -
               Da.mix.0mat.a * Denom1.a / Denom0.a^2
      d2A.a <- (Da.mix.2mat.a - (2 * Da.mix.1mat.a * Denom1.a +
                Da.mix.0mat.a * Denom2.a) / Denom0.a +
            2 * Da.mix.0mat.a * (Denom1.a / Denom0.a)^2) / Denom0.a
    }  # lalt.mix




    dl.dshape.p <- -A8.p / (1 - shape.p) + y / shape.p
    dl.dshape.p[!is.ns] <- 0  # For is.alt.mixed & is.alt.mlmed
    dl.dshape.a <-
    dl.dshape.i <- numeric(n)  # Replace some elts below
    dl.dpstr.mix <- (-1) / Numer  # \notin A, I, T
    dl.dpobs.mix <- numeric(n)  # Replace some elts below
    dl.dpobs.mix[is.ns] <- (-1) / Numer[is.ns]
    dl.dpstr.mix[is.alt.mixed] <- 0
    dl.dpstr.mix[is.alt.mlmed] <- 0
    dl.dpobs.mlm <-
    dl.dpstr.mlm <- matrix(0, n, 1)  # Might not be needed



    if (tmp3.TF[6] && lalt.mlm) {  # aka \calA_{np}
      dl.dpobs.mlm <- matrix(-1 / Numer, n, lalt.mlm)  # \notin calS
      dl.dpobs.mlm[!is.ns, ] <- 0  # For alt.mix only really
      for (jay in seq(lalt.mlm)) {
        aval <- alt.mlm[jay]
        is.alt.j.mlm <- extra$skip.mlm.a[, jay]  # Logical vector
        tmp7a <- 1 / pobs.mlm[is.alt.j.mlm, jay]
        dl.dpobs.mlm[is.alt.j.mlm, jay] <- tmp7a
      }  # jay
    }  # lalt.mlm



    dl.dshape.p[is.ns] <- dl.dshape.p[is.ns] -
                           (Denom1.p / Denom0.p)[is.ns]


    
    if (tmp3.TF[7] && linf.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, linf.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For alt.mlm and alt.mix

      for (jay in seq(linf.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7i <- Numer * d1B.pi.mlm[, jay] / Delta.mlm[, jay]
        dl.dshape.p[is.inf.j.mlm] <- tmp7i[is.inf.j.mlm]


        if (tmp3.TF[6] && lalt.mlm) {  # 20200919
          tmp9i <- (-d0B.pi.mlm[, jay] / Delta.mlm[, jay])
          dl.dpobs.mlm[is.inf.j.mlm, ] <- tmp9i[is.inf.j.mlm]
        }


        if (tmp3.TF[2] && lalt.mix) {  # 20200919 linf.mix wrong
          tmp2 <- (-d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
          dl.dpobs.mix[is.inf.j.mlm] <- tmp2[is.inf.j.mlm]
        }


        if (tmp3.TF[4] && linf.mix) {  # 20200916
          tmp2 <- (-d0B.pi.mlm[, jay] / Delta.mlm[, jay])
          dl.dpstr.mix[is.inf.j.mlm] <- tmp2[is.inf.j.mlm]
        }

          
        tmp9 <- (0 - d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
        tmp8 <- (1 - d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
        dl.dpstr.mlm[is.inf.j.mlm, ] <- tmp9[is.inf.j.mlm]
        dl.dpstr.mlm[is.inf.j.mlm, jay] <- tmp8[is.inf.j.mlm]
      }  # jay
    }  # linf.mlm > 0


    



    if (tmp3.TF[2] && lalt.mix) {  # aka \calA_{p}
      dl.dpobs.mix[is.alt.mixed] <- 1 / pobs.mix[is.alt.mixed]

      if (tmp3.TF[3] && lalt.mix > 1)
        for (jay in seq(lalt.mix)) {
          is.alt.j.mix <- extra$skip.mix.a[, jay]  # Logical vector
          tmp2 <- d1A.a[, jay] / d0A.a[, jay]
          dl.dshape.a[is.alt.j.mix] <- tmp2[is.alt.j.mix]  # ccc.
        }  # jay
    }  # lalt.mix


    if (tmp3.TF[4] && linf.mix > 0) {  # aka \calI_{p}
      for (jay in seq(linf.mix)) {
        ival <- inf.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7b <- Numer * d1B.pi.mix[, jay] / Delta.mix[, jay]
        dl.dshape.p[is.inf.j.mix] <- tmp7b[is.inf.j.mix]
        if (tmp3.TF[2] && lalt.mix) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpobs.mix[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }
        tmp8 <- (d0A.i[, jay] - d0B.pi.mix[, jay]) / Delta.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]

        if (linf.mix > 1) {
          tmp2 <- pstr.mix * d1A.i[, jay] / Delta.mix[, jay]
          dl.dshape.i[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }

        if (tmp3.TF[6] && lalt.mlm) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpobs.mlm[is.inf.j.mix, ] <- tmp2[is.inf.j.mix]
        }
        if (tmp3.TF[7] && linf.mlm) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpstr.mlm[is.inf.j.mix, ] <- tmp2[is.inf.j.mix]
        }
      }  # jay
    }  # linf.mix > 0



    tmp3.TF <- !is.na(rowSums(extra$indeta))



    if (lall.len) {  # MLM fitted
      all4.dldp <- cbind(
        if (tmp3.TF[2]) dl.dpobs.mix else NULL,
        if (tmp3.TF[4]) dl.dpstr.mix else NULL,
        if (tmp3.TF[6]) dl.dpobs.mlm else NULL,
        if (tmp3.TF[7]) dl.dpstr.mlm else NULL)

      dl.deta.mlm4 <- matrix(0, n, ncol(all4.dldp))
      ind.quad.mlm <- seq(ncol(all4.dldp))  # Contiguous & starts at 1

      for (jay in ind.quad.mlm) {
        for (sss in ind.quad.mlm) {  # Dont need all4probs, allprobs ok
          dl.deta.mlm4[, jay] <- dl.deta.mlm4[, jay] +
            allprobs[, sss] * ((sss == jay) - allprobs[, jay]) *
            all4.dldp[, sss]
        }  # sss
      }  # jay
    }  # lall.len


    

    dshape.p.deta <- dtheta.deta(shape.p, .lshape.p , .eshape.p )
    if (tmp3.TF[3])
      dshape.a.deta <- dtheta.deta(shape.a, .lshape.a , .eshape.a )
    if (tmp3.TF[5])
      dshape.i.deta <- dtheta.deta(shape.i, .lshape.i , .eshape.i )


    iptr <- 0
    ansd <- cbind(shape.p = c(dl.dshape.p * dshape.p.deta))
    ansd <- cbind(ansd,
      pobs.mix = if (tmp3.TF[2]) dl.deta.mlm4[, (iptr <- iptr+1)],
      shape.a = if (tmp3.TF[3]) dl.dshape.a * dshape.a.deta)
    ansd <- cbind(ansd,
      pstr.mix = if (tmp3.TF[4]) dl.deta.mlm4[, (iptr <- iptr+1)],
      shape.i = if (tmp3.TF[5]) dl.dshape.i * dshape.i.deta)
    if (any(tmp3.TF[6:7])) {  # The remainder
      ansd <- cbind(ansd, dl.deta.mlm4[, (iptr+1):ncol(dl.deta.mlm4)])
      ind3 <- (ncol(ansd) - lalt.mlm - linf.mlm + 1):ncol(ansd)
      colnames(ansd)[ind3] <- as.character(c(alt.mlm, inf.mlm))
    }

    c(w) * ansd
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({

    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    mean.true.p <- A8.p * shape.p / (1 - shape.p)
    cond.EY.p <- c(mean.true.p -
                   bits[["SumI1.mlm.p"]] - bits[["SumI1.mix.p"]] -
                   bits[["SumT1.p"]] -
                   bits[["SumA1.mlm.p"]] - bits[["SumA1.mix.p"]]) / c(
                   Denom0.p - bits[["SumI0.mix.p"]] -
                              bits[["SumI0.mlm.p"]])
    probns <- Numer * (1 - c(bits[["SumI0.mix.p"]] +
                             bits[["SumI0.mlm.p"]]) / Denom0.p)






    zero0n <- numeric(n)
    ned2l.dpobs.mix.shape.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpobs.mix.shape.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.p <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.shape.i <- zero0n  # mB overwritten below


      posn.pobs.mix <- as.vector(extra$indeta[2, 'launch'])
      posn.shape.a <- as.vector(extra$indeta[3, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[4, 'launch'])
      posn.shape.i <- as.vector(extra$indeta[5, 'launch'])
      posn.pobs.mlm <- as.vector(extra$indeta[6, 'launch'])
      posn.pstr.mlm <- as.vector(extra$indeta[7, 'launch'])


    ned2l.dpstr.mix2         <-  # Elt (4, 4)
    ned2l.dpobs.mlm.pstr.mix <-  # Elts (4, >=6)
    ned2l.dpobs.mix.pstr.mix <- probns / Numer^2  # ccc Elt (2, 4)
    if (all(c(lalt.mix, linf.mlm) > 0))
    ned2l.dpobs.mix.pstr.mlm <- matrix(probns / Numer^2, n, linf.mlm)
    if (all(c(linf.mix, linf.mlm) > 0))
    ned2l.dpstr.mix.pstr.mlm <- matrix(probns / Numer^2, n, linf.mlm)
    



    ned2l.dshape.p2 <- probns * (cond.EY.p / shape.p^2 +  # ccc
                                 A8.p * (1 - A8.p) / (1 - shape.p)^2 +
      Denom2.p / Denom0.p - (Denom1.p / Denom0.p)^2) + 
      (if (tmp3.TF[4] && linf.mix) Numer *
      rowSums(Numer * (d1B.pi.mix^2) / Delta.mix - d2B.pi.mix) else 0) +
      (if (tmp3.TF[7] && linf.mlm) Numer *
      rowSums(Numer * (d1B.pi.mlm^2) / Delta.mlm - d2B.pi.mlm) else 0)

    wz[, iam(1, 1, M)] <- ned2l.dshape.p2 * dshape.p.deta^2


    ned2l.dpobs.mix2 <- 1 / pobs.mix + probns / Numer^2
    if (tmp3.TF[4] && linf.mix > 0) {
      ned2l.dpobs.mix2 <-  # More just below, ccc
      ned2l.dpobs.mix2 + rowSums(d0B.pi.mix^2 / Delta.mix)
    }
    if (tmp3.TF[7] && linf.mlm > 0) {
      ned2l.dpobs.mix2 <-  # ccc.
      ned2l.dpobs.mix2 + rowSums(d0B.pi.mlm^2 / Delta.mlm)
    }
    if (tmp3.TF[2] && lalt.mix > 0)
      wz[, iam(2, 2, M)] <- ned2l.dpobs.mix2  # Link done later


    if (tmp3.TF[3] && lalt.mix > 1) {
      ned2l.dshape.a2 <- pobs.mix * (
        rowSums((Da.mix.1mat.a^2) / Da.mix.0mat.a) / Denom0.a -
        (Denom1.a / Denom0.a)^2)  # ccc.
      wz[, iam(3, 3, M)] <- ned2l.dshape.a2 * dshape.a.deta^2
    }






    if (tmp3.TF[4] && linf.mix > 0) {

      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.pi.mix)^2 / Delta.mix)


      if (tmp3.TF[2] && lalt.mix > 0)
        ned2l.dpobs.mix.shape.p <-
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ned2l.dpstr.mix.shape.p <-
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.pi.mix * (1 + Numer * (d0A.i - d0B.pi.mix) / Delta.mix))


    if (all(tmp3.TF[c(2, 4)]))
      ned2l.dpobs.mix.pstr.mix <-  # ccc
      ned2l.dpobs.mix.pstr.mix +
        rowSums(-d0B.pi.mix * (d0A.i - d0B.pi.mix) / Delta.mix)

    }  # (tmp3.TF[4] && linf.mix > 0)



    
    if (all(tmp3.TF[c(2, 4, 7)])) {  # was  lalt.mix > 0 & Delta.mix
      ned2l.dpobs.mix.pstr.mix <-  # ccc.
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.pi.mlm^2 / Delta.mlm)
    }

    if (!is.na(posn.pobs.mix) && !is.na(posn.pstr.mix))
      wz[, iam(posn.pobs.mix, posn.pstr.mix, M)] <-
        ned2l.dpobs.mix.pstr.mix  # Link done later




    if (tmp3.TF[5] && linf.mix > 1) {  # \calI_{p}, includes \theta_iota

      ned2l.dshape.p.shape.i <- pstr.mix * Numer *
        rowSums(d1A.i * d1B.pi.mix / Delta.mix)  # ccc.
      wz[, iam(1, posn.shape.i, M)] <- ned2l.dshape.p.shape.i *
        dshape.p.deta * dshape.i.deta  # All links done here

      ned2l.dshape.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i^2) / Delta.mix - d2A.i)  # ccc.
      wz[, iam(posn.shape.i, posn.shape.i, M)] <-
        ned2l.dshape.i2 * dshape.i.deta^2

      if (tmp3.TF[2]) {  # tmp3.TF[4] is TRUE, given tmp3.TF[5]
        ned2l.dpobs.mix.shape.i <-
          rowSums(-pstr.mix * d1A.i * d0B.pi.mix / Delta.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.shape.i, M)] <-  # link done later
          ned2l.dpobs.mix.shape.i  # * dshape.i.deta
      }

      if (tmp3.TF[4]) {
        ned2l.dpstr.mix.shape.i <- rowSums(  # ccc.
          d1A.i * (pstr.mix * (d0A.i - d0B.pi.mix) / Delta.mix - 1))
        wz[, iam(posn.pstr.mix, posn.shape.i, M)] <-  # link done later
          ned2l.dpstr.mix.shape.i  # * dshape.i.deta
      }

      if (tmp3.TF[6]) {
        ned2l.dpobs.mlm.shape.i <- rowSums(
          -pstr.mix * d0B.pi.mix * d1A.i / Delta.mix)  # ccc.
        for (uuu in seq(lalt.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.shape.i, M)] <-
            ned2l.dpobs.mlm.shape.i  # * dshape.i.deta done later
      }
    }  # (tmp3.TF[5] && linf.mix > 1)



        
    if (tmp3.TF[7] && linf.mlm > 0) {  # \calI_{np}, includes \phi_s

      if (lalt.mix && tmp3.TF[2])
        ned2l.dpobs.mix.shape.p <-  # ccc.
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm))

      ned2l.dpstr.mix.shape.p <-  # ccc.
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.pi.mlm^2 / Delta.mlm)
      }
    }  # tmp3.TF[7] && linf.mlm > 0


    if (!is.na(posn.pobs.mix))  # Optional (1, 2) element:
      wz[, iam(1, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.shape.p  # One link done later
    if (!is.na(posn.pstr.mix))  # Optional (1, 4) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.shape.p  # One link done later

    if (!is.na(posn.pstr.mix))  # Optional (4, 4) element
      wz[, iam(posn.pstr.mix,
               posn.pstr.mix, M)] <- ned2l.dpstr.mix2  # Link done later






    if (tmp3.TF[6] && lalt.mlm) {  # \calA_{np}, includes \omega_s
      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (uuu in seq(lalt.mlm)) {  # Diagonal elts only
        wz[, iam(ofset + uuu,
                 ofset + uuu, M)] <- 1 / pobs.mlm[, uuu]
      }  # uuu

      tmp8a <- probns / Numer^2
      if (tmp3.TF[4] && linf.mix)
        tmp8a <- tmp8a + rowSums((d0B.pi.mix^2) / Delta.mix)
      if (tmp3.TF[7] && linf.mlm)
        tmp8a <- tmp8a + rowSums((d0B.pi.mlm^2) / Delta.mlm)
      for (uuu in seq(lalt.mlm))  # All elts
        for (vvv in uuu:lalt.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
          wz[, iam(ofset + uuu, ofset + vvv, M)] + tmp8a  # All elts
    }  # lalt.mlm


 

    if (tmp3.TF[6] && lalt.mlm) {

      init0.val <- if (tmp3.TF[7] && linf.mlm) rowSums(
        d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm)) else zero0n
      ned2l.dpobs.mlm.shape.p <- init0.val  # Vector, not matrix

      if (tmp3.TF[4] && linf.mix)
        ned2l.dpobs.mlm.shape.p <-
        ned2l.dpobs.mlm.shape.p + rowSums(
          d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (vvv in seq(lalt.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.shape.p
    }  # lalt.mlm > 0





    
    if (tmp3.TF[7] && linf.mlm > 0) {  # \calI_{np}, includes \phi_s
      init0.val <- probns / Numer^2
      if (linf.mix)
        init0.val <- init0.val + rowSums((d0B.pi.mix^2) / Delta.mix)
      ned2l.dpstr.mlm2 <-
        matrix(init0.val, n, linf.mlm * (linf.mlm + 1) / 2)

      for (uuu in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mlm2[, iam(uuu, uuu, linf.mlm)] <-
          ned2l.dpstr.mlm2[, iam(uuu, uuu, linf.mlm)] +
            ((sss == uuu) - d0B.pi.mlm[, sss])^2 / Delta.mlm[, sss]
      if (linf.mlm > 1) {
        for (uuu in 1:(linf.mlm-1))
          for (vvv in (uuu+1):linf.mlm)
            for (sss in seq(linf.mlm))
              ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] <-
              ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] +
              ((sss == uuu) - d0B.pi.mlm[, sss]) *
              ((sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss]
      }  # if (linf.mlm > 1)

      ofset <- posn.pstr.mlm - 1
      for (uuu in seq(linf.mlm))
        for (vvv in uuu:linf.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] 
    }  # linf.mlm > 0



    if (tmp3.TF[7] && linf.mlm > 0) {
      ned2l.dpstr.mlm.theta.p <- matrix(0, n, linf.mlm)
      for (vvv in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mlm.theta.p[, vvv] <-
          ned2l.dpstr.mlm.theta.p[, vvv] +
          d1B.pi.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss])

      if (linf.mix && tmp3.TF[4])
        ned2l.dpstr.mlm.theta.p <-
        ned2l.dpstr.mlm.theta.p + rowSums(
          d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(linf.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta.p[, vvv]
    }  # linf.mlm > 0




    if (linf.mlm && linf.mix > 1) {

      ned2l.dpstr.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.pi.mix * d1A.i / Delta.mix)

      for (vvv in seq(linf.mlm))
        wz[, iam(posn.shape.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.i  # ccc.
    }  # linf.mlm && linf.mix > 1



    if (all(c(lalt.mlm, linf.mlm) > 0)) {
      ned2l.dpobs.mlm.pstr.mlm2 <-
        array(probns / Numer^2, c(n, lalt.mlm, linf.mlm))
      for (uuu in seq(lalt.mlm))
        for (vvv in seq(linf.mlm))
          for (sss in seq(linf.mlm))
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] - d0B.pi.mlm[, sss] *
              ((sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss]

      if (tmp3.TF[4] && linf.mix)
        ned2l.dpobs.mlm.pstr.mlm2 <-
        ned2l.dpobs.mlm.pstr.mlm2 + rowSums(d0B.pi.mix^2 / Delta.mix)

      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pstr <- posn.pstr.mlm - 1
      for (uuu in seq(lalt.mlm))
        for (vvv in seq(linf.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pstr + vvv, M)] <-
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] 
    }  # all(c(lalt.mlm, linf.mlm) > 0)






    if (all(c(lalt.mix, lalt.mlm) > 0)) {
      ned2l.dpobs.mix.pobs.mlm <- probns / Numer^2  # Initialize
      if (linf.mix)  # tmp3.TF[4]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.pi.mix^2 / Delta.mix)
      if (linf.mlm)  # tmp3.TF[7]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.pi.mlm^2 / Delta.mlm)

      for (uuu in seq(lalt.mlm))  # ccc.
        wz[, iam(posn.pobs.mix, posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pobs.mlm  # Link done later
    }



    if (all(c(lalt.mix, linf.mlm) > 0)) {  # all(tmp3.TF[c(2, 7)])
      if (linf.mix)  # tmp3.TF[4]
        ned2l.dpobs.mix.pstr.mlm <-
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.pi.mix^2 / Delta.mix)

      for (uuu in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpobs.mix.pstr.mlm[, uuu] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.pi.mlm[, sss]) *
                            d0B.pi.mlm[, sss] / Delta.mlm[, sss]

      for (uuu in seq(linf.mlm))  # ccc.
        wz[, iam(posn.pobs.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu]  # Link done later
    }



    if (all(c(linf.mix, lalt.mlm) > 0)) {  # all(tmp3.TF[c(4, 6)])
      if (linf.mlm)  # tmp3.TF[7]
        ned2l.dpobs.mlm.pstr.mix <-
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.pi.mlm^2 / Delta.mlm)
        ned2l.dpobs.mlm.pstr.mix <-  # tmp3.TF[4] && linf.mix
        ned2l.dpobs.mlm.pstr.mix -
        rowSums((d0A.i - d0B.pi.mix) * d0B.pi.mix / Delta.mix)

      for (uuu in seq(lalt.mlm))  # ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pstr.mix  # Link done later
    }



    if (all(c(linf.mix, linf.mlm) > 0)) {  # all(tmp3.TF[c(4, 7)])

      for (uuu in seq(linf.mlm))  # tmp3.TF[4]
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mix.pstr.mlm[, uuu] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.pi.mlm[, sss]) *
              d0B.pi.mlm[, sss] / Delta.mlm[, sss]
      ned2l.dpstr.mix.pstr.mlm <-
      ned2l.dpstr.mix.pstr.mlm -
      rowSums((d0A.i - d0B.pi.mix) * d0B.pi.mix / Delta.mix)

      for (uuu in seq(linf.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu]  # Link done later
    }

 



 
    if (lall.len) {
      wz.4 <- matrix(0, n, M * (M + 1) / 2)  # Or == 0 * wz
      ind.rc <- setdiff(1:M, ind.shape.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the MLM

 

      for (uuu in ind.shape.z)
        for (sss in seq(M))
          wz.4[, iam(uuu, sss, M)] <- wz[, iam(uuu, sss, M)]





      speed.up <- intercept.only && (
                  length(offset) == 1 || all(offset[1] == offset))

      ind.mlm <- iam(NA, NA, lind.rc, both = TRUE, diag = TRUE)
      n.use <- 2
      bread <- if (speed.up)
               -allprobs[1:n.use, ind.mlm$row, drop = FALSE] *
                allprobs[1:n.use, ind.mlm$col, drop = FALSE] else
               -allprobs[, ind.mlm$row, drop = FALSE] *
                allprobs[, ind.mlm$col, drop = FALSE]
      if (speed.up) {
        bread[, 1:lind.rc] <-
        bread[, 1:lind.rc] + allprobs[1:n.use, 1:lind.rc]
      } else {
        bread[, 1:lind.rc] <-
        bread[, 1:lind.rc] + allprobs[, 1:lind.rc]  # -use.refLevel
      }
      bread <- m2a(bread, M = lind.rc)  # Half wasteful really

      if (!length(extra$ind.wz.match)) {
        imat <- matrix(NA, lind.rc, lind.rc)
        for (jay in seq(lind.rc)) {
          iptr <- jay
          for (kay in (ind.rc[jay]):M) {
            if (!any(kay %in% ind.shape.z)) {
              imat[jay, iptr] <-
                which(extra$index.M$row == ind.rc[jay] &
                      extra$index.M$col == kay)
              iptr <- iptr + 1
            }  # if
          }  # kay
        }  # jay
        ind.wz.match <- imat[cbind(ind.mlm$row.ind, ind.mlm$col.ind)]
        extra$ind.wz.match <- ind.wz.match  # Assign it once
      }  # !length(extra$ind.wz.match)
      filling <- if (speed.up)
        wz[1:n.use, extra$ind.wz.match, drop = FALSE] else
        wz[, extra$ind.wz.match, drop = FALSE]
      wz.5 <- mux5(filling, bread, M = lind.rc, matrix.arg = TRUE)
      wz.4[, extra$ind.wz.match] <- if (speed.up)
        matrix(wz.5[1, ], n, ncol(wz.5), byrow = TRUE) else c(wz.5)




 


      dstar.deta <- cbind(dshape.p.deta,
                          if (tmp3.TF[3]) dshape.a.deta else NULL,
                          if (tmp3.TF[5]) dshape.i.deta else NULL)
      iptr <- 0
      if (length(ind.shape.z))
      for (uuu in ind.shape.z) {  # Could delete 3 for shape.a (orthog)
        iptr <- iptr + 1
        for (ttt in seq(lind.rc)) {
          wz.4[, iam(uuu, ind.rc[ttt], M)] <- 0  # Initialize
          for (sss in seq(lind.rc)) {
            wz.4[, iam(uuu, ind.rc[ttt], M)] <-
            wz.4[, iam(uuu, ind.rc[ttt], M)] +
              allprobs[, sss] * (max(0, sss == ttt) - allprobs[, ttt]) *
              wz[, iam(uuu, ind.rc[sss], M)] * dstar.deta[, iptr]
          }  # sss
        }  # ttt
      }  # uuu


      wz <- wz.4  # Completed
    }  # lall.len



    mytiny <- (allprobs <       sqrt(.Machine$double.eps)) |
              (allprobs > 1.0 - sqrt(.Machine$double.eps))
    atiny <- rowSums(mytiny) > 0
    if (any(atiny)) {
      ind.diags <- setdiff(1:M, ind.shape.z)  # Exclude thetas
      wz[atiny, ind.diags] <- .Machine$double.eps +
      wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
    }



    c(w) * wz
  }), list( .truncate = truncate ))))
}  # gaitlog












 gaitpoisson <-
  function(alt.mix = NULL, inf.mix = NULL,  # Unstructured probs are
           alt.mlm = NULL, inf.mlm = NULL,  # contiguous
           truncate = NULL, max.support = Inf,
           zero = c("pobs", "pstr"),  # Pruned later, handles all four
           eq.ap = FALSE,  # TRUE applies to the intercept, g_a(theta_a)
           eq.ip = FALSE,  # TRUE applies to the intercept, g_i(theta_i)
           parallel.ap = FALSE,  # TRUE applies to the intercept
           parallel.ip = FALSE,  # TRUE applies to the intercept
           llambda.p = "loglink",
           llambda.a = "loglink",
           llambda.i = "loglink",
           type.fitted = c("mean", "lambdas",
                           "pobs.mlm", "pstr.mlm",
                           "pobs.mix", "pstr.mix",
                           "Pobs.mix", "Pstr.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "ptrunc.p", "cdf.max.s"),
           gpstr.mix = ppoints(9) / 2,
           gpstr.mlm = ppoints(9) / (2 + length(inf.mlm)),
           imethod = 1,
           imux = 0.5,  # General downward multiplier for init values
           ilambda.p = NULL, ilambda.a = ilambda.p,
           ilambda.i = ilambda.p,
           ipobs.mix = NULL, ipstr.mix = NULL,  # 0.25, 
           ipobs.mlm = NULL, ipstr.mlm = NULL,  # 0.25, 
           byrow.ai = FALSE,
           ishrinkage = 0.95,
           probs.y = 0.35) {
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support, min.support = lowsup)
  lalt.mix <- length(alt.mix)
  linf.mix <- length(inf.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mlm <- length(inf.mlm)
  ltruncat <- length(truncate)
  ltrunc.use <- ltruncat > 0 || !is.infinite(max.support) 

  llambda.p <- as.list(substitute(llambda.p))
  elambda.p <- link2list(llambda.p)
  llambda.p <- attr(elambda.p, "function.name")
  llambda.p.save <- llambda.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  llambda.a <- as.list(substitute(llambda.a))
  elambda.a <- link2list(llambda.a)
  llambda.a <- attr(elambda.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  llambda.i <- as.list(substitute(llambda.i))
  elambda.i <- link2list(llambda.i)
  llambda.i <- attr(elambda.i, "function.name")


  if (is.vector(zero) && is.character(zero) && length(zero) == 2) {
    if (linf.mix + linf.mlm == 0)
      zero <- setdiff(zero, "pstr")
    if (lalt.mix + lalt.mlm == 0)
      zero <- setdiff(zero, "pobs")
  }


  lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
  if (lall.len + ltruncat == 0 && is.infinite(max.support))
    return(eval(substitute(
           poissonff(link = .llambda.p.save , zero = NULL),
           list( .llambda.p.save = llambda.p.save))))

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")
  if (!is.logical(parallel.ap) || length(parallel.ap) != 1)
    stop("argument 'parallel.ap' must be a single logical")
  if (!is.logical(parallel.ip) || length(parallel.ip) != 1)
    stop("argument 'parallel.ip' must be a single logical")


  if (lalt.mix == 1 && eq.ap)
    warning("Only one unstructured altered value (no 'lambda.a')",
            ", so setting 'eq.ap = TRUE' is meaningless")
  if (linf.mix == 1 && eq.ip)
    warning("Only one unstructured inflated value (no 'lambda.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (lalt.mlm == 1 && parallel.ap)  # Only \omega_1
    warning("Only one altered mixture probability, 'pobs", alt.mlm,
            "', so setting 'parallel.ap = TRUE' is meaningless")
  if (linf.mlm == 1 && parallel.ip)  # Only \phi_1
    warning("Only one inflated mixture probability, 'pstr", inf.mlm,
            "', so setting 'parallel.ip = TRUE' is meaningless")


  type.fitted.choices <-
            c("mean", "lambdas",
              "pobs.mlm", "pstr.mlm",
              "pobs.mix", "pstr.mix",
              "Pobs.mix", "Pstr.mix",
              "nonspecial", "Numer", "Denom.p",
              "sum.mlm.i", "ptrunc.p", "cdf.max.s")
  type.fitted <- match.arg(type.fitted[1], type.fitted.choices)[1]

  tmp7a <- if (lalt.mlm) paste0("pobs.mlm", alt.mlm) else NULL
  tmp7b <- if (linf.mlm) paste0("pstr.mlm", inf.mlm) else NULL
  tmp3 <- c(lambda.p = llambda.p,
            pobs.mix = if (lalt.mix) "multilogitlink" else NULL,
            lambda.a = if (lalt.mix > 1) llambda.a else NULL,
            pstr.mix = if (linf.mix) "multilogitlink" else NULL,
            lambda.i = if (linf.mix > 1) llambda.i else NULL,
            if (lalt.mlm) rep("multilogitlink", lalt.mlm) else NULL,
            if (linf.mlm) rep("multilogitlink", linf.mlm) else NULL)
  Ltmp3 <- length(tmp3) 
  if (lalt.mlm + linf.mlm)
    names(tmp3)[(Ltmp3 - lalt.mlm - linf.mlm + 1):Ltmp3] <-
      c(tmp7a, tmp7b)
  par1or2 <- 1  # 2
  tmp3.TF <- c(TRUE, lalt.mix > 0, lalt.mix > 1,
                     linf.mix > 0, linf.mix > 1,
                     lalt.mlm > 0, linf.mlm > 0)
  indeta.finish <- cumsum(c(par1or2, 1, par1or2,
                                     1, par1or2,
                            lalt.mlm, linf.mlm,
                            linf.mlm + 1) * c(tmp3.TF, 1)) 
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  rownames(indeta) <- c("lambda.p", "pobs.mix", "lambda.a",
                        "pstr.mix", "lambda.i", "pobs.mlm", "pstr.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.
      

  blurb1 <- "P"
  if (lalt.mlm + lalt.mix) blurb1 <- "Generally-altered P"
  if (linf.mlm + linf.mix) blurb1 <- "Generally-inflated P"
  if (ltrunc.use) blurb1 <- "Generally-truncated P"
  if ( (lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered and -inflated P"
  if ( (lalt.mlm + lalt.mix) && !(linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered and -truncated P"
  if (!(lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated and -truncated P"
  if ( (lalt.mlm + lalt.mix) &&  (linf.mlm + linf.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated and -truncated P"

      
  new("vglmff",
  blurb = c(blurb1, "oisson regression\n",
            "(GAIT-Pois(lambda.p)-Pois(lambda.a)-MLM-",
                  "Pois(lambda.i)-MLM generally)\n\n",
            "Links: ",
            namesof("lambda.p", llambda.p, earg = elambda.p, tag = FALSE),
            if (lalt.mix > 0) c(", ", "multilogit(pobs.mix)"),
            if (lalt.mix > 1) c(", ",
            namesof("lambda.a",  llambda.a, elambda.a, tag = FALSE)),
            if (lalt.mix && linf.mix) ", \n       ",
            if (linf.mix > 0) c(  if (lalt.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (linf.mix > 1) c(", ",
            namesof("lambda.i",  llambda.i, elambda.i, tag = FALSE)),
            if (lalt.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7a, collapse = "),\n"),
            ")") else NULL,
            if (linf.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7b, collapse = "),\n"),
            ")") else NULL),
  constraints = eval(substitute(expression({
    M1 <- max(extra$indeta, na.rm = TRUE)
    lalt.mix <- ( .lalt.mix )
    linf.mix <- ( .linf.mix )
    lalt.mlm <- ( .lalt.mlm )
    linf.mlm <- ( .linf.mlm )

    use.mat.mlm.a <- if (lalt.mlm) {
      if ( .parallel.ap ) matrix(1, lalt.mlm, 1) else diag(lalt.mlm)
    } else {
      NULL
    }

    use.mat.mlm.i <- if (linf.mlm) {
       if ( .parallel.ip ) matrix(1, linf.mlm, 1) else diag(linf.mlm)
    } else {
      NULL
    }

    if (lalt.mlm + linf.mlm == 0) {
      use.mat <- use.mat.mlm <- cbind(M)  # lambda.p only
    }
    if (lalt.mlm + linf.mlm) {
      nc1 <- if (length(use.mat.mlm.a)) ncol(use.mat.mlm.a) else 0
      nc2 <- if (length(use.mat.mlm.i)) ncol(use.mat.mlm.i) else 0
      use.mat.mlm <- cbind(1, matrix(0, 1, nc1 + nc2))
      if (lalt.mlm)
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, lalt.mlm, 1),
                                   use.mat.mlm.a,
                                   if (length(use.mat.mlm.i) == 0) NULL
                                   else matrix(0, lalt.mlm, nc2)))
      if (linf.mlm )
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, linf.mlm, 1 + nc1),
                                   use.mat.mlm.i))
      use.mat <- use.mat.mlm
    }  # lalt.mlm + linf.mlm






    use.mat.mix <- if ( ( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,1,0,1,  0,1,0,0,0,  0,0,0,1,0), 5, 3) else
    if (!( .eq.ap ) &&  ( .eq.ip ))
      matrix(c(1,0,0,0,1,  0,1,0,0,0,  0,0,1,0,0,
               0,0,0,1,0), 5, 4) else
    if ( ( .eq.ap ) && !( .eq.ip ))
      matrix(c(1,0,1,0,0,  0,1,0,0,0,  0,0,0,1,0,
               0,0,0,0,1), 5, 4) else
      diag(5)

    tmp3.TF <- ( .tmp3.TF )
    tmp3.TF <- tmp3.TF[-(6:7)]  # zz for now. temporary only.
    use.mat.mix <- use.mat.mix[tmp3.TF, , drop = FALSE]

    mincl <- apply(use.mat.mix, 2, min)
    maxcl <- apply(use.mat.mix, 2, max)
    use.mat.mix <- use.mat.mix[, mincl != 0 | maxcl != 0, drop = FALSE]

    if (lalt.mix + linf.mix > 0)
      use.mat <- use.mat.mix





    if (lalt.mlm + linf.mlm > 0 &&
        lalt.mix + linf.mix > 0) {
      use.mat <- rbind(use.mat.mix,
                 matrix(0, nrow(use.mat.mlm) - 1, ncol(use.mat.mix)))
      use.mat <- cbind(use.mat, matrix(0, nrow(use.mat),
                                          ncol(use.mat.mlm) - 1))
      use.mat[row(use.mat) > nrow(use.mat.mix) &
              col(use.mat) > ncol(use.mat.mix)] <- use.mat.mlm[-1, -1]
    }  # lalt.mlm + linf.mlm > 0 && lalt.mix + linf.mix > 0



    if (is.null(constraints)) {
      constraints <-
        cm.VGAM(use.mat, x = x, apply.int = TRUE,  # FALSE
                bool = .eq.ap || .eq.ip || .parallel.ap || .parallel.ip ,
                constraints = constraints)  # FALSE
    }



    if (lalt.mix + linf.mix + lalt.mlm + linf.mlm)
      constraints <-
        cm.zero.VGAM(constraints, x = x, .zero , M = M, M1 = M1,
                     predictors.names = paste0(predictors.names,
                                        names(predictors.names)))
  }), list( .zero = zero, .tmp3 = tmp3, .tmp3.TF = tmp3.TF,
            .eq.ap = eq.ap, .eq.ip = eq.ip,
            .parallel.ap = parallel.ap, .parallel.ip = parallel.ip,
            .lalt.mlm = lalt.mlm, .linf.mlm = linf.mlm,
            .lalt.mix = lalt.mix, .linf.mix = linf.mix
          ))),
  infos = eval(substitute(function(...) {
    list(M1 = .M1 ,
         Q1 = 1,
         link = .predictors.names ,  # .tmp3, as.vector strips names off
         link1parameter = as.logical( .lall.len <= 2),  # <= 1 safer
         mixture.links  = any(c( .lalt.mlm , .linf.mlm , .lalt.mix ,
                                 .linf.mix ) > 1),  # FALSE if NULL
         alt.mix = as.vector( .alt.mix ),  # Handles NULL
         alt.mlm = as.vector( .alt.mlm ),
         inf.mix = as.vector( .inf.mix ),
         inf.mlm = as.vector( .inf.mlm ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ),
         Support  = c( .lowsup , Inf, 1),  # a(b)c format as a,c,b.
         expected = TRUE,
         multipleResponses = FALSE,  # poissonff might be called if TRUE
         parameters.names = names( .predictors.names ),
         parent.name = c("poissonff", "pois"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
         baseparams.argnames  = "lambda",
         MM1 = 1,  # One parameter for 1 response (lambda)
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .llambda.p = llambda.p, .elambda.p = elambda.p,
           .llambda.a = llambda.a, .elambda.a = elambda.a,
           .llambda.i = llambda.i, .elambda.i = elambda.i,
           .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
           .alt.mix = alt.mix, .inf.mix = inf.mix,
           .lalt.mlm = lalt.mlm, .linf.mlm = linf.mlm,
           .lalt.mix = lalt.mix, .linf.mix = linf.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len
         ))),
  initialize = eval(substitute(expression({
    extra$indeta <- ( .indeta )  # Avoids recomputing it several times
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
    truncate <- as.vector( .truncate )
    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(y)
    M <- NOS * M1

    tmp3.TF <- ( .tmp3.TF )

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Since max.support = 9 is possible
              ncol.y.max = 1,
              out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gaitcombo.check(y, truncate = truncate,
                               alt.mlm = alt.mlm, alt.mix = alt.mix,
                               inf.mlm = inf.mlm, inf.mix = inf.mix,
                               max.support = .max.support )
    extra$skip.mix.a <- glist$skip.mix.a
    extra$skip.mix.i <- glist$skip.mix.i
    extra$skip.mlm.a <- glist$skip.mlm.a
    extra$skip.mlm.i <- glist$skip.mlm.i


    
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- as.vector( .type.fitted )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1
    extra$index.M <- iam(NA, NA, M, both = TRUE)  # Used in @weight


    predictors.names <- ( .predictors.names )  # Got it, named


    if (!length(etastart)) {
      lambda.a.init <-  lambda.i.init <-  # Needed
      lambda.p.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                               imu = .ilambda.p ,  # x = x,
                               ishrinkage = .ishrinkage ,
                               probs.y = .probs.y )
      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(lambda.p.init, .llambda.p , earg = .elambda.p ))


      pobs.mix.init <- numeric(n)
      if (tmp3.TF[2]) {  # lalt.mix > 0
        pobs.mix.init <- if (length( .ipobs.mix )) {
        rep_len( .ipobs.mix , n)
      } else {
        is.alt.mix <- rowSums(glist$y0.mix.a) > 0
        rep(sum(w[is.alt.mix]) / sum(w), n)
      }  
      }  # lalt.mix > 0


      if (tmp3.TF[3]) {  # Assign coln 3; lalt.mix > 1
        lambda.a.init <- if (length( .ilambda.a ))
          rep_len( .ilambda.a , n) else lambda.p.init  # A vector
        etastart[, 3] <-
          theta2eta(lambda.a.init, .llambda.a , earg = .elambda.a )
      }

      pstr.mix.init <- numeric(n)
      try.gridsearch.pstr.mix <- FALSE
      if (tmp3.TF[4]) {  # linf.mix > 0
        pstr.mix.init <- if (length( .ipstr.mix )) {
          rep_len( .ipstr.mix , n)
        } else {
          try.gridsearch.pstr.mix <- TRUE
          numeric(n)  # Overwritten by gridsearch
        }
      }  # linf.mix > 0


      if (tmp3.TF[5]) {  # linf.mix > 1
        lambda.i.init <- if (length( .ilambda.i ))
          rep_len( .ilambda.i , n) else lambda.p.init  # A vector
        etastart[, (extra$indeta[5, 'launch'])] <-
          theta2eta(lambda.i.init, .llambda.i , earg = .elambda.i )
      }  # linf.mix > 1


      if (tmp3.TF[6]) {  #  lalt.mlm
        if (length( .ipobs.mlm )) {
          pobs.mlm.init <- matrix( .ipobs.mlm , n, lalt.mlm,
                                  byrow = .byrow.ai )
        } else {
          pobs.mlm.init <- colSums(c(w) * extra$skip.mlm.a) / colSums(w)
          pobs.mlm.init <- pobs.mlm.init * as.vector( .imux )
          pobs.mlm.init <- matrix(pobs.mlm.init, n, lalt.mlm, byrow=TRUE)
        }
      } else {
        pobs.mlm.init <- matrix(0, n, 1)
      }


      try.gridsearch.pstr.mlm <- FALSE
      if (tmp3.TF[7]) {  #  linf.mlm
        try.gridsearch.pstr.mlm <- !(length( .ipstr.mlm ))
        pstr.mlm.init <- 0  # Might be overwritten by gridsearch

        if (length( .ipstr.mlm ))
          pstr.mlm.init <- as.vector( .ipstr.mlm )

        pstr.mlm.init <- matrix(pstr.mlm.init, n, linf.mlm,
                                byrow = .byrow.ai )
      } else {
        pstr.mlm.init <- matrix(0, n, 1)
      }






      gaitpois.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitpois(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  lambda.p    = extraargs$lambda.p,
                  lambda.a    = extraargs$lambda.a,
                  lambda.i    = extraargs$lambda.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }

 gaitpois.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitpois(y, pstr.mlm = pstr.mlm.val,
                  pstr.mix    = extraargs$pstr.mix,  # Differs here
                  lambda.p    = extraargs$lambda.p,
                  lambda.a    = extraargs$lambda.a,
                  lambda.i    = extraargs$lambda.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }

 gaitpois.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitpois(y, pstr.mix = pstr.mix.val, pstr.mlm = pstr.mlm.val,
                  lambda.p    = extraargs$lambda.p,
                  lambda.a    = extraargs$lambda.a,
                  lambda.i    = extraargs$lambda.i,
                  alt.mix     = extraargs$alt.mix,
                  alt.mlm     = extraargs$alt.mlm,
                  inf.mix     = extraargs$inf.mix,
                  inf.mlm     = extraargs$inf.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm, log = TRUE))
  }



      if (linf.mix + linf.mlm) {
        extraargs <- list(
              lambda.p    = lambda.p.init,
              lambda.a    = lambda.a.init,
              lambda.i    = lambda.i.init,
              alt.mix     = alt.mix,
              alt.mlm     = alt.mlm,
              inf.mix     = inf.mix,
              inf.mlm     = inf.mlm,
              truncate    = truncate,
              max.support = as.vector( .max.support ),
              pobs.mix    = pobs.mix.init ,
              pobs.mlm    = pobs.mlm.init )
        pre.warn <- options()$warn
        options(warn = -1)  # Ignore warnings during gridsearch 

        try.this <-
          if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
            grid.search2( .gpstr.mix ,  .gpstr.mlm ,
                         objfun = gaitpois.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- pstr.mlm.init
            grid.search ( .gpstr.mix ,
                         objfun = gaitpois.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- pstr.mix.init
            grid.search ( .gpstr.mlm ,
                         objfun = gaitpois.Loglikfun1.mlm,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          }

        options(warn = pre.warn)  # Restore warnings 
        if (any(is.na(try.this)))
          warning("gridsearch returned NAs. It's going to crash.",
                  immediate. = TRUE)
        if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
          pstr.mix.init <- rep_len(try.this["Value1"], n)
          pstr.mlm.init <- matrix(try.this["Value2"], n, linf.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mix = seq(5) / 100'",
                 " and/or 'gpstr.mlm = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mix) {
          pstr.mix.init <- rep_len(try.this["Value"], n)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mix = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mlm) {
          pstr.mlm.init <- matrix(try.this["Value"], n, linf.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. Try something like 'gpstr.mlm = seq(5) / 100'.")
        }
      }  # lalt.mix + lnf.mix



        


      while (any((vecTF <- pobs.mix.init + pstr.mix.init +
                           rowSums(pobs.mlm.init) +
                           rowSums(pstr.mlm.init) > 0.96875))) {
        pobs.mix.init[vecTF]   <- 0.875 * pobs.mix.init[vecTF]
        pstr.mix.init[vecTF]   <- 0.875 * pstr.mix.init[vecTF]
        pobs.mlm.init[vecTF, ] <- 0.875 * pobs.mlm.init[vecTF, ]
        pstr.mlm.init[vecTF, ] <- 0.875 * pstr.mlm.init[vecTF, ]
      }

      Numer <- 1 - rowSums(pobs.mlm.init) - rowSums(pstr.mlm.init) -
               pobs.mix.init - pstr.mix.init
        
      etastart.z <- if (lall.len == 0) NULL else
        multilogitlink(cbind(if (tmp3.TF[2]) pobs.mix.init else NULL,
                             if (tmp3.TF[4]) pstr.mix.init else NULL,
                             if (tmp3.TF[6]) pobs.mlm.init else NULL,
                             if (tmp3.TF[7]) pstr.mlm.init else NULL,
                             Numer))
      if (!is.matrix(etastart.z)) etastart.z <- cbind(etastart.z)

      nextone <- 1  # Might not be used actually
      if (tmp3.TF[2]) {
        etastart[, 2] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[4]) {  # Coln 2 or 4
        etastart[, (extra$indeta[4, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[6]) {
        ind6 <- (extra$indeta[6, 'launch']):(extra$indeta[6, 'finish'])
        etastart[, ind6] <- etastart.z[, nextone:(nextone+lalt.mlm-1)]
        nextone <- nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind7 <- (extra$indeta[7, 'launch']):(extra$indeta[7, 'finish'])
        etastart[, ind7] <- etastart.z[, nextone:ncol(etastart.z)]
      }
    }
  }), list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .ilambda.p = ilambda.p,
    .ilambda.a = ilambda.a,
    .ilambda.i = ilambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .ipstr.mix = ipstr.mix, .ipobs.mix = ipobs.mix,
    .ipstr.mlm = ipstr.mlm, .ipobs.mlm = ipobs.mlm,
    .byrow.ai = byrow.ai,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .imux = imux,
    .gpstr.mix = gpstr.mix,
    .gpstr.mlm = gpstr.mlm,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .indeta = indeta,
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"
      }
    type.fitted <-
      match.arg(type.fitted[1],
                c("mean", "lambdas",
                  "pobs.mlm", "pstr.mlm",
                  "pobs.mix", "pstr.mix",
                  "Pobs.mix", "Pstr.mix",
                  "nonspecial", "Numer", "Denom.p", "sum.mlm.i",
                  "ptrunc.p", "cdf.max.s"))[1]

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <- lambda.p  # Needed; and answer not corrupted

    if (any(tmp3.TF[c(3, 5)])) {  # At least one lambda.[ai]
      ind.lambda.z <- extra$indeta[c(1, 3, 5), 'launch']  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value
      iptr <- 1
      lambda.a <- if (!tmp3.TF[3]) lambda.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[5]) lambda.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.i , .elambda.i )
    }  # lalt.mix + linf.mix > 0







    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")


      Nextone <- 0  # Might not be used actually
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitcombo.pois(lambda.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              lambda.a = lambda.a, lambda.i = lambda.i,
              truncate = truncate, max.support = max.support)

    n.eta <- nrow(eta)
    Denom.p <- c(Bits[["cdf.max.s"]] - Bits[["SumT0.p"]] -
                 Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])

    if (any(Denom.p == 0)) {
      smallval <- min(Denom.p[Denom.p > 0])
      Denom.p[Denom.p == 0] <- 1e-09  # smallval
      warning("0s found in variable 'Denom.p'. Trying to fix it.")
    }
    Numer <- c(1 - pobs.mix - pstr.mix -
               (if (lalt.mlm) rowSums(pobs.mlm) else 0) -
               (if (linf.mlm) rowSums(pstr.mlm) else 0))

    if (!lalt.mlm && type.fitted %in% c("pobs.mlm")) {
      warning("No altered MLM values; returning an NA")
      return(NA)
    }
    if (!linf.mlm && type.fitted %in% c("sum.mlm.i", "pstr.mlm")) {
      warning("No inflated MLM values; returning an NA")
      return(NA)
    }
    if (!lalt.mix && type.fitted %in% c("Pobs.mix")) {
      warning("No altered mixture values; returning an NA")
      return(NA)
    }
    if (!linf.mix && type.fitted %in% c("sum.mix.i", "Pstr.mix")) {
      warning("No inflated mixture values; returning an NA")
      return(NA)
    }

    if (lalt.mix) {
      tmp13 <-  # dpois() does not retain the matrix format
        dpois(matrix(alt.mix,  NROW(eta), lalt.mix, byrow = TRUE),
              matrix(lambda.a, NROW(eta), lalt.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(NROW(eta), lalt.mix)
      dimnames(tmp13) <- list(rownames(eta), as.character(alt.mix))
      propn.mat.a <- tmp13
    }

    if (linf.mix) {
      tmp55 <-  # dpois() does not retain the matrix format
        dpois(matrix(inf.mix,  NROW(eta), linf.mix, byrow = TRUE),
              matrix(lambda.i, NROW(eta), linf.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(NROW(eta), linf.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(inf.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "lambdas"    = cbind(lambda.p,
                           if (tmp3.TF[3]) lambda.a else NULL,
                           if (tmp3.TF[5]) lambda.i else NULL),
      "pobs.mlm"   = pobs.mlm,  # aka omegamat, n x lalt.mlm
      "pstr.mlm"   = pstr.mlm,  # aka phimat, n x linf.mlm
      "pobs.mix"   = pobs.mix,  # n-vector
      "pstr.mix"   = pstr.mix,  # n-vector
      "Pobs.mix"   = c(pobs.mix) * propn.mat.a,  # matrix
      "Pstr.mix"   = c(pstr.mix) * propn.mat.i,
      "nonspecial" = Numer * (1 -
         (Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) / Denom.p),
      "Numer"      = Numer,
      "Denom.p"    = Denom.p,
      "sum.mlm.i"  = pstr.mlm + Numer *
             dpois(matrix(inf.mlm,  NROW(eta), linf.mlm, byrow = TRUE),
                   matrix(lambda.p, NROW(eta), linf.mlm)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)

    ynames.pobs.mlm <- as.character(alt.mlm)  # Works with NULLs
    ynames.pstr.mlm <- as.character(inf.mlm)  # Works with NULLs
    if (length(ans))
      label.cols.y(ans, NOS = NOS, colnames.y =
      switch(type.fitted,
             "lambdas"   = c("lambda.p", "lambda.a",  # Some colns NA
                     "lambda.i")[(tmp3.TF[c(1, 3, 5)])],
             "Pobs.mix"  = as.character(alt.mix),
             "Pstr.mix"  = as.character(inf.mix),
             "pobs.mlm"  = ynames.pobs.mlm,
             "sum.mlm.i" = ,  #
             "pstr.mlm"  = ynames.pstr.mlm,
             extra$colnames.y)) else ans
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  last = eval(substitute(expression({
    pred.names <- c( .predictors.names )  # Save it
    link.names <- as.vector( .predictors.names )
    parameter.names <- names(pred.names)
    predictors.names <- NULL
    for (jay in seq(M))
      predictors.names <- c(predictors.names,
        namesof(parameter.names[jay], link.names[jay], tag = FALSE,
                earg = list()))  # This line isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .elambda.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[2])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[3])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .elambda.a )
    if (tmp3.TF[4])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[5])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .elambda.i )
    if (tmp3.TF[6]) {  # lalt.mlm
      for (ii in seq(lalt.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # lalt.mlm
    if (tmp3.TF[7]) {  # linf.mlm
      for (ii in seq(linf.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # linf.mlm
  }), list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .predictors.names = predictors.names,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm 
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <- lambda.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5)])) {  # At least one lambda.[ai]
      ind.lambda.z <- extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value
      iptr <- 1
      lambda.a <- if (!tmp3.TF[3]) lambda.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[5]) lambda.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.i , .elambda.i )
    }  # lalt.mix + linf.mix > 0
    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitpois(y, lambda.p, log = TRUE,  # byrow.ai = F,
                  alt.mix = alt.mix, inf.mix = inf.mix,
                  alt.mlm = alt.mlm, inf.mlm = inf.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  lambda.a = lambda.a, lambda.i = lambda.i, 
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitpoisson"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    small. <- 1e-14
    pobs.mix <- pstr.mix <- small.
    pobs.mlm <- pstr.mlm <- matrix(small., NROW(eta), 1) # 4 rowSums()
    lambda.a <- lambda.i <- 1  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lambda.p <-
      cbind(eta2theta(eta[, 1], .llambda.p , earg = .elambda.p ))
    iptr <- 1
    ind.lambda.z <- 1  # Points to lambda.p only.
    if (lalt.mix + linf.mix > 1) {  # At least one lambda.[ai]
      ind.lambda.z <- extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value

      lambda.a <- if (lalt.mix <= 1) lambda.p else {
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.a ,
                  earg = .elambda.a )
      }
      lambda.i <- if (linf.mix <= 1) lambda.p else {
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.i ,
                  earg = .elambda.i )
      }
    }  # lalt.mix + linf.mix > 0

    
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm)
    okay.mix <-
      all(is.finite(lambda.p)) && all(0 < lambda.p) &&
      all(lambda.p < .max.support ) &&
      all(is.finite(lambda.a)) && all(0 < lambda.a) &&
      all(is.finite(lambda.i)) && all(0 < lambda.i) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(pobs.mix + pstr.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1) 4 rowSums()
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    lambda.a <- lambda.i <- lambda.p  # Needed
    ind.lambda.z <- 1  # Points to lambda.p only.

    if (lalt.mix + linf.mix > 1) {  # At least one lambda.[ai]
      ind.lambda.z <- object@extra$indeta[c(1, 3, 5), 1]  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value
      iptr <- 1
      lambda.a <- if (lalt.mix <= 1) lambda.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.a , .elambda.a )
      lambda.i <- if (linf.mix <= 1) lambda.p else
        eta2theta(eta[, (iptr <- iptr + 2)], .llambda.i , .elambda.i )
    }  # lalt.mix + linf.mix > 0

    tmp3.TF <- ( .tmp3.TF )
    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
        Nextone <- Nextone + linf.mlm  # Not needed
      }
    }  # lall.len

    rgaitpois(nsim * length(lambda.p), lambda.p,
              pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
              pobs.mix = pobs.mix, pstr.mix = pstr.mix,
              lambda.a = lambda.a, lambda.i = lambda.i, 
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              truncate = .truncate , max.support = .max.support )
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .tmp3.TF = tmp3.TF,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({

    tmp3.TF <- ( .tmp3.TF )
    calA.p  <- tmp3.TF[2]
    calI.p  <- tmp3.TF[4]
    calA.np <- tmp3.TF[6]
    calI.np <- tmp3.TF[7]

    Denom1.a <- Denom1.i <- Denom2.i <- 0  # Denom2.a is not needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lalt.mix <- length((alt.mix <- as.vector( .alt.mix )))
    linf.mix <- length((inf.mix <- as.vector( .inf.mix )))
    lalt.mlm <- length((alt.mlm <- as.vector( .alt.mlm )))
    linf.mlm <- length((inf.mlm <- as.vector( .inf.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )

    lall.len <- lalt.mix + linf.mix + lalt.mlm + linf.mlm
    pobs.mix <- pstr.mix <- 0
    pobs.mlm <- pstr.mlm <- 0  # matrix(0, NROW(eta), 1)  # 4 rowSums()
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <- lambda.p  # Needed; doesnt corrupt answer

    if (any(tmp3.TF[c(3, 5)])) {  # At least one lambda.[ai]
      ind.lambda.z <- extra$indeta[c(1, 3, 5), 'launch']  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value
      lambda.a <- if (!tmp3.TF[3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
    }  # lalt.mix + linf.mix > 0








    if (lall.len) {  # A MLM was fitted
      allprobs <- multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
                                 inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[6]) {
        ind2 <- (Nextone + 1):(Nextone + lalt.mlm)
        pobs.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(alt.mlm))
        Nextone <- Nextone + lalt.mlm
      }
      if (tmp3.TF[7]) {
        ind2 <- (Nextone + 1):(Nextone + linf.mlm)
        pstr.mlm <- allprobs[, ind2, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(inf.mlm))
      }
    }  # lall.len


    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")



    is.alt.mixed <- if (tmp3.TF[2])
      rowSums(extra$skip.mix.a) > 0 else rep(FALSE, n)
    is.inf.mixed <- if (tmp3.TF[4])
      rowSums(extra$skip.mix.i) > 0 else rep(FALSE, n)
    is.alt.mlmed <- if (tmp3.TF[6])
      rowSums(extra$skip.mlm.a) > 0 else rep(FALSE, n)
    is.inf.mlmed <- if (tmp3.TF[7])
      rowSums(extra$skip.mlm.i) > 0 else rep(FALSE, n)

    is.ns <- !is.alt.mlmed & !is.inf.mlmed  &
             !is.alt.mixed & !is.inf.mixed  # & !is.truncd

    dl.dlambda.p <- y / lambda.p - 1  # == dl.dlambda.p.usual
    dl.dlambda.p[!is.ns] <- 0  # For is.alt.mixed & is.alt.mlmed
    prob.mlm.a <- if (lalt.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (linf.mlm) rowSums(pstr.mlm) else 0  # scalar okay
    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)


    sumD.mix.1a.p <- sumD.mix.2a.p <- matrix(0, n, NOS)
    if (lalt.mix > 0) {  # \calA_p
      Da.mix.0mat.a <-  # Matches naming convention further below
      Da.mix.1mat.a <- Da.mix.2mat.a <- matrix(0, n, lalt.mix)
      for (jay in seq(lalt.mix)) {
        aval <- alt.mix[jay]
        sumD.mix.1a.p <- sumD.mix.1a.p + pmf.deriv1(aval, lambda.p)
        sumD.mix.2a.p <- sumD.mix.2a.p + pmf.deriv2(aval, lambda.p)
        pmf.a <- dpois(aval, lambda.a)
        Da.mix.0mat.a[, jay] <- pmf.a
        Da.mix.1mat.a[, jay] <- pmf.deriv1(aval, lambda.a)
      }
      Denom1.a <- rowSums(Da.mix.1mat.a)  # aka sumD.mix.1a.a
      Denom2.a <- rowSums(Da.mix.2mat.a)  # Needed; sumD.mix.2a.a <-
    }  # lalt.mix > 0




    if (linf.mix) {
      Di.mix.0mat.i <-  # wrt inflated distribution
      Di.mix.1mat.i <- Di.mix.2mat.i <- matrix(0, n, linf.mix)
      Dp.mix.0mat.i <-  # wrt parent distribution
      Dp.mix.1mat.i <- Dp.mix.2mat.i <- matrix(0, n, linf.mix)
        for (jay in seq(linf.mix)) {
          ival <- inf.mix[jay]
          pmf.i <- dpois(ival, lambda.i)
          Di.mix.0mat.i[, jay] <- pmf.i
          Di.mix.1mat.i[, jay] <- pmf.deriv1(ival, lambda.i)
          Di.mix.2mat.i[, jay] <- pmf.deriv2(ival, lambda.i)
          pmf.p <- dpois(ival, lambda.p)
          Dp.mix.0mat.i[, jay] <- pmf.p
          Dp.mix.1mat.i[, jay] <- pmf.deriv1(ival, lambda.p)
          Dp.mix.2mat.i[, jay] <- pmf.deriv2(ival, lambda.p)
        }  # jay
      Denom1.i <- rowSums(Di.mix.1mat.i)
      Denom2.i <- rowSums(Di.mix.2mat.i)
    }  # linf.mix



    bits <- moments.gaitcombo.pois(lambda.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              alt.mix = alt.mix, inf.mix = inf.mix,
              alt.mlm = alt.mlm, inf.mlm = inf.mlm,
              lambda.a = lambda.a, lambda.i = lambda.i,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p <- sumD.mlm.2a.p <- matrix(0, n, NOS)
    if (lalt.mlm)
      for (aval in alt.mlm) {
        sumD.mlm.1a.p <- sumD.mlm.1a.p + pmf.deriv1(aval, lambda.p)
        sumD.mlm.2a.p <- sumD.mlm.2a.p + pmf.deriv2(aval, lambda.p)
      }


    Denom0.p <- c(bits[["cdf.max.s"]] - bits[["SumT0.p"]] -
                  bits[["SumA0.mix.p"]] - bits[["SumA0.mlm.p"]])
    Numer <- 1 - pobs.mix - pstr.mix - prob.mlm.a - prob.mlm.i
    Denom0.a <- c(bits[["SumA0.mix.a"]])  # Not .p
    Denom0.i <- c(bits[["SumI0.mix.i"]])


 

    Dp.mlm.0Mat.i <-  # wrt parent distribution
    Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, NOS)
    if (linf.mlm > 0) {
      Dp.mlm.0Mat.i <-  # wrt parent distribution
      Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, linf.mlm)
      for (jay in seq(linf.mlm)) {
        ival <- inf.mlm[jay]
        pmf.p <- dpois(ival, lambda.p)
        Dp.mlm.0Mat.i[, jay] <- pmf.p
        Dp.mlm.1Mat.i[, jay] <- pmf.deriv1(ival, lambda.p)
        Dp.mlm.2Mat.i[, jay] <- pmf.deriv2(ival, lambda.p)
      }  # jay
    }  # linf.mlm


    

    sumD.1t.a <- sumD.2t.a <-
    sumD.1t.i <- sumD.2t.i <-
    sumD.1t.p <- sumD.2t.p <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p <- sumD.1t.p + pmf.deriv1(tval, lambda.p)
        sumD.2t.p <- sumD.2t.p + pmf.deriv2(tval, lambda.p)
        sumD.1t.a <- sumD.1t.a + pmf.deriv1(tval, lambda.a)
        sumD.2t.a <- sumD.2t.a + pmf.deriv2(tval, lambda.a)
        sumD.1t.i <- sumD.1t.i + pmf.deriv1(tval, lambda.i)
        sumD.2t.i <- sumD.2t.i + pmf.deriv2(tval, lambda.i)
      }
    sumD.1t.p <- sumD.1t.p + dpois(max.support  , lambda.p)
    sumD.2t.p <- sumD.2t.p + dpois(max.support-1, lambda.p) -
                                     dpois(max.support  , lambda.p)
    sumD.1t.a <- sumD.1t.a + dpois(max.support  , lambda.a)
    sumD.2t.a <- sumD.2t.a + dpois(max.support-1, lambda.a) -
                                     dpois(max.support  , lambda.a)
    sumD.1t.i <- sumD.1t.i + dpois(max.support  , lambda.i)
    sumD.2t.i <- sumD.2t.i + dpois(max.support-1, lambda.i) -
                                     dpois(max.support  , lambda.i)


    Denom1.p <- c(-sumD.1t.p - sumD.mlm.1a.p - sumD.mix.1a.p)
    Denom2.p <- c(-sumD.2t.p - sumD.mlm.2a.p - sumD.mix.2a.p)


    d0B.pi.mlm <- Dp.mlm.0Mat.i / Denom0.p
    d1B.pi.mlm <- Dp.mlm.1Mat.i / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.i * Denom1.p / Denom0.p^2
    d2B.pi.mlm <-     Dp.mlm.2Mat.i / Denom0.p -
                  2 * Dp.mlm.1Mat.i * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.i * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.i * (Denom1.p^2) / Denom0.p^3


    Delta.mlm <- if (linf.mlm > 0) {
      pstr.mlm + Numer * d0B.pi.mlm  # n x linf.mlm.
    } else {
      matrix(0, n, 1)  # If linf.mlm == 0, for rowSums().
    }


    if (linf.mix > 0) {
      d0A.i <- Di.mix.0mat.i / Denom0.i
      d0B.pi.mix <- Dp.mix.0mat.i / Denom0.p
      Delta.mix <- pstr.mix * d0A.i + Numer * d0B.pi.mix


      d1A.i <- (Di.mix.1mat.i - Di.mix.0mat.i *
                Denom1.i / Denom0.i) / Denom0.i
      d2A.i <- (Di.mix.2mat.i - (2 * Di.mix.1mat.i * Denom1.i +
                Di.mix.0mat.i * Denom2.i) / Denom0.i +
            2 * Di.mix.0mat.i * (Denom1.i / Denom0.i)^2) / Denom0.i

      d1B.pi.mix <- Dp.mix.1mat.i / Denom0.p -
                    Dp.mix.0mat.i * Denom1.p / Denom0.p^2
      d2B.pi.mix <-     Dp.mix.2mat.i / Denom0.p -
                    2 * Dp.mix.1mat.i * Denom1.p / Denom0.p^2 -
                        Dp.mix.0mat.i * Denom2.p / Denom0.p^2 +
                    2 * Dp.mix.0mat.i * (Denom1.p^2) / Denom0.p^3
    }  # linf.mix > 0


    if (lalt.mix) {
      d0A.a <- Da.mix.0mat.a / Denom0.a
      d1A.a <- Da.mix.1mat.a / Denom0.a -
               Da.mix.0mat.a * Denom1.a / Denom0.a^2
      d2A.a <- (Da.mix.2mat.a - (2 * Da.mix.1mat.a * Denom1.a +
                Da.mix.0mat.a * Denom2.a) / Denom0.a +
            2 * Da.mix.0mat.a * (Denom1.a / Denom0.a)^2) / Denom0.a
    }  # lalt.mix




    dl.dlambda.a <-
    dl.dlambda.i <- numeric(n)  # Replace some elts below
    dl.dpstr.mix <- (-1) / Numer  # \notin A, I, T
    dl.dpobs.mix <- numeric(n)  # Replace some elts below
    dl.dpobs.mix[is.ns] <- (-1) / Numer[is.ns]
    dl.dpstr.mix[is.alt.mixed] <- 0
    dl.dpstr.mix[is.alt.mlmed] <- 0
    dl.dpobs.mlm <-
    dl.dpstr.mlm <- matrix(0, n, 1)  # Might not be needed



    if (tmp3.TF[6] && lalt.mlm) {  # aka \calA_{np}
      dl.dpobs.mlm <- matrix(-1 / Numer, n, lalt.mlm)  # \notin calS
      dl.dpobs.mlm[!is.ns, ] <- 0  # For alt.mix only really
      for (jay in seq(lalt.mlm)) {
        aval <- alt.mlm[jay]
        is.alt.j.mlm <- extra$skip.mlm.a[, jay]  # Logical vector
        tmp7a <- 1 / pobs.mlm[is.alt.j.mlm, jay]
        dl.dpobs.mlm[is.alt.j.mlm, jay] <- tmp7a
      }  # jay
    }  # lalt.mlm



    dl.dlambda.p[is.ns] <- dl.dlambda.p[is.ns] -
                           (Denom1.p / Denom0.p)[is.ns]


    
    if (tmp3.TF[7] && linf.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, linf.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For alt.mlm and alt.mix

      for (jay in seq(linf.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7i <- Numer * d1B.pi.mlm[, jay] / Delta.mlm[, jay]
        dl.dlambda.p[is.inf.j.mlm] <- tmp7i[is.inf.j.mlm]


        if (tmp3.TF[6] && lalt.mlm) {  # 20200919
          tmp9i <- (-d0B.pi.mlm[, jay] / Delta.mlm[, jay])
          dl.dpobs.mlm[is.inf.j.mlm, ] <- tmp9i[is.inf.j.mlm]
        }


        if (tmp3.TF[2] && lalt.mix) {  # 20200919 linf.mix wrong
          tmp2 <- (-d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
          dl.dpobs.mix[is.inf.j.mlm] <- tmp2[is.inf.j.mlm]
        }


        if (tmp3.TF[4] && linf.mix) {  # 20200916
          tmp2 <- (-d0B.pi.mlm[, jay] / Delta.mlm[, jay])
          dl.dpstr.mix[is.inf.j.mlm] <- tmp2[is.inf.j.mlm]
        }

          
        tmp9 <- (0 - d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
        tmp8 <- (1 - d0B.pi.mlm[, jay]) / Delta.mlm[, jay]
        dl.dpstr.mlm[is.inf.j.mlm, ] <- tmp9[is.inf.j.mlm]
        dl.dpstr.mlm[is.inf.j.mlm, jay] <- tmp8[is.inf.j.mlm]
      }  # jay
    }  # linf.mlm > 0


    



    if (tmp3.TF[2] && lalt.mix) {  # aka \calA_{p}
      dl.dpobs.mix[is.alt.mixed] <- 1 / pobs.mix[is.alt.mixed]

      if (tmp3.TF[3] && lalt.mix > 1)
        for (jay in seq(lalt.mix)) {
          is.alt.j.mix <- extra$skip.mix.a[, jay]  # Logical vector
          tmp2 <- d1A.a[, jay] / d0A.a[, jay]
          dl.dlambda.a[is.alt.j.mix] <- tmp2[is.alt.j.mix]  # ccc.
        }  # jay
    }  # lalt.mix


    if (tmp3.TF[4] && linf.mix > 0) {  # aka \calI_{p}
      for (jay in seq(linf.mix)) {
        ival <- inf.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7b <- Numer * d1B.pi.mix[, jay] / Delta.mix[, jay]
        dl.dlambda.p[is.inf.j.mix] <- tmp7b[is.inf.j.mix]
        if (tmp3.TF[2] && lalt.mix) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpobs.mix[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }
        tmp8 <- (d0A.i[, jay] - d0B.pi.mix[, jay]) / Delta.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]

        if (linf.mix > 1) {
          tmp2 <- pstr.mix * d1A.i[, jay] / Delta.mix[, jay]
          dl.dlambda.i[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }

        if (tmp3.TF[6] && lalt.mlm) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpobs.mlm[is.inf.j.mix, ] <- tmp2[is.inf.j.mix]
        }
        if (tmp3.TF[7] && linf.mlm) {
          tmp2 <- (-d0B.pi.mix[, jay] / Delta.mix[, jay])
          dl.dpstr.mlm[is.inf.j.mix, ] <- tmp2[is.inf.j.mix]
        }
      }  # jay
    }  # linf.mix > 0



    tmp3.TF <- !is.na(rowSums(extra$indeta))



    if (lall.len) {  # MLM fitted
      all4.dldp <- cbind(
        if (tmp3.TF[2]) dl.dpobs.mix else NULL,
        if (tmp3.TF[4]) dl.dpstr.mix else NULL,
        if (tmp3.TF[6]) dl.dpobs.mlm else NULL,
        if (tmp3.TF[7]) dl.dpstr.mlm else NULL)

      dl.deta.mlm4 <- matrix(0, n, ncol(all4.dldp))
      ind.quad.mlm <- seq(ncol(all4.dldp))  # Contiguous & starts at 1

      for (jay in ind.quad.mlm) {
        for (sss in ind.quad.mlm) {  # Dont need all4probs, allprobs ok
          dl.deta.mlm4[, jay] <- dl.deta.mlm4[, jay] +
            allprobs[, sss] * ((sss == jay) - allprobs[, jay]) *
            all4.dldp[, sss]
        }  # sss
      }  # jay
    }  # lall.len


    

    dlambda.p.deta <- dtheta.deta(lambda.p, .llambda.p , .elambda.p )
    if (tmp3.TF[3])
      dlambda.a.deta <- dtheta.deta(lambda.a, .llambda.a , .elambda.a )
    if (tmp3.TF[5])
      dlambda.i.deta <- dtheta.deta(lambda.i, .llambda.i , .elambda.i )


    iptr <- 0
    ansd <- cbind(lambda.p = c(dl.dlambda.p * dlambda.p.deta))
    ansd <- cbind(ansd,
      pobs.mix = if (tmp3.TF[2]) dl.deta.mlm4[, (iptr <- iptr+1)],
      lambda.a = if (tmp3.TF[3]) dl.dlambda.a * dlambda.a.deta)
    ansd <- cbind(ansd,
      pstr.mix = if (tmp3.TF[4]) dl.deta.mlm4[, (iptr <- iptr+1)],
      lambda.i = if (tmp3.TF[5]) dl.dlambda.i * dlambda.i.deta)
    if (any(tmp3.TF[6:7])) {  # The remainder
      ansd <- cbind(ansd, dl.deta.mlm4[, (iptr+1):ncol(dl.deta.mlm4)])
      ind3 <- (ncol(ansd) - lalt.mlm - linf.mlm + 1):ncol(ansd)
      colnames(ansd)[ind3] <- as.character(c(alt.mlm, inf.mlm))
    }

    c(w) * ansd
  }), list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .alt.mix = alt.mix, .inf.mix = inf.mix,
    .alt.mlm = alt.mlm, .inf.mlm = inf.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({


    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    cond.EY.p <- c(lambda.p -
                   bits[["SumI1.mlm.p"]] - bits[["SumI1.mix.p"]] -
                   bits[["SumT1.p"]] -
                   bits[["SumA1.mlm.p"]] - bits[["SumA1.mix.p"]]) / c(
                   Denom0.p - bits[["SumI0.mix.p"]] -
                              bits[["SumI0.mlm.p"]])
    probns <- Numer * (1 - c(bits[["SumI0.mix.p"]] +
                             bits[["SumI0.mlm.p"]]) / Denom0.p)






    zero0n <- numeric(n)
    ned2l.dpobs.mix.lambda.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.lambda.a <- zero0n  # Final; nothing to do
    ned2l.dpobs.mix.lambda.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.lambda.p <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.lambda.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.lambda.i <- zero0n  # mB overwritten below


      posn.pobs.mix <- as.vector(extra$indeta[2, 'launch'])
      posn.lambda.a <- as.vector(extra$indeta[3, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[4, 'launch'])
      posn.lambda.i <- as.vector(extra$indeta[5, 'launch'])
      posn.pobs.mlm <- as.vector(extra$indeta[6, 'launch'])
      posn.pstr.mlm <- as.vector(extra$indeta[7, 'launch'])


    ned2l.dpstr.mix2         <-  # Elt (4, 4)
    ned2l.dpobs.mlm.pstr.mix <-  # Elts (4, >=6)
    ned2l.dpobs.mix.pstr.mix <- probns / Numer^2  # ccc Elt (2, 4)
    if (all(c(lalt.mix, linf.mlm) > 0))
    ned2l.dpobs.mix.pstr.mlm <- matrix(probns / Numer^2, n, linf.mlm)
    if (all(c(linf.mix, linf.mlm) > 0))
    ned2l.dpstr.mix.pstr.mlm <- matrix(probns / Numer^2, n, linf.mlm)
    



    ned2l.dlambda.p2 <- probns * (cond.EY.p / lambda.p^2 +  # ccc
      Denom2.p / Denom0.p - (Denom1.p / Denom0.p)^2) + 
      (if (tmp3.TF[4] && linf.mix) Numer *
      rowSums(Numer * (d1B.pi.mix^2) / Delta.mix - d2B.pi.mix) else 0) +
      (if (tmp3.TF[7] && linf.mlm) Numer *
      rowSums(Numer * (d1B.pi.mlm^2) / Delta.mlm - d2B.pi.mlm) else 0)


    wz[, iam(1, 1, M)] <- ned2l.dlambda.p2 * dlambda.p.deta^2


    ned2l.dpobs.mix2 <- 1 / pobs.mix + probns / Numer^2
    if (tmp3.TF[4] && linf.mix > 0) {
      ned2l.dpobs.mix2 <-  # More just below, ccc
      ned2l.dpobs.mix2 + rowSums(d0B.pi.mix^2 / Delta.mix)
    }
    if (tmp3.TF[7] && linf.mlm > 0) {
      ned2l.dpobs.mix2 <-  # ccc.
      ned2l.dpobs.mix2 + rowSums(d0B.pi.mlm^2 / Delta.mlm)
    }
    if (tmp3.TF[2] && lalt.mix > 0)
      wz[, iam(2, 2, M)] <- ned2l.dpobs.mix2  # Link done later


    if (tmp3.TF[3] && lalt.mix > 1) {
      ned2l.dlambda.a2 <- pobs.mix * (
        rowSums((Da.mix.1mat.a^2) / Da.mix.0mat.a) / Denom0.a -
        (Denom1.a / Denom0.a)^2)  # ccc.
      wz[, iam(3, 3, M)] <- ned2l.dlambda.a2 * dlambda.a.deta^2
    }






    if (tmp3.TF[4] && linf.mix > 0) {

      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.pi.mix)^2 / Delta.mix)


      if (tmp3.TF[2] && lalt.mix > 0)
        ned2l.dpobs.mix.lambda.p <-
        ned2l.dpobs.mix.lambda.p +
          rowSums(d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ned2l.dpstr.mix.lambda.p <-
      ned2l.dpstr.mix.lambda.p + rowSums(
        d1B.pi.mix * (1 + Numer * (d0A.i - d0B.pi.mix) / Delta.mix))


    if (all(tmp3.TF[c(2, 4)]))
      ned2l.dpobs.mix.pstr.mix <-  # ccc
      ned2l.dpobs.mix.pstr.mix +
        rowSums(-d0B.pi.mix * (d0A.i - d0B.pi.mix) / Delta.mix)

    }  # (tmp3.TF[4] && linf.mix > 0)



    
    if (all(tmp3.TF[c(2, 4, 7)])) {  # was  lalt.mix > 0 & Delta.mix
      ned2l.dpobs.mix.pstr.mix <-  # ccc.
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.pi.mlm^2 / Delta.mlm)
    }

    if (!is.na(posn.pobs.mix) && !is.na(posn.pstr.mix))
      wz[, iam(posn.pobs.mix, posn.pstr.mix, M)] <-
        ned2l.dpobs.mix.pstr.mix  # Link done later




    if (tmp3.TF[5] && linf.mix > 1) {  # \calI_{p}, includes \theta_iota

      ned2l.dlambda.p.lambda.i <- pstr.mix * Numer *
        rowSums(d1A.i * d1B.pi.mix / Delta.mix)  # ccc.
      wz[, iam(1, posn.lambda.i, M)] <- ned2l.dlambda.p.lambda.i *
        dlambda.p.deta * dlambda.i.deta  # All links done here

      ned2l.dlambda.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i^2) / Delta.mix - d2A.i)  # ccc.
      wz[, iam(posn.lambda.i, posn.lambda.i, M)] <-
        ned2l.dlambda.i2 * dlambda.i.deta^2

      if (tmp3.TF[2]) {  # tmp3.TF[4] is TRUE, given tmp3.TF[5]
        ned2l.dpobs.mix.lambda.i <-
          rowSums(-pstr.mix * d1A.i * d0B.pi.mix / Delta.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.lambda.i, M)] <-  # link done later
          ned2l.dpobs.mix.lambda.i  # * dlambda.i.deta
      }

      if (tmp3.TF[4]) {
        ned2l.dpstr.mix.lambda.i <- rowSums(  # ccc.
          d1A.i * (pstr.mix * (d0A.i - d0B.pi.mix) / Delta.mix - 1))
        wz[, iam(posn.pstr.mix, posn.lambda.i, M)] <-  # link done later
          ned2l.dpstr.mix.lambda.i  # * dlambda.i.deta
      }

      if (tmp3.TF[6]) {
        ned2l.dpobs.mlm.lambda.i <- rowSums(
          -pstr.mix * d0B.pi.mix * d1A.i / Delta.mix)  # ccc.
        for (uuu in seq(lalt.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.lambda.i, M)] <-
            ned2l.dpobs.mlm.lambda.i  # * dlambda.i.deta done later
      }
    }  # (tmp3.TF[5] && linf.mix > 1)



        
    if (tmp3.TF[7] && linf.mlm > 0) {  # \calI_{np}, includes \phi_s

      if (lalt.mix && tmp3.TF[2])
        ned2l.dpobs.mix.lambda.p <-  # ccc.
        ned2l.dpobs.mix.lambda.p +
          rowSums(d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm))

      ned2l.dpstr.mix.lambda.p <-  # ccc.
      ned2l.dpstr.mix.lambda.p + rowSums(
        d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.pi.mlm^2 / Delta.mlm)
      }
    }  # tmp3.TF[7] && linf.mlm > 0


    if (!is.na(posn.pobs.mix))  # Optional (1, 2) element:
      wz[, iam(1, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.lambda.p  # One link done later
    if (!is.na(posn.pstr.mix))  # Optional (1, 4) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.lambda.p  # One link done later

    if (!is.na(posn.pstr.mix))  # Optional (4, 4) element
      wz[, iam(posn.pstr.mix,
               posn.pstr.mix, M)] <- ned2l.dpstr.mix2  # Link done later






    if (tmp3.TF[6] && lalt.mlm) {  # \calA_{np}, includes \omega_s
      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (uuu in seq(lalt.mlm)) {  # Diagonal elts only
        wz[, iam(ofset + uuu,
                 ofset + uuu, M)] <- 1 / pobs.mlm[, uuu]
      }  # uuu

      tmp8a <- probns / Numer^2
      if (tmp3.TF[4] && linf.mix)
        tmp8a <- tmp8a + rowSums((d0B.pi.mix^2) / Delta.mix)
      if (tmp3.TF[7] && linf.mlm)
        tmp8a <- tmp8a + rowSums((d0B.pi.mlm^2) / Delta.mlm)
      for (uuu in seq(lalt.mlm))  # All elts
        for (vvv in uuu:lalt.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
          wz[, iam(ofset + uuu, ofset + vvv, M)] + tmp8a  # All elts
    }  # lalt.mlm


 

    if (tmp3.TF[6] && lalt.mlm) {

      init0.val <- if (tmp3.TF[7] && linf.mlm) rowSums(
        d1B.pi.mlm * (1 - Numer * d0B.pi.mlm / Delta.mlm)) else zero0n
      ned2l.dpobs.mlm.lambda.p <- init0.val  # Vector, not matrix

      if (tmp3.TF[4] && linf.mix)
        ned2l.dpobs.mlm.lambda.p <-
        ned2l.dpobs.mlm.lambda.p + rowSums(
          d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (vvv in seq(lalt.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.lambda.p
    }  # lalt.mlm > 0





    
    if (tmp3.TF[7] && linf.mlm > 0) {  # \calI_{np}, includes \phi_s
      init0.val <- probns / Numer^2
      if (linf.mix)
        init0.val <- init0.val + rowSums((d0B.pi.mix^2) / Delta.mix)
      ned2l.dpstr.mlm2 <-
        matrix(init0.val, n, linf.mlm * (linf.mlm + 1) / 2)

      for (uuu in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mlm2[, iam(uuu, uuu, linf.mlm)] <-
          ned2l.dpstr.mlm2[, iam(uuu, uuu, linf.mlm)] +
            ((sss == uuu) - d0B.pi.mlm[, sss])^2 / Delta.mlm[, sss]
      if (linf.mlm > 1) {
        for (uuu in 1:(linf.mlm-1))
          for (vvv in (uuu+1):linf.mlm)
            for (sss in seq(linf.mlm))
              ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] <-
              ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] +
              ((sss == uuu) - d0B.pi.mlm[, sss]) *
              ((sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss]
      }  # if (linf.mlm > 1)

      ofset <- posn.pstr.mlm - 1
      for (uuu in seq(linf.mlm))
        for (vvv in uuu:linf.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpstr.mlm2[, iam(uuu, vvv, linf.mlm)] 
    }  # linf.mlm > 0



    if (tmp3.TF[7] && linf.mlm > 0) {
      ned2l.dpstr.mlm.theta.p <- matrix(0, n, linf.mlm)
      for (vvv in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mlm.theta.p[, vvv] <-
          ned2l.dpstr.mlm.theta.p[, vvv] +
          d1B.pi.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss])

      if (linf.mix && tmp3.TF[4])
        ned2l.dpstr.mlm.theta.p <-
        ned2l.dpstr.mlm.theta.p + rowSums(
          d1B.pi.mix * (1 - Numer * d0B.pi.mix / Delta.mix))

      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(linf.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta.p[, vvv]
    }  # linf.mlm > 0




    if (linf.mlm && linf.mix > 1) {

      ned2l.dpstr.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.pi.mix * d1A.i / Delta.mix)

      for (vvv in seq(linf.mlm))
        wz[, iam(posn.lambda.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.i  # ccc.
    }  # linf.mlm && linf.mix > 1



    if (all(c(lalt.mlm, linf.mlm) > 0)) {
      ned2l.dpobs.mlm.pstr.mlm2 <-
        array(probns / Numer^2, c(n, lalt.mlm, linf.mlm))
      for (uuu in seq(lalt.mlm))
        for (vvv in seq(linf.mlm))
          for (sss in seq(linf.mlm))
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] - d0B.pi.mlm[, sss] *
              ((sss == vvv) - d0B.pi.mlm[, sss]) / Delta.mlm[, sss]

      if (tmp3.TF[4] && linf.mix)
        ned2l.dpobs.mlm.pstr.mlm2 <-
        ned2l.dpobs.mlm.pstr.mlm2 + rowSums(d0B.pi.mix^2 / Delta.mix)

      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pstr <- posn.pstr.mlm - 1
      for (uuu in seq(lalt.mlm))
        for (vvv in seq(linf.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pstr + vvv, M)] <-
            ned2l.dpobs.mlm.pstr.mlm2[, uuu, vvv] 
    }  # all(c(lalt.mlm, linf.mlm) > 0)






    if (all(c(lalt.mix, lalt.mlm) > 0)) {
      ned2l.dpobs.mix.pobs.mlm <- probns / Numer^2  # Initialize
      if (linf.mix)  # tmp3.TF[4]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.pi.mix^2 / Delta.mix)
      if (linf.mlm)  # tmp3.TF[7]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.pi.mlm^2 / Delta.mlm)

      for (uuu in seq(lalt.mlm))  # ccc.
        wz[, iam(posn.pobs.mix, posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pobs.mlm  # Link done later
    }



    if (all(c(lalt.mix, linf.mlm) > 0)) {  # all(tmp3.TF[c(2, 7)])
      if (linf.mix)  # tmp3.TF[4]
        ned2l.dpobs.mix.pstr.mlm <-
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.pi.mix^2 / Delta.mix)

      for (uuu in seq(linf.mlm))
        for (sss in seq(linf.mlm))
          ned2l.dpobs.mix.pstr.mlm[, uuu] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.pi.mlm[, sss]) *
                            d0B.pi.mlm[, sss] / Delta.mlm[, sss]

      for (uuu in seq(linf.mlm))  # ccc.
        wz[, iam(posn.pobs.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu]  # Link done later
    }



    if (all(c(linf.mix, lalt.mlm) > 0)) {  # all(tmp3.TF[c(4, 6)])
      if (linf.mlm)  # tmp3.TF[7]
        ned2l.dpobs.mlm.pstr.mix <-
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.pi.mlm^2 / Delta.mlm)
        ned2l.dpobs.mlm.pstr.mix <-  # tmp3.TF[4] && linf.mix
        ned2l.dpobs.mlm.pstr.mix -
        rowSums((d0A.i - d0B.pi.mix) * d0B.pi.mix / Delta.mix)

      for (uuu in seq(lalt.mlm))  # ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pstr.mix  # Link done later
    }



    if (all(c(linf.mix, linf.mlm) > 0)) {  # all(tmp3.TF[c(4, 7)])

      for (uuu in seq(linf.mlm))  # tmp3.TF[4]
        for (sss in seq(linf.mlm))
          ned2l.dpstr.mix.pstr.mlm[, uuu] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.pi.mlm[, sss]) *
              d0B.pi.mlm[, sss] / Delta.mlm[, sss]
      ned2l.dpstr.mix.pstr.mlm <-
      ned2l.dpstr.mix.pstr.mlm -
      rowSums((d0A.i - d0B.pi.mix) * d0B.pi.mix / Delta.mix)

      for (uuu in seq(linf.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu]  # Link done later
    }

 



 
    if (lall.len) {
      wz.4 <- matrix(0, n, M * (M + 1) / 2)  # Or == 0 * wz
      ind.rc <- setdiff(1:M, ind.lambda.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the MLM


 # Copy in the thetas values: the looping is overkill.
      for (uuu in ind.lambda.z)
        for (sss in seq(M))
          wz.4[, iam(uuu, sss, M)] <- wz[, iam(uuu, sss, M)]

       
  if (FALSE) {
      for (uuu in seq(lind.rc)) {
        for (vvv in uuu:lind.rc) {
          for (sss in seq(lind.rc)) {
            for (ttt in seq(lind.rc)) {
              wz.4[, iam(ind.rc[uuu], ind.rc[vvv], M)] <-
              wz.4[, iam(ind.rc[uuu], ind.rc[vvv], M)] +
                allprobs[, sss] *
                (max(0, sss == uuu) - allprobs[, uuu]) *
                wz[, iam(ind.rc[sss], ind.rc[ttt], M)] *
                allprobs[, ttt] *
                (max(0, ttt == vvv) - allprobs[, vvv])
            }  # ttt
          }  # sss
        }  # vvv
      }  # uuu
  }  # FALSE








      speed.up <- intercept.only && (
                  length(offset) == 1 || all(offset[1] == offset))



      ind.mlm <- iam(NA, NA, lind.rc, both = TRUE, diag = TRUE)
      n.use <- 2
      bread <- if (speed.up)
               -allprobs[1:n.use, ind.mlm$row, drop = FALSE] *
                allprobs[1:n.use, ind.mlm$col, drop = FALSE] else
               -allprobs[, ind.mlm$row, drop = FALSE] *
                allprobs[, ind.mlm$col, drop = FALSE]
      if (speed.up) {
        bread[, 1:lind.rc] <-
        bread[, 1:lind.rc] + allprobs[1:n.use, 1:lind.rc]
      } else {
        bread[, 1:lind.rc] <-
        bread[, 1:lind.rc] + allprobs[, 1:lind.rc]  # -use.refLevel
      }
      bread <- m2a(bread, M = lind.rc)  # Half wasteful really



      if (!length(extra$ind.wz.match)) {
        imat <- matrix(NA, lind.rc, lind.rc)
        for (jay in seq(lind.rc)) {
          iptr <- jay
          for (kay in (ind.rc[jay]):M) {
            if (!any(kay %in% ind.lambda.z)) {
              imat[jay, iptr] <-
                which(extra$index.M$row == ind.rc[jay] &
                      extra$index.M$col == kay)
              iptr <- iptr + 1
            }  # if
          }  # kay
        }  # jay
        ind.wz.match <- imat[cbind(ind.mlm$row.ind, ind.mlm$col.ind)]
        extra$ind.wz.match <- ind.wz.match  # Assign it once
      }  # !length(extra$ind.wz.match)
      filling <- if (speed.up)
        wz[1:n.use, extra$ind.wz.match, drop = FALSE] else
        wz[, extra$ind.wz.match, drop = FALSE]



      wz.5 <- mux5(filling, bread, M = lind.rc, matrix.arg = TRUE)

      wz.4[, extra$ind.wz.match] <- if (speed.up)
        matrix(wz.5[1, ], n, ncol(wz.5), byrow = TRUE) else c(wz.5)











 


      dstar.deta <- cbind(dlambda.p.deta,
                          if (tmp3.TF[3]) dlambda.a.deta else NULL,
                          if (tmp3.TF[5]) dlambda.i.deta else NULL)
      iptr <- 0
      if (length(ind.lambda.z))
      for (uuu in ind.lambda.z) {  # Could delete 3 for lambda.a (orthog)
        iptr <- iptr + 1
        for (ttt in seq(lind.rc)) {
          wz.4[, iam(uuu, ind.rc[ttt], M)] <- 0  # Initialize
          for (sss in seq(lind.rc)) {
            wz.4[, iam(uuu, ind.rc[ttt], M)] <-
            wz.4[, iam(uuu, ind.rc[ttt], M)] +
              allprobs[, sss] * (max(0, sss == ttt) - allprobs[, ttt]) *
              wz[, iam(uuu, ind.rc[sss], M)] * dstar.deta[, iptr]
          }  # sss
        }  # ttt
      }  # uuu


      wz <- wz.4  # Completed
    }  # lall.len


    if (lall.len) {  # A MLM was fitted
      mytiny <- (allprobs <       sqrt(.Machine$double.eps)) |
                (allprobs > 1.0 - sqrt(.Machine$double.eps))
      atiny <- rowSums(mytiny) > 0
      if (any(atiny)) {
        ind.diags <- setdiff(1:M, ind.lambda.z)  # Exclude thetas
        wz[atiny, ind.diags] <- .Machine$double.eps +
        wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
      }
    }  # lall.len




    c(w) * wz
  }), list( .truncate = truncate ))))
}  # gaitpoisson



















 dgaitbinom <-
  function(x, size.p, prob.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL,
           max.support = NULL,  # NA for binom?
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           byrow.ai = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           deflation = FALSE,  # Single logical
           log = FALSE) {
  log.arg <- log;  rm(log)
  if (!length(max.support))
    max.support <- max(size.p, size.a, size.i, na.rm = TRUE)  # Manually
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc   <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      max.support >= max(size.p, na.rm = TRUE))
    return(dbinom(x, size.p, prob.p, log = log.arg))


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (!deflation && any(pstr.mlm < 0, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")

  LLL <- max(length(x),        length(pobs.mix),   length(pstr.mix),
             length(size.p),   length(size.a),     length(size.i),
             length(prob.p),   length(prob.a),     length(prob.i))
  if (length(x)          < LLL) x          <- rep_len(x,          LLL)
  if (length(size.p)     < LLL) size.p     <- rep_len(size.p,     LLL)
  if (length(size.a)     < LLL) size.a     <- rep_len(size.a,     LLL)
  if (length(size.i)     < LLL) size.i     <- rep_len(size.i,     LLL)
  if (length(prob.p)     < LLL) prob.p     <- rep_len(prob.p,     LLL)
  if (length(prob.a)     < LLL) prob.a     <- rep_len(prob.a,     LLL)
  if (length(prob.i)     < LLL) prob.i     <- rep_len(prob.i,     LLL)
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dbinom(tval, size.p, prob.p)
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- pbinom(max.support, size.p, prob.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dbinom(x, size.p, prob.p) / denom.t)


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalt.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, lalt.mlm,
                          byrow = byrow.ai)
    sum.a <- .rowSums(pobs.mlm, LLL, lalt.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    for (aval in alt.mlm)
      suma <- suma + dbinom(aval, size.p, prob.p)  # Part i

    for (jay in seq(lalt.mlm)) {
      aval <- alt.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
          pmf0[vecTF] <- pobs.mlm[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # lalt.mlm



  pmf2.a <- pmf2.i <- 0
  if (lalt.mix) {
    allx.a <- lowsup:max(alt.mix)
    pmf2.a <- dgaitbinom(x,  # Outer distribution---mlm type
                         size.p = size.a,
                         prob.p = prob.a,
                         truncate = setdiff(allx.a, alt.mix),
                         max.support = max(alt.mix))
    for (aval in alt.mix) {  # Part ii added; cumulative
      suma <- suma + dbinom(aval, size.p, prob.p)
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }

  if (linf.mix) {
    allx.i <- if (length(inf.mix)) lowsup:max(inf.mix) else NULL
    pmf2.i <- dgaitbinom(x,  # Outer distribution---mlm type
                         size.i, prob.i,
                         truncate = setdiff(allx.i, inf.mix),
                         max.support = max(inf.mix))
  }




  sum.i <- 0
  if (linf.mlm) {
    pstr.mlm <-  matrix(pstr.mlm, LLL, linf.mlm,
                          byrow = byrow.ai)
    sum.i <- .rowSums(pstr.mlm, LLL, linf.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")
  }  # linf.mlm

  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix - pstr.mix
  if (linf.mlm) {
    if (deflation) {
      tmp0 <- cdf.max.s - suma - sumt
      for (jay in 1:linf.mlm) {
        vecTF <- is.finite(x) & inf.mlm[jay] == x
        pmf.i <- dbinom(inf.mlm[jay], size.p[vecTF], prob.p[vecTF])
        if (any(pstr.mlm[vecTF, jay] <
                -(tmp6[vecTF] + pstr.mlm[vecTF, jay]) * pmf.i / (
                  tmp0[vecTF] - pmf.i), na.rm = TRUE)) {
          warning("too much deflation in argument 'pstr.mlm'. ",
                  "Returning NA")
          tmp6[vecTF] <- NA
        }
      }  # for
    } else {
      if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
        warning("the sum of variables 'sum.a', 'sum.i', 'pobs.mix' ",
                "and 'pstr.mix' exceeds unity. Returning NA")
        tmp6[!skip & tmp6 < 0] <- NA
      }
    }  # deflation
  }  # linf.mlm


  pmf0[!skip] <- (tmp6 *
    dbinom(x, size.p, prob.p) / (cdf.max.s - suma - sumt))[!skip]


  if (linf.mlm) {
    for (jay in seq(linf.mlm)) {
      ival <- inf.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
        pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # linf.mlm


  pmf0 <- pmf0 + pobs.mix * pmf2.a + pstr.mix * pmf2.i



  if (log.arg) log(pmf0) else pmf0
}  # dgaitbinom






 pgaitbinom <-
  function(q, size.p, prob.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL,
           max.support = NULL,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           lower.tail = TRUE) {
  if (!length(max.support))
    max.support <- max(size.p, size.a, size.i, na.rm = TRUE)  # Manually
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc   <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      max.support >= max(size.p, na.rm = TRUE))
    return(pbinom(q, size.p, prob.p, lower.tail = lower.tail))  # log.p


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0

  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")

  LLL <- max(length(q),        length(pobs.mix),   length(pstr.mix),
             length(size.p),   length(size.a),     length(size.i),
             length(prob.p),   length(prob.a),     length(prob.i))
  if (length(q)          < LLL) q          <- rep_len(q,          LLL)
  if (length(size.p)     < LLL) size.p     <- rep_len(size.p,     LLL)
  if (length(size.a)     < LLL) size.a     <- rep_len(size.a,     LLL)
  if (length(size.i)     < LLL) size.i     <- rep_len(size.i,     LLL)
  if (length(prob.p)     < LLL) prob.p     <- rep_len(prob.p,     LLL)
  if (length(prob.a)     < LLL) prob.a     <- rep_len(prob.a,     LLL)
  if (length(prob.i)     < LLL) prob.i     <- rep_len(prob.i,     LLL)
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)

  offset.a <- offset.i <- Offset.a <- Offset.i <- numeric(LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- pbinom(max.support, size.p, prob.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dbinom(tval, size.p, prob.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (lalt.mlm) {
    pobs.mlm <- matrix(pobs.mlm, LLL, lalt.mlm, byrow = byrow.ai)
    sum.a <- .rowSums(pobs.mlm, LLL, lalt.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")

    for (jay in seq(lalt.mlm)) {
      aval <- alt.mlm[jay]
      pmf.p <- dbinom(aval, size.p, prob.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalt.mlm

  sum.i <- 0
  if (linf.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, linf.mlm, byrow = byrow.ai)
    sum.i <- .rowSums(pstr.mlm, LLL, linf.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")

    for (jay in seq(linf.mlm)) {
      ival <- inf.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # linf.mlm



  use.pobs.mix <- 0
  if (lalt.mix) {
    use.pobs.mix <- matrix(0, LLL, lalt.mix)
    for (jay in seq(lalt.mix)) {
      aval <- alt.mix[jay]
      pmf.a <- dbinom(aval, size.a, prob.a)
      pmf.p <- dbinom(aval, size.p, prob.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                      use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(lalt.mix)) {
      aval <- alt.mix[jay]
      pmf.p <- dbinom(aval, size.p, prob.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalt.mix

  use.pstr.mix <- 0
  if (linf.mix) {
    use.pstr.mix <- matrix(0, LLL, linf.mix)
    for (jay in seq(linf.mix)) {
      ival <- inf.mix[jay]
      use.pstr.mix[, jay] <- dbinom(ival, size.i, prob.i)
    }
    use.pstr.mix <- pstr.mix *
                      use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(linf.mix)) {
      ival <- inf.mix[jay]
      pmf.p <- dbinom(ival, size.p, prob.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix[vecTF, jay]
      }
    }  # jay
  }  # linf.mix

  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * (pbinom(q, size.p, prob.p) - fudge.t -
                   fudge.a) / denom1 +
         offset.i + offset.a + Offset.i + Offset.a
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitbinom






 qgaitbinom <-
  function(p, size.p, prob.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL,
           max.support = NULL,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p) {
  if (!length(max.support))
    max.support <- max(size.p, size.a, size.i, na.rm = TRUE)  # Manually
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc   <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      max.support >= max(size.p, na.rm = TRUE))
    return(qbinom(p, size.p, prob.p ))  # lower.tail, log.p = FALSE


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")


  LLL <- max(length(p),        length(pobs.mix),   length(pstr.mix),
             length(size.p),   length(size.a),     length(size.i),
             length(prob.p),   length(prob.a),     length(prob.i))
  if (length(p)          < LLL) p          <- rep_len(p,          LLL)
  if (length(size.p)     < LLL) size.p     <- rep_len(size.p,     LLL)
  if (length(size.a)     < LLL) size.a     <- rep_len(size.a,     LLL)
  if (length(size.i)     < LLL) size.i     <- rep_len(size.i,     LLL)
  if (length(prob.p)     < LLL) prob.p     <- rep_len(prob.p,     LLL)
  if (length(prob.a)     < LLL) prob.a     <- rep_len(prob.a,     LLL)
  if (length(prob.i)     < LLL) prob.i     <- rep_len(prob.i,     LLL)
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)

  pobs.mlm <- matrix(pobs.mlm, LLL, max(lalt.mlm, 1),
                       byrow = byrow.ai)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(linf.mlm, 1),
                       byrow = byrow.ai)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + size.p + size.a + size.i +
             prob.p + prob.a + prob.i

  bad0.p <- !is.finite(size.p) | size.p <= 0 |
            !is.finite(prob.p) | prob.p <= 0 | 1 <= prob.p
  bad0.a <- !is.finite(size.a) | size.a <= 0 |
            !is.finite(prob.a) | prob.a <= 0 | 1 <= prob.a
  bad0.i <- !is.finite(size.i) | size.i <= 0 |
            !is.finite(prob.i) | prob.i <= 0 | 1 <= prob.i
  bad0 <- bad0.p | bad0.a | bad0.i
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- rep_len(max.support + 0.5, LLL)
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitbinom(Hi, size.p, prob.p,
                    alt.mix = alt.mix, alt.mlm = alt.mlm,
                    inf.mix = inf.mix, inf.mlm = inf.mlm,
                    truncate = truncate, max.support = max.support,
                    pstr.mix = pstr.mix, pobs.mix = pobs.mix,
                    pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
                    size.a = size.p, size.i = size.p,
                    prob.a = prob.p, prob.i = prob.p,
                    byrow.ai = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitbinom(Hi[!done],
                       size.p[!done], prob.p[!done],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       size.a = size.a[!done],
                       size.i = size.i[!done],
                       prob.a = prob.a[!done],
                       prob.i = prob.i[!done],
                       byrow.ai = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, size.p, prob.p,
                      alt.mix = NULL, alt.mlm = NULL,
                      inf.mix = NULL, inf.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix = 0, pstr.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0,
                      size.a = size.p, size.i = size.p,
                      prob.a = prob.p, prob.i = prob.p,
                      byrow.ai = FALSE, p)
      pgaitbinom(q, size.p = size.p, prob.p = prob.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix, inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       size.a = size.a, prob.a = prob.a,
                       size.i = size.i, prob.i = prob.i,
                       byrow.ai = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitbinom(min.support.use,
                        size.p = size.p, prob.p = prob.p,
                        alt.mix = alt.mix, alt.mlm = alt.mlm,
                        inf.mix = inf.mix, inf.mlm = inf.mlm,
                        truncate = truncate, max.support = max.support,
                        pobs.mix = pobs.mix,
                        pstr.mix = pstr.mix,
                        pobs.mlm = pobs.mlm,
                        pstr.mlm = pstr.mlm,
                        size.a = size.a, prob.a = prob.a,
                        size.i = size.i, prob.i = prob.i,
                        byrow.ai = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      size.p = size.p[!lhs],
                      prob.p = prob.p[!lhs],
                      alt.mix = alt.mix, alt.mlm = alt.mlm,
                      inf.mix = inf.mix,
                      inf.mlm = inf.mlm,
                      truncate = truncate, max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      size.a = size.a[!lhs],
                      prob.a = prob.a[!lhs],
                      size.i = size.i[!lhs],
                      prob.i = prob.i[!lhs],
                      byrow.ai = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitbinom(faa,
                       size.p = size.p[!lhs],
                       prob.p = prob.p[!lhs],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs],
                       prob.a = prob.a[!lhs],
                       size.i = size.i[!lhs],
                       prob.i = prob.i[!lhs],
                       byrow.ai = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitbinom(faa + 1,
                       size.p = size.p[!lhs],
                       prob.p = prob.p[!lhs],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix, inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs],
                       prob.a = prob.a[!lhs],
                       size.i = size.i[!lhs],
                       prob.i = prob.i[!lhs],
                       byrow.ai = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitbinom(min.support.use, size.p, prob.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix, inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       size.a = size.a, size.i = size.i,
                       prob.a = prob.a, prob.i = prob.i,
                       byrow.ai = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitbinom






 rgaitbinom <-
  function(n, size.p, prob.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = NULL,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p) {
    qgaitbinom(runif(n), size.p, prob.p,
              alt.mix = alt.mix,
              alt.mlm = alt.mlm,
              inf.mix = inf.mix,
              inf.mlm = inf.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              size.a = size.a, size.i = size.i,
              prob.a = prob.a, prob.i = prob.i,
              byrow.ai = byrow.ai)
}  # rgaitbinom









 gait.errorcheck <-
  function(alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL,
           truncate = NULL,
           max.support = Inf,
           min.support = 0) {
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc <- length(truncate)

  if (!is.numeric(max.support) || is.na(max.support) ||
      length(max.support) != 1 || max.support < min.support ||
      round(max.support) != max.support ||
      (length(truncate) && (
          min(truncate, na.rm = TRUE) < min.support ||
          max.support <= max(truncate, na.rm = TRUE))))
    stop("bad input for argument 'max.support' and/or ",
         "'truncate'")

  allargs <- c(alt.mix, alt.mlm, inf.mix, inf.mlm)
  allargs <- c(allargs, truncate)  # No NA, NaN, -Inf or Inf allowed
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm)
    if (!is.Numeric(allargs, integer.valued = TRUE) ||
        any(allargs < min.support) ||
        any(max.support < allargs))
      stop("bad input for arguments 'alt.mix', 'alt.mlm',",
           " 'inf.mix' and/or 'inf.mlm'")
  if (length(unique(allargs)) < lalt.mix + lalt.mlm +
                                linf.mix + linf.mlm + ltrunc)
      stop("duplicate values found in arguments 'alt.mix', ",
           "'alt.mlm', 'inf.mix', 'inf.mlm' and 'truncate'")


}  # gait.errorcheck









 dgaitnbinom <-
  function(x, size.p, prob.p = NULL, munb.p = NULL,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           byrow.ai = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p,
           deflation = FALSE,  # Single logical
           log = FALSE) {
  log.arg <- log;  rm(log)
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)


  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")


  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc     <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(if (is.prob)
      dnbinom(x, size = size.p, prob = prob.p, log = log.arg) else
      dnbinom(x, size = size.p, mu   = munb.p, log = log.arg))


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (!deflation && any(pstr.mlm < 0, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")

  LLL <- max(length(x),        length(pobs.mix),   length(pstr.mix),
             length(size.p),   length(size.a),     length(size.i),
             length(munb.p),   length(munb.a),     length(munb.i),
             length(prob.p),   length(prob.a),     length(prob.i))
  if (length(x)          < LLL) x          <- rep_len(x,          LLL)
  if (length(size.p)     < LLL) size.p     <- rep_len(size.p,     LLL)
  if (length(size.a)     < LLL) size.a     <- rep_len(size.a,     LLL)
  if (length(size.i)     < LLL) size.i     <- rep_len(size.i,     LLL)
  if (is.prob) {
  if (length(prob.p)     < LLL) prob.p     <- rep_len(prob.p,     LLL)
  if (length(prob.a)     < LLL) prob.a     <- rep_len(prob.a,     LLL)
  if (length(prob.i)     < LLL) prob.i     <- rep_len(prob.i,     LLL)
  } else {
  if (length(munb.p)     < LLL) munb.p     <- rep_len(munb.p,     LLL)
  if (length(munb.a)     < LLL) munb.a     <- rep_len(munb.a,     LLL)
  if (length(munb.i)     < LLL) munb.i     <- rep_len(munb.i,     LLL)
  }
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc) {
    if (is.prob) {  # Need tval <= max.support
      for (tval in truncate)
        sumt <- sumt + dnbinom(tval, size.p, prob = prob.p)
    } else {
      for (tval in truncate)
        sumt <- sumt + dnbinom(tval, size.p, mu   = munb.p)
    }
  }
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- if (is.prob)
                pnbinom(max.support, size.p, prob = prob.p) else
                pnbinom(max.support, size.p, mu   = munb.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- if (is.prob)  # dgtnbinom
              ifelse(vecTF.t, 0, dnbinom(max.support, size.p,
                                         prob = prob.p) / denom.t) else
              ifelse(vecTF.t, 0, dnbinom(max.support, size.p,
                                         mu   = munb.p) / denom.t)


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalt.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, lalt.mlm,
                          byrow = byrow.ai)
    sum.a <- .rowSums(pobs.mlm, LLL, lalt.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    if (is.prob) {  # Part i
      for (aval in alt.mlm)
        suma <- suma + dnbinom(aval, size.p, prob = prob.p)
    } else {
      for (aval in alt.mlm)
        suma <- suma + dnbinom(aval, size.p, mu   = munb.p)
    }

      
    for (jay in seq(lalt.mlm)) {
      aval <- alt.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
          pmf0[vecTF] <- pobs.mlm[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # lalt.mlm



  pmf2.a <- pmf2.i <- 0
  if (lalt.mix) {
    allx.a <- lowsup:max(alt.mix)
    pmf2.a <- dgaitnbinom(x,  # Outer distribution---mlm type
                          size.p = size.a,
                          prob.p = prob.a,
                          munb.p = munb.a,
                          truncate = setdiff(allx.a, alt.mix),
                          max.support = max(alt.mix))
    for (aval in alt.mix) {
      suma <- suma + (if (is.prob)  # Part ii added; cumulative
                      dnbinom(aval, size = size.p, prob = prob.p) else
                      dnbinom(aval, size = size.p, mu   = munb.p))
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }


  if (linf.mix) {
    allx.i <- if (length(inf.mix)) lowsup:max(inf.mix) else NULL
    pmf2.i <- dgaitnbinom(x,  # Outer distribution---mlm type
                          size.p = size.i,
                          prob.p = prob.i,
                          munb.p = munb.i,
                          truncate = setdiff(allx.i, inf.mix),
                          max.support = max(inf.mix))
  }




  sum.i <- 0
  if (linf.mlm) {
    pstr.mlm <-  matrix(pstr.mlm, LLL, linf.mlm,
                          byrow = byrow.ai)
    sum.i <- .rowSums(pstr.mlm, LLL, linf.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")
  }  # linf.mlm

  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix - pstr.mix
  if (linf.mlm) {
    if (deflation) {
      tmp0 <- cdf.max.s - suma - sumt
      for (jay in 1:linf.mlm) {
        vecTF <- is.finite(x) & inf.mlm[jay] == x
        pmf.i <- if (is.prob)
    dnbinom(inf.mlm[jay], size.p[vecTF], prob = prob.p[vecTF]) else
    dnbinom(inf.mlm[jay], size.p[vecTF], mu   = munb.p[vecTF])
        if (any(pstr.mlm[vecTF, jay] <
                -(tmp6[vecTF] + pstr.mlm[vecTF, jay]) * pmf.i / (
                  tmp0[vecTF] - pmf.i), na.rm = TRUE)) {
          warning("too much deflation in argument 'pstr.mlm'. ",
                  "Returning NA")
          tmp6[vecTF] <- NA
        }
      }  # for
    } else {
      if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
        warning("the sum of variables 'sum.a', 'sum.i', 'pobs.mix' ",
                "and 'pstr.mix' exceeds unity. Returning NA")
        tmp6[!skip & tmp6 < 0] <- NA
      }
    }  # deflation
  }  # linf.mlm


  pmf0[!skip] <- (tmp6 * (if (is.prob) 
    dnbinom(x, size.p, prob = prob.p) else
    dnbinom(x, size.p, mu   = munb.p)) / (
    cdf.max.s - suma - sumt))[!skip]  # added


  if (linf.mlm) {
    for (jay in seq(linf.mlm)) {
      ival <- inf.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
        pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # linf.mlm


  pmf0 <- pmf0 + pobs.mix * pmf2.a + pstr.mix * pmf2.i

  if (log.arg) log(pmf0) else pmf0
}  # dgaitnbinom







 pgaitnbinom <-
  function(q, size.p, prob.p = NULL, munb.p = NULL,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p,
           lower.tail = TRUE) {
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)

  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")

  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc     <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(if (is.prob)
           pnbinom(q, size = size.p, prob = prob.p,
                   lower.tail = lower.tail) else
           pnbinom(q, size = size.p, mu   = munb.p,
                   lower.tail = lower.tail))  # log.p

  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0

  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")

  LLL <- max(length(q),        length(pobs.mix),   length(pstr.mix),
             length(size.p),   length(size.a),     length(size.i),
             length(munb.p),   length(munb.a),     length(munb.i),
             length(prob.p),   length(prob.a),     length(prob.i))
  offset.a <- offset.i <- Offset.a <- Offset.i <- numeric(LLL)
  if (length(q)          < LLL) q          <- rep_len(q,          LLL)
  if (length(size.p)     < LLL) size.p     <- rep_len(size.p,     LLL)
  if (length(size.a)     < LLL) size.a     <- rep_len(size.a,     LLL)
  if (length(size.i)     < LLL) size.i     <- rep_len(size.i,     LLL)
  if (is.prob) {
  if (length(prob.p)     < LLL) prob.p     <- rep_len(prob.p,     LLL)
  if (length(prob.a)     < LLL) prob.a     <- rep_len(prob.a,     LLL)
  if (length(prob.i)     < LLL) prob.i     <- rep_len(prob.i,     LLL)
  } else {
  if (length(munb.p)     < LLL) munb.p     <- rep_len(munb.p,     LLL)
  if (length(munb.a)     < LLL) munb.a     <- rep_len(munb.a,     LLL)
  if (length(munb.i)     < LLL) munb.i     <- rep_len(munb.i,     LLL)
  }
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)



  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- if (is.prob)
                pnbinom(max.support, size.p, prob = prob.p) else
                pnbinom(max.support, size.p, mu   = munb.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- if (is.prob) dnbinom(tval, size.p, prob = prob.p) else
                            dnbinom(tval, size.p, mu   = munb.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (lalt.mlm) {
    pobs.mlm <- matrix(pobs.mlm, LLL, lalt.mlm, byrow = byrow.ai)
    sum.a <- .rowSums(pobs.mlm, LLL, lalt.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")

    for (jay in seq(lalt.mlm)) {
      aval <- alt.mlm[jay]
      pmf.p <- if (is.prob) dnbinom(aval, size.p, prob = prob.p) else
                            dnbinom(aval, size.p, mu   = munb.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalt.mlm

  sum.i <- 0
  if (linf.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, linf.mlm, byrow = byrow.ai)
    sum.i <- .rowSums(pstr.mlm, LLL, linf.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")

    for (jay in seq(linf.mlm)) {
      ival <- inf.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # linf.mlm



  use.pobs.mix <- 0
  if (lalt.mix) {
    use.pobs.mix <- matrix(0, LLL, lalt.mix)
    for (jay in seq(lalt.mix)) {
      aval <- alt.mix[jay]
      pmf.a <- if (is.prob) dnbinom(aval, size.a, prob = prob.a) else
                            dnbinom(aval, size.a, mu   = munb.a)
      pmf.p <- if (is.prob) dnbinom(aval, size.p, prob = prob.p) else
                            dnbinom(aval, size.p, mu   = munb.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                      use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(lalt.mix)) {
      aval <- alt.mix[jay]
      pmf.p <- if (is.prob) dnbinom(aval, size.p, prob = prob.p) else
                            dnbinom(aval, size.p, mu   = munb.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalt.mix

  use.pstr.mix <- 0
  if (linf.mix) {
    use.pstr.mix <- matrix(0, LLL, linf.mix)
    for (jay in seq(linf.mix)) {
      ival <- inf.mix[jay]
      use.pstr.mix[, jay] <- if (is.prob)
        dnbinom(ival, size.i, prob = prob.i) else
        dnbinom(ival, size.i, mu   = munb.i)
    }
    use.pstr.mix <- pstr.mix *
                      use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(linf.mix)) {
      ival <- inf.mix[jay]
      pmf.p <- if (is.prob) dnbinom(ival, size.p, prob = prob.p) else
                            dnbinom(ival, size.p, mu   = munb.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix[vecTF, jay]
      }
    }  # jay
  }  # linf.mix

  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * ((if (is.prob)
         pnbinom(q, size.p, prob = prob.p) else
         pnbinom(q, size.p, mu   = munb.p)) -
         fudge.t - fudge.a) / denom1 +
         offset.i + offset.a + Offset.i + Offset.a
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitnbinom








 qgaitnbinom <-
  function(p, size.p, prob.p = NULL, munb.p = NULL,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p) {
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)

  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")

  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc     <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(if (is.prob)
      qnbinom(p, size = size.p, prob = prob.p) else
      qnbinom(p, size = size.p, mu   = munb.p))  # lower.tail, log.p

  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")


  LLL <- max(length(p),        length(pobs.mix),   length(pstr.mix),
             length(size.p),   length(size.a),     length(size.i),
             length(munb.p),   length(munb.a),     length(munb.i),
             length(prob.p),   length(prob.a),     length(prob.i))
  if (length(p)          < LLL) p          <- rep_len(p,          LLL)
  if (length(size.p)     < LLL) size.p     <- rep_len(size.p,     LLL)
  if (length(size.a)     < LLL) size.a     <- rep_len(size.a,     LLL)
  if (length(size.i)     < LLL) size.i     <- rep_len(size.i,     LLL)
  if (is.prob) {
  if (length(prob.p)     < LLL) prob.p     <- rep_len(prob.p,     LLL)
  if (length(prob.a)     < LLL) prob.a     <- rep_len(prob.a,     LLL)
  if (length(prob.i)     < LLL) prob.i     <- rep_len(prob.i,     LLL)
  } else {
  if (length(munb.p)     < LLL) munb.p     <- rep_len(munb.p,     LLL)
  if (length(munb.a)     < LLL) munb.a     <- rep_len(munb.a,     LLL)
  if (length(munb.i)     < LLL) munb.i     <- rep_len(munb.i,     LLL)
  }
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)


  pobs.mlm <- matrix(pobs.mlm, LLL, max(lalt.mlm, 1),
                       byrow = byrow.ai)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(linf.mlm, 1),
                       byrow = byrow.ai)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + size.p + size.a + size.i + (if (is.prob)
             prob.p + prob.a + prob.i else 
             munb.p + munb.a + munb.i)

  bad0.p <- !is.finite(size.p) | size.p <= 0 |
            (if (is.prob)
              !is.finite(prob.p) | prob.p <= 0 | 1 <= prob.p else
              !is.finite(munb.p) | munb.p <= 0)
  bad0.a <- !is.finite(size.a) | size.a <= 0 |
            (if (is.prob)
              !is.finite(prob.a) | prob.a <= 0 | 1 <= prob.a else
              !is.finite(munb.a) | munb.a <= 0)
  bad0.i <- !is.finite(size.i) | size.i <= 0 |
            (if (is.prob)
              !is.finite(prob.i) | prob.i <= 0 | 1 <= prob.i else
              !is.finite(munb.i) | munb.i <= 0)
  bad0 <- bad0.p | bad0.a | bad0.i
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitnbinom(Hi, size.p, prob.p = prob.p, munb.p = munb.p,
               alt.mix = alt.mix, alt.mlm = alt.mlm,
               inf.mix = inf.mix, inf.mlm = inf.mlm,
               truncate = truncate, max.support = max.support,
               pstr.mix = pstr.mix, pobs.mix = pobs.mix,
               pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
               prob.a = prob.a, prob.i = prob.i,
               munb.a = munb.a, munb.i = munb.i,
               byrow.ai = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitnbinom(Hi[!done], size.p = size.p[!done],
                       prob.p = prob.p[!done], munb.p = munb.p[!done],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       size.a = size.a[!done], size.i = size.i[!done],
                       prob.a = prob.a[!done], prob.i = prob.i[!done],
                       munb.a = munb.a[!done], munb.i = munb.i[!done],
                       byrow.ai = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, size.p, prob.p = NULL, munb.p = NULL,
                      alt.mix = NULL, alt.mlm = NULL,
                      inf.mix = NULL, inf.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix = 0, pstr.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0,
                      size.a = size.p, size.i = size.p,
                      prob.a = prob.p, prob.i = prob.p,
                      munb.a = munb.p, munb.i = munb.p,
                      byrow.ai = FALSE, p)
      pgaitnbinom(q, size.p = size.p, prob.p = prob.p, munb.p = munb.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       size.a = size.a, size.i = size.i,
                       prob.a = prob.a, prob.i = prob.i,
                       munb.a = munb.a, munb.i = munb.i,
                       byrow.ai = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitnbinom(min.support.use, size.p = size.p,
                       prob.p = prob.p, munb.p = munb.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       size.a = size.a, size.i = size.i,
                       prob.a = prob.a, prob.i = prob.i,
                       munb.a = munb.a, munb.i = munb.i,
                       byrow.ai = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      size.p = size.p[!lhs],
                      prob.p = prob.p[!lhs], munb.p = munb.p[!lhs],
                      alt.mix = alt.mix, alt.mlm = alt.mlm,
                      inf.mix = inf.mix,
                      inf.mlm = inf.mlm,
                      truncate = truncate, max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      size.a = size.a[!lhs], size.i = size.i[!lhs],
                      prob.a = prob.a[!lhs], prob.i = prob.i[!lhs],
                      munb.a = munb.a[!lhs], munb.i = munb.i[!lhs],
                      byrow.ai = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitnbinom(faa, size.p = size.p[!lhs],
                       prob.p = prob.p[!lhs], munb.p = munb.p[!lhs],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs], size.i = size.i[!lhs],
                       prob.a = prob.a[!lhs], prob.i = prob.i[!lhs],
                       munb.a = munb.a[!lhs], munb.i = munb.i[!lhs],
                       byrow.ai = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitnbinom(faa + 1, size.p = size.p[!lhs],
                       prob.p = prob.p[!lhs], munb.p = munb.p[!lhs],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs], size.i = size.i[!lhs],
                       prob.a = prob.a[!lhs], prob.i = prob.i[!lhs],
                       munb.a = munb.a[!lhs], munb.i = munb.i[!lhs],
                       byrow.ai = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitnbinom(min.support.use, size.p = size.p,
                       prob.p = prob.p, munb.p = munb.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       size.a = size.a, size.i = size.i,
                       prob.a = prob.a, prob.i = prob.i,
                       munb.a = munb.a, munb.i = munb.i,
                       byrow.ai = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitnbinom





 rgaitnbinom <-
  function(n, size.p, prob.p = NULL, munb.p = NULL,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           byrow.ai = FALSE,
           size.a = size.p, size.i = size.p,
           prob.a = prob.p, prob.i = prob.p,
           munb.a = munb.p, munb.i = munb.p) {
    qgaitnbinom(runif(n),
              size.p = size.p, prob.p = prob.p, munb.p = munb.p,
              alt.mix = alt.mix,
              alt.mlm = alt.mlm,
              inf.mix = inf.mix,
              inf.mlm = inf.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              size.a = size.a, prob.a = prob.a, munb.a = munb.a,
              size.i = size.i, prob.i = prob.i, munb.i = munb.i,
              byrow.ai = byrow.ai)
}  # rgaitnbinom













 moments.gaitcombo <-
  function(theta.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # Vector
           pobs.mlm = 0,  # Matrix
           pstr.mix = 0,
           pstr.mlm = 0,  # Ditto
           byrow.ai = FALSE,  # For pobs.mlm and pstr.mlm
           theta.a = theta.p, theta.i = theta.p,
           moments2 = FALSE,  # Use this for variances.
           rmlife1 = 0, rmlife2 = 0,
           dfun = "dpois") {
      

  NOS <- 1
  nnn <- length(theta.p)
  pfun <- dfun
  substring(pfun, 1) <- "p"  # Replace the "d" by a "p"
  cdf.max.s <- do.call(pfun, list(max.support, theta.p))
  LALT.MIX <- length(alt.mix)
  LALT.MLM <- length(alt.mlm)
  LINF.MIX <- length(inf.mix)
  LINF.MLM <- length(inf.mlm)
  LTRUNCAT <- length(truncate)

  if (LALT.MLM == 0) {
    if (!all(pobs.mlm == 0))
      warning("ignoring argument 'pobs.mlm'")
    pobs.mlm <- 0
  }
  if (LINF.MLM == 0) {
    if (!all(pstr.mlm == 0))
      warning("ignoring argument 'pstr.mlm'")
    pstr.mlm <- 0
  }
  if (LALT.MIX == 0) {
    if (!all(pobs.mix == 0))
      warning("ignoring argument 'pobs.mix'")
    pobs.mix <- 0
  }
  if (LINF.MIX == 0) {
    if (!all(pstr.mix == 0))
      warning("ignoring argument 'pstr.mix'")
    pstr.mix <- 0
  }

  SumT0.p <- matrix(0, nnn, NOS)  # Does not include upper RHS tail
  SumT1.p <- matrix(rmlife1, nnn, NOS)  # Includes RHS tail
  SumT2.p <- matrix(rmlife2, nnn, NOS)  # Includes RHS tail
  if (LTRUNCAT)
    for (tval in truncate) {
      pmf.p <- do.call(dfun, list(tval, theta.p))
      SumT0.p <- SumT0.p + pmf.p  # Need tval<=max.support
      SumT1.p <- SumT1.p + pmf.p * tval
      if (moments2)
        SumT2.p <- SumT2.p + pmf.p * tval^2
    }

  use.pobs.mix <- use.pobs.mlm <-  # So rowSums() works below.
  use.pstr.mix <- use.pstr.mlm <- matrix(0, nnn, 1)
  aprd1.mix <- aprd1.mlm <-    # aprd1.m?? is an innerprod
  aprd2.mix <- aprd2.mlm <- 0  # aprd2.m?? is an innerprod
  SumA0.mix.p <- SumA0.mlm.p <-
  SumA0.mix.a <- SumA0.mlm.a <-
  SumA1.mix.p <- SumA1.mlm.p <-
  SumA1.mix.a <- SumA1.mlm.a <-
  SumA2.mix.p <- SumA2.mlm.p <-
  SumA2.mix.a <- SumA2.mlm.a <-
  SumA0.mix.x <- SumA0.mlm.x <-
  SumA1.mix.x <- SumA1.mlm.x <-
  SumA2.mix.x <- SumA2.mlm.x <- matrix(0, nnn, NOS)
  if (LALT.MLM)
    use.pobs.mlm <- matrix(pobs.mlm, nnn, LALT.MLM, byrow = byrow.ai)
  if (LINF.MLM)
    use.pstr.mlm <- matrix(pstr.mlm, nnn, LINF.MLM, byrow = byrow.ai)
  if (LALT.MIX)
    use.pobs.mix <- matrix(pobs.mix, nnn, 1)
  if (LINF.MIX)
    use.pstr.mix <- matrix(pstr.mix, nnn, 1)


  if (LALT.MIX) {
    for (jay in seq_len(LALT.MIX)) {
      aval <- alt.mix[jay]
      pmf.p <- do.call(dfun, list(aval, theta.p))
      pmf.a <- do.call(dfun, list(aval, theta.a))
      SumA0.mix.p <- SumA0.mix.p + pmf.p
      SumA0.mix.a <- SumA0.mix.a + pmf.a
      SumA1.mix.p <- SumA1.mix.p + pmf.p * aval
      SumA1.mix.a <- SumA1.mix.a + pmf.a * aval
      if (moments2) {
        SumA2.mix.p <- SumA2.mix.p + pmf.p * aval^2
        SumA2.mix.a <- SumA2.mix.a + pmf.a * aval^2
      }
    }  # for jay
    aprd1.mix <- use.pobs.mix * SumA1.mix.a / SumA0.mix.a
    if (moments2)
      aprd2.mix <- use.pobs.mix * SumA2.mix.a / SumA0.mix.a
  }  # LALT.MIX


  if (LALT.MLM) {
    for (jay in seq_len(LALT.MLM)) {
      aval <- alt.mlm[jay]
      pmf.x <- use.pobs.mlm[, jay]
      pmf.p <- do.call(dfun, list(aval, theta.p))
      pmf.a <- do.call(dfun, list(aval, theta.a))
      SumA0.mlm.p <- SumA0.mlm.p + pmf.p
      SumA0.mlm.a <- SumA0.mlm.a + pmf.a
      SumA1.mlm.p <- SumA1.mlm.p + pmf.p * aval
      SumA1.mlm.a <- SumA1.mlm.a + pmf.a * aval
      SumA1.mlm.x <- SumA1.mlm.x + pmf.x * aval
      if (moments2) {
        SumA2.mlm.p <- SumA2.mlm.p + pmf.p * aval^2
        SumA2.mlm.a <- SumA2.mlm.a + pmf.a * aval^2
        SumA2.mlm.x <- SumA2.mlm.x + pmf.x * aval^2
      }
    }  # for jay
    aprd1.mlm <- SumA1.mlm.x
    if (moments2)
      aprd2.mlm <- SumA2.mlm.x
  }  # LALT.MLM


  iprd1.mix <- iprd1.mlm <-    # iprd1.m?? is an innerprod
  iprd2.mix <- iprd2.mlm <- 0  # iprd2.m?? is an innerprod
  SumI0.mix.p <- SumI0.mlm.p <-
  SumI0.mix.i <- SumI0.mlm.i <-
  SumI1.mix.p <- SumI1.mlm.p <-
  SumI1.mix.i <- SumI1.mlm.i <-
  SumI2.mix.p <- SumI2.mlm.p <-
  SumI2.mix.i <- SumI2.mlm.i <-
  SumI0.mix.x <- SumI0.mlm.x <-
  SumI1.mix.x <- SumI1.mlm.x <-
  SumI2.mix.x <- SumI2.mlm.x <- matrix(0, nnn, NOS)
  if (LINF.MIX) {
    for (jay in seq_len(LINF.MIX)) {
      ival <- inf.mix[jay]
      pmf.p <- do.call(dfun, list(ival, theta.p))
      pmf.i <- do.call(dfun, list(ival, theta.i))
      SumI0.mix.p <- SumI0.mix.p + pmf.p
      SumI0.mix.i <- SumI0.mix.i + pmf.i
      SumI1.mix.p <- SumI1.mix.p + pmf.p * ival
      SumI1.mix.i <- SumI1.mix.i + pmf.i * ival
      if (moments2) {
        SumI2.mix.p <- SumI2.mix.p + pmf.p * ival^2
        SumI2.mix.i <- SumI2.mix.i + pmf.i * ival^2
      }
    }  # for jay
    iprd1.mix <- use.pstr.mix * SumI1.mix.i / SumI0.mix.i
    if (moments2)
      iprd2.mix <- use.pstr.mix * SumI2.mix.i / SumI0.mix.i
  }  # LINF.MIX


  if (LINF.MLM) {
    for (jay in seq_len(LINF.MLM)) {
      ival <- inf.mlm[jay]
      pmf.x <- use.pstr.mlm[, jay]
      pmf.p <- do.call(dfun, list(ival, theta.p))
      pmf.i <- do.call(dfun, list(ival, theta.i))
      SumI0.mlm.p <- SumI0.mlm.p + pmf.p
      SumI0.mlm.i <- SumI0.mlm.i + pmf.i
      SumI1.mlm.p <- SumI1.mlm.p + pmf.p * ival
      SumI1.mlm.i <- SumI1.mlm.i + pmf.i * ival
      SumI1.mlm.x <- SumI1.mlm.x + pmf.x * ival
      if (moments2) {
        SumI2.mlm.p <- SumI2.mlm.p + pmf.p * ival^2
        SumI2.mlm.i <- SumI2.mlm.i + pmf.i * ival^2
        SumI2.mlm.x <- SumI2.mlm.x + pmf.x * ival^2
      }
    }  # for jay
    iprd1.mlm <- SumI1.mlm.x
    if (moments2)
      iprd2.mlm <- SumI2.mlm.x
  }  # LINF.MLM


  use.this <- 1 - rowSums(use.pobs.mlm) - rowSums(use.pstr.mlm) -
              use.pobs.mix - use.pstr.mix
  ans <- list('cdf.max.s'   = cdf.max.s,
              'SumT0.p'     = SumT0.p,
              'SumT1.p'     = SumT1.p,
              'SumA0.mix.a' = SumA0.mix.a,
              'SumA0.mix.p' = SumA0.mix.p,
              'SumA1.mix.a' = SumA1.mix.a,
              'SumA1.mix.p' = SumA1.mix.p,
              'SumA0.mlm.a' = SumA0.mlm.a,
              'SumA0.mlm.p' = SumA0.mlm.p,
              'SumA1.mlm.a' = SumA1.mlm.a,
              'SumA1.mlm.p' = SumA1.mlm.p,
              'SumI0.mix.i' = SumI0.mix.i,
              'SumI0.mix.p' = SumI0.mix.p,
              'SumI1.mix.i' = SumI1.mix.i,
              'SumI1.mix.p' = SumI1.mix.p,
              'SumI0.mlm.i' = SumI0.mlm.i,
              'SumI0.mlm.p' = SumI0.mlm.p,
              'SumI1.mlm.i' = SumI1.mlm.i,
              'SumI1.mlm.p' = SumI1.mlm.p,
              'aprd1.mix'   = aprd1.mix,
              'aprd1.mlm'   = aprd1.mlm,
              'iprd1.mix'   = iprd1.mix,
              'iprd1.mlm'   = iprd1.mlm,
              'use.this'    = use.this)

  if (moments2) {  # Add more info
  ans <- c(ans,
         list(  # 'rmlife2'     = rmlife2,   # May be scalar
              'aprd2.mix'   = aprd2.mix,
              'aprd2.mlm'   = aprd2.mlm,
              'iprd2.mix'   = iprd2.mix,
              'iprd2.mlm'   = iprd2.mlm,
              'SumT2.p'     = SumT2.p,
              'SumA2.mix.p' = SumA2.mix.p,
              'SumA2.mix.a' = SumA2.mix.a,
              'SumI2.mix.p' = SumI2.mix.p,
              'SumI2.mix.i' = SumI2.mix.i,
              'SumA2.mlm.p' = SumA2.mlm.p,
              'SumA2.mlm.a' = SumA2.mlm.a,
              'SumI2.mlm.p' = SumI2.mlm.p,
              'SumI2.mlm.i' = SumI2.mlm.i))
  }
  ans
}  # moments.gaitcombo 






 moments.gaitcombo.pois <-
  function(lambda.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # Vector
           pobs.mlm = 0,  # Matrix
           pstr.mix = 0,
           pstr.mlm = 0,  # Ditto
           byrow.ai = FALSE,  # For pobs.mlm and pstr.mlm
           lambda.a = lambda.p, lambda.i = lambda.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  rmlife1 <- ppois(max.support - 1, lambda.p, lower.tail = FALSE) *
             lambda.p
  rmlife2 <- ppois(max.support - 2, lambda.p, lower.tail = FALSE) * 
             lambda.p^2 + rmlife1

  mylist1 <- moments.gaitcombo(theta.p = lambda.p,
           alt.mix = alt.mix, alt.mlm = alt.mlm,
           inf.mix = inf.mix, inf.mlm = inf.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           byrow.ai = byrow.ai,  # type.fitted = type.fitted,
           theta.a = lambda.a, theta.i = lambda.i,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dpois")


  themean <- with(mylist1,
                  aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm +
                  use.this * (lambda.p - SumA1.mix.p -
                              SumA1.mlm.p - SumT1.p) / (
                  cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))

  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitcombo.pois






 moments.gaitcombo.log <-
  function(shape.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0, pobs.mlm = 0,  # Vector and matrix resp.
           pstr.mix = 0, pstr.mlm = 0,  # Ditto
           byrow.ai = FALSE,  # For pobs.mlm and pstr.mlm
           shape.a = shape.p, shape.i = shape.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  A8.p <- -1 / log1p(-shape.p)
  rmlife1 <- A8.p * (shape.p^(max.support + 1)) / (1 - shape.p)
  rmlife2 <- A8.p * ((shape.p^(max.support + 1)) *
                     (max.support + 1 / (1 - shape.p))
                     / (1 - shape.p))

  mylist1 <- moments.gaitcombo(theta.p = shape.p,
           alt.mix = alt.mix, alt.mlm = alt.mlm,
           inf.mix = inf.mix, inf.mlm = inf.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           byrow.ai = byrow.ai,  # type.fitted = type.fitted,
           theta.a = shape.a, theta.i = shape.i,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dlog")


  themean <- with(mylist1,
                  aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm +
                  use.this *
                  (-shape.p / (log1p(-shape.p) * (1 - shape.p)) -
                  SumA1.mix.p - SumA1.mlm.p - SumT1.p) / (
                  cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))

  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitcombo.log






 moments.gaitcombo.zeta <-
  function(shape.p,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0, pobs.mlm = 0,  # Vector and matrix resp.
           pstr.mix = 0, pstr.mlm = 0,  # Ditto
           byrow.ai = FALSE,  # For pobs.mlm and pstr.mlm
           shape.a = shape.p, shape.i = shape.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  rmlife1 <- if (is.finite(max.support)) zeta(shape.p) * (1 -
    pzeta(max.support, shape.p - 1)) / zeta(shape.p + 1) else
    numeric(length(shape.p))
  rmlife1[shape.p <= 1] <- NA  # NA or Inf, not sure
  rmlife2 <- if (is.finite(max.support)) zeta(shape.p - 1) * (1 -
    pzeta(max.support, shape.p - 2)) / zeta(shape.p + 1) else
    numeric(length(shape.p))
  rmlife2[shape.p <= 2] <- NA  # NA or Inf, not sure

  mylist1 <- moments.gaitcombo(theta.p = shape.p,
           alt.mix = alt.mix, alt.mlm = alt.mlm,
           inf.mix = inf.mix, inf.mlm = inf.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           byrow.ai = byrow.ai,  # type.fitted = type.fitted,
           theta.a = shape.a, theta.i = shape.i,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dzeta")

  themean <-
    with(mylist1,
         aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm + use.this *
         (ifelse(shape.p > 1, zeta(shape.p) / zeta(shape.p + 1), NA) -
          SumA1.mix.p - SumA1.mlm.p - SumT1.p) / (
          cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))
  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitcombo.zeta







 dgaitpois <-
  function(x, lambda.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           byrow.ai = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           lambda.a = lambda.p, lambda.i = lambda.p,
           deflation = FALSE,  # Single logical
           log = FALSE) {
  log.arg <- log;  rm(log)
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc     <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(dpois(x, lambda.p, log = log.arg))


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (!deflation && any(pstr.mlm < 0, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")

  LLL <- max(length(x),        length(pobs.mix),   length(pstr.mix),
             length(lambda.p), length(lambda.a),   length(lambda.i))
  if (length(x)          < LLL) x          <- rep_len(x,          LLL)
  if (length(lambda.p)   < LLL) lambda.p   <- rep_len(lambda.p,   LLL)
  if (length(lambda.a)   < LLL) lambda.a   <- rep_len(lambda.a,   LLL)
  if (length(lambda.i)   < LLL) lambda.i   <- rep_len(lambda.i,   LLL)
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dpois(tval, lambda.p)  # Need tval <= max.support
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- ppois(max.support, lambda.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dpois(x, lambda.p) / denom.t)  # dgtpois


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (lalt.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, lalt.mlm,
                          byrow = byrow.ai)
    sum.a <- .rowSums(pobs.mlm, LLL, lalt.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    for (aval in alt.mlm)
      suma <- suma + dpois(aval, lambda.p)  # Part i

    for (jay in seq(lalt.mlm)) {
      aval <- alt.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
          pmf0[vecTF] <- pobs.mlm[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # lalt.mlm



  pmf2.a <- pmf2.i <- 0
  if (lalt.mix) {
    allx.a <- lowsup:max(alt.mix)
    pmf2.a <- dgaitpois(x, lambda.a,  # Outer distribution---mlm type
                        truncate = setdiff(allx.a, alt.mix),
                        max.support = max(alt.mix))
    for (aval in alt.mix) {
      suma <- suma + dpois(aval, lambda.p)  # Part ii added; cumulative
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }

  if (linf.mix) {
    allx.i <- if (length(inf.mix)) lowsup:max(inf.mix) else NULL
    pmf2.i <- dgaitpois(x, lambda.i,  # Outer distribution---mlm type
                        truncate = setdiff(allx.i, inf.mix),
                        max.support = max(inf.mix))
  }




  sum.i <- 0
  if (linf.mlm) {
    pstr.mlm <-  matrix(pstr.mlm, LLL, linf.mlm,
                          byrow = byrow.ai)
    sum.i <- .rowSums(pstr.mlm, LLL, linf.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")
  }  # linf.mlm

  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix - pstr.mix
  if (linf.mlm) {
    if (deflation) {
      tmp0 <- cdf.max.s - suma - sumt
      for (jay in 1:linf.mlm) {
        vecTF <- is.finite(x) & inf.mlm[jay] == x
        pmf.i <- dpois(inf.mlm[jay], lambda.p[vecTF])
        if (any(pstr.mlm[vecTF, jay] <
                -(tmp6[vecTF] + pstr.mlm[vecTF, jay]) * pmf.i / (
                  tmp0[vecTF] - pmf.i), na.rm = TRUE)) {
          warning("too much deflation in argument 'pstr.mlm'. ",
                  "Returning NA")
          tmp6[vecTF] <- NA
        }
      }  # for
    } else {
      if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
        warning("the sum of variables 'sum.a', 'sum.i', 'pobs.mix' ",
                "and 'pstr.mix' exceeds unity. Returning NA")
        tmp6[!skip & tmp6 < 0] <- NA
      }
    }  # deflation
  }  # linf.mlm


  pmf0[!skip] <- (tmp6 *
    dpois(x, lambda.p) / (cdf.max.s - suma - sumt))[!skip]  # added


  if (linf.mlm) {
    for (jay in seq(linf.mlm)) {
      ival <- inf.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
          pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # linf.mlm


  pmf0 <- pmf0 + pobs.mix * pmf2.a + pstr.mix * pmf2.i




if (FALSE) {
  allx.a <- if (length(alt.mix))   lowsup:max(alt.mix  ) else NULL
  allx.i <- if (length(inf.mix)) lowsup:max(inf.mix) else NULL
  use.ms <- if (is.finite(max.support)) {
    whatsleft <- setdiff(lowsup:max.support, alt.mix)
    if (length(whatsleft) > 0)
      max(whatsleft) else max.support
  } else {
    max.support
  }
  pmf1 <- dgaitpois(x, lambda.p,  # Inner distribution
                    truncate = c(alt.mix[alt.mix < use.ms],
                                 alt.mlm[alt.mlm < use.ms],  # added
                                 truncate),  # No ivec
                    max.support = use.ms)
  if (linf.mix)
    pmf2.i <- dgaitpois(x, lambda.i,  # Outer distribution---mlm type
                        truncate = setdiff(allx.i, inf.mix),
                        max.support = max(inf.mix))
  pmf0.check <- pobs.mix * pmf2.a +
                pstr.mix * pmf2.i +  # mixprob * pmf2ai +
                (1 - sum.a - sum.i - pobs.mix - pstr.mix) * pmf1
}  # TRUE or FALSE

  if (log.arg) log(pmf0) else pmf0
}  # dgaitpois






 pgaitpois <-
  function(q, lambda.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           byrow.ai = FALSE,
           lambda.a = lambda.p, lambda.i = lambda.p,
           lower.tail = TRUE) {
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc     <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(ppois(q, lambda.p, lower.tail = lower.tail))  # log.p


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0

  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")

  LLL <- max(length(q),        length(pobs.mix), length(pstr.mix),
             length(lambda.p), length(lambda.a),   length(lambda.i))
  offset.a <- offset.i <- Offset.a <- Offset.i <- numeric(LLL)
  if (length(q)          < LLL) q          <- rep_len(q,          LLL)
  if (length(lambda.p)   < LLL) lambda.p   <- rep_len(lambda.p,   LLL)
  if (length(lambda.a)   < LLL) lambda.a   <- rep_len(lambda.a,   LLL)
  if (length(lambda.i)   < LLL) lambda.i   <- rep_len(lambda.i,   LLL)
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- ppois(max.support, lambda.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dpois(tval, lambda.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (lalt.mlm) {
    pobs.mlm <- matrix(pobs.mlm, LLL, lalt.mlm, byrow = byrow.ai)
    sum.a <- .rowSums(pobs.mlm, LLL, lalt.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")

    for (jay in seq(lalt.mlm)) {
      aval <- alt.mlm[jay]
      pmf.p <- dpois(aval, lambda.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalt.mlm

  sum.i <- 0
  if (linf.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, linf.mlm, byrow = byrow.ai)
    sum.i <- .rowSums(pstr.mlm, LLL, linf.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")

    for (jay in seq(linf.mlm)) {
      ival <- inf.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # linf.mlm



  use.pobs.mix <- 0
  if (lalt.mix) {
    use.pobs.mix <- matrix(0, LLL, lalt.mix)
    for (jay in seq(lalt.mix)) {
      aval <- alt.mix[jay]
      pmf.a <- dpois(aval, lambda.a)
      pmf.p <- dpois(aval, lambda.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                      use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(lalt.mix)) {
      aval <- alt.mix[jay]
      pmf.p <- dpois(aval, lambda.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # lalt.mix

  use.pstr.mix <- 0
  if (linf.mix) {
    use.pstr.mix <- matrix(0, LLL, linf.mix)
    for (jay in seq(linf.mix)) {
      ival <- inf.mix[jay]
      use.pstr.mix[, jay] <- dpois(ival, lambda.i)
    }
    use.pstr.mix <- pstr.mix *
                      use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(linf.mix)) {
      ival <- inf.mix[jay]
      pmf.p <- dpois(ival, lambda.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix[vecTF, jay]
      }
    }  # jay
  }  # linf.mix

  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * (ppois(q, lambda.p) - fudge.t - fudge.a) / denom1 +
         offset.i + offset.a + Offset.i + Offset.a
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitpois






 qgaitpois <-
  function(p, lambda.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           byrow.ai = FALSE,
           lambda.a = lambda.p, lambda.i = lambda.p) {
  lowsup <- 0
  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)
  ltrunc     <- length(truncate)
  if (lalt.mix + lalt.mlm + linf.mix + linf.mlm + ltrunc == 0 &&
      is.infinite(max.support))
    return(qpois(p, lambda.p))  # lower.tail = TRUE, log.p = FALSE


  if (lalt.mix == 0) pobs.mix <- 0
  if (lalt.mlm == 0) pobs.mlm <- 0
  if (linf.mix == 0) pstr.mix <- 0
  if (linf.mlm == 0) pstr.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")


  LLL <- max(length(p),        length(pobs.mix),   length(pstr.mix),
             length(lambda.p), length(lambda.a),   length(lambda.i))
  if (length(p)          < LLL) p          <- rep_len(p,          LLL)
  if (length(lambda.p)   < LLL) lambda.p   <- rep_len(lambda.p,   LLL)
  if (length(lambda.a)   < LLL) lambda.a   <- rep_len(lambda.a,   LLL)
  if (length(lambda.i)   < LLL) lambda.i   <- rep_len(lambda.i,   LLL)
  if (length(pobs.mix)   < LLL) pobs.mix   <- rep_len(pobs.mix,   LLL)
  if (length(pstr.mix)   < LLL) pstr.mix   <- rep_len(pstr.mix,   LLL)

  pobs.mlm <- matrix(pobs.mlm, LLL, max(lalt.mlm, 1), byrow = byrow.ai)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(linf.mlm, 1), byrow = byrow.ai)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + lambda.p

  bad0 <- !is.finite(lambda.p) | lambda.p <= 0 |
          !is.finite(lambda.a) | lambda.a <= 0 |
          !is.finite(lambda.i) | lambda.i <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitpois(Hi, lambda.p,
                   alt.mix = alt.mix, alt.mlm = alt.mlm,
                   inf.mix = inf.mix, inf.mlm = inf.mlm,
                   truncate = truncate, max.support = max.support,
                   pstr.mix = pstr.mix, pobs.mix = pobs.mix,
                   pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
                   lambda.a = lambda.a, lambda.i = lambda.i,
                   byrow.ai = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitpois(Hi[!done], lambda.p[!done],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       lambda.a = lambda.a[!done],
                       lambda.i = lambda.i[!done],
                       byrow.ai = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, lambda.p,
                      alt.mix = NULL, alt.mlm = NULL,
                      inf.mix = NULL, inf.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix = 0, pstr.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0,
                      lambda.a = lambda.p, lambda.i = lambda.p,
                      byrow.ai = FALSE, p)
      pgaitpois(q, lambda.p = lambda.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       byrow.ai = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitpois(min.support.use, lambda.p = lambda.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       byrow.ai = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      lambda.p = lambda.p[!lhs],
                      alt.mix = alt.mix, alt.mlm = alt.mlm,
                      inf.mix = inf.mix,
                      inf.mlm = inf.mlm,
                      truncate = truncate, max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      lambda.a = lambda.a[!lhs],
                      lambda.i = lambda.i[!lhs],
                      byrow.ai = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitpois(faa, lambda.p[!lhs],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       lambda.a = lambda.a[!lhs],
                       lambda.i = lambda.i[!lhs],
                       byrow.ai = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitpois(faa + 1, lambda.p[!lhs],
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       lambda.a = lambda.a[!lhs],
                       lambda.i = lambda.i[!lhs],
                       byrow.ai = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitpois(min.support.use, lambda.p,
                       alt.mix = alt.mix, alt.mlm = alt.mlm,
                       inf.mix = inf.mix,
                       inf.mlm = inf.mlm,
                       truncate = truncate, max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       byrow.ai = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitpois





 rgaitpois <-
  function(n, lambda.p,
           alt.mix = NULL,
           alt.mlm = NULL,
           inf.mix = NULL,
           inf.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           byrow.ai = FALSE,
           lambda.a = lambda.p, lambda.i = lambda.p) {
    qgaitpois(runif(n), lambda.p,
              alt.mix = alt.mix,
              alt.mlm = alt.mlm,
              inf.mix = inf.mix,
              inf.mlm = inf.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              lambda.a = lambda.a, lambda.i = lambda.i,
              byrow.ai = byrow.ai)
}  # rgaitpois













specialsvglm <-
  function(object, ...) {
  infos <- object@family@infos()
  ans <- list(alt.mix  = infos$alt.mix,
              alt.mlm  = infos$alt.mlm,
              inf.mix  = infos$inf.mix,
              inf.mlm  = infos$inf.mlm,
              truncate = infos$truncate)
  if (is.numeric(tmp7e <- infos$max.support))
    ans <- c(ans, max.support = tmp7e)
  ans
}  # specialsvglm


if (!isGeneric("specials"))
  setGeneric("specials", function(object, ...)
             standardGeneric("specials"),
             package = "VGAM")
setMethod("specials", "vglm",
          function(object, ...)
          specialsvglm(object, ...))




if (!isGeneric("altered"))
  setGeneric("altered", function(object, ...)
             standardGeneric("altered"),
             package = "VGAM")

setMethod("altered", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          c(tmp$alt.mix, tmp$alt.mlm)})




if (!isGeneric("inflated"))
  setGeneric("inflated", function(object, ...)
             standardGeneric("inflated"),
             package = "VGAM")

setMethod("inflated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          c(tmp$inf.mix, tmp$inf.mlm)})


if (!isGeneric("truncated"))
  setGeneric("truncated", function(object, ...)
             standardGeneric("truncated"),
             package = "VGAM")

setMethod("truncated", "vglm",
          function(object, ...) {
          ans <- specialsvglm(object, ...)
          if (any(names(ans) == "max.support"))
            ans[c("truncate", "max.support")] else
            ans[["truncate"]]
          })









  setGeneric("is.altered", function(object, ...)
             standardGeneric("is.altered"),
             package = "VGAM")
setMethod("is.altered", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(c(tmp$alt.mix, tmp$alt.mlm)))})



  setGeneric("is.inflated", function(object, ...)
             standardGeneric("is.inflated"),
             package = "VGAM")
setMethod("is.inflated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(c(tmp$inf.mix, tmp$inf.mlm)))})



  setGeneric("is.truncated", function(object, ...)
             standardGeneric("is.truncated"),
             package = "VGAM")
setMethod("is.truncated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(tmp$truncated)) ||
          (length(tmp$max.support) > 0 && is.finite(tmp$max.support))
          })







 y.gaitcombo.check <-
  function(y, truncate = NULL,
           alt.mix = NULL, alt.mlm = NULL,
           inf.mix = NULL, inf.mlm = NULL,
           max.support = Inf, min.support = 0) {
  lalt.mix <- length(alt.mix)
  lalt.mlm <- length(alt.mlm)
  linf.mix <- length(inf.mix)
  linf.mlm <- length(inf.mlm)




  n <- length(y)
  css.mix.a <- css.mix.i <- skip.mix.a <- skip.mix.i <-  # Default
  css.mlm.a <- css.mlm.i <- skip.mlm.a <- skip.mlm.i <- NULL

  if (length(truncate) && any(y %in% truncate))
    stop("some response values == values in argument 'truncate'")
  if (max.support < max(y))
    stop("some response values are greater than the ",
         "'max.support' argument")


  y0.mix.a <- y0.mlm.a <- y0.mix.i <- y0.mlm.i <- NULL
  if (lalt.mix > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mix.a <- matrix(0, n, lalt.mix)
    for (jay in seq(lalt.mix))
      y0.mix.a[, jay] <- as.numeric(y == alt.mix[jay])
    skip.mix.a <- matrix(as.logical(y0.mix.a), n, lalt.mix)  # dim lost
    if (any((css.mix.a <- colSums(skip.mix.a)) == 0))
      stop("some 'alt.mix' argument values have no response values: ",
           paste(alt.mix[css.mix.a == 0], collapse = ", "))          
  }  # lalt.mix

  if (lalt.mlm > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mlm.a <- matrix(0, n, lalt.mlm)
    for (jay in seq(lalt.mlm))
      y0.mlm.a[, jay] <- as.numeric(y == alt.mlm[jay])
    skip.mlm.a <- matrix(as.logical(y0.mlm.a), n, lalt.mlm)  # dim lost
    if (any((css.mlm.a <- colSums(skip.mlm.a)) == 0))
      stop("some 'alt.mlm' argument values have no response values: ",
           paste(alt.mlm[css.mlm.a == 0], collapse = ", "))          
  }  # lalt.mlm

      
  if (linf.mix > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mix.i <- matrix(0, n, linf.mix)
    for (jay in seq(linf.mix))
      y0.mix.i[, jay] <- as.numeric(y == inf.mix[jay])
    skip.mix.i <- matrix(as.logical(y0.mix.i), n, linf.mix)  # dim lost
      if (any((css.mix.i <- colSums(skip.mix.i)) == 0))
      stop("some 'inf.mix' argument values have no response values: ",
           paste(inf.mix[css.mix.i == 0], collapse = ", "))          
  }  # linf.mix

      
  if (linf.mlm > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mlm.i <- matrix(0, n, linf.mlm)
    for (jay in seq(linf.mlm))
      y0.mlm.i[, jay] <- as.numeric(y == inf.mlm[jay])
    skip.mlm.i <- matrix(as.logical(y0.mlm.i), n, linf.mlm)  # dim lost
      if (any((css.mlm.i <- colSums(skip.mlm.i)) == 0))
      stop("some 'inf.mlm' argument values have no response values: ",
           paste(inf.mlm[css.mlm.i == 0], collapse = ", "))          
  }  # linf.mlm


  list(css.mix.a = css.mix.a, skip.mix.a = skip.mix.a,
       css.mix.i = css.mix.i, skip.mix.i = skip.mix.i,
       css.mlm.a = css.mlm.a, skip.mlm.a = skip.mlm.a,
       css.mlm.i = css.mlm.i, skip.mlm.i = skip.mlm.i,
        y0.mix.a =  y0.mix.a,   y0.mlm.a  =  y0.mlm.a, 
        y0.mix.i =  y0.mix.i,   y0.mlm.i  =  y0.mlm.i) 
}  # y.gaitcombo.check











 dgaitplot <-
  function(   # xx, pmf,
           theta.p,  # scalar, else a 2-vector
           fam = "pois",  # "zeta", "log", "genpois0",
           alt.mix = NULL, inf.mix = NULL,  # Unstructured probs are
           alt.mlm = NULL, inf.mlm = NULL,  # contiguous
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # scalar
           pobs.mlm = 0,  # vector of length alt.mlm
           pstr.mix = 0,  # scalar
           pstr.mlm = 0,  # vector of length alt.mlm
           byrow.ai = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           theta.a = theta.p,  # scalar, else a 2-vector
           theta.i = theta.p,  # scalar, else a 2-vector
           deflation = FALSE,  # Single logical
           plot.it = TRUE, new.plot = TRUE,
           offset.x = 0,  # Allows multiple side-by-side plots
           type.plot = "h",  # Matches 'type' argument
           xlim = c(0, min(100, max.support + 2)),
           ylim = NULL,
           xlab = "",  # Was "y" prior to using oma
           ylab = "Probability",
           main = "",
           cex.main = 1.2, posn.main = NULL,
           lty.p     = "solid",  # longdash, dashed, twodash, solid
           lty.a.mix = "longdash",
           lty.a.mlm = "longdash",
           lty.i.mix = "dashed",
           lty.i.mlm = "dashed",
           col.p = "pink2",  # "salmon"  # salmon1,..., salmon4
           col.a.mix = "#007FFF",  #  Wiki azure
           col.a.mlm = "blue",
           col.i.mix = "#3F00FF",  # Indigo, diffs a bit from "purple"
           col.i.mlm = "purple",   # 
           col.t = "tan",
           cex.p = 1,
           lwd.p = NULL, lwd.a = NULL, lwd.i = NULL,  # Default: par()$lwd
           iontop = TRUE,
           las = 0, lend = "round",  # "round", "butt", "square", 0:2
           axes.x = TRUE, axes.y = TRUE,
           Plot.trunc = TRUE, cex.t = 1, pch.t = 1,
           baseparams.argnames = NULL,   #,  # Optional safety
           ...
           ) {  # ... ignored currently.


  if (!length(lwd.p)) lwd.p <- par()$lwd
  if (!length(lwd.a)) lwd.a <- par()$lwd
  if (!length(lwd.i)) lwd.i <- par()$lwd


  gait.errorcheck(alt.mix, alt.mlm, inf.mix, inf.mlm,
                  truncate, max.support)  #, min.support = lowsup

  MM <- length(theta.p)
  if (MM != 1 && MM !=  2)
    stop("can only handle 1 or 2 parameters")
  
  xx <- seq(xlim[1], xlim[2])
  ind.a.mlm <- ind.a.mix <- ind.trunc <-
  ind.i.mlm <- ind.i.mix <- FALSE

  if (length(truncate))
    ind.trunc <- (xx %in% truncate | max.support < xx)
  if (length(alt.mix)) {
    ind.a.mix <- xx %in% alt.mix
  } else {
    pobs.mix <- 0  # Make sure
  }
  if (length(alt.mlm)) {
    ind.a.mlm <- xx %in% alt.mlm
  } else {
    pobs.mlm <- 0  # Make sure
  }
  if (length(inf.mix)) {
    ind.i.mix <- xx %in% inf.mix
  } else {
    pstr.mix <- 0  # Make sure
  }
  if (length(inf.mlm)) {
    ind.i.mlm <- xx %in% inf.mlm
  } else {
    pstr.mlm <- 0  # Make sure
  }

  special.xx <- ind.a.mix | ind.a.mlm | ind.i.mix | ind.i.mlm |
                ind.trunc
  if (length(pobs.mix) != 1) stop("bad input for argument 'pobs.mix'")
  if (length(pstr.mix) != 1) stop("bad input for argument 'pstr.mix'")
  if (length(alt.mlm) && length(pobs.mlm) > length(alt.mlm))
    warning("bad input for argument 'pobs.mlm'?")
  if (length(inf.mlm) && length(pstr.mlm) > length(inf.mlm))
    warning("bad input for argument 'pstr.mlm'?")

      
  if (any(ind.a.mlm))
    pobs.mlm <- matrix(pobs.mlm, 1,  # length(xx),
                       length(alt.mlm), byrow = byrow.ai)
  if (any(ind.i.mlm))
    pstr.mlm <- matrix(pstr.mlm, 1,  # length(xx),
                       length(inf.mlm), byrow = byrow.ai)

      
  dfun <- paste0("dgait", fam)

  pmf.p <- if (MM == 1) do.call(dfun, list(x =  xx, theta.p)) else
           do.call(dfun, list(x =  xx, theta.p[1], theta.p[2]))


  alist <- list(  # x = xx,  # theta.p,
            alt.mix = alt.mix, alt.mlm = alt.mlm,
            inf.mix = inf.mix, inf.mlm = inf.mlm,
            truncate = truncate, max.support = max.support,
            pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
            pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
            byrow.ai = byrow.ai)

  if (length(baseparams.argnames)) {
    alist[[paste0(baseparams.argnames[1], ".p")]] <- theta.p[1]
    alist[[paste0(baseparams.argnames[1], ".a")]] <- theta.a[1]
    alist[[paste0(baseparams.argnames[1], ".i")]] <- theta.i[1]
    if (MM == 2) {
      alist[[paste0(baseparams.argnames[2], ".p")]] <- theta.p[2]
      alist[[paste0(baseparams.argnames[2], ".a")]] <- theta.a[2]
      alist[[paste0(baseparams.argnames[2], ".i")]] <- theta.i[2]
    }
  } else {
    if (MM == 1) {
      alist <- c(alist, list(theta.p, theta.a, theta.i))
    } else {  # MM == 2
      alist <- c(alist,  # Unnamed, for lambda.p, etc.:
                 list(theta.p[1], theta.p[2],
                      theta.a[1], theta.i[1],  # Order is crucial.
                      theta.a[2], theta.i[2]))
    }
  }


  dlist <- alist
  dlist$x <- xx
  dlist$deflation <- deflation
  dlist$log <- FALSE
  pmf.z <- do.call(dfun, dlist)


  mlist <- alist
  mlist$type.fitted <- "All"
  mlist$moments2 <- TRUE
  mom.fun <- paste0("moments.gaitcombo.", fam)
  Bits <- do.call(mom.fun, mlist)

  myylim <- if (is.null(ylim))
    c(0, max(0, pmf.z, na.rm = TRUE) * 1.04) else ylim

  if (plot.it) {
    if (new.plot)
      plot(xx[!special.xx] + offset.x,
           pmf.z[!special.xx], type = "n",
           las = las, axes = axes.x && axes.y,  # axes,
           xlim = xlim,
           ylim = myylim,  # c(0, 0.11),  # Fix ylim here

         
           xlab = xlab, ylab = ylab)  # Parent distribution

    lines(xx[!special.xx] + offset.x,
          pmf.z[!special.xx], type = type.plot, col = col.p,
          lwd = lwd.p,  # Ordinary points
          cex = cex.p, lend = lend, lty = lty.p)

    if (length(posn.main)) {
      posn.main <- rep(posn.main, 2)
      text(posn.main[1], posn.main[2], labels = main,
           cex = cex.main)
    } else {
      title(main = main)  # , cex = cex.main
    }
  }


      
  if (plot.it && !(axes.x && axes.y)) {
    box()
    axis(1, labels = FALSE, tick = TRUE)
    axis(2, labels = FALSE, tick = TRUE)
  }
  if (plot.it && !axes.x && axes.y) {  # Avoid clutter
    axis(1, labels = FALSE)
    axis(2, las = 1)
  }
  if (plot.it && axes.x && !axes.y) {  # Avoid clutter
    axis(1)
    axis(2, labels = FALSE)
  }


  if (plot.it && length(alt.mlm))  # Altered mlm
    lines(xx[ ind.a.mlm] + offset.x, pmf.z[ind.a.mlm],
          lwd = lwd.a, type = type.plot,
          col = col.a.mlm, lty = lty.a.mlm, lend = lend)


  if (plot.it && length(alt.mix))  # Altered mix
    lines(xx[ ind.a.mix] + offset.x, pmf.z[ind.a.mix],
          lwd = lwd.a, type = type.plot,
          col = col.a.mix, lty = lty.a.mix, lend = lend)




  if (any(c(ind.i.mix, ind.i.mlm))) {  # Inflated mix or mlm

    Denom.p <- as.vector(Bits[["cdf.max.s"]] - Bits[["SumT0.p"]] -
                         Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    if (any(Denom.p == 0))
      stop("0s found in the denominator (variable 'Denom.p')")
    Numer <- as.vector(1 -
               (if (length(alt.mix)) pobs.mix else 0) -
               (if (length(inf.mix)) pstr.mix else 0) -
               (if (length(alt.mlm)) rowSums(rbind(pobs.mlm)) else 0) -
               (if (length(inf.mlm)) rowSums(rbind(pstr.mlm)) else 0))
    Numer <- Numer[1]
    if (Numer < 0 && !deflation)
      warning("variable 'Numer' is negative and 'deflation = FALSE'")
  }


      

  if (any(ind.i.mix)) {  # Inflated mix

    start.pt.mix <- Numer * pmf.p[ind.i.mix] / Denom.p
    spikes.mix <- pmf.z[ind.i.mix]  # Top of the spike

  if (plot.it)
    segments(xx[ind.i.mix] + offset.x,
             if (iontop) start.pt.mix else spikes.mix - start.pt.mix,
             xx[ind.i.mix] + offset.x,
             spikes.mix,  # This value is unchanged
             lwd = if (iontop) lwd.i else lwd.p,
             lty = if (iontop) lty.i.mix else lty.p,
             col = if (iontop) col.i.mix else col.p, lend = lend)


  if (plot.it)
    lines(xx[ind.i.mix] + offset.x,
          if (iontop) start.pt.mix else spikes.mix - start.pt.mix,
          lend = lend,
          lwd = if (iontop) lwd.p else lwd.i,
          lty = if (iontop) lty.p else lty.i.mix,
          col = if (iontop) col.p else col.i.mix,
          type = type.plot)  # Blend in
  }  # ind.i.mix


      

  if (any(ind.i.mlm)) {  # Inflated mlm
    start.pt.mlm <- Numer * pmf.p[ind.i.mlm] / Denom.p
    spikes.mlm <- pmf.z[ind.i.mlm]  # Top of the spike

  if (plot.it)
    segments(xx[ind.i.mlm] + offset.x,
             if (iontop) start.pt.mlm else spikes.mlm - start.pt.mlm,
             xx[ind.i.mlm] + offset.x,
             spikes.mlm,  # This value is unchanged
             lwd = if (iontop) lwd.i else lwd.p,
             col = if (iontop) col.i.mlm else col.p,
             lend = lend,
             lty = if (iontop) lty.i.mlm else lty.p)


  if (plot.it)
    lines(xx[ind.i.mlm] + offset.x,
          if (iontop) start.pt.mlm else spikes.mlm - start.pt.mlm,
          lend = lend,
          lwd = if (iontop) lwd.p else lwd.i,
          lty = if (iontop) lty.p else lty.i.mlm,
          col = if (iontop) col.p else col.i.mlm,
          type = type.plot)  # Blend in
  }  # ind.i.mlm

      

  if (Plot.trunc && plot.it && any(ind.trunc)) {
    lhs.tvec <- xx[ind.trunc] + offset.x * 0
    lhs.tvec <- lhs.tvec[ceiling(par()$usr[1]) + 0.5 < lhs.tvec]
    if (length(lhs.tvec))
      points(lhs.tvec, numeric(length(lhs.tvec)),
             col = col.t, cex = cex.t, pch = pch.t)
  }
  if (Plot.trunc && plot.it) {
    rhs.tvec <- if (floor(par()$usr[2]) > max.support + 2)
                (max.support + 1):(floor(par()$usr[2]) - 1) else NULL
    if (length(rhs.tvec))
      points(rhs.tvec + offset.x * 0,
             numeric(length(rhs.tvec)),
             col = col.t, cex = cex.t, pch = pch.t)
  }


  ans <- pmf.p
  if (any(ind.i.mix)) ans[ind.i.mix] <- start.pt.mix
  if (any(ind.i.mlm)) ans[ind.i.mlm] <- start.pt.mlm
  invisible(list(x = xx,
                 pmf.z = pmf.z,
                 mid.pmf = ans))
}  # dgaitplot







plotdgait.vglm <-
  function(object, ...) {

  infos.list <- object@family@infos()
  specvals <- specials(object)


  Inside <- sapply(specvals, is.null)
  if (length(Inside) == 5 && all(Inside))
    stop("'object' has no special values. ",
         "Is it a GAIT regression object?")
  if (length(Inside) == 6 && all(Inside[1:5]) &&
      infos.list$max.support == infos.list$Support[2])
    stop("'object' has no special values. ",
         "Is it really a GAIT regression object?")

  if (!is.numeric(MM1 <- infos.list$MM1))
    MM1 <- 1  # Default really
  if (MM1 > 2)
    stop("Can only handle 1- or 2-parameter distributions")

  etamat <- predict(object)  # n x M
  eta.p <- etamat[, 1:MM1, drop = FALSE]  # n x MM1
  theta.p1 <- as.vector(eta2theta(eta.p[, 1], linkfun(object)[1]))
  theta.p2 <- if (MM1 == 2)
    as.vector(eta2theta(eta.p[, 2], linkfun(object)[2])) else NULL
  theta.p <- cbind(theta.p1, theta.p2)
  if (!is.logical(intercept.only <- object@misc$intercept.only))
    stop("cannot determine whether 'object' is intercept-only")
  if (!intercept.only)
    warning("argument 'object' is not intercept-only")



  Pobs.mix <- if (length(specvals$alt.mix))
    fitted(object, type.fitted = "Pobs.mix") else cbind(0, 0)
  Pobs.mlm <- if (length(specvals$alt.mlm))
    fitted(object, type.fitted = "pobs.mlm") else cbind(0, 0)
  Pstr.mix <- if (length(specvals$inf.mix))
    fitted(object, type.fitted = "Pstr.mix") else cbind(0, 0)
  Pstr.mlm <- if (length(specvals$inf.mlm))
    fitted(object, type.fitted = "pstr.mlm") else cbind(0, 0)
  tfitted  <- paste0(infos.list$baseparams.argnames, "s")
  Thetas   <- fitted(object, type.fitted = tfitted)


 indeta <- object@extra$indeta
 if (MM1 == 1) {
   theta.a <- if (any(is.na(indeta[3, ]))) theta.p else
                as.vector(eta2theta(etamat[, (indeta[3, 1])],
                          linkfun(object)[(indeta[3, 1])]))
   theta.i <- if (any(is.na(indeta[5, ]))) theta.p else
                as.vector(eta2theta(etamat[, (indeta[5, 1])],
                          linkfun(object)[(indeta[5, 1])]))
 } else {
   stop("cannot handle MM1 != 1 just now")
 }

      

      
  sum.mlm.i <- if (length(specvals$inf.mlm))
    fitted(object, type.fitted = "sum.mlm.i") else cbind(0, 0)


  dgaitplot(theta.p[1, ], fam = infos.list$parent.name[2],
            alt.mix = specvals$alt.mix, inf.mix = specvals$inf.mix, 
            alt.mlm = specvals$alt.mlm, inf.mlm = specvals$inf.mlm,
            truncate = specvals$truncate,
            theta.a = theta.a, theta.i = theta.i,
            max.support = specvals$max.support,
            pobs.mix = sum(Pobs.mix[1, ]),
            pobs.mlm = Pobs.mlm[1, ],
            pstr.mix = sum(Pstr.mix[1, ]),
            pstr.mlm = Pstr.mlm[1, ],
            byrow.ai = TRUE,  # Important really here 20201008
            baseparams.argnames = infos.list$baseparams.argnames,  # safe
            ...)
}  # plotdgait.vglm



 if (!isGeneric("plotdgait"))
   setGeneric("plotdgait",
              function(object, ...) standardGeneric("plotdgait"))

  setMethod("plotdgait", signature(object = "vglm"),
            function(object, ...)
            invisible(plotdgait.vglm(object, ...)))

















