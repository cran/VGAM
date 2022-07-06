# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.












 dgaitdzeta <-
  function(x, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
       byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
       shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
           log = FALSE) {


  log.arg <- log;  rm(log)
  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(dzeta(x, shape.p, log = log.arg))

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")


  LLL <- max(length(x),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(shape.p),  length(shape.a),
             length(shape.i),  length(shape.d))
  if (length(x)        < LLL) x         <- rep_len(x,         LLL)
  if (length(shape.p)  < LLL) shape.p   <- rep_len(shape.p,   LLL)
  if (length(shape.a)  < LLL) shape.a   <- rep_len(shape.a,   LLL)
  if (length(shape.i)  < LLL) shape.i   <- rep_len(shape.i,   LLL)
  if (length(shape.d)  < LLL) shape.d   <- rep_len(shape.d,   LLL)
  if (length(pobs.mix) < LLL) pobs.mix  <- rep_len(pobs.mix,  LLL)
  if (length(pstr.mix) < LLL) pstr.mix  <- rep_len(pstr.mix,  LLL)
  if (length(pdip.mix) < LLL) pdip.mix  <- rep_len(pdip.mix,  LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dzeta(tval, shape.p)  # Need tval <= max.support
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- pzeta(max.support, shape.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dzeta(x, shape.p) / denom.t)  # dgtlog


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (la.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <- .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    for (aval in a.mlm)
      suma <- suma + dzeta(aval, shape.p)  # Part i

    for (jay in seq(la.mlm)) {
      aval <- a.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
          pmf0[vecTF] <- pobs.mlm[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # la.mlm



  pmf2.a <- pmf2.i <- pmf2.d <- 0
  if (la.mix) {
    allx.a <- lowsup:max(a.mix)
    pmf2.a <- dgaitdzeta(x, shape.a,  # Outer distribution---mlm type
                         truncate = setdiff(allx.a, a.mix),
                         max.support = max(a.mix))
    for (aval in a.mix) {
      suma <- suma + dzeta(aval, shape.p)  # Part ii added; cumulative
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }

  if (li.mix) {
    allx.i <- if (length(i.mix)) lowsup:max(i.mix) else NULL
    pmf2.i <- dgaitdzeta(x, shape.i,  # Outer distribution---mlm type
                         truncate = setdiff(allx.i, i.mix),
                         max.support = max(i.mix))
  }


      
  sum.d <- 0  # numeric(LLL)
  if (ld.mlm) {
    pdip.mlm <-  matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <- .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")
  }  # ld.mlm


  if (ld.mix) {
    allx.d <- lowsup:max(d.mix)
    pmf2.d <- dgaitdzeta(x, shape.d,  # Outer distn---mlm type
                         truncate = setdiff(allx.d, d.mix),
                         max.support = max(d.mix))
  }  # ld.mix






  sum.i <- 0
  if (li.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, li.mlm, byrow = byrow.aid)
    sum.i <- .rowSums(pstr.mlm, LLL, li.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")
  }


  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix - pstr.mix + sum.d + pdip.mix
  if (li.mlm) {
    if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
      warning("the vector of normalizing constants contains ",
              "some negative values. Replacing them with NAs")
      tmp6[!skip & tmp6 < 0] <- NA
    }
  }  # li.mlm


  pmf0[!skip] <-
    (tmp6 * dzeta(x, shape.p) / (cdf.max.s - suma - sumt))[!skip]



  if (li.mlm) {
    for (jay in seq(li.mlm)) {
      ival <- i.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
          pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # li.mlm




  if (ld.mlm) {
    for (jay in seq(ld.mlm)) {
      dval <- d.mlm[jay]
      if (any(vecTF <- is.finite(x) & dval == x)) {
          pmf0[vecTF] <- pmf0[vecTF] - pdip.mlm[vecTF, jay]
      }
    }  # jay
  }  # ld.mlm




  pmf0 <- pmf0 + pobs.mix * pmf2.a + pstr.mix * pmf2.i -
                 pdip.mix * pmf2.d


  if (log.arg) log(pmf0) else pmf0
}  # dgaitdzeta






 pgaitdzeta <-
  function(q, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
           lower.tail = TRUE) {


  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(pzeta(q, shape.p, lower.tail = lower.tail))  # log.p

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0

  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")

  LLL <- max(length(q),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(shape.p),  length(shape.a),
             length(shape.i),  length(shape.d))
  offset.a <- offset.i <- offset.d <-
  Offset.a <- Offset.i <- Offset.d <- numeric(LLL)
  if (length(q)        < LLL) q        <- rep_len(q,         LLL)
  if (length(shape.p)  < LLL) shape.p  <- rep_len(shape.p,   LLL)
  if (length(shape.a)  < LLL) shape.a  <- rep_len(shape.a,   LLL)
  if (length(shape.i)  < LLL) shape.i  <- rep_len(shape.i,   LLL)
  if (length(shape.d)  < LLL) shape.d  <- rep_len(shape.d,   LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix,  LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix,  LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix,  LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- pzeta(max.support, shape.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dzeta(tval, shape.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (la.mlm) {
    pobs.mlm <- matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <- .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")

    for (jay in seq(la.mlm)) {
      aval <- a.mlm[jay]
      pmf.p <- dzeta(aval, shape.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # la.mlm

  sum.i <- 0
  if (li.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, li.mlm, byrow = byrow.aid)
    sum.i <- .rowSums(pstr.mlm, LLL, li.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")

    for (jay in seq(li.mlm)) {
      ival <- i.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # li.mlm



  use.pobs.mix <- 0
  if (la.mix) {
    use.pobs.mix <- matrix(0, LLL, la.mix)
    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.a <- dzeta(aval, shape.a)
      pmf.p <- dzeta(aval, shape.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                    use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.p <- dzeta(aval, shape.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # la.mix

  use.pstr.mix <- 0
  if (li.mix) {
    use.pstr.mix <- matrix(0, LLL, li.mix)
    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      use.pstr.mix[, jay] <- dzeta(ival, shape.i)
    }
    use.pstr.mix <- pstr.mix *
                    use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      pmf.p <- dzeta(ival, shape.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix[vecTF, jay]
      }
    }  # jay
  }  # li.mix






  sum.d <- 0
  if (ld.mlm) {
    pdip.mlm <- matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <-  .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")

    for (jay in seq(ld.mlm)) {
      dval <- d.mlm[jay]
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        offset.d[vecTF] <- offset.d[vecTF] + pdip.mlm[vecTF, jay]
      }
    }  # jay
  }  # ld.mlm



  use.pdip.mix <- 0
  if (ld.mix) {
    use.pdip.mix <- matrix(0, LLL, ld.mix)
    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      use.pdip.mix[, jay] <- dzeta(dval, shape.d)
    }
    use.pdip.mix <- pdip.mix *
                    use.pdip.mix / rowSums(use.pdip.mix)

    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      pmf.p <- dzeta(dval, shape.p)
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        Offset.d[vecTF] <- Offset.d[vecTF] + use.pdip.mix[vecTF, jay]
      }
    }  # jay
  }  # ld.mix






  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix + sum.d + pdip.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * (pzeta(q, shape.p) - fudge.t - fudge.a) / denom1 +
         offset.a + offset.i - offset.d +
         Offset.a + Offset.i - Offset.d
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitdzeta






 qgaitdzeta <-
  function(p, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           shape.a = shape.p, shape.i = shape.p,
           shape.d = shape.p) {


  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc     <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(qzeta(p, shape.p))  # lower.tail = TRUE, log.p = FALSE

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")


  LLL <- max(length(p),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(shape.p),  length(shape.a),
             length(shape.i),  length(shape.d))
  if (length(p)        < LLL) p        <- rep_len(p,        LLL)
  if (length(shape.p)  < LLL) shape.p  <- rep_len(shape.p,  LLL)
  if (length(shape.a)  < LLL) shape.a  <- rep_len(shape.a,  LLL)
  if (length(shape.i)  < LLL) shape.i  <- rep_len(shape.i,  LLL)
  if (length(shape.d)  < LLL) shape.d  <- rep_len(shape.d,  LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)

  pobs.mlm <- matrix(pobs.mlm, LLL, max(la.mlm, 1),
                     byrow = byrow.aid)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(li.mlm, 1),
                     byrow = byrow.aid)
  pdip.mlm <- matrix(pdip.mlm, LLL, max(ld.mlm, 1),
                     byrow = byrow.aid)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + shape.p

  bad0 <- !is.finite(shape.p) | shape.p <= 0 |
          !is.finite(shape.a) | shape.a <= 0 |
          !is.finite(shape.i) | shape.i <= 0 |
          !is.finite(shape.d) | shape.d <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  use.base <- 2  # bbbb
  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else use.base * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitdzeta(Hi, shape.p,
                   a.mix = a.mix, a.mlm = a.mlm,
                   i.mix = i.mix, i.mlm = i.mlm,
                   d.mix = d.mix, d.mlm = d.mlm,
                   truncate = truncate, max.support = max.support,
                   pstr.mix = pstr.mix, pobs.mix = pobs.mix,
                   pdip.mix = pdip.mix,
                   pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
                   pdip.mlm = pdip.mlm,
                   shape.a = shape.a, shape.i = shape.i,
                   shape.d = shape.d,
                   byrow.aid = FALSE)

  iter <- 0
  max.iter <- round(log(.Machine$double.xmax, base = use.base)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- use.base * Hi[!done] + 10.5  # Bug fixed bbbb
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitdzeta(Hi[!done], shape.p[!done],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pdip.mix = pdip.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       pdip.mlm = pdip.mlm[!done, , drop = FALSE],
                       shape.a = shape.a[!done],
                       shape.i = shape.i[!done],
                       shape.d = shape.d[!done],
                       byrow.aid = FALSE))
    iter <- iter + 1
  }  # while

      foo <- function(q, shape.p,
                      a.mix = NULL, a.mlm = NULL,
                      i.mix = NULL, i.mlm = NULL,
                      d.mix = NULL, d.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix = 0, pstr.mix = 0, pdip.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0, pdip.mlm = 0,
                      shape.a = shape.p, shape.i = shape.p,
                      shape.d = shape.p,
                      byrow.aid = FALSE, p)
      pgaitdzeta(q, shape.p = shape.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       shape.a = shape.a, shape.i = shape.i,
                       shape.d = shape.p,
                       byrow.aid = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitdzeta(min.support.use, shape.p = shape.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       shape.a = shape.a, shape.i = shape.i,
                       shape.d = shape.d,
                       byrow.aid = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      shape.p = shape.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      shape.a  = shape.a[!lhs],
                      shape.i  = shape.i[!lhs],
                      shape.d  = shape.d[!lhs],
                      byrow.aid = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitdzeta(faa, shape.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       shape.a  = shape.a[!lhs],
                       shape.i  = shape.i[!lhs],
                       shape.d  = shape.d[!lhs],
                       byrow.aid = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitdzeta(faa + 1, shape.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      shape.a  = shape.a[!lhs],
                      shape.i  = shape.i[!lhs],
                      shape.d  = shape.d[!lhs],
                      byrow.aid = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitdzeta(min.support.use, shape.p,
                          a.mix = a.mix, a.mlm = a.mlm,
                          i.mix = i.mix, i.mlm = i.mlm,
                          d.mix = d.mix, d.mlm = d.mlm,
                          truncate = truncate,
                          max.support = max.support,
                          pobs.mix = pobs.mix,
                          pstr.mix = pstr.mix,
                          pdip.mix = pdip.mix,
                          pobs.mlm = pobs.mlm,
                          pstr.mlm = pstr.mlm,
                          pdip.mlm = pdip.mlm,
                          shape.a = shape.a, shape.i = shape.i,
                          shape.d = shape.d,
                          byrow.aid = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitdzeta





 rgaitdzeta <-
  function(n, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,
           shape.a = shape.p, shape.i = shape.p,
           shape.d = shape.p) {
    qgaitdzeta(runif(n), shape.p,
              a.mix = a.mix,
              a.mlm = a.mlm,
              i.mix = i.mix,
              i.mlm = i.mlm,
              d.mix = d.mix,
              d.mlm = d.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              shape.a = shape.a, shape.i = shape.i,
              shape.d = shape.d,
              byrow.aid = byrow.aid)
}  # rgaitdzeta







 dgaitdlog <-
  function(x, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
           log = FALSE) {

  log.arg <- log;  rm(log)
  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
   ltrunc     <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(dlog(x, shape.p, log = log.arg))

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")

  LLL <- max(length(x),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(shape.p),  length(shape.a),
             length(shape.i),  length(shape.d))
  if (length(x)         < LLL) x        <- rep_len(x,        LLL)
  if (length(shape.p)   < LLL) shape.p  <- rep_len(shape.p,  LLL)
  if (length(shape.a)   < LLL) shape.a  <- rep_len(shape.a,  LLL)
  if (length(shape.i)   < LLL) shape.i  <- rep_len(shape.i,  LLL)
  if (length(shape.d)   < LLL) shape.d  <- rep_len(shape.d,  LLL)
  if (length(pobs.mix)  < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix)  < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix)  < LLL) pdip.mix <- rep_len(pdip.mix, LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dlog(tval, shape.p)  # Need tval <= max.support
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- plog(max.support, shape.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dlog(x, shape.p) / denom.t)  # dgtlog


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (la.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <- .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    for (aval in a.mlm)
      suma <- suma + dlog(aval, shape.p)  # Part i

    for (jay in seq(la.mlm)) {
      aval <- a.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
          pmf0[vecTF] <- pobs.mlm[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # la.mlm



  pmf2.a <- pmf2.i <- pmf2.d <- 0
  if (la.mix) {
    allx.a <- lowsup:max(a.mix)
    pmf2.a <- dgaitdlog(x, shape.a,  # Outer distribution---mlm type
                        truncate = setdiff(allx.a, a.mix),
                        max.support = max(a.mix))
    for (aval in a.mix) {
      suma <- suma + dlog(aval, shape.p)  # Part ii added; cumulative
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }

  if (li.mix) {
    allx.i <- if (length(i.mix)) lowsup:max(i.mix) else NULL
    pmf2.i <- dgaitdlog(x, shape.i,  # Outer distribution---mlm type
                        truncate = setdiff(allx.i, i.mix),
                        max.support = max(i.mix))
  }





  sum.d <- 0  # numeric(LLL)
  if (ld.mlm) {
    pdip.mlm <-  matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <- .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")
  }  # ld.mlm


  if (ld.mix) {
    allx.d <- lowsup:max(d.mix)
    pmf2.d <- dgaitdlog(x, shape.d,  # Outer distn---mlm type
                        truncate = setdiff(allx.d, d.mix),
                        max.support = max(d.mix))
  }  # ld.mix







  sum.i <- 0
  if (li.mlm) {
    pstr.mlm <-  matrix(pstr.mlm, LLL, li.mlm, byrow = byrow.aid)
    sum.i <- .rowSums(pstr.mlm, LLL, li.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")
  }

  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix - pstr.mix + sum.d + pdip.mix
  if (li.mlm) {
    if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
      warning("the vector of normalizing constants contains ",
              "some negative values. Replacing them with NAs")
      tmp6[!skip & tmp6 < 0] <- NA
    }
  }  # li.mlm



  pmf0[!skip] <- (tmp6 *
    dlog(x, shape.p) / (cdf.max.s - suma - sumt))[!skip]



  if (li.mlm) {
    for (jay in seq(li.mlm)) {
      ival <- i.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
          pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # li.mlm



  if (ld.mlm) {
    for (jay in seq(ld.mlm)) {
      dval <- d.mlm[jay]
      if (any(vecTF <- is.finite(x) & dval == x)) {
          pmf0[vecTF] <- pmf0[vecTF] - pdip.mlm[vecTF, jay]
      }
    }  # jay
  }  # ld.mlm


  pmf0 <- pmf0 + pobs.mix * pmf2.a + pstr.mix * pmf2.i -
                 pdip.mix * pmf2.d

  if (log.arg) log(pmf0) else pmf0
}  # dgaitdlog






 pgaitdlog <-
  function(q, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
           lower.tail = TRUE) {


  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
   ltrunc     <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(plog(q, shape.p, lower.tail = lower.tail))  # log.p

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0

  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")

  LLL <- max(length(q),
             length(pobs.mix),length(pstr.mix), length(pdip.mix),
             length(shape.p), length(shape.a),
             length(shape.i), length(shape.d))
  offset.a <- offset.i <- offset.d <-
  Offset.a <- Offset.i <- Offset.d <- numeric(LLL)
  if (length(q)        < LLL) q        <- rep_len(q,        LLL)
  if (length(shape.p)  < LLL) shape.p  <- rep_len(shape.p,  LLL)
  if (length(shape.a)  < LLL) shape.a  <- rep_len(shape.a,  LLL)
  if (length(shape.i)  < LLL) shape.i  <- rep_len(shape.i,  LLL)
  if (length(shape.d)  < LLL) shape.d  <- rep_len(shape.d,  LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)


  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- plog(max.support, shape.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dlog(tval, shape.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (la.mlm) {
    pobs.mlm <- matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <- .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")

    for (jay in seq(la.mlm)) {
      aval <- a.mlm[jay]
      pmf.p <- dlog(aval, shape.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # la.mlm

  sum.i <- 0
  if (li.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, li.mlm, byrow = byrow.aid)
    sum.i <- .rowSums(pstr.mlm, LLL, li.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")

    for (jay in seq(li.mlm)) {
      ival <- i.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # li.mlm



  use.pobs.mix <- 0
  if (la.mix) {
    use.pobs.mix <- matrix(0, LLL, la.mix)
    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.a <- dlog(aval, shape.a)
      pmf.p <- dlog(aval, shape.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                    use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.p <- dlog(aval, shape.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # la.mix

  use.pstr.mix <- 0
  if (li.mix) {
    use.pstr.mix <- matrix(0, LLL, li.mix)
    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      use.pstr.mix[, jay] <- dlog(ival, shape.i)
    }
    use.pstr.mix <- pstr.mix *
                    use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      pmf.p <- dlog(ival, shape.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix[vecTF, jay]
      }
    }  # jay
  }  # li.mix




  sum.d <- 0
  if (ld.mlm) {
    pdip.mlm <- matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <- .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")

    for (jay in seq(ld.mlm)) {
      dval <- d.mlm[jay]
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        offset.d[vecTF] <- offset.d[vecTF] + pdip.mlm[vecTF, jay]
      }
    }  # jay
  }  # ld.mlm



  use.pdip.mix <- 0
  if (ld.mix) {
    use.pdip.mix <- matrix(0, LLL, ld.mix)
    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      use.pdip.mix[, jay] <- dlog(dval, shape.d)
    }
    use.pdip.mix <- pdip.mix *
                    use.pdip.mix / rowSums(use.pdip.mix)

    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      pmf.p <- dlog(dval, shape.p)
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        Offset.d[vecTF] <- Offset.d[vecTF] + use.pdip.mix[vecTF, jay]
      }
    }  # jay
  }  # ld.mix





  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix + sum.d + pdip.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * (plog(q, shape.p) - fudge.t - fudge.a) / denom1 +
         offset.a + offset.i - offset.d +
         Offset.a + Offset.i - Offset.d
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitdlog






 qgaitdlog <-
  function(p, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           shape.a = shape.p, shape.i = shape.p,
           shape.d = shape.p) {


  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc     <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(qlog(p, shape.p))  # lower.tail = TRUE, log.p = FALSE


  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")


  LLL <- max(length(p),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(shape.p),  length(shape.a),
              length(shape.i), length(shape.d))
  if (length(p)        < LLL) p        <- rep_len(p,        LLL)
  if (length(shape.p)  < LLL) shape.p  <- rep_len(shape.p,  LLL)
  if (length(shape.a)  < LLL) shape.a  <- rep_len(shape.a,  LLL)
  if (length(shape.i)  < LLL) shape.i  <- rep_len(shape.i,  LLL)
  if (length(shape.d)  < LLL) shape.d  <- rep_len(shape.d,  LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)

  pobs.mlm <- matrix(pobs.mlm, LLL, max(la.mlm, 1), byrow = byrow.aid)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(li.mlm, 1), byrow = byrow.aid)
  pdip.mlm <- matrix(pdip.mlm, LLL, max(ld.mlm, 1), byrow = byrow.aid)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + shape.p

  bad0 <- !is.finite(shape.p) | shape.p <= 0 |
          !is.finite(shape.a) | shape.a <= 0 |
          !is.finite(shape.i) | shape.i <= 0 |
          !is.finite(shape.d) | shape.d <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitdlog(Hi, shape.p,
                  a.mix = a.mix, a.mlm = a.mlm,
                  i.mix = i.mix, i.mlm = i.mlm,
                  d.mix = d.mix, d.mlm = d.mlm,
                  truncate = truncate, max.support = max.support,
                  pstr.mix = pstr.mix, pobs.mix = pobs.mix,
                  pdip.mix = pdip.mix,
                  pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
                  pdip.mlm = pdip.mlm,
                  shape.a = shape.a, shape.i = shape.i,
                  shape.d = shape.d,
                  byrow.aid = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitdlog(Hi[!done], shape.p[!done],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pdip.mix = pdip.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       pdip.mlm = pdip.mlm[!done, , drop = FALSE],
                       shape.a = shape.a[!done],
                       shape.i = shape.i[!done],
                       shape.d = shape.d[!done],
                       byrow.aid = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, shape.p,
                      a.mix = NULL, a.mlm = NULL,
                      i.mix = NULL, i.mlm = NULL,
                      d.mix = NULL, d.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix = 0, pstr.mix = 0, pdip.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0, pdip.mlm = 0,
                      shape.a = shape.p, shape.i = shape.p,
                      shape.d = shape.p,
                      byrow.aid = FALSE, p)
      pgaitdlog(q, shape.p = shape.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       shape.a = shape.a, shape.i = shape.i,
                       shape.d = shape.d,
                       byrow.aid = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitdlog(min.support.use, shape.p = shape.p,
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pobs.mix = pobs.mix,
                      pstr.mix = pstr.mix,
                      pdip.mix = pdip.mix,
                      pobs.mlm = pobs.mlm,
                      pstr.mlm = pstr.mlm,
                      pdip.mlm = pdip.mlm,
                      shape.a = shape.a, shape.i = shape.i,
                      shape.d = shape.d,
                      byrow.aid = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      shape.p = shape.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      shape.a = shape.a[!lhs],
                      shape.i = shape.i[!lhs],
                      shape.d = shape.d[!lhs],
                      byrow.aid = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitdlog(faa, shape.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      shape.a  = shape.a[!lhs],
                      shape.i  = shape.i[!lhs],
                      shape.d  = shape.d[!lhs],
                      byrow.aid = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitdlog(faa + 1, shape.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      shape.a  = shape.a[!lhs],
                      shape.i  = shape.i[!lhs],
                      shape.d  = shape.d[!lhs],
                      byrow.aid = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitdlog(min.support.use, shape.p,
                         a.mix = a.mix, a.mlm = a.mlm,
                         i.mix = i.mix, i.mlm = i.mlm,
                         d.mix = d.mix, d.mlm = d.mlm,
                         truncate = truncate,
                         max.support = max.support,
                         pobs.mix = pobs.mix,
                         pstr.mix = pstr.mix,
                         pdip.mix = pdip.mix,
                         pobs.mlm = pobs.mlm,
                         pstr.mlm = pstr.mlm,
                         pdip.mlm = pdip.mlm,
                         shape.a = shape.a, shape.i = shape.i,
                         shape.d = shape.d,
                         byrow.aid = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitdlog





 rgaitdlog <-
  function(n, shape.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,
           shape.a = shape.p, shape.i = shape.p,
           shape.d = shape.p) {
    qgaitdlog(runif(n), shape.p,
              a.mix = a.mix,
              a.mlm = a.mlm,
              i.mix = i.mix,
              i.mlm = i.mlm,
              d.mix = d.mix,
              d.mlm = d.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix,
              pdip.mlm = pdip.mlm,
              shape.a = shape.a, shape.i = shape.i,
              shape.d = shape.d,
              byrow.aid = byrow.aid)
}  # rgaitdlog



















dlog <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(shape))
  if (length(x)     != N) x     <- rep_len(x,     N)
  if (length(shape) != N) shape <- rep_len(shape, N)
  ox <- !is.finite(x)
  zero <- ox | round(x) != x | x < 1
  ans <- rep_len(0.0, length(x))
  if (log.arg) {
    ans[ zero] <- log(0.0)
    ans[!zero] <- x[!zero] * log(shape[!zero]) - log(x[!zero]) -
                  log(-log1p(-shape[!zero]))
    ans[ox] <- log(0)  # 20141212 KaiH
  } else {
    ans[!zero] <- -(shape[!zero]^(x[!zero])) / (x[!zero] *
                   log1p(-shape[!zero]))
    ans[ox] <- 0.0  # 20141212 KaiH
  }
  ans[shape < 0 | 1 < shape] <- NaN
  ans
}



 plog <- function(q, shape, lower.tail = TRUE, log.p = FALSE) {

  if (any(is.na(q))) stop("NAs not allowed for argument 'q'")
  if (any(is.na(shape)))
    stop("NAs not allowed for argument 'shape'")


  N <- max(length(q), length(shape))
  if (length(q)     != N) q     <- rep_len(q,     N)
  if (length(shape) != N) shape <- rep_len(shape, N)





  bigno <- 10
  owen1965 <- (q * (1 - shape) > bigno)
  if (specialCase <- any(owen1965)) {
    qqq <- q[owen1965]
    ppp <- shape[owen1965]
    pqp <- qqq * (1 - ppp)
    bigans <- (ppp^(1+qqq) / (1-ppp)) * (1/qqq -
              1 / (            pqp * (qqq-1)) +
              2 / ((1-ppp)   * pqp * (qqq-1) * (qqq-2)) -
              6 / ((1-ppp)^2 * pqp * (qqq-1) * (qqq-2) * (qqq-3)) +
        24 / ((1-ppp)^3 * pqp * (qqq-1) * (qqq-2) * (qqq-3) * (qqq-4)))
      bigans <- 1 + bigans / log1p(-ppp)
  }

  floorq <- pmax(1, floor(q))  # Ensures at least 1 element per q value
  floorq[owen1965] <- 1
  seqq <- sequence(floorq)
  seqp <- rep(shape, floorq)
  onevector <- (seqp^seqq / seqq) / (-log1p(-seqp))
  rlist <-  .C("tyee_C_cum8sum",
               as.double(onevector), answer = double(N),
               as.integer(N), as.double(seqq),
               as.integer(length(onevector)), notok = integer(1))
  if (rlist$notok != 0)
    stop("error in C function 'cum8sum'")
  ans <- if (log.p) log(rlist$answer) else rlist$answer
  if (specialCase)
    ans[owen1965] <- if (log.p) log(bigans) else bigans
  ans[q < 1] <- if (log.p) log(0.0) else 0.0
  ans[shape < 0 | 1 < shape] <- NaN
  if (lower.tail) ans else 1 - ans
}



 qlog <- function(p, shape) {

  LLL <- max(length(p), length(shape))
  if (length(p)     < LLL) p     <- rep_len(p,     LLL)
  if (length(shape) < LLL) shape <- rep_len(shape, LLL)
  ans <- rep_len(0, LLL)

  lowsup <- 1
  lo <- rep_len(lowsup - 0.5, LLL)
  approx.ans <- lo  # True at lhs
  hi <- 2 * lo + 10.5
  dont.iterate <- p == 1 | shape <= 0 | 1 < shape
  done <- p <= plog(hi, shape) | dont.iterate
  while (!all(done)) {  # 20200307; bug fixed
    lo[!done] <- hi[!done]
    hi[!done] <- 2 * hi[!done] + 10.5
    done[!done] <- (p[!done] <= plog(hi[!done], shape[!done]))
  }

  foo <- function(q, shape, p)
    plog(q, shape) - p

  lhs <- (p <= dlog(1, shape)) | dont.iterate
  approx.ans[!lhs] <-
    bisection.basic(foo, lo[!lhs], hi[!lhs], tol = 1/16,
                    shape = shape[!lhs], p = p[!lhs])
  faa <- floor(approx.ans)
  ans <- ifelse(plog(faa, shape) < p & p <= plog(faa+1, shape),
                faa+1, faa)

  ans[p == 1] <- Inf
  ans[shape <= 0] <- NaN
  ans[1 < shape] <- NaN

  ans
}  # qlog



rlog <- function(n, shape) {
  qlog(runif(n), shape)
}






 logff <-
  function(lshape = "logitlink",
           gshape = -expm1(-7 * ppoints(4)),
           zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  new("vglmff",
  blurb = c("Logarithmic distribution f(y) = a * shape^y / y, ",
             "y = 1, 2, 3,...,\n",
             "            0 < shape < 1, a = -1 / log(1-shape)  \n\n",
             "Link:    ", namesof("shape", lshape, earg = eshape),
             "\n", "\n",
             "Mean:    a * shape / (1 - shape)", "\n"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 1
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         dpqrfun = "log",  # dlog, qlog, etc.
         link2name = "logff",  # logffMlink, logffMdlink,  logffQlink,
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <- namesof(mynames1, .lshape , earg = .eshape ,
                                tag = FALSE)


    if (!length(etastart)) {
      logff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
        sum(c(w) * dlog(x = y, shape = shapeval, log = TRUE))
      }
      Init.shape <- matrix(0, n, M)
      shape.grid <- .gshape

      for (ilocal in 1:ncoly) {
        Init.shape[, ilocal] <- grid.search(shape.grid,
                                            objfun = logff.Loglikfun,
                                            y = y[, ilocal],  # x = x,
                                            w = w[, ilocal])
      }  # for
      etastart <- theta2eta(Init.shape, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    A8 <- -1 / log1p(-shape)
    A8 * shape / (1 - shape)
  }, list( .lshape = lshape, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <- c(rep_len( .lshape , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }
  }), list( .lshape = lshape, .eshape = eshape ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlog(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("logff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    okay0 <- if ( .lshape == "logfflink") all(0 < eta) else TRUE
    okay1 <- if (okay0) {
      shape <- eta2theta(eta, .lshape , earg = .eshape )
      all(is.finite(shape)) && all(0 < shape & shape < 1)
    } else {
      FALSE
    }
    okay0 && okay1
  }, list( .lshape = lshape, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    rlog(nsim * length(shape), shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),


  deriv = eval(substitute(expression({
    M1 <- 1
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    A8 <- -1 / log1p(-shape)
    dl.dshape <- -A8 / (1 - shape) + y / shape
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- A8 * (1 - A8 * shape) / (shape * (1-shape)^2)
    wz <- c(w) * ned2l.dshape2 * dshape.deta^2
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # logff





















