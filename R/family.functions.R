# These functions are
# Copyright (C) 1998-2025 T.W. Yee, University of Auckland.
# All rights reserved.








fill1 <- fill2 <- fill3 <- fill4 <-
  function(x, values = 0, ncolx = ncol(x)) {
  x <- as.matrix(x)
  matrix(values, nrow(x), ncolx, byrow = TRUE)
}





extract.arg <- function(a) {
  s <- substitute(a)
  as.character(s)
}






remove.arg <- function(string) {

  nc <- nchar(string)
  bits <- substring(string, 1:nc, 1:nc)
  b1 <- (1:nc)[bits == "("]
  b1 <- if (length(b1)) b1[1] - 1 else nc
  if (b1 == 0)
    return("")
  string <- paste(bits[1:b1], collapse = "")
  string
}






add.arg <- function(string, arg.string) {

  if (arg.string == "")
    return(string)
  nc <- nchar(string)
  lastc <- substring(string, nc, nc)
  if (lastc == ")") {
    if (substring(string, nc - 1, nc - 1) == "(") {
      paste0(substring(string, 1, nc - 2),
             "(", arg.string, ")")
    } else {
      paste0(substring(string, 1, nc - 1),
             ", ", arg.string, ")")
    }
  } else {
    paste0(string, "(", arg.string, ")")
  }
}  # add.arg






get.arg <- function(string) {

  nc <- nchar(string)
  bits <- substring(string, 1:nc, 1:nc)
  b1 <- (1:nc)[bits == "("]
  b2 <- (1:nc)[bits == ")"]
  b1 <- if (length(b1)) min(b1) else return("")
  b2 <- if (length(b2)) max(b2) else return("")
  if (b2 - b1 == 1) "" else
    paste(bits[(1+b1):(b2-1)], collapse = "")
}  # get.arg







 eifun <- function(i, n)
    cbind(as.numeric((1:n) == i))



 eifun <-
 I.col <- function(i, n)
    diag(n)[, i, drop = FALSE]



 eijfun <- function(i, n) {
  temp <- matrix(0, n, 1)
  if (length(i))
    temp[i, ] <- 1
  temp
}






tapplymat1 <-
  function(mat, function.arg = c("cumsum", "diff", "cumprod")) {


  if (!missing(function.arg))
    function.arg <- as.character(substitute(function.arg))
  function.arg <- match.arg(function.arg,
                            c("cumsum", "diff", "cumprod"))[1]

  type <- switch(function.arg, cumsum = 1, diff = 2, cumprod = 3,
           stop("argument 'function.arg' not matched"))

  if (!is.matrix(mat))
    mat <- as.matrix(mat)
  NR <- nrow(mat)
  NC <- ncol(mat)
  fred <- .C("tapply_mat1", mat = as.double(mat), as.integer(NR),
             as.integer(NC), as.integer(type))  # , PACKAGE = "VGAM"
  dim(fred$mat) <- c(NR, NC)
  dimnames(fred$mat) <- dimnames(mat)
  switch(function.arg,
         cumsum  = fred$mat,
         diff    = fred$mat[, -1, drop = FALSE],
         cumprod = fred$mat)
}






matrix.power <- function(wz, M, power, fast = TRUE) {




  n <- nrow(wz)
  index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
  dimm.value <- if (is.matrix(wz)) ncol(wz) else 1
  if (dimm.value > M*(M+1)/2)
    stop("too many columns")


  if (M == 1 || dimm.value == M) {
    WW <- wz^power  # May contain NAs
    return(t(WW))
  }

  if (fast) {
    k <- veigen(t(wz), M = M)  # matrix.arg)
    evals  <- k$values   # M x n
    evects <- k$vectors  # M x M x n
  } else {
    stop("sorry, cannot handle matrix-band form yet")
    k <- unlist(apply(wz, 3, eigen), use.names = FALSE)
    dim(k) <- c(M, M+1, n)
    evals  <- k[,  1, , drop = TRUE]  # M x n
    evects <- k[, -1, , drop = TRUE]  # M x M x n
  }

  temp <- evals^power    # Some values may be NAs


  index <- as.vector( matrix(1, 1, M) %*% is.na(temp) )


  index <- (index == 0)
  if (!all(index)) {
    warning("Some weight matrices have negative ",
            "eigenvalues. They will be assigned NAs")
    temp[,!index] <- 1
  }

  WW <- mux55(evects, temp, M = M)
  WW[,!index] <- NA
  WW
}  # matrix.power






ResSS.vgam <- function(z, wz, M) {


  if (M == 1)
    return(sum(c(wz) * c(z^2)))

  wz.z <- mux22(t(wz), z, M = M, as.matrix = TRUE)
  sum(wz.z * z)
}  # ResSS.vgam






wweighted.mean <-
    function(y, w = NULL, matrix.arg = TRUE) {
  if (!matrix.arg)
    stop("currently, 'matrix.arg' must be TRUE")
  y <- as.matrix(y)
  M <- ncol(y)
  n <- nrow(y)
  if (M == 1) {
    if (missing(w))
      mean(y) else sum(w * y) / sum(w)
  } else {
    if (missing(w)) y %*% rep(1, n) else {
      numer <- mux22(t(w), y, M, as.matrix = TRUE)
      numer <- t(numer) %*% rep(1, n)
      denom <- t(w) %*% rep(1, n)
      denom <- matrix(denom, 1, length(denom))
      if (matrix.arg)
        denom <- m2a(denom, M = M)[, , 1]
      c(solve(denom, numer))
    }
  }
}  # wweighted.mean






veigen <- function(x, M) {


  n <- ncol(x)
  index <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
  dimm.value <- nrow(x)  # usually M or M(M+1)/2

  z <- .Fortran("veigenf",
      as.integer(M),
      as.integer(n),
      as.double(x),
      values = double(M * n),
      as.integer(1),
      vectors = double(M*M*n),
      double(M),
      double(M),
      wk = double(M*M),
      as.integer(index$row), as.integer(index$col),
      as.integer(dimm.value),
      error.code = integer(1))

  if (z$error.code)
    stop("eigen algorithm (rs) returned error code ", z$error.code)
  ord <- M:1
  dim(z$values) <- c(M, n)
  z$values <- z$values[ord, , drop = FALSE]
  dim(z$vectors) <- c(M, M, n)
  z$vectors <- z$vectors[, ord, , drop = FALSE]
  return(list(values  = z$values,
              vectors = z$vectors))
}  # veigen






ima <- function(j, k, M) {
  if (length(M) > 1 || M <= 0 || j <= 0 ||
      k <= 0 || j > M || k > M)
    stop("input wrong in ima()")
  m <- diag(M)
  m[col(m) <= row(m)] <- 1:(M*(M+1)/2)
  if (j >= k) m[j, k] else m[k, j]
}






checkwz <-
  function(wz, M, trace = FALSE,
           wzepsilon = .Machine$double.eps^0.75) {
  if (wzepsilon > 0.5)
    warning("argument 'wzepsilon' is probably too large")
  if (!is.matrix(wz))
    wz <- as.matrix(wz)
  wzsubset <- wz[, 1:M, drop = FALSE]
  if (any(is.na(wzsubset)))
    stop("NAs in the working weights var 'wz'")
  if (any(!is.finite(wzsubset)))
    stop("Some elements in the working weights",
         " variable 'wz' are not finite")

  if ((temp <- sum(wzsubset < wzepsilon)))
    warning(temp, " diagonal elements of the ",
      "working weights variable 'wz' have ",
       "been replaced by ", signif(wzepsilon, 5))
  wz[, 1:M] <- pmax(wzepsilon, wzsubset)
  wz
}  # checkwz








label.cols.y <-
  function(answer,
           colnames.y = NULL,
           NOS = 1,
           percentiles = c(25, 50, 75),
           one.on.one = TRUE,
           byy = TRUE) {
  if (!is.matrix(answer))
    answer <- as.matrix(answer)

  if (one.on.one) {
    colnames(answer) <-
      if (length(colnames.y) == ncol(answer))
        colnames.y else NULL
    return(answer)
  }



  if (is.null(percentiles))
    percentiles <- c(25, 50, 75)  # Restore to the default

  if (!is.Numeric(percentiles) ||
      min(percentiles) <= 0 ||
      max(percentiles) >= 100)
    stop("values of 'percentiles' should be in [0, 100]")

  percentiles <- signif(percentiles, digits = 5)

  ab1 <- rep(as.character(percentiles), length = ncol(answer))
  ab1 <- paste0(ab1, "%")
  if (NOS > 1) {
    suffix.char <- if (length(colnames.y) == NOS)
      colnames.y else as.character(1:NOS)
    ab1 <- paste0(ab1,
                  rep(suffix.char,
                      each = length(percentiles)))
  }
  colnames(answer) <- ab1


  if (byy) {
    answer <-
      answer[, interleave.VGAM(.M = NCOL(answer),
                               M1 = NOS),   # length(percentiles)),
             drop = FALSE]
  }
  answer
}  # label.cols.y






b3hfun <- 
  function(wz, dmat, w,
           b3h.wz = 0.5   # == 0.5 by default
           ) {
  if (!is.matrix(wz))   wz   <- cbind(wz)
  if (!is.matrix(w))    w    <- cbind(w)
  if (!is.matrix(dmat)) dmat <- cbind(dmat)
  
  M <- ncol(dmat)
  i5 <- iam(NA, NA, M, both = TRUE, diag = TRUE)
  ncwz <- ncol(wz)
  if (ncol(w) != 1 && ncol(w) != M)
    stop("ncol(w) must be 1 or M")
  ncwz <- ncol(wz)  # Prune i5
  ind.w <- rep(1:ncol(w), length = ncwz)
  i5$row.index <- i5$row.index[1:ncwz]
  i5$col.index <- i5$col.index[1:ncwz]
  dmat <- as.matrix(dmat) / c(w)  # Undo 'c(w) *'
  bhhhmat <-  # Unweighted xprod matrix
      dmat[, i5$row, drop = FALSE] *
      dmat[, i5$col, drop = FALSE]
  ans4 <- matrix(0, nrow(dmat), ncwz)
  for (jay in 1:ncwz)
    ans4[, jay] <-
      weighted.mean(bhhhmat[, jay],
                    w[, ind.w[jay]])
  ans4 <- c(w) * ans4

    if (length(b3h.wz)) {
      ans4 <- if (b3h.wz == 0) wz else
        if (b3h.wz == 1) ans4 else
        ans4 * b3h.wz + (1 - b3h.wz) * wz
    }  # length(oim.bhhh) && bhhh

  ans4
}  # b3hfun























