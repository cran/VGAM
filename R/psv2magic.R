# These functions are
# Copyright (C) 1998-2016 T.W. Yee, University of Auckland.
# All rights reserved.




 psv2magic <-
    function(x.VLM, constraints, lambda.vlm, ps.list) {




  colperm <- function(x, from, to) {
    ncx <- ncol(x)
    if (length(from) != length(to) ||
        any(from != round(from)) ||
        any(from < 1 | from > ncx) ||
        any(duplicated(from)) ||
        any(sort(from) != sort(to)))
      stop("invalid column permutation indices")
    perm <- seq(length = ncx)
    perm[to] <- perm[from]
    x[, perm]
  }



  assignx <- ps.list$assignx
  nassignx <- names(assignx)
  indexterms <- ps.list$indexterms
  which.X.ps <- ps.list$which.X.ps
  term.labels <- ps.list$term.labels
  ncol.X.ps <- sapply(which.X.ps, length)
  ncolHlist.model <- unlist(lapply(constraints, ncol))


  ncolHlist.new <- ncolHlist.model
  if (names(constraints)[[1]] == "(Intercept)") {
    ncolHlist.new <- ncolHlist.new[-1]
    nassignx <- nassignx[-1]
  }


  ncol.H.ps <- ncolHlist.new[indexterms]
  num.ps.terms <- length(which.X.ps)


  allterms <- length(term.labels)
  ncol.allterms <- sapply(assignx, length)

  ncol.model <- if (names(constraints)[[1]] == "(Intercept)")
                ncol.allterms[-1] else ncol.allterms
  jay <- 0
  jjoffset <- if (names(constraints)[[1]] == "(Intercept)")
              ncolHlist.model[1] else 0
  perm.list <- list()
  for (ii in seq_along(term.labels)) {
    if (indexterms[ii]) {
      jay <- jay + 1
      perm.list[[jay]] <-
        matrix(jjoffset + 1:(ncol.X.ps[jay] * ncol.H.ps[jay]),
               ncol = ncol.H.ps[jay], byrow = TRUE)
      jjoffset <- jjoffset +  ncol.H.ps[[jay]] * ncol.X.ps[[jay]]
    } else {
      jjoffset <- jjoffset + ncolHlist.new[ii] * ncol.model[ii]
    }
  }
  vindex.min <- sapply(perm.list, min)  # function(x) min(x)
  vindex.max <- sapply(perm.list, max)  # function(x) max(x)
  o1 <- vector("list", length(ncol.H.ps))  # list()
  for (ii in seq_along(ncol.H.ps)) {
    o1[[ii]] <- vindex.min[ii]:vindex.max[ii]
  }
  ooo <- unlist(o1)  # do.call("c", o1)
  ppp <- unlist(perm.list)  # do.call("c", perm.list)


  off.list <- vector("list", num.ps.terms)  # list()
  for (ii in 1:num.ps.terms) {
    index <- 0
    off.list[[ii]] <- numeric()
    for (jay in 1: ncol.H.ps[ii]) {
      off.list[[ii]][jay] <- vindex.min[ii] + index
      index <- ncol.X.ps[ii] * jay
    }
  }

  rl <-
    list(x.VLM.new = colperm(x.VLM, ppp, ooo),
         sp = unlist(lambda.vlm),
         S.arg = rep(ps.list$S.arg, ncol.H.ps),  # Argument 'S' of magic()
         off = unlist(off.list))
  rl
}


