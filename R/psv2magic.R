# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.




 psv2magic <-
    function(x.VLM, constraints, spar.vlm, sm.osps.list) {




  colperm <- function(x, from, to) {

    ncx <- ncol(x)
    if (length(from) != length(to) ||
        any(from != round(from)) ||
        any(from < 1 | ncx < from) ||
        any(duplicated(from)) ||
        any(sort(from) != sort(to)))
      stop("invalid column permutation indices")
    perm <- seq_len(ncx)
    perm[to] <- perm[from]
    x[, perm]
  }



  assignx <- sm.osps.list$assignx
  nassignx <- names(assignx)
  indexterms <- sm.osps.list$indexterms
  which.X.sm.osps <- sm.osps.list$which.X.sm.osps
  term.labels <- sm.osps.list$term.labels
  ncol.X.sm.osps <- sapply(which.X.sm.osps, length)
  ncolHlist.model <- unlist(lapply(constraints, ncol))


  ncolHlist.new <- ncolHlist.model
  if (names(constraints)[[1]] == "(Intercept)") {
    ncolHlist.new <- ncolHlist.new[-1]
    nassignx <- nassignx[-1]
  }


  ncol.H.ps <- ncolHlist.new[indexterms]
  num.osps.terms <- length(which.X.sm.osps)


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
        matrix(jjoffset + 1:(ncol.X.sm.osps[jay] * ncol.H.ps[jay]),
               nrow = ncol.X.sm.osps[jay],  # Redundant really
               ncol = ncol.H.ps[jay], byrow = TRUE)
      jjoffset <- jjoffset + ncol.H.ps[[jay]] * ncol.X.sm.osps[[jay]]
    } else {
      jjoffset <- jjoffset + ncolHlist.new[ii] * ncol.model[ii]
    }
  }  # for ii
  vindex.min <- sapply(perm.list, min)  # function(x) min(x)
  vindex.max <- sapply(perm.list, max)  # function(x) max(x)
  oo1 <- vector("list", length(ncol.H.ps))  # list()
  for (ii in seq_along(ncol.H.ps)) {
    oo1[[ii]] <- seq.int(vindex.min[ii], vindex.max[ii])
  }
  ooo <- unlist(oo1, use.names = FALSE)  # do.call("c", oo1)
  ppp <- unlist(perm.list, use.names = FALSE)  # do.call("c", perm.list)


  OFF.list <- vector("list", num.osps.terms)  # list()
  for (ii in 1:num.osps.terms) {
    index <- 0
    OFF.list[[ii]] <- numeric()
    for (jay in 1:(ncol.H.ps[ii])) {
      OFF.list[[ii]][jay] <- vindex.min[ii] + index
      index <- ncol.X.sm.osps[ii] * jay
    }
  }


  list(x.VLM.new = if (identical(ppp, ooo)) x.VLM else
                   colperm(x.VLM, ppp, ooo),
       sp = unlist(spar.vlm),
       S.arg = rep(sm.osps.list$S.arg, ncol.H.ps),  # Argument 'S' of magic()
       OFF = unlist(OFF.list))
}  # psv2magic


