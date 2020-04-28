# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.

















 dtheta.deta <-
  function(theta,
           link = "identitylink",
           earg = list(theta = theta,  # Needed
                       inverse = TRUE,  # 20150711: big change!!!!
                       deriv = 1,
                       short = TRUE,
                       tag = FALSE)) {

  function.name  <- link

  function.name2 <- attr(earg, "function.name")
  if (length(function.name2) && function.name != function.name2) {
    warning("apparent conflict in name of link function")
  }

  earg[["theta"]] <- theta  # New data

  if (length(earg$inverse))
    earg[["inverse"]] <- TRUE else
    earg$inverse <- TRUE

  earg[["deriv"]] <- 1  # New


  do.call(function.name, earg)
}  # dtheta.deta 






 d2theta.deta2 <-
  function(theta,
           link = "identitylink",
           earg = list(theta = theta,  # Needed
                       inverse = TRUE,  # 20150711: big change!!!!
                       deriv = 2,
                       short = TRUE,
                       tag = FALSE)) {


  function.name  <- link

  function.name2 <- attr(earg, "function.name")
  if (length(function.name2) && function.name != function.name2)
    warning("apparent conflict in name of link function in ",
            "D2theta.deta2()")

  earg[["theta"]] <- theta  # New data


  if (length(earg$inverse))
    earg[["inverse"]] <- TRUE else
    earg$inverse <- TRUE

  earg[["deriv"]] <- 2  # New

  do.call(function.name, earg)
}  # d2theta.deta2






 d3theta.deta3 <-
  function(theta,
           link = "identitylink",
           earg = list(theta = theta,
                       inverse = TRUE,
                       deriv = 3,
                       short = TRUE,
                       tag = FALSE)) {

  function.name  <- link
  earg[["theta"]] <- theta  # New data

  if (length(earg$inverse))
    earg[["inverse"]] <- TRUE else
    earg$inverse <- TRUE

  earg[["deriv"]] <- 3  # New
  do.call(function.name, earg)
}  # d3theta.deta3






 theta2eta <-
  function(theta,
           link = "identitylink",
           earg = list(theta = NULL)) {

  function.name  <- link

  function.name2 <- attr(earg, "function.name")
  if (length(function.name2) && function.name != function.name2)
    warning("apparent conflict in name of link function")

  earg[["theta"]] <- theta  # New data

  do.call(function.name, earg)
}  # theta2eta






 eta2theta <-
  function(theta,  # This is really eta.
           link = "identitylink",
           earg = list(theta = NULL),
           special.fun = "multilogitlink",
           delete.coln = TRUE  # Only for "multilogitlink"
           ) {


  orig.earg <- earg
  if (!is.list(earg))
    stop("argument 'earg' is not a list")

  level1 <- length(earg) > 3 &&
            length(intersect(names(earg),
              c("theta", "inverse", "deriv", "short", "tag"))) > 3

  if (level1)
    earg <- list(oneOnly = earg)






  llink <- length(link)

  if (llink != length(earg))
    stop("length of argument 'link' differs from ",
         "length of argument 'earg'")
  if (llink == 0)
    stop("length(earg) == 0 not allowed")
  if (llink == 1) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,


    if (is.list(earg[[1]]))
      earg <- earg[[1]]

    function.name  <- link  # First chance

    function.name2 <- attr(earg, "function.name")  # May be, e.g., NULL
    if (length(function.name2) && function.name != function.name2)
      warning("apparent conflict in name of link function")

    earg[["theta"]] <- theta  # New data

    earg[["inverse"]] <- TRUE  # New
    return(do.call(function.name, earg))
  }  # llink == 1





  if (!is.matrix(theta) &&
      length(theta) == length(earg))
    theta <- rbind(theta)

  vecTF <- link == special.fun 
  Ans <- NULL
  iii <- 1
  while (iii <= llink) {
    first.index <- last.index <- iii  # Ordinary case
    special.case <- vecTF[iii]  # && sum(vecTF) < length(vecTF)

    if (special.case) {
      next.i <- iii+1
      while (next.i <= llink) {
        if (vecTF[next.i]) {
          last.index <- next.i
          next.i <- next.i + 1
        } else {
          break
        }
      }  # while
    }  # special.case

    iii <- iii + last.index - first.index + 1  # For next time
    use.earg <- earg[[first.index]]
    use.earg[["inverse"]] <- TRUE  # New
    use.earg[["theta"]] <- theta[, first.index:last.index,
                                 drop = FALSE]  # New
    use.function.name <- link[first.index]  # "multilogitlink"

    if (first.index != last.index && special.case) {
      adjusted.M <- last.index - first.index + 1
      use.earg$M <- adjusted.M
    }
 
    Ans2 <- do.call(use.function.name, use.earg)
    if (special.case && special.fun == "multilogitlink" &&
        delete.coln)
      Ans2 <- Ans2[, -use.earg$refLevel, drop = FALSE]

    Ans <- cbind(Ans, Ans2)
  }  # while (iii <= llink)

  if (length(orig.earg) == ncol(Ans) &&
      length(names(orig.earg)) > 0 &&
      ncol(Ans) > 0)
    colnames(Ans) <- names(orig.earg)

  Ans
}  # eta2theta






 namesof <- function(theta,
                     link = "identitylink",
                     earg = list(tag = tag, short = short),
                     tag = FALSE,
                     short = TRUE) {

  funname.only <- strsplit(as.character(link), "(", fixed = TRUE)
  funname.only <- (funname.only[[1]])[1]
  link <- funname.only

  earg[["theta"]] <- as.character(theta)

  earg[["tag"]] <- tag
  earg[["short"]] <- short


  do.call(link, earg)
}  # namesof





link2list <- function(link
                      ) {

  ans <- link

  fun.name <- as.character(ans[[1]])


  big.list <- as.list(as.function(get(fun.name)))


  big.list[[length(big.list)]] <- NULL  # Kill the body of code





  t.index <- pmatch(names(ans[-1]), names(big.list))
  t.index
  if (anyNA(t.index))
    stop("in '", fun.name, "' could not match argument(s) ",
         paste('"', names(ans[-1])[is.na(t.index)], '"', sep = "",
               collapse = ", "))


  Big.list <- big.list
  trivial.call <- (length(t.index) == 0)
  if (!trivial.call) {
    Big.list[t.index] <- ans[-1]
  }


  attr(Big.list, "function.name") <- fun.name


  Big.list
}  # link2list




