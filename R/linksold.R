# These functions are
# Copyright (C) 1998-2020 T.W. Yee, University of Auckland.
# All rights reserved.








ToString <- function(x)
  paste(x, collapse = ", ")
















 # loglink <-
 loge <-
  function(theta,
           bvalue = NULL,  # .Machine$double.xmin is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {


  .Deprecated("loglink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("loglink(",  theta, ")", sep = "") else
        paste("log(",  theta, ")", sep = "")
    if (tag)
      string <- paste("Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(theta),
           "1" = theta,
           "2" = theta,
           "3" = theta,
           "4" = theta,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
       "0" = log(theta),
       "1" =  1 / theta,
       "2" = -1 / theta^2,
       "3" =  2 / theta^3,
       "4" = -6 / theta^4,
       stop("argument 'deriv' unmatched"))
  }
}  # loglink






 logneg <-
  function(theta,
           bvalue = NULL,  # .Machine$double.xmin is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {


      .Deprecated("logneglink")

      

  if (is.character(theta)) {
    string <- if (short)
        paste("logneglink(",  theta, ")",  sep = "") else
        paste( "log(-(",  theta, "))", sep = "")
    if (tag)
      string <- paste("Log negative:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = -exp(theta),
           "1" = theta,
           "2" = theta,
           "3" = theta,
           "4" = theta,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(-theta),
           "1" =  1 / theta,
           "2" = -1 / theta^2,
           "3" =  2 / theta^3,
           "4" = -6 / theta^4,
           stop("argument 'deriv' unmatched"))
  }
}  # logneglink






 logoff <-
  function(theta,
           offset = 0,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {

  .Deprecated("logofflink")

      
  if (!is.Numeric(offset))
    stop("bad input for argument 'offset'")

  if (is.character(theta)) {
    string <- if (short)
      paste("logofflink(", theta,
            ", offset = ", as.character(offset),
            ")", sep = "") else
      paste("log(",
            as.character(offset),
            "+",
            as.char.expression(theta),
            ")", sep = "")
    if (tag)
      string <- paste("Log with offset:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(theta) - offset,
           "1" = theta + offset,
           "2" = theta + offset,
           "3" = theta + offset,
           "4" = theta + offset,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(theta + offset),
           "1" =  1 / (theta + offset),
           "2" = -1 / (theta + offset)^2,
           "3" =  2 / (theta + offset)^3,
           "4" = -6 / (theta + offset)^4,
           stop("argument 'deriv' unmatched"))
  }
}  # logofflink






if (FALSE)
 identitylink <-
  function(theta,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- theta
    if (tag)
      string <- paste("Identity:", string)
    return(string)
  }

  switch(as.character(deriv),
         "0" = theta,
         "1" = theta * 0 + 1,
         "2" = theta * 0,  # zz Does not handle Inf and -Inf
         "3" = theta * 0,  # zz Does not handle Inf and -Inf
         "4" = theta * 0,  # zz Does not handle Inf and -Inf
         "5" = theta * 0,  # zz Does not handle Inf and -Inf
         stop("argument 'deriv' unmatched"))
}  # identitylink






 negidentity <-
  function(theta,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("negidentitylink")

      
  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste("-", theta, sep = "")
    if (tag)
      string <- paste("Negative-identity:", string)
    return(string)
  }

  switch(as.character(deriv),
         "0" = -theta,
         "1" = theta * 0 - 1,
         "2" = theta * 0,  # zz Does not handle Inf and -Inf
         "3" = theta * 0,  # zz Does not handle Inf and -Inf
         "4" = theta * 0,  # zz Does not handle Inf and -Inf
         "5" = theta * 0,  # zz Does not handle Inf and -Inf
         stop("argument 'deriv' unmatched"))
}  # negidentitylink






 logit <-
  function(theta,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("logitlink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("logitlink(",  # "logit(",
               theta,
              ")", sep = "") else
        paste("log(",
              as.char.expression(theta),
              "/(1-",
              as.char.expression(theta),
              "))", sep = "")
    if (tag)
      string <- paste("Logit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }
  if (inverse) {
    switch(as.character(deriv),
        "0" = plogis(theta),
        "1" =    1 / Recall(theta = theta,
                      bvalue = bvalue,
                      inverse = FALSE, deriv = deriv),
        "2" = theta * (1 - theta) * (1 - 2 * theta),
        "3" = (1 - 6 * theta * (1 - theta)) *
              theta * (1 - theta),
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
       "0" = qlogis(theta),
       "1" = 1 / (theta * (1 - theta)),
       "2" = (2 * theta - 1) / (theta * (1 - theta))^2,
       "3" = 2 * (1 - 3 * theta * (1 - theta)) / (theta * (1 - theta))^3,
       stop("argument 'deriv' unmatched"))
  }
}  # logitlink






 loglog <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {



  .Deprecated("logloglink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("logloglink(",  theta, ")",  sep = "") else
        paste("log(log(", theta, "))", sep = "")
    if (tag)
      string <- paste("Log-Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 1.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(exp(theta)),
           "1" = (theta * log(theta)),
           "2" = { junk <- log(theta)
                   theta  * junk * (1 + junk) },
           "3" = { Junk <- theta * log(theta)
                   Junk * ((1 + log(theta))^2 + Junk / theta)
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(log(theta)),
           "1" = 1 / (theta * log(theta)),
           "2" = { junk <- log(theta)
                   -(1 + junk) / (theta * junk)^2
           },
           "3" = { Junk <- theta * log(theta)
                   (2 * (1 + log(theta))^2 / Junk - 1 / theta) / Junk^2
                 },
           stop("argument 'deriv' unmatched"))
  }
}  # logloglink





if (FALSE)
 loglogloglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste("loglogloglink(",  theta, ")",  sep = "") else
        paste("log(log(log(", theta, ")))", sep = "")
    if (tag)
      string <- paste("Log-Log-Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= exp(1.0)] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(exp(exp(theta))),
           "1" = theta * log(theta) * log(log(theta)),
           "2" = { junk <- log(theta)
                   logjunk <- log(junk)
                   theta * junk * logjunk * (1 + logjunk * (1 + junk))
                 },
           "3" = { junk <- log(theta)
                   logjunk <- log(junk)
                   theta * junk^2 * logjunk^3  * (
                   3 + junk + 1 / junk + 3 / logjunk +
                   3 / (junk * logjunk) +
                   1 / (junk * logjunk^2))
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(log(log(theta))),
           "1" = 1 / (theta * log(theta) * log(log(theta))),
           "2" = { junk <- log(theta)
                   logjunk <- log(junk)
                   (-1 / (theta^2 * junk * logjunk)) *
                   (1 + (1 / junk) * (1 + 1 / logjunk))
           },
           "3" = { junk <- log(theta)
                   logjunk <- log(junk)
                   (3 + 2 * junk + 2 / junk +
                   3 / logjunk +
                   3 / (junk * logjunk) +
                   2 / (junk * logjunk^2)) / (
                   theta^3 * junk^2 * logjunk)
                 },
           stop("argument 'deriv' unmatched"))
  }
}  # loglogloglink









 cloglog <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("clogloglink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("clogloglink(",    theta, ")",  sep = "") else
        paste("log(-log(1-",
              as.char.expression(theta),
              "))", sep = "")
    if (tag)
      string <- paste("Complementary log-log:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = { -expm1(-exp(theta)) },
           "1" = -((1 - theta) * log1p(-theta)),
           "2" = { junk <- log1p(-theta)
                   -(1 - theta) * (1 + junk) * junk },
           "3" = {
             junk <- log1p(-theta)
             Junk <- (1 - theta) * junk
             -Junk * (Junk / (1 - theta) + (1 + junk)^2)
           },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(-log1p(-theta)),
           "1" = -1 / ((1 - theta) * log1p(-theta)),
           "2" = {  junk <- log1p(-theta)
               -(1 + junk) / ((1 - theta) * junk)^2
           },
           "3" = {
             junk <- log1p(-theta)
             Junk <- (1 - theta) * junk
             (1 / (1 - theta) - 2 * (1 + junk)^2 / Junk) / Junk^2
           },
           stop("argument 'deriv' unmatched"))
  }
}  # clogloglink






 probit <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("probitlink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("probitlink(", theta, ")", sep = "") else
        paste("qnorm(",  theta, ")", sep = "")
    if (tag)
      string <- paste("Probit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1 - bvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
                   ans <- pnorm(theta)
                   if (is.matrix(theta))
                     dim(ans) <- dim(theta)
                   ans
                  },
           "1" = {  # 1st deriv
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
         },
           "2" = {  # 2nd deriv
        Junk <- qnorm(theta)
        ans <- -Junk * dnorm(Junk)
        if (is.vector(theta)) ans else
        if (is.matrix(theta)) {
          dim(ans) <- dim(theta)
          ans
        } else {
          warning("can only handle vectors and matrices;",
                  " converting to vector")
          ans
        }
      },
           "3" = {
             Junk <- qnorm(theta)
             junk <- dnorm(Junk)
             junk * (Junk^2 - 1)
           },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
        ans <- qnorm(theta)
        if (is.matrix(theta))
            dim(ans) <- dim(theta)
        ans
     },
           "1" = {  # 1st deriv
       if (is.matrix(theta)) {
         ans <- 1 / dnorm(qnorm(theta))
         dim(ans) <- dim(theta)
         ans
       } else {
         1 / dnorm(qnorm(as.vector(theta)))
       }
      },
           "2" = {  # 2nd deriv
        Junk <- qnorm(theta)

        ans <- Junk / (dnorm(Junk))^2

        if (is.vector(theta)) ans else
        if (is.matrix(theta)) {
          dim(ans) <- dim(theta)
          ans
        } else {
          warning("can only handle vectors and matrices;",
                  " converting to vector")
          ans
        }
      },
           "3" = {
             Junk <- qnorm(theta)
             junk <- dnorm(Junk)
             (1 + 2 * Junk^2) / junk^3
           },
           stop("argument 'deriv' unmatched"))
  }
}  # probitlink






if (FALSE)
 explink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste("explink(", theta, ")", sep = "") else
        paste("exp(", theta, ")", sep = "")
    if (tag)
      string <- paste("Exp:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (inverse) {
    switch(as.character(deriv),
           "0" = log(theta),
           "1" =      exp(    -theta),
           "2" = -    exp(-2 * theta),  # 20170610 Fixes up a bug
           "3" =  2 * exp(-3 * theta),
           "4" = -6 * exp(-4 * theta),
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = exp(theta),
           "1" = exp(theta),
           "2" = exp(theta), 
           "3" = exp(theta),
           "4" = exp(theta),
           stop("argument 'deriv' unmatched"))
  }
}  # explink






 reciprocal <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("reciprocallink")

      
  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste("1/", theta, sep = "")
    if (tag)
      string <- paste("Reciprocal:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta == 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = 1 / theta,
           "1" = -   theta^2,
           "2" =  2 * theta^3,
           "3" = -6 * theta^4,
           "4" = 24 * theta^5,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" =  1 / theta,
           "1" = -1 / theta^2,
           "2" =  2 / theta^3,
           "3" = -6 / theta^4,
           "4" = 24 / theta^5,
           stop("argument 'deriv' unmatched"))
  }
}  # reciprocallink






 negloge <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("negloglink")

      
  if (is.character(theta)) {
      string <- if (short)
          paste("negloglink(", theta, ")", sep = "") else
          paste("-log(",  theta, ")", sep = "")
      if (tag)
        string <- paste("Negative log:", string)
      return(string)
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(-theta),
           "1" = -theta,
           "2" =  theta,
           "3" = -theta,
           "4" =  theta,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = -log(theta),
           "1" = -1/theta,
           "2" =  1/theta^2,
           "3" = -2/theta^3,
           "4" =  6/theta^4,
           stop("argument 'deriv' unmatched"))
  }
}  # negloglink






 negreciprocal <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE,
           deriv = 0, short = TRUE, tag = FALSE) {
  .Deprecated("negreciprocallink")

      
  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste("-1/", theta, sep = "")
    if (tag)
      string <- paste("Negative reciprocal:", string)
    return(string)
  }


  if (!inverse && length(bvalue))
    theta[theta == 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = -1 / theta,
           "1" =      theta^2,
           "2" =  2 * theta^3,
           "3" =  6 * theta^4,
           "4" = 24 * theta^5,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" =  -1 / theta,
           "1" =   1 / theta^2,
           "2" =  -2 / theta^3,
           "3" =   6 / theta^4,
           "4" = -24 / theta^5,
           stop("argument 'deriv' unmatched"))
  }
}  # negreciprocallink






if (FALSE)
  igcanlink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {

  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste("-1/", theta, sep = "")
    if (tag)
      string <- paste("Negative inverse:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" =  1 / sqrt(-2*theta),
           "1" =      theta^3,
           "2" =  3 * theta^5, 
           "3" = 15 * theta^7, 
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" =  -1 / (2 * theta^2),
           "1" =   1 / theta^3,
           "2" =  -3 / theta^4,
           "3" =  12 / theta^5,
           "4" = -60 / theta^6,
           stop("argument 'deriv' unmatched"))
  }
}  # igcanlink






 rhobit <-
  function(theta,
           bminvalue = NULL,
           bmaxvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("rhobitlink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("rhobitlink(", theta, ")", sep = "") else
        paste("log((1+",
              as.char.expression(theta),
              ")/(1-",
              as.char.expression(theta),
              "))", sep = "")
    if (tag)
      string <- paste("Rhobit:", string)
    return(string)
  }

  if (!inverse) {
   if (length(bminvalue)) theta[theta <= -1.0] <- bminvalue
   if (length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = { junk <- exp(theta)
                   expm1(theta) / (junk + 1.0) },
           "1" = (1 - theta^2) / 2,
           "2" = (-theta / 2) * (1 - theta^2),
           "3" = (3 * theta^2 - 1) * (1 - theta^2) / 4,
             stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
                 log1p(theta) - log1p(-theta) },
           "1" = 2 / (1 - theta^2),
           "2" = (4*theta) / (1 - theta^2)^2,
           "3" = 4 * (1 + 3 * theta^2) / (1 - theta^2)^3,
           stop("argument 'deriv' unmatched"))
  }
}  # rhobitlink






 fisherz <-
  function(theta,
           bminvalue = NULL,
           bmaxvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("fisherzlink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("fisherzlink(", theta, ")", sep = "") else
        paste("(1/2) * log((1+",
              as.char.expression(theta),
              ")/(1-",
              as.char.expression(theta),
              "))", sep = "")
    if (tag)
      string <- paste("Fisher's Z transformation:", string)
    return(string)
  }

  if (!inverse) {
     if (length(bminvalue)) theta[theta <= -1.0] <- bminvalue
     if (length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = tanh(theta),
           "1" = 1 - theta^2,
           "2" = 2 * (-theta) * (1 - theta^2),
           "3" = (3 * theta^2 - 1) * (1 - theta^2) * 2,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = atanh(theta),
           "1" = 1 / (1 - theta^2),
           "2" = (2*theta) / (1 - theta^2)^2,
           "3" = 2 * (1 + 3 * theta^2) / (1 - theta^2)^3,
           stop("argument 'deriv' unmatched"))
    }
}  # fisherzlink






 multilogit <-
  function(theta,
           refLevel = "(Last)",
           M = NULL,  # stop("argument 'M' not specified"),
           whitespace = FALSE,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           all.derivs = FALSE,
           short = TRUE, tag = FALSE) {

  .Deprecated("multilogitlink")

      

  fillerChar <- ifelse(whitespace, " ", "")

  if (length(refLevel) != 1)
    stop("the length of argument 'refLevel' must be one")

  if (is.character(refLevel)) {
    if (refLevel != "(Last)")
      stop('if a character, refLevel must be "(Last)"')
    refLevel <- -1
  } else
  if (is.factor(refLevel)) {
    if (is.ordered(refLevel))
      warning("argument 'refLevel' is from an ordered factor")
    refLevel <- as.character(refLevel) == levels(refLevel)
    refLevel <- (seq_along(refLevel))[refLevel]
    if (!is.Numeric(refLevel, length.arg = 1,
                    integer.valued = TRUE, positive = TRUE))
      stop("could not coerce 'refLevel' into a single positive integer")
  } else
  if (!is.Numeric(refLevel, length.arg = 1,
                  integer.valued = TRUE))
    stop("'refLevel' must be a single (positive?) integer")




  if (is.character(theta)) {
    is.M <- is.finite(M) && is.numeric(M)
    string <- if (short) {
        paste("multilogitlink(", theta, ")", sep = "")
    } else {
        theta <- as.char.expression(theta)

         if (refLevel < 0) {
           ifelse(whitespace,
             paste("log(", theta, "[,j] / ",
                   theta, "[,",
                   ifelse(is.M, M+1, "M+1"),
                   "]), j = 1:",
                   ifelse(is.M, M, "M"), sep = ""),
             paste("log(", theta, "[,j]/",
                   theta, "[,",
                   ifelse(is.M, M+1, "M+1"),
                   "]), j=1:",
                   ifelse(is.M, M, "M"), sep = ""))
         } else {
             if (refLevel == 1) {
               paste("log(", theta, "[,", "j]",
                   fillerChar, "/", fillerChar,
                     "", theta, "[,", refLevel, "]), j",
                     fillerChar, "=", fillerChar, "2:",
                   ifelse(is.M, (M+1), "(M+1)"),
                     sep = "")
             } else {
               paste("log(", theta, "[,", "j]", fillerChar, "/",
                     "", theta, "[,", refLevel, "]), j",
                     fillerChar, "=", fillerChar,
                     "c(1:", refLevel-1, ",",
                     fillerChar,
                     refLevel+1, ":",
                     ifelse(is.M, (M+1), "(M+1)"),
                     ")", sep = "")
             }
         }
    }
    if (tag)
      string <- paste("Multinomial logit link:", string)
    return(string)
  }



  M.orig <- M
  M <- NCOL(theta) - !(inverse && deriv == 0)
  if (M < 1)
    ifelse(inverse,
           stop("argument 'eta' should have at least one column"),
           stop("argument 'theta' should have at least two columns"))
  if (is.numeric(M.orig) && M != M.orig) {
    warning("argument 'M' does not seem right but using it")
    M <- M.orig
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (!inverse && length(bvalue))
    theta[theta >= 1.0] <- 1 - bvalue



  foo <- function(eta, refLevel = -1, M) {
    phat <- if ((refLevel < 0) || (refLevel == M+1)) {
      care.exp2(cbind(eta, 0.0))
    } else if ( refLevel == 1) {
      care.exp2(cbind(0.0, eta))
    } else {
      use.refLevel <- if ( refLevel < 0) M+1 else refLevel
      etamat <- cbind(eta[, 1:( refLevel - 1), drop = FALSE],
                      0.0,
                      eta[, ( refLevel ):M, drop = FALSE])
      care.exp2(etamat)
    }
    ans <- phat / rowSums(phat)
    colnames(ans) <- NULL
    ans
  }  # foo


  if (inverse) {
    use.refLevel <- if (refLevel < 0) ncol(theta) else refLevel
    switch(as.character(deriv),
      "0" = {
      foo(theta, refLevel, M = M)  # log(theta[, -jay] / theta[, jay])
            },
      "1" = if (all.derivs) {
              index <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
              theta <- theta[, -use.refLevel, drop = FALSE]
              wz <- -theta[, index$row, drop = FALSE] *
                     theta[, index$col, drop = FALSE]
              wz[, 1:M] <- wz[, 1:M] + theta
              wz
            } else {
              theta[, -use.refLevel,  drop = FALSE] *
              theta[,  use.refLevel] / (
              theta[, -use.refLevel,  drop = FALSE] +
              theta[,  use.refLevel])
            },
      "2" = (theta*(1-theta)*(1-2*theta))[, -use.refLevel,  drop = FALSE],
      "3" = {
              temp1 <- theta * (1 - theta)
              (temp1 * (1 - 6 * temp1))[, -use.refLevel,  drop = FALSE]
            },
      stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
      ans <- if (refLevel < 0) {
        log(theta[, -ncol(theta)] / theta[, ncol(theta)])
      } else {
        use.refLevel <- if (refLevel < 0) ncol(theta) else refLevel
        log(theta[, -( use.refLevel )] / theta[, use.refLevel ])
      }
      colnames(ans) <- NULL
      ans
      },
      "1" = care.exp(-log(theta) - log1p(-theta)),
      "2" = (2 * theta - 1) / care.exp(2*log(theta) + 2*log1p(-theta)),
      "3" = {
        temp1 <- care.exp(log(theta) + log1p(-theta))
        2 * (1 - 3 * temp1) / temp1^3
      },
      stop("argument 'deriv' unmatched"))
  }
}  # multilogitlink










 foldsqrt <-
  function(theta,  #  = NA  , = NULL,
           min = 0, max = 1, mux = sqrt(2),
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("foldsqrtlink")

      
  if (!is.Numeric(min, length.arg = 1))
    stop("bad input for 'min' component")
  if (!is.Numeric(max, length.arg = 1))
    stop("bad input for 'max' component")
  if (!is.Numeric(mux, length.arg = 1, positive = TRUE))
    stop("bad input for 'mux' component")
  if (min >= max)
    stop("'min' >= 'max' is not allowed")

  if (is.character(theta)) {
    string <- if (short)
      paste("foldsqrtlink(", theta, ")", sep = "") else {
    theta <- as.char.expression(theta)
      if (abs(mux-sqrt(2)) < 1.0e-10)
        paste("sqrt(2*", theta, ") - sqrt(2*(1-", theta, "))",
              sep = "") else
      paste(as.character(mux),
            " * (sqrt(", theta, "-", min, ") - sqrt(",
            max, "-", theta, "))",
            sep = "")
    }
    if (tag)
      string <- paste("Folded square root:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
      mid <- (min + max) / 2
      boundary <- mux * sqrt(max - min)
      temp <- pmax(0, (theta/mux)^2 * (2*(max-min) - (theta/mux)^2))
      ans <- theta
      if (any(ind5 <- theta <  0))
        ans[ind5] <- mid - 0.5 * sqrt(temp[ind5])
      if (any(ind5 <- theta >= 0))
        ans[ind5] <- mid + 0.5 * sqrt(temp[ind5])
      ans[theta < -boundary] <- NA
      ans[theta >  boundary] <- NA
      ans
        },
       "1" = (2 / mux ) / (1/sqrt(theta-min) + 1/sqrt(max-theta)),
       "2" = stop("use the chain rule formula to obtain this"),
       "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = mux * (sqrt(theta-min) - sqrt(max-theta)),
           "1" = (1/sqrt(theta-min) + 1/sqrt(max-theta)) * mux / 2,
           "2" = -(mux / 4) * ((theta-min)^(-3/2) - (max-theta)^(-3/2)),
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
           },
           stop("argument 'deriv' unmatched"))
  }
}  # foldsqrtlink





if (FALSE)
 powerlink <-
  function(theta,
           power = 1,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
    exponent <- power
    if (exponent == 0)
      stop("use the 'loge' link")

  if (is.character(theta)) {
    string <- if (short)
        paste("powerlink(", theta, ", power = ",
              as.character(exponent), ")",
              sep = "") else
        paste(as.char.expression(theta),
              "^(", as.character(exponent), ")", sep = "")
    if (tag)
      string <- paste("Power link:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = theta^(1/exponent),
           "1" = (theta^(1-exponent)) / exponent,
           "2" = ((1-exponent) / exponent^2) * (theta^(1 - 2*exponent)),
      "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = theta^exponent,
           "1" = exponent / (theta^(1-exponent)),
           "2" = exponent * (exponent-1) * (theta^(exponent-2)),
      "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  }
}  # powerlink







 extlogit <-
  function(theta,
           min = 0, max = 1,
           bminvalue = NULL,
           bmaxvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("extlogitlink")

      

    A <- min
    B <- max
   if (!inverse && length(bminvalue)) theta[theta <= A] <- bminvalue
   if (!inverse && length(bmaxvalue)) theta[theta >= B] <- bmaxvalue



  if (is.character(theta)) {
    string <- if (short) {
      if (A != 0 || B != 1)
        paste("extlogitlink(", theta,
              ", min = ", A,
              ", max = ", B, ")", sep = "") else
        paste("extlogitlink(", theta, ")", sep = "")
    } else {
      paste("log((",
            as.char.expression(theta),
            "-min)/(max-",
            as.char.expression(theta),
            "))", sep = "")
    }
    if (tag)
      string <- paste("Extended logit:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
           junk <- care.exp(theta)
           (A + B * junk) / (1.0 + junk) },
           "1" = ((theta - A) * (B - theta)) / (B-A),
           "2" = (A + B - 2 * theta) * (theta - A) *
                 (B - theta) / (B-A)^2,
      "3" = { #3rd deriv
           (theta - A) * (B - theta) * ((2 * theta - A - B)^2 -
       2 * (theta - A) * (B - theta)) / (B - A)^3
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
           log((theta - A)/(B - theta))},
           "1" = (B-A) / ((theta - A) * (B - theta)),
           "2" = ((2 * theta - A - B) * (B-A)) / ((theta - A) *
                 (B - theta))^2,
      "3" = { #3rd deriv
           (B - A) * (2 / ((theta - A) * (B - theta))^2) *
           (1 + (2 * theta - A - B)^2 / ((theta - A) * (B - theta)))
      },
           stop("argument 'deriv' unmatched"))
  }
}  # extlogitlink






 logc <-
  function(theta,
           bvalue = NULL,  # .Machine$double.xmin is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("logclink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("logclink(", theta, ")", sep = "") else {
        theta <- as.char.expression(theta)
        paste("log(1-", theta, ")", sep = "")
    }
    if (tag)
      string <- paste("Log Complementary:", string)
    return(string)
  }


  if (!inverse && length(bvalue)) {
    theta[theta >= 1.0] <- bvalue;
  }
  if (inverse) {
    switch(as.character(deriv),
           "0" = -expm1(theta),
           "1" = theta - 1,
           "2" = theta - 1,
           "3" = theta - 1,
           "4" = theta - 1,
           "5" = theta - 1,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log1p(-theta),
           "1" =  -1 / (1 - theta),
           "2" =  -1 / (1 - theta)^2,
           "3" =  -2 / (1 - theta)^3,
           "4" =  -6 / (1 - theta)^4,
           "5" = -24 / (1 - theta)^5,
           stop("argument 'deriv' unmatched"))
  }
}  # logclink






 cauchit <-
  function(theta,
           bvalue = .Machine$double.eps,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("cauchitlink")

      
  if (is.character(theta)) {
    string <- if (short)
        paste("cauchitlink(", theta, ")", sep = "") else {
        theta <- as.char.expression(theta)
        paste("tan(pi*(", theta, "-0.5))", sep = "")
    }
    if (tag)
      string <- paste("Cauchit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }
  if (inverse) {
    switch(as.character(deriv),
           "0" = 0.5 + atan(theta) / pi,
           "1" = (cos(pi * (theta-0.5)))^2  / pi,
           "2" = {
             temp2 <- cos(pi * (theta-0.5))
             temp4 <- sin(pi * (theta-0.5))
             -2 * temp4 * temp2^3 / pi
           },
           "3" = {
             temp2 <- cos(pi * (theta-0.5))
             temp5 <- tan(pi * (theta-0.5))
             2 * temp2^6 * (3 * temp5^2 - 1) / pi
           },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = tan(pi * (theta-0.5)),
           "1" = pi / (cos(pi * (theta-0.5)))^2,
           "2" =  {
           temp2 <- cos(pi * (theta-0.5))
           temp3 <- tan(pi * (theta-0.5))
           (temp3 * 2 * pi^2) / temp2^2
         },
           "3" =  {
           temp2 <- cos(pi * (theta-0.5))
           temp3 <- tan(pi * (theta-0.5))
           2 * pi^3 * (1 + 3 * temp3^2) / temp2^2
         },
           stop("argument 'deriv' unmatched"))
  }
}  # cauchitlink







 golf <-
  function(theta,
           lambda = 1,
           cutpoint = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {



  .Deprecated("gordlink")

      


  if (!is.Numeric(lambda, positive = TRUE))
    stop('could not determine lambda or lambda has negative values')
  if (is.Numeric(cutpoint))
    if (any(cutpoint < 0) ||
        !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should contain ",
            "non-negative integer values")

  if (is.character(theta)) {
    string <- if (short) {
      lenl <- length(lambda) > 1
      lenc <- length(cutpoint) > 1
      paste("gordlink(", theta,
            ", lambda = ",
            if (lenl) "c(" else "",
            ToString(lambda),
            if (lenl) ")" else "",
            if (is.Numeric(cutpoint))
      paste(", cutpoint = ",
            if (lenc) "c(" else "",
      ToString(cutpoint),
            if (lenc) ")" else "",
      sep = "") else "",
                  ")", sep = "")
    } else {
      theta <- as.char.expression(theta)
      if (is.Numeric(cutpoint)) {
        paste("-3*log(1-qnorm(", theta,
              ")/(3*sqrt(lambda)))",
              " + log(cutpoint)", sep = "")
      } else {
        paste("-3*log(1-qnorm(", theta,
              ")/(3*sqrt(lambda)))", sep = "")
      }
    }
    if (tag)
      string <- paste("Gamma-ordinal link function:", string)
    return(string)
  }


  thmat <- cbind(theta)
  lambda <- rep_len(lambda, ncol(thmat))  # Allow recycling for lambda
  if (is.Numeric(cutpoint))
    cutpoint <- rep_len(cutpoint, ncol(thmat))
  if (ncol(thmat) > 1) {
    answer <- thmat
    for (ii in 1:ncol(thmat))
      answer[, ii] <- Recall(theta = thmat[, ii],
                             lambda = lambda[ii],
                             cutpoint = if (is.Numeric(cutpoint))
                                        cutpoint[ii] else NULL,
                            inverse = inverse, deriv = deriv)
    return(answer)
  }


  answer <- if (inverse) {
    switch(as.character(deriv),
           "0" = {
      if (is.Numeric(cutpoint)) {
        pnorm((1-care.exp(-(theta-log(cutpoint))/3)) * 3 * sqrt(lambda))
      } else {
        pnorm((1-care.exp(-theta/3)) * 3 * sqrt(lambda))
      }
    },

      "1" = 1 / Recall(theta = theta,
                 lambda = lambda,
                 cutpoint = cutpoint,
                 inverse = FALSE, deriv = deriv),
      "2" = stop('cannot currently handle deriv = 2',
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
    )
  } else {
    smallno <- 1 * .Machine$double.eps
    Theta <- theta
    Theta <- pmin(Theta, 1 - smallno)  # Since theta==1 is a possibility
    Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
    Ql <- qnorm(Theta)
    switch(as.character(deriv),
           "0" = {
        temp <- Ql / (3*sqrt(lambda))
        temp <- pmin(temp, 1.0 - smallno)  # 100 / .Machine$double.eps
        origans <- -3*log1p(-temp) +
        if (is.Numeric(cutpoint)) log(cutpoint) else 0
        1 / origans
      },
        "1" = {
        origans <- (1 - Ql / (3*sqrt(lambda))) * sqrt(lambda) * dnorm(Ql)
        1 / origans
      },
        "2" = {  stop('cannot currently handle deriv = 2') },
        "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  }
  if (!is.Numeric(answer))
    warning("the answer contains some NAs")
  answer
}  # gordlink, aka golf







 polf <-
  function(theta,  # = 1,
           cutpoint = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("pordlink")

      


  if (!is.Numeric(cutpoint))
    stop("could not determine the cutpoint")
  if (any(cutpoint < 0) ||
      !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should",
            " contain non-negative integer values")


  if (is.character(theta)) {
    string <- if (short) {
      lenc <- length(cutpoint) > 1
      paste("pordlink(", theta,
             ", cutpoint = ",
            if (lenc) "c(" else "",
            ToString(cutpoint),
            if (lenc) ")" else "",
            ")", sep = "")
    } else {
      theta <- as.char.expression(theta)
      paste("2*log(0.5*qnorm(", theta,
            ") + sqrt(cutpoint+7/8))", sep = "")
    }
    if (tag)
      string <- paste("Poisson-ordinal link function:", string)
    return(string)
  }



    thmat <- cbind(theta)
    if (ncol(thmat) > 1) {
        answer <- thmat
        cutpoint <- rep_len(cutpoint, ncol(thmat))
        for (ii in 1:ncol(thmat))
            answer[, ii] <- Recall(theta = thmat[, ii],
                                 cutpoint = cutpoint,
                                 inverse = inverse, deriv = deriv)
        return(answer)
    }

  answer <-
  if (inverse) {
    switch(as.character(deriv),
           "0" = {
 # deriv == 0
          origans <-
          if (any(cp.index <- cutpoint == 0)) {
              tmp <- theta
              tmp[cp.index] <-
              clogloglink(theta = theta[cp.index],
                      inverse = inverse, deriv = deriv)
              tmp[!cp.index] <-
                pnorm(2 * exp(theta[!cp.index]/2) -
                      2 * sqrt(cutpoint[!cp.index] + 7/8))
              tmp
          } else {
            pnorm(2 * exp(theta/2) - 2 * sqrt(cutpoint + 7/8))
          }
        1 / origans
      },
      "1" =           1 / Recall(theta = theta,
                     cutpoint = cutpoint,
                     inverse = FALSE, deriv = deriv),
      "2" =        stop('cannot currently handle deriv = 2'),
      "3" =        { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
             stop("argument 'deriv' unmatched"))

  } else {
    if (any(cp.index <- cutpoint == 0)) {
        clogloglink(theta = theta,
                inverse = inverse, deriv = deriv)


    } else {
      smallno <- 1 * .Machine$double.eps
      SMALLNO <- 1 * .Machine$double.xmin
      Theta <- theta
      Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is possible
      Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
      Ql <- qnorm(Theta)


    switch(as.character(deriv),
           "0" = {
      temp <- 0.5 * Ql + sqrt(cutpoint + 7/8)
      temp <- pmax(temp, SMALLNO)
      origans <- 2 * log(temp)
      1 / origans
    },
      "1" =  {
      origans <- (Ql/2 + sqrt(cutpoint + 7/8)) * dnorm(Ql)
      1 / origans
      },

      "2" = {  stop('cannot currently handle deriv = 2') },
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
      stop("argument 'deriv' unmatched"))
    }
  }
  if (!is.Numeric(answer))
    warning("the answer contains some NAs")
  answer
}  # pordlink, aka polf







 nbolf <-
  function(theta,
           cutpoint = NULL,
           k = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("nbordlink")

      





  kay <- k
  if (!is.Numeric(kay, positive = TRUE))
    stop("could not determine 'k' or it is not positive-valued")
  if (!is.Numeric(cutpoint))
    stop("could not determine the cutpoint")
  if (any(cutpoint < 0) ||
      !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should",
            " contain non-negative integer values")

  if (is.character(theta)) {
    string <- if (short) {
        lenc <- length(cutpoint) > 1
        lenk <- length(kay) > 1
        paste("nbordlink(", theta,
              ", cutpoint = ",
              if (lenc) "c(" else "",
              ToString(cutpoint),
              if (lenc) ")" else "",
              ", k = ",
              if (lenk) "c(" else "",
              ToString(kay),
              if (lenk) ")" else "",
              ")", sep = "")
      } else {
        theta <- as.char.expression(theta)
        paste("2*log(sqrt(k) * sinh(qnorm(", theta,
              ")/(2*sqrt(k)) + ",
              "asinh(sqrt(cutpoint/k))))", sep = "")
      }
      if (tag)
        string <- paste("Negative binomial-ordinal link function:",
                        string)
      return(string)
  }


  thmat <- cbind(theta)
  kay <- rep_len(kay, ncol(thmat))  # Allow recycling for kay
  cutpoint <- rep_len(cutpoint, ncol(thmat)) # Allow recycling 4 cutpoint
  if (ncol(thmat) > 1) {
    answer <- thmat
    for (ii in 1:ncol(thmat))
        answer[, ii] <- Recall(theta = thmat[, ii],
                             cutpoint = cutpoint[ii],
                             k = kay[ii],
                             inverse = inverse, deriv = deriv)
    return(answer)
  }

  answer <-
  if (inverse) {
    switch(as.character(deriv),
           "0" = {
      if (cutpoint == 0) {
        1.0 - (kay / (kay + care.exp(theta)))^kay
      } else {
          pnorm((asinh(exp(theta/2)/sqrt(kay)) -
                 asinh(sqrt(cutpoint/kay))) * 2 * sqrt(kay))
      }
       }, "0" =  {
      1 / Recall(theta = theta,
                 cutpoint = cutpoint,
                 k = kay,
                 inverse = FALSE, deriv = deriv)
    }, "0" = {
     stop('cannot currently handle deriv = 2')
   },
           "0" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
             stop("argument 'deriv' unmatched"))




  } else {
    smallno <- 1 * .Machine$double.eps
    SMALLNO <- 1 * .Machine$double.xmin
    Theta <- theta
    Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is possible
    Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
    if (cutpoint == 0) {
    switch(as.character(deriv),
           "0" = {
      temp <- (1 - Theta)^(-1/kay) - 1
      temp <- pmax(temp, SMALLNO)
      origans <- log(kay) + log(temp)
      1 / origans
    },
           "1" = {
      origans <- (kay / (1 - Theta)^(1/kay) - kay) *
          (1 - Theta)^(kay+1/kay)
      1 / origans
      },
      "2" = {  stop('cannot handle deriv = 2') },
      "3" = {  stop('cannot handle deriv = 2') },
      stop("argument 'deriv' unmatched"))
    } else {
      Ql <- qnorm(Theta)
    switch(as.character(deriv),
           "0" = {
            temp <- sqrt(kay) * sinh(Ql/(2*sqrt(kay)) +
                   asinh(sqrt(cutpoint/kay)))
            temp <- pmax(temp, SMALLNO)
            origans <- 2 * log(temp)
            1 / origans
          },
           "1" = {
            arg1 <- (Ql/(2*sqrt(kay)) + asinh(sqrt(cutpoint/kay)))
            origans <- sqrt(kay) * tanh(arg1) * dnorm(Ql)
            1 / origans
          },
           "2" = {  stop('cannot currently handle deriv = 2') },
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
            stop("argument 'deriv' unmatched"))
    }
  }
  if (!is.Numeric(answer))
    warning("the answer contains some NAs")
  answer
}  # nbordlink, aka nbolf






 nbolf2 <-
  function(theta,
           cutpoint = NULL,
           k = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  .Deprecated("nbord2link")

      

warning("20150711; this function has not been updated")




  kay <- k
  if (!is.Numeric(kay, positive = TRUE))
    stop("could not determine argument 'k' or ",
         "it is not positive-valued")
  if (!is.Numeric(cutpoint))
    stop("could not determine the cutpoint")
  if (any(cutpoint < 0) ||
      !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should ",
            "contain non-negative integer values")

  if (is.character(theta)) {
    string <- if (short) {
      lenc <- length(cutpoint) > 1
      lenk <- length(kay) > 1
      paste("nbord2link(", theta,
            ", earg = list(cutpoint = ",
            if (lenc) "c(" else "",
            ToString(cutpoint),
            if (lenc) ")" else "",
            ", k = ",
            if (lenk) "c(" else "",
            ToString(kay),
            if (lenk) ")" else "",
            "))", sep = "")
  } else {
    theta <- as.char.expression(theta)
    paste("3*log(<a complicated expression>)", sep = "")
  }
  if (tag)
    string <- paste("Negative binomial-ordinal link function 2:",
                   string)
  return(string)
  }


    thmat <- cbind(theta)
    kay <- rep_len(kay, ncol(thmat))  # Allow recycling for kay
    if (ncol(thmat) > 1) {
        answer <- thmat
        for (ii in 1:ncol(thmat))
            answer[, ii] <- Recall(theta = thmat[, ii],
                                 cutpoint = cutpoint[ii],
                                 k = kay[ii],
                                 inverse = inverse, deriv = deriv)
        return(answer)
    }

    answer <-
    if (inverse) {
        if (deriv > 0) {
            1 / Recall(theta = theta,
                       cutpoint = cutpoint,
                       k = kay,
                       inverse = FALSE, deriv = deriv)
        } else {
            if (cutpoint == 0) {
                1.0 - (kay / (kay + care.exp(theta)))^kay
            } else {

            a1 <- -(9*cutpoint+8) / (cutpoint+1)
            a2 <- (9*kay-1) / (kay * (cutpoint+1)^(1/3))
            a3 <- 9 / (kay * (cutpoint+1)^(2/3))
            a4 <- 9 / (cutpoint+1)
            B <- exp(theta/3)
            mymat <- rbind(a1^2*a2^2 + 2*a1*a2^3*B + B^2*a2^4, 0,
                    -2*a1*a2*a3*B - 2*a2^2*a3*B^2 - a1^2*a3 - a2^2*a4,
                           0,
                    B^2 * a3^2 + a3 * a4)
            ans <- Re(t(apply(mymat, 2, polyroot)))
            theta2 <- invfun <- pnorm(-ans)  # pnorm(-x) = 1-pnorm(x)
            for (ii in 1:4) {
              theta2[, ii] <-
                Recall(theta = theta2[, ii],
                       cutpoint = cutpoint,
                       k = kay,
                       inverse = FALSE, deriv = deriv)
            }
            rankmat <- t(apply(abs(theta2 - theta), 1, rank))
            for (ii in 2:4) {
              if (any(index4 <- (rankmat[, ii] == 1))) {
                invfun[index4, 1] <- invfun[index4, ii]
              }
            }
            invfun[, 1]
            }
        }
    } else {
        smallno <- 1 * .Machine$double.eps
        SMALLNO <- 1 * .Machine$double.xmin
        Theta <- theta
        Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is possible
        Theta <- pmax(Theta, smallno)  # Since theta == 0 is possible
        if (cutpoint == 0) {
    switch(as.character(deriv),
           "0" = {
            temp <- (1 - Theta)^(-1/kay) - 1
            temp <- pmax(temp, SMALLNO)
            log(kay) + log(temp)},
            "0" = (kay / (1 - Theta)^(1/kay) - kay) *
                  (1 - Theta)^(kay+1/kay),
            "0" = {  stop("cannot handle 'deriv = 2'") },
           "0" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
            stop("argument 'deriv' unmatched"))
        } else {
            Ql <- qnorm(Theta)
            a1 <- -(9*cutpoint+8) / (cutpoint+1)
            a2 <- (9*kay-1) / (kay * (cutpoint+1)^(1/3))
            a3 <- 9 / (kay * (cutpoint+1)^(2/3))
            a4 <- 9 / (cutpoint+1)
            discrim <- a1^2 * a3 + a2^2 * a4 - Ql^2 * a3 * a4
            denomin <- Ql^2 * a3 - a2^2
            numerat <- (a1*a2 - Ql * sqrt(discrim))
            argmax1 <- numerat / denomin
    switch(as.character(deriv),
           "0" = {
                argmax2 <- (a1*a2 + Ql * sqrt(discrim)) / denomin
                temp <- ifelse(argmax1 > 0, argmax1, argmax2)
                temp <- pmax(temp, SMALLNO)
                3 * log(temp)},
           "1" = {
                 BB <- (sqrt(discrim) - Ql^2 * a3 *
                       a4 / sqrt(discrim)) / dnorm(Ql)
                 CC <- 2 * Ql * a3 / dnorm(Ql)
                 dA.dtheta <- (-denomin * BB - numerat * CC) / denomin^2
                 argmax1 / (3 * dA.dtheta)
                },
                "2" = {  stop('cannot currently handle deriv = 2') },
                   "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
                stop("argument 'deriv' unmatched"))
      }
  }
  if (!is.Numeric(answer))
    warning("the answer contains some NAs")
  answer
}  #  nbord2link, aka nbolf2






 Cut <-
  function(y, breaks = c(-Inf, quantile(c(y), prob = (1:4)/4))) {
  y <- as.matrix(y)


  temp <- cut(y, breaks = breaks, labels = FALSE)
  temp <- c(temp)  # integer vector of integers
  if (anyNA(temp))
    warning("there are NAs")
  answer <- if (ncol(y) > 1) matrix(temp, nrow(y), ncol(y)) else temp
  if (ncol(y) > 1) {
    ynames <- dimnames(y)[[2]]
    if (!length(ynames))
      ynames <- paste("Y", 1:ncol(y), sep = "")
    xnames <- dimnames(y)[[1]]
    if (!length(xnames)) xnames = as.character(1:nrow(y))
    dimnames(answer) <- list(xnames, ynames)
  }
  attr(answer, "breaks") <- breaks
  answer
}  # Cut



 checkCut <- function(y) {
  if (!is.Numeric(y, positive = TRUE, integer.valued = TRUE))
    stop("argument 'y' must contain positive integers only")
  uy <- unique(y)
  L <- max(uy)
  oklevels <- 1:L
  if (L == 1)
    stop("only one unique value")
  for (ii in oklevels) {
    if (all(ii != uy))
      stop("there is no ", ii, " value")
  }
  TRUE
}  # checkCut






if (FALSE)
 nbcanlink <-
  function(theta,
           size = NULL,
           wrt.param = NULL,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    lastchars1 <- substr(theta, nchar(theta), nchar(theta))
    lastchars2 <- ifelse(nchar(theta) > 1,
                    substr(theta, nchar(theta) - 1, nchar(theta) - 1),
                    rep("", length(theta)))

    size.names <- rep("size", length(theta))
    dig1 <- lastchars1 %in% as.character(0:9)
    dig2 <- lastchars2 %in% as.character(0:9)
    size.names <- ifelse(dig1,
                         paste("size",            lastchars1, sep = ""),
                         size.names)
    size.names <- ifelse(dig2,
                         paste("size", lastchars2, lastchars1, sep = ""),
                         size.names)

    string <- if (short)
      paste("nbcanlink(", theta,
             ", ", theta, "(", size.names, ")",  # Added 20180803
             ")", sep = "") else {
      theta <- as.char.expression(theta)
      paste("log(", theta, " / (", theta, " + ", size.names, "))",
            sep = "")
    }
    if (tag)
      string <- paste("Nbcanlink:", string)
    return(string)
  }



  kmatrix <- size
  theta <- cbind(theta)
  kmatrix <- cbind(kmatrix)
  if (ncol(kmatrix) != ncol(theta))
    stop("arguments 'theta' and 'size' do not have ",
         "an equal number of cols")
  if (nrow(kmatrix) != nrow(theta))
    stop("arguments 'theta' and 'size' do not have ",
         "an equal number of rows")


  if (deriv > 0) {
    if (!(wrt.param %in% 1:2))
      stop("argument 'wrt.param' should be 1 or 2")
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
       ans <- (kmatrix / expm1(-theta))
       if (is.matrix(ans))
         dimnames(ans) <- NULL else
         names(ans) <- NULL
       ans
       },

        "1" = if (wrt.param == 1)
           (theta * (theta + kmatrix)) / kmatrix else
          -(theta + kmatrix),

       "2" = if (wrt.param == 1)
       (2 * theta + kmatrix) * theta *
           (theta + kmatrix) / kmatrix^2 else
        theta + kmatrix,
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      })
  } else {
    ans <-
    switch(as.character(deriv),
           "0" = log(theta / (theta + kmatrix)),

           "1" = if (wrt.param == 1) kmatrix / (theta * (theta +
                                     kmatrix)) else
                 -1 / (theta + kmatrix),

       "2" = if (wrt.param == 1)
             (2 * theta + kmatrix) *
             (-kmatrix) / (theta * (theta + kmatrix))^2 else
             1 / (theta + kmatrix)^2,
       "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      })

     if (is.matrix(ans))
       dimnames(ans) <- NULL else
       names(ans) <- NULL
     ans
  }
}  # nbcanlink








linkfun.vglm <- function(object, earg = FALSE, ...) {
  if (!any(slotNames(object) == "extra"))
    stop("cannot access the 'extra' slot of the object")
  if (!any(slotNames(object) == "misc"))
    stop("cannot access the 'misc' slot of the object")

  M <- npred(object)

  misc <- object@misc
  LINKS1 <- misc$link
  EARGS1 <- misc$earg

  extra <- object@extra
  LINKS2 <- extra$link
  EARGS2 <- extra$earg

  if (length(LINKS1) != M && length(LINKS2) != M) {
    if (LINKS1 != "multilogitlink" && LINKS2 != "multilogitlink")
      warning("the length of the 'links' component is not ", M)
  }

  if (length(LINKS1)) {
    if (earg) list(link = LINKS1, earg = EARGS1) else LINKS1
  } else {
    if (earg) list(link = LINKS2, earg = EARGS2) else LINKS2
  }
}  # linkfun.vglm



if (!isGeneric("linkfun"))
  setGeneric("linkfun", function(object, ...)
             standardGeneric("linkfun"))



setMethod("linkfun", "vglm", function(object, ...)
    linkfun.vglm(object, ...))







 if (FALSE)
 logitoffsetlink <-
  function(theta,
           offset = 0,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste("logitoffsetlink(",
               theta,
              ", ", offset[1],
              ")", sep = "") else
        paste("log(",
              as.char.expression(theta),
              "/(1-",
              as.char.expression(theta),
              ")",
              " - ", offset[1],
              ")", sep = "")
    if (tag)
      string <- paste("Logit-with-offset:", string)
    return(string)
  }




  if (inverse) {
    switch(as.character(deriv),
           "0" = {
           exp.eta <- exp(theta)
           (exp.eta + offset) / (1 + exp.eta + offset)
           },
           "1" = 1 / Recall(theta = theta,
                      offset = offset,
                      inverse = FALSE, deriv = deriv),
           "2" = theta * (1 - theta) * (1 - 2 * theta),
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
       temp2 <- log(theta / (1 - theta) - offset)
       temp2
       },
       "1" = 1 / ((1 - theta) * (theta - (1-theta) * offset)),
       "2" = (2 * (theta - offset * (1-theta)) - 1) / (
       (theta - (1-theta)*offset) * (1-theta))^2,
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
       stop("argument 'deriv' unmatched"))
  }
}  # logitoffsetlink








