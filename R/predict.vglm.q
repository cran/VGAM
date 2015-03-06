# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.





predictvglm <-
  function(object,
           newdata = NULL,
           type = c("link", "response", "terms"),  # "parameters",
           se.fit = FALSE,
           deriv = 0,
           dispersion = NULL,
           untransform = FALSE,
           extra = object@extra, ...) {
  na.act <- object@na.action
  object@na.action <- list()

  if (missing(extra)) {
  }

  if (deriv != 0)
    stop("'deriv' must be 0 for predictvglm()")

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("link", "response", "terms"))[1]

  if (untransform &&
     (type == "response" || type == "terms" || se.fit || deriv != 0))
    stop("argument 'untransform=TRUE' only if 'type=\"link\", ",
         "se.fit = FALSE, deriv=0'")




  pred <-
    if (se.fit) {
      switch(type,
             response = {
               warning("'type=\"response\"' and 'se.fit=TRUE' not valid ",
                       "together; setting 'se.fit = FALSE'")
               se.fit <- FALSE
               predictor <- predict.vlm(object, newdata = newdata,
                                        type = type, se.fit = se.fit,
                                        deriv = deriv, 
                                        dispersion = dispersion, ...) 
               fv <- object@family@linkinv(predictor, extra)


               fv <- as.matrix(fv)
               dn1 <- dimnames(fv)[[1]]
               dn2 <- dimnames(object@fitted.values)[[2]]
               if (nrow(fv) == length(dn1) &&
                   ncol(fv) == length(dn2))
                 dimnames(fv) <- list(dn1, dn2)


               fv
             },
             link = {
               predict.vlm(object, newdata = newdata,
                           type = "response", se.fit = se.fit,
                           deriv = deriv, dispersion = dispersion, ...) 
             },
             terms = {
               predict.vlm(object, newdata = newdata,
                           type = type, se.fit = se.fit,
                           deriv = deriv, dispersion = dispersion, ...) 
             })  # End of switch
  } else {
    if (is.null(newdata)) {
      switch(type, 
             link = object@predictors, 
             response = object@fitted.values,
             terms = {
                 predict.vlm(object, newdata = newdata,
                             type = type, se.fit = se.fit,
                             deriv = deriv, dispersion = dispersion, ...) 
             })
    } else {
      if (!(length(object@offset) == 1 && object@offset == 0))
        warning("zero offset used") 
      switch(type, 
             response = {




               predictor <- predict.vlm(object, newdata = newdata,
                                        type = type, se.fit = se.fit,
                                        deriv = deriv, 
                                        dispersion = dispersion, ...)



               M <- object@misc$M

               fv <- object@family@linkinv(predictor, extra)
               if (M > 1 && is.matrix(fv)) {

               fv <- as.matrix(fv)
               dn1 <- dimnames(fv)[[1]]
               dn2 <- dimnames(object@fitted.values)[[2]]
               if (nrow(fv) == length(dn1) &&
                   ncol(fv) == length(dn2))
                 dimnames(fv) <- list(dn1, dn2)
               } else {
               }
                 fv
               },
               link = {
                 predict.vlm(object, newdata = newdata,
                             type = "response", se.fit = se.fit,
                             deriv = deriv, dispersion = dispersion, ...)
               },
               terms = {
                 predict.vlm(object, newdata = newdata,
                             type = type, se.fit = se.fit,
                             deriv = deriv, dispersion = dispersion, ...) 
               })  # End of switch
        }
  }

  if (!length(newdata) && length(na.act)) {
    if (se.fit) {
      pred$fitted.values <- napredict(na.act[[1]], pred$fitted.values)
      pred$se.fit <- napredict(na.act[[1]], pred$se.fit)
    } else {
      pred <- napredict(na.act[[1]], pred)
    }
  }
  
  if (untransform) untransformVGAM(object, pred) else pred
} # predictvglm




setMethod("predict", "vglm", function(object, ...) 
  predictvglm(object, ...))






predict.rrvglm <- function(object, 
                          newdata = NULL, 
                          type = c("link", "response", "terms"),
                          se.fit = FALSE, 
                          deriv = 0,
                          dispersion = NULL, 
                          extra = object@extra, ...) {

  if (se.fit) {
    stop("20030811; predict.rrvglm(..., se.fit=TRUE) not complete yet") 
    pred <- 
    switch(type,
           response = {
             warning("'type=\"response\"' and 'se.fit=TRUE' not valid ",
                     "together; setting 'se.fit = FALSE'")
             se.fit <- FALSE
               predictor <- predict.vlm(object, newdata = newdata,
                                        type = type, se.fit = se.fit,
                                        deriv = deriv, 
                                        dispersion = dispersion, ...) 
             fv <- object@family@linkinv(predictor, extra)


             fv <- as.matrix(fv)
             dn1 <- dimnames(fv)[[1]]
             dn2 <- dimnames(object@fitted.values)[[2]]
             if (nrow(fv) == length(dn1) &&
                 ncol(fv) == length(dn2))
               dimnames(fv) <- list(dn1, dn2)


             fv
           },
           link = {
             type <- "response"
             predict.vlm(object, newdata = newdata,
                         type = type, se.fit = se.fit,
                         deriv = deriv, dispersion = dispersion, ...) 
           },
           terms = {
             predict.vlm(object, newdata = newdata,
                         type = type, se.fit = se.fit,
                         deriv = deriv, dispersion = dispersion, ...) 
           }
          )
  } else {
    return(predictvglm(object, newdata = newdata,
                       type = type, se.fit = se.fit,
                       deriv = deriv, 
                       dispersion = dispersion, ...))
  }

  na.act <- object@na.action

  if (!length(newdata) && length(na.act)) {
    if (se.fit) {
      pred$fitted.values <- napredict(na.act[[1]], pred$fitted.values)
      pred$se.fit <- napredict(na.act[[1]], pred$se.fit)
    } else {
      pred <- napredict(na.act[[1]], pred)
    }
  }

  pred
}


setMethod("predict", "rrvglm", function(object, ...) 
  predict.rrvglm(object, ...))



untransformVGAM <- function(object, pred) {
 

  M <- object@misc$M
  Links <- object@misc$link
  if (length(Links) != M && length(Links) != 1)
   stop("cannot obtain the link functions to untransform the object")

  upred <- pred
  earg <- object@misc$earg





  LINK <- object@misc$link  # link.names # This should be a character vector.
  EARG <- object@misc$earg  # This could be a NULL
  if (is.null(EARG))
    EARG <- list(theta = NULL)
  if (!is.list(EARG))
    stop("the 'earg' component of 'object@misc' must be a list")

  if (length(LINK) != M &&
      length(LINK) != 1)
    stop("cannot obtain the link functions to untransform 'object'")



  if (!is.character(LINK))
    stop("the 'link' component of 'object@misc' should ",
         "be a character vector")

  learg <- length(EARG)
  llink <- length(LINK)
  if (llink != learg)
    stop("the 'earg' component of 'object@misc' should ",
         "be a list of length ", learg)


  level1 <- length(EARG) > 3 &&
            length(intersect(names(EARG),
              c("theta", "inverse", "deriv", "short", "tag"))) > 3
  if (level1)
    EARG <- list(oneOnly = EARG)



  learg <- length(EARG)





  for (ii in 1:M) {
    TTheta <- pred[, ii]  # Transformed theta


    use.earg <- if (llink == 1) EARG[[1]] else EARG[[ii]]
    function.name <- if (llink == 1) LINK else LINK[ii]


    use.earg[["inverse"]] <- TRUE  # New
    use.earg[["theta"]] <- TTheta  # New
    Theta <- do.call(function.name, use.earg)





    upred[, ii] <- Theta
  }

  dmn2 <- if (length(names(object@misc$link)))
    names(object@misc$link) else {
      if (length(object@misc$parameters)) object@misc$parameters else NULL
  }
  dimnames(upred) <- list(dimnames(upred)[[1]], dmn2)
  upred
}






