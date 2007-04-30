# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.





summary.lms <- function(object, ...)
{
    ans <- NextMethod("summary")

    ans$testing <- object$lin[1:4,]

    class(ans) <- c("summary.lms", class(ans))


    ans
}



summary.rc.exponential <- function(object, ...)
{

    ans <- NextMethod("summary")

    ans$num.censored <- attr(object$terms, "num.censored")

    class(ans) <- c("summary.rc.exponential", class(ans))


    ans
}


