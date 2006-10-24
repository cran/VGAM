# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.



family.vglm <- function(object, ...) 
    object$vfamily

print.vfamily <- function(x, ...)
{
    f <- x$vfamily
    if(is.null(f))
        stop("not a VGAM family function")

    nn <- x$blurb
    if(is.null(nn))
        invisible(return(x))

    cat("Family: ", f[1], "\n") 
    if(length(f)>1) cat("Classes:", paste(f, collapse=", "), "\n")
    cat("\n")

    for(i in 1:length(nn))
        cat(nn[i])
    cat("\n")
    invisible(return(x))
}



