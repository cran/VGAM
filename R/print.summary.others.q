# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.





printsummary.lms <- function(x, digits = NULL, quote = TRUE, prefix = "")
{

    printsummary.vglm(x, digits = NULL, quote = TRUE, prefix = "")

    cat("\nLMS method!!!\n")

    cat("\nfirst 4 rows of $lin are\n")
    print.matrix(x$testing)

    invisible(NULL)
}



printsummary.rc.exponential <- function(x, digits = NULL, quote = TRUE, prefix = "")
{

    printsummary.vglm(x, digits = NULL, quote = TRUE, prefix = "")


    cat("\nNumber of censored observations: ", x$num.censored, "\n")

    invisible(NULL)
}

