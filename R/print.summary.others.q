# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.




if(FALSE)
printsummary.lms <- function(x, digits = NULL, quote = TRUE, prefix = "") {

    printsummary.vglm(x, digits = NULL, quote = TRUE, prefix = "")



    invisible(NULL)
}



printsummary.rc.exponential <- function(x, digits = NULL, quote = TRUE,
        prefix = "") {

    printsummary.vglm(x, digits = NULL, quote = TRUE, prefix = "")


    cat("\nNumber of censored observations: ", x$num.censored, "\n")

    invisible(NULL)
}

