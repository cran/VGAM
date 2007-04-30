# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.



coefvlm <- function(object, matrix.out=FALSE, label=TRUE, compress=TRUE)
{

    ans <- object@coefficients
    if(!label)
        names(ans) <- NULL
    if(!matrix.out && compress)
        return(ans)

 
    ncolx <- object@misc$p   # = length(object@constraints)
    M <- object@misc$M

    xij <- object@control$xij
    Blist <- object@constraints
    if(!length(xij) && all(trivial.constraints(Blist))) {
        B <- matrix(ans, nrow=ncolx, ncol=M, byrow=TRUE)
    } else {
        B <- matrix(as.numeric(NA), nrow=ncolx, ncol=M)

        if(length(xij)) {

            Xmat = object@x # model.matrix(object)
            atx = attr(Xmat, "assign")
           tmp9=attributes(lm2vlm.model.matrix(Xmat,object@constraints,xij=xij))
            nasgn = tmp9$assign
            vasgn = tmp9$vassign
            if(!length(atx) || !length(nasgn) || !length(vasgn))
                stop("can't get atx, nasgn and/or vasgn")

            if(inherits(xij, "formula"))
                xij = list(xij)

            for(i in 1:length(xij)) {
                tform = terms(xij[[i]])
                atform = attr(tform, "term.labels")
                for(k in 1:length(atform)) {
                    for(s in atx[[(atform[[k]])]]) {
                        if(length(Blist[[s]])) {
                            Blist[[s]] = ei(k, M)   # changed 
                        }
                    }
                }
            }


            ncolBlist <- unlist(lapply(Blist, ncol))   # Modified 
            ans.save = ans    # small
            ans = rep(as.numeric(NA), len=sum(ncolBlist))   # big
            ans.copied = rep(FALSE, len=length(ans))
            ans.save.copied = rep(FALSE, len=length(ans.save))

            ptr0 = rep(as.numeric(NA), len=length(xij))
            for(i in 1:length(xij)) {
                tform = terms(xij[[i]])
                atform = attr(tform, "term.labels")
                response.name = (dimnames(attr(tform, "factors"))[[1]])[1]

                for(s in 1:M)
                    if(length(nasgn[[atform[s]]])) {
                        ptr0[i] = s
                        break
                    }
                dest.col.index = nasgn[[atform[ptr0[i]]]]
                rindex = vlabel(response.name, 
                                length(dest.col.index), M=M, separator="")
                dummy = ans.save*0 + 1:length(ans.save)  # names retained
                dest.col.index = dummy[rindex]

                for(k in ((1:M))) {
                    from.col.index = nasgn[[atform[k]]]   # May be NULL
                    if(length(from.col.index)) {
                        ans[from.col.index] = ans.save[dest.col.index]
                        ans.copied[from.col.index] = TRUE
                        ans.save.copied[dest.col.index] = TRUE
                    }
                }
            }

            if(any(!ans.copied)) {
                ans[!ans.copied] = ans.save[!ans.save.copied]
            }
            names(ans) = vlabel(names(ncolBlist), ncolBlist, 
                                M=M, separator="")

        }

        if(!matrix.out && !compress)
            return(ans) 

        ncolBlist <- unlist(lapply(Blist, ncol)) 
        nasgn <- names(Blist)
        temp <- c(0, cumsum(ncolBlist))
        for(i in 1:length(nasgn)) {
            index <- (temp[i]+1):temp[i+1]
            cm <- Blist[[nasgn[i]]]
            B[i,] <- cm %*% ans[index]
        }
    }

    if(label) {
        d1 <- object@misc$colnames.x
        d2 = object@misc$predictors.names # Could be NULL
        dimnames(B) <- list(d1, d2)
    }

    if(compress && length(xij)) {
        ci2 = NULL
        for(i in 1:length(xij)) {
            tform = terms(xij[[i]])
            atform = attr(tform, "term.labels")
            response.name = (dimnames(attr(tform, "factors"))[[1]])[1]

            dest.col.index = atx[[atform[ptr0[i]]]]
 

            for(k in ((1:M)[-ptr0[i]])) {
                from.col.index = atx[[atform[k]]]   # May be NULL
                if(length(from.col.index)) {
                    B[dest.col.index,] = B[dest.col.index,] + B[from.col.index,]
                    tmp5 = dimnames(B)[[1]]
                    tmp5[dest.col.index] = vlabel(response.name, 
                               length(dest.col.index), M=M, separator="")
                    dimnames(B) = list(tmp5, dimnames(B)[[2]])
                    ci2 = c(ci2, from.col.index)
                }
            }
        }
        B = B[-ci2,,drop=FALSE]   # Delete rows not wanted 
    }

    B
} # end of coefvlm




    setMethod("coefficients", "vlm", function(object, ...)
               coefvlm(object, ...))
    setMethod("coef", "vlm", function(object, ...)
               coefvlm(object, ...))
    setMethod("coefficients", "vglm", function(object, ...)
               coefvlm(object, ...))
    setMethod("coef", "vglm", function(object, ...)
               coefvlm(object, ...))





Coef.vlm <- function(object, ...) {
    LL <- length(object@family@vfamily)
    funname = paste("Coef.", object@family@vfamily[LL], sep="")
    if(exists(funname)) {
        newcall = paste("Coef.", object@family@vfamily[LL],
                        "(object, ...)", sep="")
        newcall = parse(text=newcall)[[1]]
        eval(newcall)
    } else
    if(length(tmp2 <- object@misc$link) &&
       object@misc$intercept.only &&
       trivial.constraints(object@constraints)) {

        answer = eta2theta(rbind(coef(object)),
                           link=object@misc$link,
                           earg=object@misc$earg)
        answer = c(answer)
        if(length(ntmp2 <- names(tmp2)) == object@misc$M)
            names(answer) = ntmp2
        answer
    } else {
        coef(object, ... )
    }
}

setMethod("Coefficients", "vlm", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coef", "vlm", function(object, ...)
               Coef.vlm(object, ...))


if(!is.R()) {
setMethod("Coefficients", "vglm", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coef", "vglm", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coefficients", "vgam", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coef", "vgam", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coefficients", "rrvglm", function(object, ...)
               Coef.vlm(object, ...))
setMethod("Coef", "rrvglm", function(object, ...)
               Coef.vlm(object, ...))
}



if(FALSE) 
coef.rrvglm = function(object, 
                       type = c("all", "vlm"), 
                       matrix.out=FALSE, label=TRUE) {

}


