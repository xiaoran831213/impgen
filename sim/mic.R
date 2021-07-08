library(mice)
mic <- function(N=1e3, M=1e1, NAF=.1, out=NULL, ...)
{
    arg <- get.arg()
    out <- arg$out

    if(!is.null(out) && file.exists(out))
    {
        cat("XT: ouput exists, ", out, "\n", sep="")
        ret <- readRDS(out)
    }
    else
    {
        set.seed(arg$seed)
        ## get genotype
        ret <- get.gmx(N, M, NAF=NAF, psd=-1e-6, ...)
        gmx <- ret$gmx
        obs <- ret$obs
        ## impute with mice
        ret <- acs(gmx, obs, imp.mic, id="mic", ...)
        set.seed(NULL)
        ret <- cbind(arg, ret)
        
        if(!is.null(out))
            saveRDS(ret, out)
    }
    ret
}
