imp.mic <- function(obs, ...)
{
    obs <-data.frame(obs)
    for(i in seq_along(obs))
        obs[[i]] <- factor(obs[[i]])
    ret <- mice(obs, method="polr", ...)
    hat <- as.integer(as.matrix(complete(ret)))
    dim(hat) <- dim(obs)
    cat(sprintf("MICE: %d NA --> %d\n", sum(is.na(obs)), sum(is.na(hat))))
    hat
}

acs <- function(gmx, obs, mtd, id=NULL, ...)
{
    tic <- proc.time()
    hat <- mtd(obs, ...)
    toc <- proc.time()
    tcs <- unname((toc - tic)['elapsed'])

    if(is.null(id))
    {
        mtd <- as.character(substitute(mtd))
        id <- sub("^imp[.]", "", mtd)
    }
    crh <- cor(hat, use="pair")
    crg <- cor(gmx)

    mce <- mce(crh, crg)
    mae <- mae(hat, gmx)
    nna <- sum(is.na(hat))
    egn <- eigen(crh, TRUE, TRUE)$values[ncol(gmx)]

    rpt <- c(nna=nna, mce=mce, mae=mae, egn=egn, tcs=tcs)
    ret <- .d(mtd=id, key=names(rpt), val=unname(rpt))
    ret
}


sim <- function(N=1e3, M=1e1, NAF=.1, out=NULL, ...)
{
    arg <- get.arg()
    set.seed(arg$seed)
    arg$out <- NULL

    if(!is.null(out) && file.exists(out))
    {
        cat("XT: ouput exists, ", out, "\n", sep="")
        return(readRDS(out))
    }

    ## gmx: true genotype; obs: observed
    flood(get.gmx(N, M, NAF=NAF, ...))

    ## complete correlation matrix
    cmx <- cor(gmx)

    ## pairwise complete correlation
    pcc <- cor(obs, use='pair') 

    ## emc: EM estimated correlation
    emc <- cov2cor(emc(obs, 1e2, 1e-5, msg=0)$v)

    ## compete
    ret <- list()
    ## ret[[length(ret) + 1]] <- acs(gmx, obs, imp.rnd, id="drw", ...)
    ret[[length(ret) + 1]] <- acs(gmx, obs, imp.mod, id="mod", ...)
    ## ret[[length(ret) + 1]] <- acs(gmx, obs, imp.reg, id="zs1", wgt=1, ...)
    ret[[length(ret) + 1]] <- acs(gmx, obs, imp.reg, id="zs2", wgt=2, ...)
    ret[[length(ret) + 1]] <- acs(gmx, obs, imp.box, id="b02", hrd=0, wgt=2, ...)
    ret[[length(ret) + 1]] <- acs(gmx, obs, imp.box, id="b12", hrd=1, wgt=2, ...)
    ret[[length(ret) + 1]] <- .d(mtd='pcc', key='mce', val=mce(pcc, cmx))
    ret[[length(ret) + 1]] <- .d(mtd='emc', key='mce', val=mce(emc, cmx))
    ## ret[[length(ret) + 1]] <- acs(gmx, obs, imp.box, id="b22", hrd=2, wgt=2, ...)
    ret[[length(ret) + 1]] <- acs(gmx, obs, imp.mic, id="mic", ...)
    set.seed(NULL)

    ret <- within(cbind(arg, do.call(rbind, ret)), val <- round(val, 4))
    if(!is.null(out))
        saveRDS(ret, out)
    ret
}
## r <- sim(1500, 10, NAF=.05, seed=1146)
