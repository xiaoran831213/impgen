#' Cached genotype data
#'
#' A piece of genotype from 1000 genome project.
#'
#' Taken from chromosome 17, region q12. The MAF is no less than 0.05.
if(!exists("c17"))
    c17 <- readRDS("17q12.rds")
if(!exists("f17"))
    f17 <- as.vector(table(as.vector(c17), useNA="no") / sum(!is.na(c17)))

## non-colinear variables
ncv <- function(ldm, lo=0, up=1, ...)
{
    r <- abs(ldm)
    r[upper.tri(r, TRUE)] <- lo + (up - lo) / 2
    b <- which(lo <= r & r <= up, arr.ind=TRUE)
    apply(r, 2, function(.) all(lo <= . & . <= up))
}

maf <- function(gmx)
{
    af <- colMeans(gmx, na.rm=TRUE) / 2
    pmin(af, 1 - af)
}

#' Draw a genotype matrix
#'
#' @param N: samples to draw (rows)
#' @param M: variants to craw (columns)
#' @param nlr: non-linearity thereshold
#' @param psd: postive semi-definite threshold
#' @param bat: the genotype data, if NULL, use beta distribution
#' @param a: shape parameter alpha
#' @param b: shape parameter beta
#'
#' The between variants correlation must not exceed \code{nlr}.
#'
#' The smallest eigenvalue of the LD matrix must not below (-\code{psd}) * (the
#' largest eigenvalue).
get.gmx <- function(N=5e2, M=5, NAF=.1, MAF=.04, psd=NULL, ...)
{
    if(is.null(psd))
        psd <- sqrt(.Machine$double.eps)

    ## M variants is in demand, reserve 3 times more
    P <- min(M * 4, ncol(c17))
    while(TRUE)
    {
        i <- sample.int(nrow(c17), N)                 # N
        j <- seq(sample(ncol(c17) - P, 1) + 1, l=P)   # P
        gmx <- c17[i, j]                              # N x P
        gmx <- gmx[, maf(gmx) > MAF]
        ldm <- cor(gmx, use="pair")                   # P x P

        ## drop variants to enforce non-linearity; if there are less than M
        ## remaining, try again with a bigger reserve.
        kpp <- ncv(ldm, ...)
        gmx <- gmx[, kpp]
        if(NCOL(gmx) < M)
        {
            P <- min(P + M - NCOL(gmx), ncol(c17))
            ## cat("P = ", P, "\n", sep="")
            next
        }
        
        ## select M variants now
        j <- seq(sample(ncol(gmx) - M, 1) + 1, l=M)
        gmx <- gmx[, j, drop=FALSE]
        ldm <- cor(gmx, use="pair")
        
        ## if min(eigenvalue) < (threshold) * max(eigenvalue), try again
        egv <- eigen(ldm, TRUE, TRUE)$values
        if (egv[M] < psd * egv[1L])
        {
            cat("Non-PSD!\n")
            next
        }
        break
    }

    idx <- which(is.na(gmx))
    gmx[idx] <- sample(0:2, length(idx), TRUE, f17)
    obs <- set.nan(gmx, NAF)

    ## remove rows that are entirely NA.
    i <- rowSums(is.na(obs)) < M
    gmx <- gmx[i, ]
    obs <- obs[i, ]
    
    list(gmx=gmx, obs=obs)
}
