#' Impute by Regresion
#'
#' @param g genotype of N row samples and P column variants
#'
#' @return a list of imputed values and dosage genotype
imp.reg <- function(g, brk=brk.max, wgt=1, ...)
{
    P <- ncol(g)
    N <- nrow(g)

    ## NA mask, and pairwise completion count
    Z <- is.na(g)
    C <- 1 * !Z

    ## center and scale X_j, for j = 1 ... P (the SNPs)
    x <- as.matrix(scale(g, TRUE, FALSE))

    ## setting NA to 0 after centering, amounts to opt out incomplete pairs.
    x[Z] <- 0

    ## sum_i(x_ir * x_is * i_rs), sum of product between of x_r and x_s,
    ## complete pairs only
    xy <- crossprod(x)

    ## sum_i(x_ir * x_ir * i_rs), sum of square of x_r in the complete pairs of
    ## x_r and x_s
    xx <- crossprod(x^2, C)
    
    ## sum_i(x_ir * x_is * i_rs) / sum_i(x_ir * x_ir * i_rs) = xy / xx_rs, the
    ## regression coeficient
    rc <- xy / xx

    ## sum_i(x_is * x_is * i_rs) - 2 * rc sum_i(x_ir * x_is * i_rs) + rc^2 *
    ## sum_i(x_ir * x_ir * i_rs), squared residual
    e2 <- t(xx) - xy * xy / xx

    ## denominator
    ## sum_i(x_ir * i_rs) / sum_i(i_rs), mean of x_r in the complete pairs of
    ## x_r and x_s
    su <- crossprod(x, C) # the sum
    nn <- crossprod(C)    # the non-NA count
    mu <- su / nn         # the mean
    d2 <- xx - 2 * su * mu + mu^2
    
    ## squared standard error
    s2 <- e2 / d2 / (nn - 2)
    ## diag(s2) <- 1 / (colSums(C) - 2)
    ## diag(s2) <- (N - diag(nn)) / ((diag(nn) - 2) * (diag(nn) - 1))

    ## z-scores
    ## zs <- rc / sqrt(s2)
    
    ## sum of prediction
    ## h_is = sum_r(x_ir * rc * i_rs) / sum_r(i_rs)
    cc <- rowSums(C) - C + 1
    if(wgt == 0)
        x <- x %*% rc / cc
    else if(wgt == 1)
    {
        diag(s2) <- 1 / (colSums(C) - 2)
        s2[s2 <= 0 ] <- min(s2[s2 > 0]) # avoid s2=0 caused by perfect fit
        zs <- rc / sqrt(s2)
        x <- x %*% zs / cc
    }
    else
    {
        diag(s2) <- (N - diag(nn)) / ((diag(nn) - 2) * (diag(nn) - 1))
        s2[s2 <= 0 ] <- min(s2[s2 > 0]) # avoid s2=0 caused by perfect fit
        zs <- rc / sqrt(s2)
        x <- x %*% zs / cc
    }
    brk(g, x)
}
