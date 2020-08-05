#' Gaussian Kernel
gau <- function(x, y=NULL, sigma=1)
{
    x2 <- colSums(x^2)
    y2 <- if(is.null(y)) x2 else colSums(y^2)
    xy <- crossprod(x, y)

    ds <- (outer(x2, y2, `+`) - 2 * xy) / nrow(x)
    ## exp(-ds / (1 * sigma))
    max(ds) - ds
}

imp.gau <- function(g, brk=brk.max, ...)
{
    P <- ncol(g)
    N <- nrow(g)

    ## NA mask, and pairwise completion count
    Z <- is.na(g)

    ## center and scale X_j, for j = 1 ... P (the SNPs)
    x <- drop(scale(g, TRUE, FALSE))
    x[Z] <- 0

    ## d_ij = x_i' x_j, with pairwise complete samples between variant i and j;
    ## d_ij resembles a distance measure.
    D <- gau(x)
    
    ## e_ij = x_i' x_i, with pairwise complete samples between variant i and j;
    E <- crossprod(x^2, !Z)
    
    ## b_ij = d_ij / e_ij, regression coeficient
    B <- D # / E
    
    ## sum of prediction
    ## h_i = sum_{j=1}^M x_j b_ji / (non-na x_j count)
    x <- x %*% B / (P - rowSums(Z))

    brk(g, x)
}

#' Impute by Regresion
#'
#' @param g genotype of N row samples and P column variants
#'
#' @return a list of imputed values and dosage genotype
imp.reg <- function(g, brk=brk.max, wgt=0, ...)
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
    q <- x^2

    ## p_rs =  sum_i(x_ir *  x_is * i_rs)  = x_r' x_s:  sum of  product between
    ## complete samples of variant r and s, it resembles a distance measure.
    xy_rs <- crossprod(x)

    ## v_rs = sum_i(x_ir * x_ir * i_rs) = (x_r^2)' i_s, with pairwise complete
    ## samples between variant i and j;
    xx_rs <- crossprod(x * x, C)
    stopifnot(all.equal(xx_rs, crossprod(x * x, C)))

    ## regression coeficient  rc = sum_i(x_ir *  x_is * i_rs) /  sum_i(x_ir *
    ## x_ir * i_rs) = xy_rs / xx_rs
    rc_rs <- xy_rs / xx_rs

    ## mu_rs = sum_i(x_ir * i_rs) / sum_i(i_rs)
    sm_rs <- crossprod(x, !Z)
    sz_rs <- crossprod(!Z)
    mu_rs <- sm_rs / sz_rs

    ## square residual:
    ## e2_rs = sum_i(x_is * x_is * i_rs) - 2 * rc_rs sum_i(x_ir * x_is * i_rs)
    ## + rc_rs^2 * sum_i(x_ir * x_ir * i_rs)
    e2_rs <- t(xx_rs) - xy_rs * xy_rs / xx_rs
    d2_rs <- xx_rs - 2 * sm_rs * mu_rs + mu_rs^2
    s2_rs <- e2_rs / d2_rs / (sz_rs - 2)
    diag(s2_rs) <- 1 / (colSums(!Z) - 2)
    zs_rs <- rc_rs / sqrt(s2_rs)

    ## sum of prediction
    ## h_is = sum_r(x_ir * rc_rs * i_rs) / sum_r(i_rs)
    if(wgt == 0)
        x <- x %*% rc_rs / (P - rowSums(Z) + 1)
    else
        x <- x %*% zs_rs / (P - rowSums(Z) + 1)

    brk(g, x)
}

