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
