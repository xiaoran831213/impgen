#' Impute by conditional normal
#'
#' @param g genotype of N row samples and P column variants
#' 
#' @return a list of imputed values and dosage genotype
imp.cdn <- function(g)
{
    P <- ncol(g)
    N <- nrow(g)

    ## NA mask, and pairwise completion count
    Z <- is.na(g)
    C <- N - crossprod(Z)
    a <- colMeans(g, na.rm=TRUE)
    
    ## center and scale X_j, for j = 1 ... P (the SNPs)
    x <- sweep(g, 2, a, `-`)
    ## x <- scale(g)
    ## x <- g
    x[Z] <- 0
    
    ## d_ij =  cov(x_i, x_j), covariance by pairwise completion
    D <- crossprod(x) / (C - 1)
    A <- chol2inv(chol(D)) # inverse

    for(i in seq(nrow(g)))
    {
        u <- which(is.na(g[i, ]))
        if(length(u) < 1)
            next
        x[i, +u] <- D[+u, -u, drop=FALSE] %*% A[-u, -u] %*% x[i, -u]
    }
    imp.brk(g, x)
}


imp.loo <- function(g)
{
    P <- ncol(g)
    N <- nrow(g)

    ## NA mask, and pairwise completion count
    Z <- is.na(g)
    C <- N - crossprod(Z)
    a <- colMeans(g, na.rm=TRUE)
    
    ## center and scale X_j, for j = 1 ... P (the SNPs)
    x <- sweep(g, 2, a, `-`)
    ## x <- scale(g)
    x <- g
    x[Z] <- 0
    
    ## d_ij =  cov(x_i, x_j), covariance by pairwise completion
    D <- crossprod(x) / (C - 1)
    A <- chol2inv(chol(D)) # inverse

    x <- t(t(x) - A %*% t(x) / diag(A))
    imp.brk(g, x)
}
