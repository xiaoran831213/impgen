#' Impute by Regresion
#'
#' @param g genotype of N row samples and P column variants
#' 
#' @return a list of imputed values and dosage genotype
imp.reg <- function(g)
{
    P <- ncol(g)
    N <- nrow(g)

    ## NA mask, and pairwise completion count
    Z <- is.na(g)
    C <- N - crossprod(Z)
    a <- colMeans(g, na.rm=TRUE)
    
    ## center and scale X_j, for j = 1 ... P (the SNPs)
    x <- sweep(g, 2, a, `-`)
    x[Z] <- 0
    
    ## d_ij =  x_i' x_j *  (n_ij / c_ij),  where n_ij =  N, and c_ij  counts the
    ## complete pairs;
    ## d_ij resembles a distance measure.
    D <- crossprod(x) * (N / C)

    ## simple regression, use x_i to preduct x_j
    ## hat{x_j}(x_i) = b_ij x_i, needs regresion coeficient b_ij
    ## b_ij = (x_i' x_j) / (x_i' x_i) = d_ij / d_ii, i,j = 1 ... P
    B <- D / diag(D)
    
    ## for a sample h, the ith SNP
    ## h_i = sum_{j=1, j!=i} x_j b_ji / (non-na x_j count)
    diag(B) <- 0
    x <- x %*% B / (P - 1)

    sweep(x, 2, a, `+`)
}
