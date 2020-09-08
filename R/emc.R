#' mean and covariance of incomplete multivariate normal data
#'
#' Assuming the  row samples in an  imcomplete \code{N * M}  matrix \code{x}
#' follows multivate  normal, using expectation conditional  maximization (ECM)
#' to estimate its \code{N * 1} mean vector and \code{M * M} covariance matrix.
#'
#' @param x \code{N * M} matrix with \code{N} samples of a \code{M}-dimensional
#' random vector;
#' @param itr Maximum number of iterations for ECM algorithm (def=50);
#' @param tol Convergence tolerance for ECM (def=1e-8)
#' @param msg print messages on screen?
#'
#' @return the miximum likelihood 
emc <- function(x, itr=NULL, tol=NULL, msg=1)
{
    if(is.null(itr))
        itr <- 50
    if(is.null(tol))
        tol <- sqrt(.Machine$double.eps)

    ## step 2 - initialization
    N <- nrow(x)
    M <- ncol(x)
    Nan <- is.na(x)

    DOF <- N - M - sum(rowSums(Nan) == M)
    stopifnot(DOF > 0)

    NanCols <- colSums(Nan)
    stopifnot(all(NanCols <= N - 2))

    os <- list()
    x0 <- x
    m0 <- colMeans(x0, na.rm=TRUE)
    v0 <- cov(x0, use="pair")
    ## o0 <- obj(x0, m0, v0)
    ## step 3 - main loop
    while(itr > 0)
    {
        ## Step 4 - mean expectation and conditional maximization
        m1 <- rep(0, M)
        c1 <- 0
        for(n in seq(N))
        {
            d <- x0[n, ]
            a <- is.na(d)
            b <- !a
            if(any(a))
            {
                .s <- try(solve(v0[b, b], (d - m0)[b]), TRUE)
                if(inherits(.s, "try-error"))
                    next
                d[a] <- m0[a] + v0[a, b, drop=FALSE] %*% .s
            }
            m1 <- m1 + d
            c1 <- c1 + 1
        }
        m1 <- m1 / c1

        ## step 5 covariance expectation and conditional maximization
        v1 <- matrix(0, M, M)
        c1 <- 0
        for(n in seq(N))
        {
            d <- x0[n, ]
            a <- is.na(d)
            b <- !a
            e <- matrix(0, M, M)
            if(any(a))
            {
                .s <- try(solve(v0[b, b, drop=FALSE], (d - m1)[b]), TRUE)
                .t <- try(solve(v0[b, b, drop=FALSE], v0[b, a]), TRUE)
                if(inherits(.s, "try-error") || inherits(.t, "try-error"))
                    next
                if(inherits(try(v0[a, b] %*% .s), "try-error"))
                {
                    print("aaa")
                }
                d[a] <- m0[a] + v0[a, b] %*% .s
                e[a, a] <- v0[a, a] - v0[a, b] %*% .t
            }
            v1 <- v1 + tcrossprod(d - m1) + e
            c1 <- c1 + 1
        }
        v1 <- v1 / c1
        
        ## Step 6 - evaluate objective and test for convergence
        dif <- mce(v0, v1)
        if(msg) # print out?
            cat(sprintf("%04d: %-7.1e", itr, dif), "\n", sep="")
        if(abs(dif) < tol) # stop?
            break
        itr <- itr - 1 # the next iteration
        m0 <- m1
        v0 <- v1
    }
    list(m=m1, v=v1)
}

obj <- function(x, m=NULL, v=NULL)
{
    n <- nrow(x)
    p <- ncol(x)

    if(is.null(m)) # mean
        m <- colMeans(x, na.rm=TRUE)

    if(is.null(v)) # covariance
        v <- cov(x, use="pair")

    x <- sweep(x, 2, m) # centered
    x[is.na(x)] <- 0    # naive imputation

    ## u <- chol(v)
    ## a <- chol2inv(u)
    ## l <- -.5 * sum(x %*% a * x) - n * sum(log(diag(u))) - .5 * n * p * log(2*pi)
    a <- with(svd(v),
    {
        i <- d > d[1L] * sqrt(.Machine$double.eps)
        if(all(i))
            v %*% (t(u) / d)
        else
            v[, i, drop=FALSE] %*% (t(u[, i, drop=FALSE]) / d[i])
    })
    d <- as.numeric(determinant(v)$modulus)

    l <- -.5 * sum(x %*% a * x) - n * d / 2 # - .5 * n * p * log(2*pi)
    l / n
}


try1 <- function(N=1e3, M=10, f=.1)
{
    ## gmx: true genotype
    ## obs: incomplete observation
    flood(get.gmx(N, M, f, psd=1e-6))
    ## cmx: true covariance
    cmx <- cor(gmx)

    ## emc covariance
    emc <- cov2cor(emc(obs, 50, 1e-4)$v)
    
    ## pairwise complete
    pwc <- cor(obs, use="pair")

    ## imputation
    imp <- cor(imp.box(obs, wgt=2, hrd=0))
    
    err <- c(
        emc=mce(emc, cmx),
        pwc=mce(pwc, cmx),
        imp=mce(imp, cmx))
    err
}
