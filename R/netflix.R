#' Ridge Regression
#'
#' @param x matrix of predictors
#' @param y vector of outcomes
#' @param la regularization hyper-parameter lambda, NULL for auto-decide
#' @param nl number of lambda(s) to try for auto-decision.
#'
#' @return a glmnet object for fitted model.
ridge <- function(x, y, la=NULL, nl=10, a=0, ...)
{
    ## ret <- glmnet(x, y, alpha=a, lambda=la, nlambda=nl, standardize=FALSE, ...)
    with(cv.glmnet(x, y, alpha=a, lambda=la, nlambda=nl, standardize=FALSE, ...),
    {
        beta <- unname(glmnet.fit$beta[, lambda.min == lambda])
        mu <- unname(glmnet.fit$a0[lambda.min == lambda])
        list(lambda=lambda.min, beta=beta, mu=mu)
    })
}


#' Imputation by Matrix Fctorization
#'
#' A imitation of Netflix's film recommendation algorithm in genotype.
#'
#' @param g matrix of the N by P genotype
#' @param f the dimensionality latent factor space
#' @param la the regularization parameters
#'
#' \code{la} can be  a vector of regularization weights in  decending order, an
#' integer to specify  the number of weights to automatically  try, or NULL for
#' the default (100 weights).
imp.gmx <- function(g, f=2, la=0, nl=1, itr=100, thd=1e-7, ...)
{
    N <- nrow(g); P <- ncol(g); L <- N * P

    ru <- as.vector(t(g)) # user contains item
    ri <- as.vector(g)    # item contains user
    ku <- which(!is.na(ru))
    ki <- which(!is.na(ri))
    
    xu <- kronecker(Diagonal(N), rep(1, P))
    bu <- rowMeans(g, TRUE)
    pu <- matrix(rnorm(N * f, 0.2), N, f)

    xi <- kronecker(Diagonal(P), rep(1, N))
    bi <- colMeans(g, TRUE)
    qi <- Matrix(rnorm(P * f, 0.2), P, f)

    mu <- -mean(g, na.rm=TRUE)

    ## gather paramters
    s1 <- list(mu=mu, bu=bu, bi=bi, pu=pu, qi=qi)
    while(itr > 0)
    {
        s0 <- s1
        ## find bu and pu while fixing bi and qi
        fu <- kronecker(Diagonal(N), qi)
        re <- ridge(cbind(xu, fu)[ku, ], (ru - xi %*% bi)[ku, ], la=la, nl=nl, ...)
        bu <- unname(re$beta[+seq(N)])
        pu <- unname(re$beta[-seq(N)])
        pu <- matrix(pu, N, f, TRUE)

        ## find bi and qi while fixing bu and pu
        fi <- kronecker(Diagonal(P), pu)
        re <- ridge(cbind(xi, fi)[ki, ], (ri - xu %*% bu)[ki, ], la=la, nl=nl, ...)
        bi <- unname(re$beta[+seq(P)])
        qi <- unname(re$beta[-seq(P)])
        qi <- matrix(qi, P, f, TRUE)

        ## the grand mean
        mu <- re$mu

        ## convergence and prediction check
        s1 <- list(mu=mu, bu=bu, pu=pu, bi=bi, qi=qi)
        sq <- sum(sapply(names(s1), function(.) mean((s0[[.]] - s1[[.]])^2)))
        hi <- as.vector(mu + xi %*% bi + fi %*% as.vector(t(qi)))
        cr <- cor(hi, ri, use="p")
        er <- mean((hi - ri)^2, na.rm=TRUE)
        cat(sprintf("sq=%02e, er=%1e, cr=%.2f\n", sq, er, cr))
        if(sq < thd)
            break
        itr <- itr - 1
    }
    dim(hi) <- dim(g)
    hi
}

im1.gmx <- function(g, f=2, la=0, nl=1, itr=100, thd=1e-7, ...)
{
    N <- nrow(g); P <- ncol(g); L <- N * P

    ru <- as.vector(t(g)) # user contains item
    ri <- as.vector(g)    # item contains user
    ku <- which(!is.na(ru))
    ki <- which(!is.na(ri))
    
    pu <- matrix(rnorm(N * f, 0.2), N, f)

    xi <- kronecker(Diagonal(P), rep(1, N))
    bi <- colMeans(g, TRUE)
    qi <- Matrix(rnorm(P * f, 0.2), P, f)
    
    mu <- -mean(g, na.rm=TRUE)

    ## gather paramters
    s1 <- list(mu=mu, bi=bi, pu=pu, qi=qi)
    while(itr > 0)
    {
        s0 <- s1
        ## find pu while fixing bi and qi
        fu <- kronecker(Diagonal(N), qi)
        re <- ridge(fu[ku, ], (ru - xi %*% bi)[ku], la=la, nl=nl, thresh=thd, ...)
        pu <- matrix(re$beta, N, f, TRUE)
        l1 <- re$lambda

        ## find bi and qi while fixing bu and pu
        fi <- kronecker(Diagonal(P), pu)
        re <- ridge(cbind(xi, fi)[ki, ], ri[ki], la=la, nl=nl, thresh=thd, ...)
        bi <- re$beta[+seq(P)]
        qi <- matrix(re$beta[-seq(P)], P, f, TRUE)
        l2 <- re$lambda

        ## the grand mean
        mu <- re$mu

        ## convergence check
        s1 <- list(mu=mu, pu=pu, bi=bi, qi=qi)
        sq <- 0
        for(. in names(s1))
            sq <- sq + mean((s0[[.]] - s1[[.]])^2)
        hi <- as.vector(mu + xi %*% bi + fi %*% as.vector(t(qi)))
        cr <- cor(hi, ri, use="p")
        er <- mean((hi - ri)^2, na.rm=TRUE)
        ms <- sprintf("sq=%02e, er=%1e, cr=%.2f, la=%1e, %1e\n", sq, er, cr, l1, l2)
        cat(ms)
        if(sq < thd)
            break
        itr <- itr - 1
    }
    dim(hi) <- dim(g)
    hi
}


test0 <- function(N=c(5e2, 5e2), P=c(5, 5), noise=.2, r=0.02)
{
    R <- rep(seq(1, l=length(N), b=1),         N)
    C <- rep(seq(0, l=length(P), b=length(N)), P)
    x <- matrix(seq(-2, 2, l=length(N) * length(P))[outer(R, C, `+`)], sum(N), sum(P))
    
    ## x <- x + outer(rnorm(sum(N), 0, .5), rnorm(sum(P), 0, .5), `+`)

    x <- x +  rnorm(sum(N) * sum(P), 0, noise)

    d <- x
    d[sample(length(d), length(d) * r)] <- NA
    list(x=x, d=d)
}
## t1 <- within(test(c(500, 500), c(10, 10), noise=.1, r=.1), imp.gmx(d, 10))
