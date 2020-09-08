imp.box <- function(g, wgt=1, hrd=1, itr=1, thd=1e-12, ...)
{
    N <- nrow(g); M <- ncol(g); Z <- is.na(g); C <- !Z

    .zs <- function(g)
    {
        ## column center; setting NA to 0 after opt out incomplete pairs.
        x <- as.matrix(scale(g, TRUE, FALSE))
        x[Z] <- 0

        ## sum_i(x_ir * x_is * i_rs), sum of product between of x_r and x_s,
        ## complete pairs only
        xy <- crossprod(x)

        ## sum_i(x_ir * x_ir * i_rs), sum of square of x_r in the complete
        ## pairs of x_r and x_s
        xx <- crossprod(x^2, C)
        
        ## sum_i(x_ir * x_is * i_rs) / sum_i(x_ir  * x_ir * i_rs) = xy / xx_rs,
        ## the regression coeficient
        rc <- xy / xx

        ## sum_i(x_ir * i_rs) / sum_i(i_rs), mean  of x_r in the complete pairs
        ## of x_r and x_s
        rs <- crossprod(x, C) # the sum
        nn <- crossprod(C)    # the non-NA count
        mu <- rs / nn         # the mean

        ## sum_i(x_is * x_is * i_rs) - 2 * rc sum_i(x_ir * x_is * i_rs) + rc^2 *
        ## sum_i(x_ir * x_ir * i_rs), squared residual
        e2 <- t(xx) - xy * xy / xx

        ## denominator
        d2 <- xx - 2 * rs * mu + mu^2
        
        ## squared standard error
        s2 <- e2 / d2 / (nn - 2)
        if(wgt == 1)
            diag(s2) <- 1 / (colSums(C) - 2)
        if(wgt == 2)
            diag(s2) <- (N - diag(nn)) / ((diag(nn) - 2) * (diag(nn) - 1))
        ## z-scores
        s2[s2 <= 0 ] <- min(s2[s2 > 0]) # avoid s2=0 caused by perfect fit
        rc / sqrt(s2)
    }
    
    p <- g
    x <- array(c(g == 0 & C, g == 1 & C, g == 2 & C), c(N, M, 3)) * 1
    
    ## g <- imp.mod(g)
    c_0 <- Inf # consistancy
    w <- 1
    while(itr > 0)
    {
        ## complete pairs, for each genotype {0, 1, 2}
        n <- array(apply(x, 3, crossprod, C), c(M, M, 3))
        if(wgt > 0)
        {
            w <- .zs(p)^2
            w <- w / mean(diag(w))
        }
        ## transition
        y <- array(0, c(N, M, 3))
        for(i in 1:3)
        {
            n_i <- crossprod(x[, , i], C) # pairwise complete count for x == 0, 1, or 2
            r_i <- 1/ n_i
            r_i[is.infinite(r_i)] <- 0
            for(j in 1:3)
            {
                y[, , j] <- y[, , j] + x[, , i] %*% (w * crossprod(x[, , i], x[, , j]) * r_i)
            }
        }

        ## balance the contribution of predictors.
        y <- y / rowSums(C)
        y <- y / array(rowSums(y, dims=2), c(N, M, 3))
        
        ## imputation tentative 
        q <- 0 * y[, , 1] + 1* y[, , 2] + 2 * y[, , 3]

        ## consistancy check
        c_1 <- mean((round(q) - g)^2, na.rm=TRUE)
        if(c_1 - c_0 > thd)
            break
        itr <- itr - 1

        ## prepare for the next iteration
        c_0 <- c_1
        x <- y
        p <- q
        C[] <- TRUE
    }

    if(hrd == 1)
        p <- round(p)
    if(hrd == 2)
        p <- brk.kmn(g, p)
    g[Z] <- p[Z]
    g
}

