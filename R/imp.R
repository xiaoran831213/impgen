## imputation tools

#' Break imputed values into descrete dosage
#'
#' @param g matrix of original genotyped in allele dosage format;
#' @param x matrix of imputed values;
#' 
#' @return a matrix of {0, 1, 2} allele dosage.
brk.max <- function(g, x)
{
    stopifnot(all(dim(g) == dim(x)))
    h <- matrix(0L, nrow(x), ncol(x))
    for(j in seq(ncol(g)))
    {
        h[x[, j] > max(x[g[, j] == 0, j], -Inf, na.rm=TRUE), j] <- 1L
        h[x[, j] > max(x[g[, j] == 1, j], -Inf, na.rm=TRUE), j] <- 2L
    }
    g[is.na(g)] <- h[is.na(g)]
    g
}

brk.kmn <- function(g, x)
{
    stopifnot(all(dim(g) == dim(x)))
    h <- matrix(0L, nrow(x), ncol(x))
    for(j in seq(ncol(g)))
    {
        m <- unlist(tapply(x[, j], g[, j], mean, simplify = FALSE))
        r <- try(kmeans(x[, j], m)$cluster - 1)
        if(inherits(r, "try-error"))
        {
            print("bad")
        }
        h[, j] <- r
    }
    g[is.na(g)] <- h[is.na(g)]
    g
}

brk.pol <- function(g, x)
{
    stopifnot(all(dim(g) == dim(x)))
    h <- matrix(0L, nrow(x), ncol(x))
    for(j in seq(ncol(g)))
    {
        f <- factor(as.integer(g[, j]), c(0, 1, 2))
        m <- MASS::polr(f ~ x[, j])
        v <- predict(m, x[, j])
        h[, j] < as.integer(v) - 1
    }
    g[is.na(g)] <- h[is.na(g)]
    g
}


#' Imputation by random sample
imp.rnd <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- sample(x[-i], length(i))
        x
    })
}

#' Imputation by Mode
imp.mod <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        u <- unique(x[-i])
        a <- tabulate(match(x[-i], u))
        x[i] <- u[a == max(a)]
        x
    })
}
