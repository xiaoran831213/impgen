## imputation tools

#' Break imputed values into descrete dosage
#'
#' @param g matrix of original genotyped in allele dosage format;
#' @param x matrix of imputed values;
#' 
#' @return a matrix of {0, 1, 2} allele dosage.
imp.brk <- function(g, x)
{
    stopifnot(all(dim(g) == dim(x)))
    h <- matrix(0L, nrow(x), ncol(x))
    for(j in seq(ncol(g)))
    {
        h[x[, j] > max(x[g[, j] == 0, j], na.rm=TRUE), j] <- 1L
        h[x[, j] > max(x[g[, j] == 1, j], na.rm=TRUE), j] <- 2L
    }
    g[is.na(g)] <- h[is.na(g)]
    g
}

#' Assign allele dosage by closest center values
#'
#' @param g matrix of original genotyped in allele dosage format;
#' @param x matrix of imputed values;
#' 
#' @return a matrix of {0, 1, 2} allele dosage.
imp.ccv <- function(g, x)
{
    ## center values
    for(j in seq(ncol(g)))
    {
        c0 <- mean(x[g[, j] == 0, j], na.rm=TRUE)
        c1 <- mean(x[g[, j] == 1, j], na.rm=TRUE)
        c2 <- mean(x[g[, j] == 2, j], na.rm=TRUE)
        d0 <- abs(x[, j] - c0)
        d1 <- abs(x[, j] - c1)
        d2 <- abs(x[, j] - c2)
        x[d0 <= d1 && d0 <= d2, j] <- 0
        x[d1 <= d0 && d1 <= d2, j] <- 1
        x[d2 <= d0 && d2 <= d1, j] <- 2
    }
    x
}

#' Imputation by average
imp.avg <- function(g)
{
    apply(g, 2L, function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); x})
}

#' Imputation by random sample
imp.rnd <- function(g)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- sample(x[-i], length(i))
        x
    })
}

#' Imputation by Mode
imp.mod <- function(g)
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
