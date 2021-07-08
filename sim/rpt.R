library(ggplot2)

.th <- theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(), 
    strip.text.x = element_text(size=14, face="bold"),
    strip.text.y = element_text(size=14, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    legend.title=element_blank(),
    legend.text=element_text(size=14, face="bold"),
    legend.position = "bottom",
    legend.box = "horizontal")

## cap the values
.cp <- function(dat, grp, val='val', cap=0.01, mtd=c('both', 'upper', 'lower'))
{
    grp <- split(dat, dat[, grp])
    mtd <- match.arg(mtd, c('both', 'upper', 'lower'))
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        if(mtd == 'upper')
            v <- pmin(v, quantile(v, 1-cap, na.rm=TRUE))
        else if(mtd == 'lower')
            v <- pmax(v, quantile(v, 0+cap, na.rm=TRUE))
        else
        {
            v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
            v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        }
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

.lb <- function (labels, multi_line = TRUE) 
{
    labels <- label_value(labels, multi_line = multi_line)
    dc <- c(
        `F05`="'% 5' ~~ missing",
        `F15`="'%15' ~~ missing",
        `F25`="'%25' ~~ missing",
        `mae`="mean ~~ absolute ~~ error",
        `mce`="mean ~~ correlation ~~ error",
        `nna`="missings ~~ remained",
        `tcs`="running ~~ time ~~ (sec)")
    if (multi_line)
    {
        lapply(unname(labels), lapply, function(values)
        {
            if(values %in% names(dc))
                values <- dc[[values]]
            c(parse(text = as.character(values)))
        })
    }
    else
    {
        lapply(labels, function(values)
        {
            values <- paste0("list(", values, ")")
            lapply(values, function(expr) c(parse(text = expr)))
        })
    }
}

get.sim <- function(sim, cache=TRUE)
{
    rds <- paste0(sim, '.rds')
    if(file.exists(rds) && cache)
        sim <- readRDS(rds)
    else
    {
        sim <- lapply(dir(sim, '^[0-9A-F]+[.]...$', full=TRUE), function(f)
        {
            print(f)
            r <- readRDS(f)
            r$out <- NULL
            r
        })
        sim <- do.call(rbind, sim)
        saveRDS(sim, rds)
    }
    sim$out <- NULL
    invisible(sim)
}

get.agg <- function(sim)
{
    sim <- get.sim(sim)
    sim <- subset(sim, se=-seed)
    grp <- sim[, c("M", "tag", "mtd", "key")]
    agg <- by(sim, grp, function(r)
    {
        ret <- r[1, ]
        ret$val <- mean(as.numeric(r$val))
        ret
    })
    agg <- do.call(rbind, agg)

    agg <- within(agg, val[key == 'nna'] <- (val / (N * M * NAF))[key == 'nna'])
    agg
}

plt.sim <- function(sim, out = paste0(sim, '.pdf'))
{
    rpt <- get.agg(sim)
    key.sel <- c("mce", "mae", "nna")
    mtd.sel <- c("mic", "b01", "b02", "b11", "b12", "pcc", "emc")
    rpt <- subset(rpt, mtd %in% mtd.sel & key %in% key.sel)

    ## method factors and lablels
    mtd.lvl <- c("mic", "b01", "b02", "b11", "b12", "pcc", "emc")
    mtd.lbl <- mtd.lvl
    names(mtd.lbl) <- mtd.lvl
    mtd.lbl["mic"]="MICE"
    mtd.lbl["drw"]="allele sample"
    mtd.lbl["mod"]="major allele"
    mtd.lbl["nfx"]="Netflix"
    mtd.lbl["zs0"]="beta(s)"
    mtd.lbl["b02"]="Soft"
    mtd.lbl["b12"]="Hard"
    mtd.lbl["pcc"]="Pairwise Complete"
    mtd.lbl["emc"]="EM"
    mtd.lbl <- mtd.lbl[mtd.lvl]
    mtd.lbl <- paste(mtd.lbl, "    ", sep="")
    rpt <- within(rpt, mtd <- factor(mtd, mtd.lvl, mtd.lbl))
    
    g <- ggplot(rpt, aes(x=M, y=val))
    g <- g + geom_line(aes(color=mtd), size=1.5)
    g <- g + facet_grid(key ~ tag, scales="free", labeller=.lb)
    g <- g + .th

    nfy <- length(unique(rpt$key))
    ufy <- 10 / nfy
    nfx <- length(unique(rpt$tag))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    options(bitmapType = 'cairo', device = 'pdf')
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy), scale = .85)

    invisible(rpt)
}
