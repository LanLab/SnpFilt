args=(commandArgs(TRUE))

for (i in 1:6)
{
    cat(i, ":", args[[i]], "\n")
    eval(parse(text=args[[i]]))
}
miss.arg <- setdiff(c("datfile", "outpdf", "snpfile", "wsize", "index"), ls())
if (length(miss.arg) > 0)
    stop("Missing arguments:", miss.arg, "\n")

libpaths=c(RLIB, .libPaths())    # RLIB is defined in snpfilt.pl
cat("libpaths =", libpaths, "\n")
require("caTools", lib.loc=libpaths)

dat <- read.csv(gzfile(datfile), stringsAsFactors=F, header=T)
dp.med <- median(dat$dptot)
dp.mad <- mad(dat$dptot)
qq.mean <- mean(dat$qual)
qq.sd <- sd(dat$qual)
qq.thres <- qq.mean - 3*qq.sd

snp <- read.csv(snpfile, header=T, stringsAsFactors=F, comment="#")
if (nrow(snp) == 0)
    stop ("Empty SNPfile\n")

if (nchar(gsub("\\s", "", index))==0)
{
    i.list <- 1:nrow(snp)
} else {
    i.list <- as.numeric(unlist(strsplit(index, ",")))
}

i.list <- i.list[!is.na(i.list)]
if (length(i.list) == 0)
    stop("No valid SNPs requested\n")
cat("i =", i.list, "\n")

pdf(outpdf, 6, 8)
par(mfrow=c(3,1))
par(mar=c(3,4,1,1))

old.ctg=""
old.cpos=0
snp.col=c("black", rgb(0.5,0.5,0.5,0.5))

for (i in i.list)
{
    if (is.na(i) || i > nrow(snp))
        stop("Error: invalid i\n")
    if (snp$ctg[i]==old.ctg && abs(snp$cpos[i]-old.cpos) < wsize)
        next

    ind <- which(dat[,1]==snp$ctg[i] & dat[,2] > (snp$cpos[i]-wsize) & dat[,2] < (snp$cpos[i]+wsize))
    snp.ind <- which(snp$ctg==snp$ctg[i])
    snp.type <- (snp$flt[snp.ind] > 0) + 1
    
    plot(dat[ind,2], dat$dptot[ind], type="l", xaxs="i", ylim=c(0,2.2*dp.med),
        xlab=snp$ctg[snp.ind[1]], ylab="Coverage")
    ldp <- runmean(dat$dptot[ind], 201, alg="C", endrule="mean", align="center")
    lines(dat[ind,2], ldp, col=rgb(0,1,0,0.7))
    lines(dat[ind,2], dat$dpFR[ind], type="l", col="blue")
    
    abline(h=10, lty=2, col="grey70")
    abline(h=20, lty=2, col="black")
    abline(h=dp.med + 2*dp.mad, lty=2, col=rgb(0,1,0,0.7))

    err <- dat$dptot[ind] - (dat$dpFR[ind]+dat$dpRR[ind])
    lines(dat[ind,2], err, col=rgb(1,0,0,0.6))
    legend("topright", c("DPtot", "Running mean", "DP-ForwardRef", "DP-nonConsensus"),
        bty="n", cex=0.6, lty=1, col=c("black", rgb(0,1,0,0.7), "blue", rgb(1,0,0,0.6)))

    abline(v=snp$cpos[snp.ind], col=snp.col[snp.type])

    plot(dat[ind,2], dat$qual[ind], type="l", xaxs="i", ylim=c(0,300),
        xlab="Contig position", ylab="Base quality")
    abline(h=qq.thres, col="grey70", lty=2)
    abline(v=snp$cpos[snp.ind], col=snp.col[snp.type])
    qq.prop <- mean(dat$qual[ind] < qq.thres)
    qq.str <- sprintf("Proportion Q < qthres = %.3e", qq.prop)
    legend("bottomright", qq.str, bty="n", cex=0.6)

    plot(dat[ind,2], dat$mq[ind], type="n", xaxs="i", ylim=c(0,65),
        xlab="Contig position", ylab="Mapping quality")
    abline(h=58, col="grey70", lty=2)
    lines(dat[ind,2], dat$mq[ind])
    abline(v=snp$cpos[snp.ind], col=snp.col[snp.type])
    legend("bottomleft", snp$ctg[i], bty="n", cex=0.6)

    old.ctg  <- snp$ctg[i]
    old.cpos <- snp$cpos[i]
}
dev.off()
    
