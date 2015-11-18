args=(commandArgs(TRUE))

for (i in 1:3)
{
    cat(i, ":", args[[i]], "\n")
    eval(parse(text=args[[i]]))
}
miss.arg <- setdiff(c("datfile", "outfile"), ls())
if (length(miss.arg) > 0)
    stop("Missing arguments:", miss.arg, "\n")

libpaths=c(RLIB, .libPaths())    # RLIB is defined in snpfilt.pl
cat("libpaths =", libpaths, "\n")
require("caTools", lib.loc=libpaths)

dat <- read.table(gzfile(datfile), stringsAsFactors=F, sep=",", header=T)
dp.med <- median(dat$dptot)
dp.mad <- mad(dat$dptot)
qq.med <- median(dat$qual)
qq.mad <- mad(dat$qual)
qq.mean <- mean(dat$qual)
qq.sd <- sd(dat$qual)
qq.thres <- qq.mean - 3*qq.sd

cat(sprintf("Median DP=%f, MAD(DP)=%f\nMedian Qual=%f, Mean Qual=%f\nSD(QQ)=%f, MAD(QQ)=%f\n", 
    dp.med, dp.mad, qq.med, qq.mean, qq.sd, qq.mad))

n <- nrow(dat)
flt <- data.frame(ctg=dat[,1], cpos=dat[,2], flt=rep(10,n), seq=dat[,3])
ctg.list <- unique(dat[,1])

for (ctg in ctg.list)
{
    a.ind <- which(dat[,1] == ctg)
    ldp <- runmean(dat$dptot[a.ind], 201, alg="C", endrule="mean", align="center")
    ldp <- (ldp - dp.med)/dp.mad
    het <- 1 - (dat$dpFR[a.ind]+dat$dpRR[a.ind])/dat$dptot[a.ind]
    het.ind <- which(het > 0.3)
    dpH.ind <- which(ldp > 2)

    lqq <- runmean(dat$qual[a.ind]<qq.thres, 2001, alg="C", endrule="mean", align="center")
    qq.ind <- which(lqq > 0.025)    # more than 50/2000 sites have qual>qq.thres
    dpL.ind <- which(dat$dpFR[a.ind] < 1 | dat$dpRR[a.ind] < 1 | dat$dptot[a.ind] < 20)
    dp1.ind <- which(dat$dpFR[a.ind] < 10)
    mq.ind <- which(dat$mq[a.ind] < 58)

    tmp1 <- rep(0, length(a.ind))
    tmp1[dpH.ind] <- 1
    tmp2 <- runmax(tmp1, 201, endrule="max", align="center")
    rm(tmp1)

    tmp2[mq.ind] <-  2 
    tmp2[dpL.ind] <- 3
    tmp3 <- runmax(tmp2, 381, endrule="max", align="center")
    rm(tmp2)
     
    tmp3[dp1.ind] <- 4
    tmp3[het.ind] <- 5 
    tmp4 <- runmax(tmp3, 21, endrule="max", align="center")
    tmp4[qq.ind] <- 6

    flt[a.ind,3] <- tmp4
    rm(tmp3, tmp4)
}
outgz <- gzfile(outfile)
write.csv(flt, outgz, quote=F, row.names=F)

print(table(flt$flt))

