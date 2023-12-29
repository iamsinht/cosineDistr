library(dplyr)
library(gridExtra)
library(CMAPToolkit)
library(pracma)
library(MASS)


# This script is used to illustrate a number of the properties of cosine similarity in high dimensional datasets

cosoutdir <- "~/Work/bhk/analysis/cosine"

if (!file.exists(file.path(cosoutdir, "cosineVsDim_n=5k.rds"))){
  a <- cosineVsDim(dims=c(seq(1,4), seq(5, 250, 5), seq(275, 1000, 25), seq(1050,5000,50)))
  saveRDS(a, file=file.path(cosoutdir, "cosineVsDim_n=5k.rds"))
} else {
  a <- readRDS(file=file.path(cosoutdir, "cosineVsDim_n=5k.rds"))
}

# Figure 2b:
pdf(file.path(cosoutdir, "cosineVsDim_variance_n=2k.pdf"), width=7, height=5)
ggplot(a[a$dim < 2000,], aes(x=dim, y=sd, color="blue")) + geom_line(size=2, color="blue") +
  geom_line(aes(y=1/sqrt(dim)), size=0.5, color="red") + 
  scale_color_manual(name="Data", breaks=c("Empirical", "Theory"), values=c("Empirical"="blue", "Fit"="red")) + 
  theme_minimal() + xlab("Dimension") + ylab("Stdev of Cosine similarity") + 
  ggtitle("Sdev of cosine of dimension N") + theme(text=element_text(size=14)) + ylim(c(0, 0.5))
dev.off()

# Log
pdf(file.path(cosoutdir, "cosineVsDim_variance_n=2k_Log.pdf"), width=8, height=6)
ggplot(a[a$dim < 2000,], aes(x=dim, y=sd)) + geom_line(linewidth=3, aes(color="Empirical")) +
  geom_line(aes(y=1/sqrt(dim), color="Theory"), linewidth=1) + theme_minimal() +
  scale_color_manual(name="Data", breaks=c("Empirical", "Theory"), values=c("Empirical"="blue", "Theory"="red")) + 
  xlab("Dimension") + ylab("Stdev of Cosine similarity") + theme(legend.position="right") + 
  ggtitle("Standard deviation of cosine of unit normals of dimension N") + theme(text=element_text(size=14)) + ylim(c(0, 0.5)) + 
  scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
dev.off()

pdf(file.path(cosoutdir, "cosineVsDim_variance_inset_n=1k.pdf"), width=10, height=6)
ggplot(a[a$dim <= 1000,], aes(x=dim)) + geom_line(aes(y=sd), size=2, color="blue") +
  geom_line(aes(y=1/sqrt(dim)), size=0.5, color="red") + theme_minimal() + xlab("Dimension") + ylab("Sdev of Cosine similarity") + 
  ggtitle("Sdev of cosine distribution of vectors of standard normals of dimension N") + theme(text=element_text(size=14)) 
dev.off()

#### Illustration of rescaling an axis ####
ii <- 100
x0 <- t(mvrnorm(n=1000, mu=numeric(ii), Sigma=diag(ii)))
xcos0 <- CMAPToolkit::cosine(x0, x0)

sigmat <- diag(ii)
sigmat[1,1] <- 25
x1 <- t(mvrnorm(n=1000, mu=numeric(ii), Sigma=sigmat))
xcos1 <- CMAPToolkit::cosine(x1, x1)

plot(density(xcos0[upper.tri(xcos0)], bw=0.01), col="orange", xlim=c(-1,1), lwd=3, xlab="Cosine Similarity", ylab="Density",
     main="Distribution of cosine similarities, N=100 dimensions")
lines(density(xcos1[upper.tri(xcos1)], bw=0.01), col="darkblue", lwd=3)
legend(x="topleft", legend=c("Standard normal", "Rescaled dim1 x 5"), lwd=c(2,2), col=c("orange", "darkblue"))


#### Rescaled sampling: ####
# Supp Fig for minimum
if (!file.exists(file=file.path(cosoutdir, "rescaleVarianceDf.rds"))) {
  multfact <- seq(log(0.01), log(100), log(100)/50)
  scaledf <- data.frame(scaleFactor=numeric(), mean=numeric(), sd=numeric())
  
  ii <- 100
  sigmat <- diag(ii)
  for (myfact in multfact){
    print(which(multfact == myfact))
    sigmat[1,1] <- exp(myfact)
    x <- t(mvrnorm(n=1000, mu=numeric(ii), Sigma=sigmat))
    xcos <- CMAPToolkit::cosine(x, x)
    
    scaledf <- rbind(scaledf, data.frame(scaleFactor=exp(myfact), mean=mean(xcos[upper.tri(xcos)]), sd=sd(xcos[upper.tri(xcos)])))
  }
  saveRDS(object = scaledf, file=file.path(cosoutdir, "rescaleVarianceDf.rds"))
} else {
  scaledf <- readRDS(file=file.path(cosoutdir, "rescaleVarianceDf.rds"))
}

pdf(file.path(cosoutdir, "cosineStdevRescaling_N=100.pdf"), width=8, height=6)
semilogx(scaledf$scaleFactor, scaledf$sd, lwd=2, col="blueviolet", xlab="Scale Factor", 
         ylab="Stdev", main="Stdev of Cosine Similarity with rescaling 1 of N=100 dimensions", pch=16)
dev.off()

pdf(file.path(cosoutdir, "rescaleNormSchematic.pdf"), width=6, height=6)
plot(rnorm(1000), rnorm(1000), pch=16, xlim=c(-5, 5), ylim=c(-5,5), xlab="Dim1", ylab="Dim2")
dev.off()
pdf(file.path(cosoutdir, "rescaleNormSchematic2.pdf"), width=6, height=6)
plot(2*rnorm(1000), rnorm(1000), pch=16, xlim=c(-5, 5), ylim=c(-5,5), xlab="Dim1", ylab="Dim2")
dev.off()


#### Compare distributions ####
x10 <- getCosineDist(mydim=10, N=1000, mkfig=0)
x100 <- getCosineDist(mydim=100, N=1000, mkfig=0)
x1000 <- getCosineDist(mydim=1000, N=1000, mkfig=0)

pdf(file.path(cosoutdir, "cosineVsDim_dists.pdf"), width=8, height=6)
plot(density(x1000[upper.tri(x1000)], bw=0.01), col="forestgreen", lwd=3, xlim=c(-1,1), 
     main="Distribution of cosine similarities for standard normals of dimension N")
lines(density(x10[upper.tri(x10)], bw=0.01), col="red", lwd=3)
lines(density(x100[upper.tri(x100)], bw=0.01), col="blue", lwd=3)
legend(x="topleft", legend=c("N=10", "N=100", "N=1000"), lwd=c(4,4,4), col=c("red", "darkblue", "forestgreen"))
dev.off()



#### Characterize L1000 distribution: ####
datapath <- "~/Work/bhk/data/l1k/2020/"
l1kmeta <- read_l1k_meta(datapath, version=2020)
attach(l1kmeta)

mycells <- c("A375", "A549", "HELA", "HEPG2", "MCF7", "NPC", "PC3", "U2OS")
l1kCellDf <- data.frame(cellid = character(), mean = numeric(), sd = numeric())


pdf(file.path(cosoutdir, "L1K_landmark_cos.pdf"), width=10, height=8)
cos978 <- getCosineDist(mydim=978, N=1000, mkfig=0)
plot(density(cos978[upper.tri(cos978)], bw=0.01), col="black", lwd=2, xlim=c(-0.5, 0.5), xlab="Cosine Similarity", ylab="Density", main="L1000 Cosine Similarities, landmark")

for (cella in mycells){
  print(cella)
  ds <- parse_gctx(get_level5_ds(datapath), rid = geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], 
                   cid=siginfo$sig_id[siginfo$cell_id == cella & siginfo$pert_type == "trt_cp"])
  ix <- sample(length(ds@cid), 2000)
  dscos <- CMAPToolkit::cosine(ds@mat[, ix], ds@mat[, ix])
  
  l1kCellDf <- rbind(l1kCellDf, data.frame(cellid = cella, mean = mean(dscos[upper.tri(dscos)]), sd = sd(dscos[upper.tri(dscos)])))
  lines(density(dscos[upper.tri(dscos)], bw=0.01), col=rainbow(8)[match(cella, mycells)], lwd=2)
}
legend(x="topright", legend=c("Null", mycells), col=c("black", rainbow(8)), lwd=3)
dev.off()

l1kCellDf$effectDim <- a$dim[sapply(l1kCellDf$sd, FUN=function(x) max(which(a$sd > x)))+1]


printdf <- rbind(data.frame(cellid="Std Norm 978", mean=mean(cos978[upper.tri(cos978)]), sd=sd(cos978[upper.tri(cos978)]), effectDim=978), l1kCellDf[1:8,])
pdf(file.path(cosoutdir, "L1K_landmark_table.pdf"), width=4, height=6)
grid.table(printdf %>% mutate(across(where(is.numeric), ~ round(., 5))), rows=NULL)
dev.off()


#### Explore dependence of additional genes

ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$cell_id == "HEPG2" & siginfo$pert_type == "trt_cp"])
ds <- subset_gct(ds, cid=sample(ds@cid, 1000))

l1kdf <- data.frame(ngenes=numeric(), mean=numeric(), sd=numeric())

for (k in seq(25)){
  for (ii in c(50, 100, 200, 400, 600, 800, 978)){
    print(ii)
    gsub <- sample(geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], ii)
    ix <- match(gsub, ds@rid)
    tempcos <- cosineDistr::cosine(ds@mat[ix, ], ds@mat[ix, ])
    l1kdf <- rbind(l1kdf, data.frame(ngenes =ii, mean=mean(tempcos[upper.tri(tempcos)]), sd=sd(tempcos[upper.tri(tempcos)])))
  }
  
  for (jj in c(500, 1000, 2000, 5000, 11350)){
    gsub <- union(geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], sample(geneinfo$pr_gene_id[geneinfo$pr_is_lm == 0], jj))
    ix <- match(gsub, ds@rid)
    print(sprintf("jj = %d, length(ix) == %d", jj, length(ix)))
    tempcos <- cosineDistr::cosine(ds@mat[ix, ], ds@mat[ix, ])
    l1kdf <- rbind(l1kdf, data.frame(ngenes = 978+jj, mean=mean(tempcos[upper.tri(tempcos)]), sd=sd(tempcos[upper.tri(tempcos)])))
  }
}

l1kdf$effectDim <- (1/l1kdf$sd)^2

saveRDS(l1kdf, file=file.path(cosoutdir, "L1K_HEPG2_geneSubsample_MeanVsDim.rds"))

pdf(file.path(cosoutdir, "L1K_HEPG2_geneSubsample_sdcos.pdf"), width=8, height=6)
semilogx(l1kdf$ngenes, l1kdf$sd, pch=16, col="chartreuse4", cex=2, xlab="Number of Genes", ylab="SD Cos", main="HEPG2 N=1000 signatures SD Cos vs N genes", 
         ylim=c(0, 0.25))
dev.off()
pdf(file.path(cosoutdir, "L1K_HEPG2_geneSubsample_effectDim.pdf"), width=8, height=6)
semilogx(l1kdf$ngenes, l1kdf$effectDim, pch=16, col="purple", cex=2, xlab="Number of Genes", ylab="Effective Dimension", main="HEPG2 N=1000 signatures EffectDim vs N genes")
dev.off()