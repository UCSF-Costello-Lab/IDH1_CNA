######
## Set up
######
setwd('~/Documents/Costello_Lab/LG3/Infinium/methylumi/2016_12_processed/COPIED_FROM_CLUSTER/IDH1_plus_TCGA/')
# Load data
load("data/metadata.merged.filtered.Rdata")
load("big_data/data.clean.Rdata")
# Load code
source("code/samples_setup.R")
source("code/heatmap.3.R")
source("code/myplclust.R")
# Load libraries
library(limma)
library(gplots)
# Set new color scheme
#IDH_status_factor <- factor(IDH_status_tmp, levels=c("mut", "wt", "amp_mut", "del_mut", "gain_mut", "gain_wt", "pre_change","del_wt"))
IDH_status_factor <- factor(IDH_status_tmp, levels=c("gain_wt", "del_wt", "amp_mut", "del_mut", "gain_mut", "pre_change", "mut","wt"))
color_IDH <- mut_cols[IDH_status_factor]


######
## GCIMP
######
## Load data
load("big_data/data.betas.Rdata")

## Pull relevant data
# these are the closest probes to the 8 sites defined in Noushmehr et al Cancer Cell 2010 
gcimp450probes <- c('cg16849041', 'cg26399201', 'cg19320816', 'cg21245652', 'cg17403875', 'cg17120764', 'cg16257983', 'cg09088508')
gcimp450data <- data.betas[gcimp450probes, rev(c(idx.IDH, idx.notIDH, idx.normals, IDH.luchman.oneper61.pre, IDH.luchman.oneper.post))]
write.table(gcimp450data, file="data/gcimp.8sites.data.txt", sep="\t", quote=FALSE)

## Plot 8 sites
xlab <- rownames(gcimp450data)
ylab <- colnames(gcimp450data)
pdf("paperplots/gcimp.8sites.pdf", width=4, height=12)
par(cex=0.5)
data.factor <- apply(gcimp450data, 2, function(x) { as.numeric(cut(x, breaks=c(0,0.2,1))) } )
toplot <- list(x=1:length(xlab), y=1:length(ylab), z=as.matrix(data.factor))
col <- c("blue", "red")
image(toplot, axes=F, col=col)
axis(3, at=1:length(xlab), labels=xlab, tick=F, las=2)
axis(2, at=1:length(ylab), labels=ylab, tick=F, las=1)
dev.off()

## Pull relevant data
samp.toselect <- c(idx.LGG, idx.GBM) ## select most variable CpG sites using only TCGA tumors
samp.tocluster <- c(idx.LGG, idx.GBM, idx.costello, IDH.luchman.oneper61.pre, IDH.luchman.oneper.post); samp.header <- "TCGA.Costello.Luchman"
ridx <- which(rowSums(is.na(data.clean[ , samp.tocluster])) > 0)
toselect <- data.clean[-ridx, samp.toselect]
tocluster <- data.clean[-ridx, samp.tocluster]
sd.all <- apply(toselect, 1, sd); m0 <- "sd"
cutoff <- 0.005; count <- round(cutoff * nrow(toselect)); print(count)
probes <- names(sort(sd.all, decreasing=TRUE)[1:count])
sd.cutoff.data <- tocluster[probes, ]
plotme <- as.matrix(sd.cutoff.data) 

## Plot heatmap
m1 <- "euclidean"
m2 <- "ward.D"
title <- paste0("heatmap.", samp.header, ".", m0, cutoff, ".", count, "probes." , m1, ".", m2)
dist.m1 <- function(x) dist(x, method=m1)
hclust.m2 <- function(x) hclust(x, method=m2)
heatmap.annotations <- as.matrix(cbind(color_IDH[samp.tocluster], color_group[samp.tocluster]))
colnames(heatmap.annotations) <- c("IDH_status", "group")
pal <- colorRampPalette(c("blue","white","red"))
cols <- pal(20)
pdf(file=paste0("paperplots/", title, ".pdf"), height=8, width=16)
hm <- heatmap.3(plotme, trace="none", dist=dist.m1, hclustfun=hclust.m2, col=cols, density.info='none', labRow=NA, main=title, keysize=1, ColSideColors=heatmap.annotations, scale="none")
legend('left', col=c(mut_cols, group_cols), legend=c(levels(IDH_status_factor), levels(group_factor)), pch=15)
dev.off()

## Cleanup
rm(data.betas, toselect, tocluster, sd.cutoff.data, plotme)


######
## PCA of TCGA applied to our samples
######
# ## Samples selection
s.tobuild <- c(idx.TCGA.mut, idx.TCGA.wt)
s.toview <- c(idx.costello)
s.all <- c(s.tobuild, s.toview)
r.nona <- apply( data.clean[, s.all], 1, function(x) { all(!is.na(x)) })  ## we will use only probes without NA values

## Build on TCGA
tobuild.subdat.1 <- data.clean[which(r.nona), s.tobuild]
tobuild.mean <- apply(tobuild.subdat.1, 1, mean)
tobuild.sd <- apply(tobuild.subdat.1, 1, sd)
tobuild.subdat.2 <- tobuild.subdat.1 - tobuild.mean  ## mean center each probe
tobuild.subdat.3 <- tobuild.subdat.2 / tobuild.sd  ## scale sd per probe
fit <- prcomp(t(tobuild.subdat.3), center = FALSE, scale. = FALSE)
var_prop <- round((fit$sdev^2/sum(fit$sdev^2)), 4) * 100
save(fit, file="big_data/pca.fit.Rdata")

## Apply to Costello samples
toview.subdat.1 <- data.clean[which(r.nona), s.toview]
toview.subdat.2 <- toview.subdat.1 - tobuild.mean ## mean probe center based on TCGA
toview.subdat.3 <- toview.subdat.2 / tobuild.sd ## scale probe sd based on TCGA
testing <- predict(fit, t(toview.subdat.3))
save(testing, file="big_data/pca.testing.Rdata")

## Load
load("big_data/pca.fit.Rdata", verbose=TRUE)
load("big_data/pca.testing.Rdata", verbose=TRUE)
var_prop <- round((fit$sdev^2/sum(fit$sdev^2)), 4) * 100
s.tobuild <- c(idx.TCGA.mut, idx.TCGA.wt)
s.toview <- c(idx.costello)
s.toplot <- s.toview
s.toplot.names <- meta.all$sample_ID[idx.costello]
s.toplot <- which(rownames(testing) %in% s.toplot.names)

## Plot
flag <- IDH_status_factor[s.toview] %in% c("wt","mut")
flag.toplot <- intersect(which(!flag), s.toplot)
bg_cex <- .6; legend_cex <- 0.5; main_cex <- 0.5; text_cex <- 0.5; points_cex <- c(1,0.8)
pdf("paperplots/PCA.TCGAvsLG3.pdf", useDingbats=FALSE, width=5, height=5)
  plot(fit$x[,1], fit$x[,2], pch=c(6, 4)[(meta.all$IDH_status[s.tobuild] == "wt") +1], cex=bg_cex, col="grey50", xlab="", ylab="")
  title(xlab=paste0('PC 1 (', var_prop[1], '%)'), ylab=paste0('PC 2 (', var_prop[2], '%)') )
  legend("topleft", pch=c(6, 4), legend=c("mut","wt"), cex=legend_cex, title="TCGA data", col="grey50")
  # make flag for IDH change vs not
  points(testing[s.toplot,1], testing[s.toplot,2], pch=16, cex=points_cex[flag+1], col=color_IDH[s.toplot])
  #text(testing[which(!flag),1], testing[which(!flag),2], rownames(testing)[which(!flag)], pos=4, col=color_IDH[which(!flag)], cex=text_cex)
  #text(testing[flag.toplot,1], testing[flag.toplot,2], rownames(testing)[flag.toplot], pos=4, col=color_IDH[flag.toplot], cex=text_cex)
  legend('bottomright', col = mut_cols, legend = levels(IDH_status_factor), pch = 16, cex = legend_cex, title = 'IDH Status')
dev.off()

## Clean up
rm(fit, testing)


######
## Plot average methylation per sample
######
## Plot
avg.meth <- colMeans(data.clean[ , 1:nrow(meta.all)], na.rm=TRUE)
main=paste0("all mean ", nrow(data.clean))
pdf("paperplots/meth.average.pdf", useDingbats=FALSE, height=6, width=6)
boxplot(avg.meth[idx.adult], avg.meth[idx.notIDH], avg.meth[idx.IDH], names=c("normal", "LGG", "IDH"), pch=16, main=main, ylim=c(.4,.6))
points(rep(3, length(idx.IDH)), avg.meth[idx.IDH], col=color_IDH[idx.IDH], pch=16)
text(rep(3, length(idx.IDH)), avg.meth[idx.IDH], col=color_IDH[idx.IDH], labels=meta.all$sample_ID2[idx.IDH], pos=4, cex=0.5)
dev.off()



######
## Change in methylation from initial to recurrence
######
## Subtract
pre <- c(pri.v1.pairs.r1g2.notIDH, pri.v1.pairs.r1g3.notIDH, pri.v1.pairs.r1g4.notIDH, IDH.pre.all)
post <- c(rec1.v1.pairs.r1g2.notIDH, rec1.v1.pairs.r1g3.notIDH, rec1.v1.pairs.r1g4.notIDH, IDH.post.all)
dd <- data.clean[ , post] - data.clean[ , pre]
names(dd) <- sapply(names(dd), function(x) { strsplit(x, "_")[[1]][1] } )

## Subset probes on CGI status
cgi <- which(data.clean$Relation_to_UCSC_CpG_Island == "Island")
shore <- which(data.clean$Relation_to_UCSC_CpG_Island %in% c("N_Shore", "S_Shore"))
shelf <- which(data.clean$Relation_to_UCSC_CpG_Island %in% c("N_Shelf", "S_Shelf"))
shsh <- c(shore, shelf)
ncgi <- which(data.clean$Relation_to_UCSC_CpG_Island == "")

## Plot change by CGI status
col.list <- color_IDH[post]
dd1 <- apply(dd[ncgi,], 2, mean, na.rm=TRUE); yl.mean <- c(-max(abs(dd1)), max(abs(dd1)))  ## this one has the most extreme values, so set all y-axes to be the same
pdf("paperplots/differences.pdf", height=5, width=12)
par(mfrow=c(1,4), cex=0.5)
dd1 <- apply(dd, 2, mean, na.rm=TRUE); barplot(dd1, las=2, col=col.list, main=paste0("all mean ", nrow(dd)), ylim=yl.mean)
dd1 <- apply(dd[cgi,], 2, mean, na.rm=TRUE); barplot(dd1, las=2, col=col.list, main=paste0("CGI mean ", length(cgi)), ylim=yl.mean)
dd1 <- apply(dd[shsh,], 2, mean, na.rm=TRUE); barplot(dd1, las=2, col=col.list, main=paste0("shore&shelf mean ", length(shsh)), ylim=yl.mean)
dd1 <- apply(dd[ncgi,], 2, mean, na.rm=TRUE); barplot(dd1, las=2, col=col.list, main=paste0("notCGI mean ", length(ncgi)), ylim=yl.mean)
dev.off()

## Cluster
cos.pre <- c(pri.v1.pairs.r1g2.notIDH, pri.v1.pairs.r1g3.notIDH, pri.v1.pairs.r1g4.notIDH, IDH.pre.all)
cos.post <- c(rec1.v1.pairs.r1g2.notIDH, rec1.v1.pairs.r1g3.notIDH, rec1.v1.pairs.r1g4.notIDH, IDH.post.all)
dd.cos <- data.clean[, cos.post] - data.clean[,cos.pre]
names(dd.cos) <- sapply(names(dd.cos), function(x) { strsplit(x, "_")[[1]][1] } )
col.list <- color_IDH[cos.post]

samp.tocluster <- 1:ncol(dd.cos)
ridx <- which(rowSums(is.na(dd[ , samp.tocluster])) > 0)  
tocluster <- dd[-ridx, samp.tocluster]

## Select probes
sd.all <- apply(tocluster, 1, sd); m0 <- "sd"
cutoff <- 0.005; count <- round(cutoff * nrow(tocluster)); print(count)
sd.cutoff.data <- tocluster[names(sort(sd.all, decreasing=TRUE)[1:count]), ]
plotme <- as.matrix(sd.cutoff.data) 

## Plot heatmap
m1 <- "euclidean"
m2 <- "ward.D"
title <- paste0("differences.", m0, cutoff, ".", count, "probes." , m1, ".", m2, ".heatmap")
dist.m1 <- function(x) dist(x, method=m1)
hclust.m2 <- function(x) hclust(x, method=m2)
pal <- colorRampPalette(c("blue","white","red"))
cols <- pal(20)
pdf(paste0("paperplots/",title,".pdf"))
heatmap.2(as.matrix(sd.cutoff.data), trace="none", dist=dist.m1, hclustfun=hclust.m2, keysize=1, density.info="none", col=cols, labRow=NA, main=title, scale="none", ColSideColors=col.list)
dev.off()

## PCA of methylation change
r.nona <- apply(dd, 1, function(x) { all(!is.na(x)) })
dd.pca1 <- dd[r.nona, ]
dd.pca2 <- dd.pca1 - apply(dd.pca1, 1, mean)
dd.pca3 <- dd.pca2 / apply(dd.pca1, 1, sd)
fit <- prcomp(t(dd.pca3), center = FALSE, scale. = FALSE)
var_prop <- round((fit$sdev^2/sum(fit$sdev^2)), 4) * 100
pdf("paperplots/methdiff.PCA.pdf", useDingbats=FALSE, width=5, height=5)
plot(fit$x[,1], fit$x[,2], pch = 20, xlab = paste0('PC 1 (', var_prop[1], '%)'),ylab = paste0('PC 2 (', var_prop[2], '%)'), col=col.list)
text(fit$x[,1], fit$x[,2], names(fit$x[,1]), pos=4, col=col.list)
dev.off()



####
# Correlations between differences
####
## From: http://www.sthda.com/english/wiki/print.php?id=78
## removed "pmat" from the code
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
    )
}
cor.diffs <- cor(dd, use="complete.obs", method="spearman")
flat.diffs <- flattenCorrMatrix(cor.diffs)
range(flat.diffs$cor)
pdf("paperplots/differences.correlation.hist.pdf", width=5, height=4)
hist(flat.diffs$cor, breaks=20, xlim=c(-0.6,0.6))
text(flat.diffs$cor, 10, paste0(flat.diffs$row,"_",flat.diffs$column), srt=90, pos=3, cex=0.5)
dev.off()


######
## Plot methylation density plots for each patient
######
## Plot
patients <- pat.IDH
pdf(file="paperplots/density_by_patient.pdf", width=10, height=3)
par(mfrow=c(1,5))
for(i in 1:length(patients)) {
	k <- which(meta.all$patient_ID==patients[i])
	plot(0,0,type="n",xlim=c(-0.05,1.05),ylim=c(0,5),xaxt="n", xlab="", ylab="")
	axis(1, at=c(0,1))
	for(l in 1:length(k)) {	lines(density(data.clean[ , k[l]], na.rm=T), col=color_IDH[k[l]]) }
	title(main=patients[i], line=-1)
}
dev.off()


####
# Cluster cultures/xenografts with tumors
####
## Select samples
samp.tocluster <- c(pri.v1.pairs.notIDH, rec1.v1.pairs.notIDH, IDH.pre.all, IDH.post.all, IDH.luchman.oneper61.pre, IDH.luchman.oneper.post)
header <- "cluster.single.LG3.IDH.LuchOnePer61"

## Remove NA probes
ridx <- which(rowSums(is.na(data.clean[ , samp.tocluster])) > 0)  
tocluster <- data.clean[-ridx, samp.tocluster]

## Set clustering parameters
m0 <- "sd"
m1 <- "euclidean"
m2 <- "ward.D"
cutoff <- 0.005
count <- round(cutoff * nrow(tocluster)); print(count)

## Select probes
sd.all <- apply(tocluster, 1, sd)
sd.cutoff.data <- tocluster[names(sort(sd.all, decreasing=TRUE)[1:count]), ]

## Plot heatmap
title <- paste0(header, ".", m0, cutoff, ".", count, "probes." , m1, ".", m2, ".heatmap")
dist.m1 <- function(x) dist(x, method=m1)
hclust.m2 <- function(x) hclust(x, method=m2)
pal <- colorRampPalette(c("blue","white","red"))
cols <- pal(20)
pdf(paste0("paperplots/",title,".pdf"))
heatmap.2(as.matrix(sd.cutoff.data), trace="none", dist=dist.m1, hclustfun=hclust.m2, keysize=1, density.info="none", col=cols, labRow=NA, main=title, scale="none", ColSideColors= color_IDH[samp.tocluster])
dev.off()



######
## Supervised differential methylation
######
## Calculate methylation change
# pull beta values
het_p <- data.clean[ , c(pri.v1.pairs.r1g2.notIDH, pri.v1.pairs.r1g3.notIDH, pri.v1.pairs.r1g4.notIDH)]
het_r <- data.clean[ , c(rec1.v1.pairs.r1g2.notIDH, rec1.v1.pairs.r1g3.notIDH, rec1.v1.pairs.r1g4.notIDH)]
idh_p <- data.clean[ , IDH.pre.sub]
idh_r <- data.clean[ , IDH.post.sub]
het_d <- het_r - het_p
idh_d <- idh_r - idh_p
dd <- cbind.data.frame(idh_d, het_d)

# logit transform
trans <- function(x) log(x/(1-x))
t_het_p <- trans(het_p)
t_het_r <- trans(het_r)
t_idh_p <- trans(idh_p)
t_idh_r <- trans(idh_r)
t_het_d <- t_het_r - t_het_p
t_idh_d <- t_idh_r - t_idh_p
t_dd <- cbind.data.frame(t_idh_d, t_het_d)

# calculate averages for plotting
het_d_avg <- rowMeans(het_d, na.rm=TRUE)
idh_d_avg <- rowMeans(idh_d, na.rm=TRUE)

## Run limma
design <- model.matrix(~factor(c(rep(1,ncol(t_idh_d)), rep(2,ncol(t_het_d)))))
colnames(design) <- c("intercept", "slope")
eset <- t_dd
fit <- lmFit(eset, design)
fit2 <- ebayes(fit)
fit2.pvalues <- fit2$p[,2]
pval <- fit2.pvalues
pval[is.na(pval)] <- 1
pval.bh <- p.adjust(pval, method="BH")

## Set significance
pcut <- 1E-2 # 29018
idx <- unname(which((pval.bh < pcut) & (idh_d_avg - het_d_avg < -0.2) & (idh_d_avg < -0.2) )); length(idx)
sig.color <- "green"

## Plots
pdf(paste0("paperplots/limma.volcano.bh",pcut,".pdf"))
plot(idh_d_avg - het_d_avg, -log10(pval.bh), pch=20, cex=0.7, ylab=expression(-log[10]~pBH), xlab= expression(Delta~bar(Delta~beta)), col="grey", main="IDH vs het", xlim=c(-1,1), yaxt='n', type="n")
#points((idh_d_avg - het_d_avg)[idx], -log10(pval.bh)[idx], pch=20, cex=1, col=sig.color)
axis(2, las=1)
abline(v=c(-.2, .2), lty=2, lwd=2, col="grey")
abline(h=-log10(pcut), lty=2, lwd=2, col="grey")
dev.off()
png(paste0("paperplots/limma.volcano.bh",pcut,".png"), units="in", width=8, height=8, res=300)
plot(idh_d_avg - het_d_avg, -log10(pval.bh), pch=20, cex=0.7, ylab=expression(-log[10]~pBH), xlab= expression(Delta~bar(Delta~beta)), col="grey", main="IDH vs het", xlim=c(-1,1), yaxt='n')
points((idh_d_avg - het_d_avg)[idx], -log10(pval.bh)[idx], pch=20, cex=1, col=sig.color)
axis(2, las=1)
dev.off()

l_idh_d_avg <- rowMeans(data.clean[, IDH.luchman.oneper.post] - data.clean[, IDH.luchman.oneper61.pre], na.rm=TRUE)
l_idh_d_avg.sig <- l_idh_d_avg[idx]; l_idh_d_avg.sig.med <- median(l_idh_d_avg.sig, na.rm=TRUE)
l_idh_d_avg.bkg <- l_idh_d_avg[-idx]; l_idh_d_avg.bkg.med <- median(l_idh_d_avg.bkg, na.rm=TRUE)
pdf(paste0("paperplots/limma.density.bh",pcut,".pdf"))
plot(density(l_idh_d_avg.bkg, na.rm=TRUE), xlim=c(-1,1))
lines(density(l_idh_d_avg.sig), col=sig.color)
legend("topleft", legend=paste0(c(pcut, "background"), "; num probes ", c(length(l_idh_d_avg.sig), length(l_idh_d_avg.bkg)), "; median ", round(c(l_idh_d_avg.sig.med, l_idh_d_avg.bkg.med), digits=6)), col=c(sig.color, "black"), lwd=1)
text(-1, 20, labels=paste0("delta median: ", round(l_idh_d_avg.sig.med-l_idh_d_avg.bkg.med, digits=6)), pos=4 )
dev.off()


load("big_data/annotations_with_gencode_1500_TSS_1000_2.RData", verbose=TRUE)
cgi <- which(data.clean$Relation_to_UCSC_CpG_Island == "Island")
shore <- which(data.clean$Relation_to_UCSC_CpG_Island %in% c("N_Shore", "S_Shore"))
shelf <- which(data.clean$Relation_to_UCSC_CpG_Island %in% c("N_Shelf", "S_Shelf"))
shsh <- c(shore, shelf)
ncgi <- which(data.clean$Relation_to_UCSC_CpG_Island == "")
prom.probes <- unique(anno$IlmnID[which(!is.na(anno$gene_name))])
prom <- which(data.clean$IlmnID %in% prom.probes)
prom.cgi <- intersect(prom,cgi)
prom.not.cgi <- setdiff(prom,cgi)

ll <-     list(cgi,  shsh,  ncgi,  prom,  prom.not.cgi )
ll.names <- c("cgi","shsh","ncgi","prom","prom.not.cgi")
boot.summary <- data.frame(row.names=ll.names)
nrep <- 10000
for(i in 1:length(ll)) {
	print(i)
    l.idx <- ll[[i]]
    inter <- intersect(l.idx, idx)
    perc.in <- length(inter)/length(idx)*100
    perc.bg <- length(l.idx)/length(pval.bh)*100
    boot.l <- replicate(nrep, sum(sample(data.clean$IlmnID, length(idx)) %in% data.clean$IlmnID[l.idx]))
    boot.l.high <- sum(length(inter) > boot.l)
    boot.l.low <- sum(length(inter) < boot.l)
    boot.summary[i, "perc_in"] <- perc.in
    boot.summary[i, "perc_bg"] <- perc.bg
    boot.summary[i, "boot_high"] <- boot.l.high
    boot.summary[i, "boot_low"] <- boot.l.low
}
pdf("paperplots/limma.genomic.space.sig.pdf")
barplot(t(as.matrix(boot.summary[, c("perc_in", "perc_bg")])), beside=TRUE, col=c("green","black"), legend=TRUE, las=2)
dev.off()
write.table(boot.summary, file="data/limma.genomic.space.sig.txt", sep="\t", quote=FALSE)
