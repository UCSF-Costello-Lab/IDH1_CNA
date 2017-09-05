## Load config file
source("code/methylumi.setup.R")
setwd(savedir)

## Read in relevant package
library(methylumi)
print(packageVersion("methylumi"))
print(sessionInfo())

### Load metadata, confirm sample names & file locations
load(paste0(metadir, metadata))
print("*** confirm sample IDs can be column names -- output should be TRUE")
identical(meta.all$sample_ID2, make.names(meta.all$sample_ID2, unique=TRUE))

filenames <- list.files(path=idatdir)
sample.green <- paste0(meta.all$barcode,"_Grn.idat")
sample.red <- paste0(meta.all$barcode,"_Red.idat")
sample.combined <- matrix(rbind(sample.green,sample.red))
match.samples <- match(sample.combined,filenames)
print("*** confirm all files were found -- output should be 0:")
sum(is.na(match.samples))  ## should be 0
filenames <- filenames[match.samples]


### Data pre-processing & background correction
## Run main processing step
dat1 <- methylumIDAT(meta.all$barcode,idatPath=idatdir)

## Look at rate of probes to throw out that have high p-values
dat1cm <- colMeans(pvals(dat1)>0.05)*100
summary(dat1cm)
save(dat1cm, file="dat1cm.prefilter.Rdata")
badsamples.count <- sum(dat1cm > 0.5) ## an arbitrary but reasonable cutoff
print("The following is the count of samples with >0.5% of probes with high p-values (>0.05):")
print(badsamples.count)
badsamples.idx <- which(dat1cm > 0.5)
badsamples.drop <- setdiff(badsamples.idx, which(meta.all$source %in% c("Costello", "Luchman"))) ## keep all of our samples
badsamples.keep <- setdiff(badsamples.idx, badsamples.drop)
if(length(badsamples.keep) > 0) {
	print("These samples have >0.5% of probes with high p-values (>0.05) but will be retained:")
	print(meta.all$sample_ID[badsamples.keep])
}
if(length(badsamples.drop) > 0) {
	print("These samples have >0.5% of probes with high p-values (>0.05) and will be excluded:")
	print(meta.all$sample_ID[badsamples.drop])

	## remove bad samples
	meta.all <- meta.all[-badsamples.drop, ]

	## Re-run main processing without those samples
	dat1 <- methylumIDAT(meta.all$barcode,idatPath=idatdir)
	dat1cm <- colMeans(pvals(dat1)>0.05)*100
	summary(dat1cm)
}
save(meta.all, file=paste0(metadir, "metadata.merged.filtered.Rdata"))
save(dat1cm, file="dat1cm.postfilter.Rdata")


## Run the background correction
dat1.noob <- methylumi.bgcorr(dat1)
## Strip unnecessary features
dat1.noob <- stripOOB(dat1.noob)
## Normalize probes
dat1.nn <- normalizeMethyLumiSet(dat1.noob)
## Pull betas out
dat1.betas <- betas(dat1.nn)
## Pull p-values out
pvals1 <- pvals(dat1)
## Make NA betas with high p-values
dat1.betas[pvals1>0.05] <- NA

## Fix column names to be patient_sample rather than barcodes
print("*** confirm that samples are in the correct order before renaming -- output should be TRUE:")
identical(colnames(dat1.betas), meta.all$barcode)   ## should return TRUE
colnames(dat1.betas) <- meta.all$sample_ID2

# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

### Save data to file
## All betas
data.betas <- cbind.data.frame(dat1.betas)
save(data.betas, file="data.betas.Rdata")

## Annotate
annotations <- read.delim("/home/mazort/LG3/infinium/HumanMethylation450_15017482_v.1.2.csv", header=TRUE, sep=",", as.is=TRUE)
match.annotations <- match(rownames(data.betas),annotations$IlmnID)
annotation.table <- annotations[match.annotations,c("IlmnID", "CHR", "MAPINFO", "Probe_SNPs", "Probe_SNPs_10", "Random_Loci", "Methyl27_Loci", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island", "Phantom", "DMR", "Enhancer", "HMM_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DHS")]
data.betas.anno <- cbind(data.betas, annotation.table)
save(data.betas.anno, file="data.betas.anno.Rdata")

## Filter out probes
# Probes on chrX/chrY
chrX <- which(data.betas.anno$CHR == "X")
chrY <- which(data.betas.anno$CHR == "Y")
# SNPs: supplementary SNP list
ill.SNPsupp <- read.delim("/home/mazort/LG3/infinium/annotations/humanmethylation450_15017482_v.1.2.snpupdate.table.v3.txt", header=T, sep="\t", as.is=T)
snps.supp <- which(data.betas.anno$IlmnID %in% ill.SNPsupp$TargetID)  ## includes TARGET SNPs with Distanace = 0  --> 153,133 probes
#snps.supp <- which(data.betas.anno$IlmnID %in% ill.SNPsupp$TargetID[ill.SNPsupp$MAF >= 0.1])
## could have MAF cutoff (ie 80,266 with MAF > 0.1)
# Multimappers
price <- read.table("/home/mazort/LG3/infinium/annotations/GPL16304-47833.txt", header=T, sep="\t", as.is=T)
multi <- which(data.betas.anno$IlmnID %in% price$ID[price$XY_Hits=="XY_YES" | price$Autosomal_Hits=="A_YES"])
# probes that don't map to genome
nomap <- which(data.betas.anno$CHR == "") ## 65 probes
# Make new data frame that EXCLUDES all sex chromosomes AND all SNP probes AND all multimappers
toexclude <- unique(c(chrX, chrY, snps.supp, multi, nomap))  ## 188,235 probes to remove
data.clean <- data.betas.anno[-toexclude, ]  ## 297,342 probes
save(data.clean, file="data.clean.Rdata")
