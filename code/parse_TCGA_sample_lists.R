### Find which TCGA samples have 450k + somatic mutations + copy number

## Make sample list for each of the 3 LGG/GBM data types
# Methylation downloaded from TCGA Data Matrix 10/29/15
file1 <- "/costellolab/mazort/TCGA_data/LGG_450k/METADATA/JHU_USC__HumanMethylation450/jhu-usc.edu_LGG.HumanMethylation450.1.12.0.sdrf.txt"
dat1 <- read.delim(file1, as.is=TRUE)
LGG_meth_pre <- unique(sapply(dat1$Comment..TCGA.Barcode., function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ))
LGG_meth <- unique(unname(sapply(LGG_meth_pre, function(x) { substr(x,1,nchar(x)-1) })))

file1 <- "/costellolab/mazort/TCGA_data/GBM_450k/METADATA/JHU_USC__HumanMethylation450/jhu-usc.edu_GBM.HumanMethylation450.1.6.0.sdrf.txt"
dat1 <- read.delim(file1, as.is=TRUE)
GBM_meth_pre <- unique(sapply(dat1$Comment..TCGA.Barcode., function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ))
GBM_meth <- unique(unname(sapply(GBM_meth_pre, function(x) { substr(x,1,nchar(x)-1) })))


# SNP Copy Number downloaded from CBio Portal 11/4/15
file1 <- "downloaded_data/CBIO_LGG_log2.txt"
file2 <- "downloaded_data/CBIO_LGG_GISTIC.txt"
dat1 <- read.delim(file1, as.is=TRUE, skip=1)
dat2 <- read.delim(file2, as.is=TRUE, skip=1)
identical(dat1$COMMON[which(!(is.na(dat1$IDH1)))], dat2$COMMON[which(!(is.na(dat2$IDH1)))])
LGG_SNP <- dat1$COMMON[which(!(is.na(dat1$IDH1)))]

file1 <- "downloaded_data/CBIO_GBM_log2.txt"
file2 <- "downloaded_data/CBIO_GBM_GISTIC.txt"
dat1 <- read.delim(file1, as.is=TRUE, skip=1)
dat2 <- read.delim(file2, as.is=TRUE, skip=1)
identical(dat1$COMMON[which(!(is.na(dat1$IDH1)))], dat2$COMMON[which(!(is.na(dat2$IDH1)))])
GBM_SNP <- dat1$COMMON[which(!(is.na(dat1$IDH1)))]



# Somatic Mutations downloaded from TCGA File Download 11/4/15
# Use UCSC mutation calls because it (a) has the most samples and (b) includes RNA allele counts
file1 <- "/costellolab/mazort/TCGA_data/Mutations_LGG_GBM/LGG/Somatic_Mutations/UCSC__Automated_Mutation_Calling/Level_2/ucsc.edu_LGG.IlluminaGA_DNASeq_automated.Level_2.1.4.0.somatic.maf"
dat1 <- read.delim(file1, as.is=TRUE, comment.char="#")
LGG_muts_pre <- unique(sapply(dat1$Tumor_Sample_Barcode, function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ))
LGG_muts <- unique(unname(sapply(LGG_muts_pre, function(x) { substr(x,1,nchar(x)-1) })))

file1 <- "/costellolab/mazort/TCGA_data/Mutations_LGG_GBM/GBM/Somatic_Mutations/UCSC__Automated_Mutation_Calling/Level_2/ucsc.edu_GBM.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.maf"
dat1 <- read.delim(file1, as.is=TRUE, comment.char="#")
GBM_muts_pre <- unique(sapply(dat1$Tumor_Sample_Barcode, function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ))
GBM_muts <- unique(unname(sapply(GBM_muts_pre, function(x) { substr(x,1,nchar(x)-1) })))


## How many samples do we have and how many overlap?
length(LGG_meth) ##  532
length(LGG_SNP)  ##  513
length(LGG_muts) ##  516
length(intersect(LGG_meth, intersect(LGG_SNP, LGG_muts)))  # 513

length(GBM_meth) ##  156
length(GBM_SNP)  ##  577
length(GBM_muts) ##  315
length(intersect(GBM_meth, intersect(GBM_SNP, GBM_muts)))  # 134


## Save TCGA IDs
LGG_IDs <- intersect(LGG_meth, intersect(LGG_SNP, LGG_muts))
GBM_IDs <- intersect(GBM_meth, intersect(GBM_SNP, GBM_muts))
save(LGG_IDs, GBM_IDs, file="data/TCGA_IDs_with_450k_SNParray_mutations.Rdata")
