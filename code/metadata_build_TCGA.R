
load("data/TCGA_IDs_with_450k_SNParray_mutations.Rdata")

options(stringsAsFactors = FALSE)
LGG.meta <- data.frame(sample_ID = LGG_IDs); rownames(LGG.meta) <- LGG.meta$sample_ID
GBM.meta <- data.frame(sample_ID = GBM_IDs); rownames(GBM.meta) <- GBM.meta$sample_ID

## Add 450k file names
# LGG
file1 <- "/costellolab/mazort/TCGA_data/LGG_450k/METADATA/JHU_USC__HumanMethylation450/jhu-usc.edu_LGG.HumanMethylation450.1.12.0.sdrf.txt"
dat1 <- read.delim(file1, as.is=TRUE)
dat1$sample_ID <- sapply(sapply(dat1$Comment..TCGA.Barcode., function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ), function(x) { substr(x,1,nchar(x)-1) })
sum(!(LGG.meta$sample_ID %in% dat1$sample_ID)) ## sanity check
tmp1 <- dat1[ which(dat1$Label == "Cy3" & dat1$sample_ID %in% LGG.meta$sample_ID), c("sample_ID","Array.Data.File")]
grn.idat <- tmp1$Array.Data.File[match(tmp1$sample_ID, LGG.meta$sample_ID)]
barcodes_grn <- unname(sapply(grn.idat, function(x) { paste(strsplit(x, "_")[[1]][1:2], collapse="_") } ))
tmp1 <- dat1[ which(dat1$Label == "Cy5" & dat1$sample_ID %in% LGG.meta$sample_ID), c("sample_ID","Array.Data.File")]
red.idat <- tmp1$Array.Data.File[match(tmp1$sample_ID, LGG.meta$sample_ID)]
barcodes_red <- unname(sapply(red.idat, function(x) { paste(strsplit(x, "_")[[1]][1:2], collapse="_") } ))
identical(barcodes_grn, barcodes_red)
LGG.meta$barcode <- barcodes_grn

# GBM
file1 <- "/costellolab/mazort/TCGA_data/GBM_450k/METADATA/JHU_USC__HumanMethylation450/jhu-usc.edu_GBM.HumanMethylation450.1.6.0.sdrf.txt"
dat1 <- read.delim(file1, as.is=TRUE)
dat1$sample_ID <- sapply(sapply(dat1$Comment..TCGA.Barcode., function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ), function(x) { substr(x,1,nchar(x)-1) })
sum(!(GBM.meta$sample_ID %in% dat1$sample_ID)) ## sanity check
tmp1 <- dat1[ which(dat1$Label == "Cy3" & dat1$sample_ID %in% GBM.meta$sample_ID), c("sample_ID","Array.Data.File")]
grn.idat <- tmp1$Array.Data.File[match(tmp1$sample_ID, GBM.meta$sample_ID)]
barcodes_grn <- unname(sapply(grn.idat, function(x) { paste(strsplit(x, "_")[[1]][1:2], collapse="_") } ))
tmp1 <- dat1[ which(dat1$Label == "Cy5" & dat1$sample_ID %in% GBM.meta$sample_ID), c("sample_ID","Array.Data.File")]
red.idat <- tmp1$Array.Data.File[match(tmp1$sample_ID, GBM.meta$sample_ID)]
barcodes_red <- unname(sapply(red.idat, function(x) { paste(strsplit(x, "_")[[1]][1:2], collapse="_") } ))
identical(barcodes_grn, barcodes_red)
GBM.meta$barcode <- barcodes_grn


## Add IDH1/2 mutation status
# LGG
file1 <- "/costellolab/mazort/TCGA_data/Mutations_LGG_GBM/LGG/Somatic_Mutations/UCSC__Automated_Mutation_Calling/Level_2/ucsc.edu_LGG.IlluminaGA_DNASeq_automated.Level_2.1.4.0.somatic.maf"
dat1 <- read.delim(file1, as.is=TRUE, comment.char="#")
dat1$sample_ID <- sapply(sapply(dat1$Tumor_Sample_Barcode, function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ), function(x) { substr(x,1,nchar(x)-1) })
sum(!(LGG.meta$sample_ID %in% dat1$sample_ID)) ## sanity check
tmp1 <- dat1[which(dat1$Hugo_Symbol %in% c("IDH1", "IDH2") & dat1$sample_ID %in% LGG.meta$sample_ID), ]
tmp1$IDH_status <- paste(tmp1$Hugo_Symbol, tmp1$Chrom, tmp1$Start_Position, tmp1$Variant_Classification, sep="_")
tmp1$IDH_exome_VF <- tmp1$Tumor_Alt_Count / (tmp1$Tumor_Ref_Count + tmp1$Tumor_Alt_Count)
tmp1$IDH_RNA_VF <- tmp1$RNA_Tumor_Alt_Count / (tmp1$RNA_Tumor_Ref_Count + tmp1$RNA_Tumor_Alt_Count)
LGG.meta$IDH_status <- "wt"; LGG.meta[tmp1$sample_ID, "IDH_status"] <- tmp1$IDH_status
LGG.meta$IDH_exome_VF <- NA; LGG.meta[tmp1$sample_ID, "IDH_exome_VF"] <- tmp1$IDH_exome_VF
LGG.meta$IDH_RNA_VF <- NA; LGG.meta[tmp1$sample_ID, "IDH_RNA_VF"] <- tmp1$IDH_RNA_VF

# the first time I ran this, I found several false negative IDH calls in the LGG dataset -- the above mutation list doesn't call an IDH1/2 mutation, but one of the other centers does
ss <- c("TCGA-DU-5851-01","TCGA-QH-A65X-01","TCGA-DB-A4XF-01","TCGA-DU-A7TI-01","TCGA-FG-A87N-01","TCGA-DU-6399-01","TCGA-DU-6407-01","TCGA-FG-6690-01","TCGA-HT-7472-01","TCGA-QH-A86X-01","TCGA-QH-A6X3-01")
file1 <- "/costellolab/mazort/TCGA_data/Mutations_LGG_GBM/LGG/Somatic_Mutations/MDA__Mutation_Calling/Level_2/mdanderson.org_LGG.IlluminaGA_DNASeq.Level_2.1.6.somatic.maf"
dat1 <- read.delim(file1, as.is=TRUE, comment.char="#")
dat1$sample_ID <- sapply(sapply(dat1$Tumor_Sample_Barcode, function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ), function(x) { substr(x,1,nchar(x)-1) })
for(mm in ss) {
	#print(mm)
	tmp1 <- dat1[which(dat1$Hugo_Symbol %in% c("IDH1", "IDH2") & dat1$sample_ID == mm), ]
	tmp1$IDH_status <- paste(tmp1$Hugo_Symbol, tmp1$Chrom, tmp1$Start_Position, tmp1$Variant_Classification, sep="_")
	#print(tmp1$IDH_status)
	LGG.meta[mm,"IDH_status"] <- tmp1$IDH_status
	LGG.meta[mm,"IDH_exome_VF"] <- tmp1$t_alt_count/tmp1$t_depth
}

ss <- c("TCGA-HT-8106-01","TCGA-HT-A618-01","TCGA-DU-A7TG-01","TCGA-FG-8189-01")
file1 <- "/costellolab/mazort/TCGA_data/Mutations_LGG_GBM/LGG/Somatic_Mutations/BI__Automated_Mutation_Calling/Level_2/PR_TCGA_LGG_PAIR_Capture_All_Pairs_QCPASS_v7.aggregated.capture.tcga.uuid.automated.somatic.maf"
dat1 <- read.delim(file1, as.is=TRUE, comment.char="#")
dat1$sample_ID <- sapply(sapply(dat1$Tumor_Sample_Barcode, function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ), function(x) { substr(x,1,nchar(x)-1) })
for(mm in ss) {
	#print(mm)
	tmp1 <- dat1[which(dat1$Hugo_Symbol %in% c("IDH1", "IDH2") & dat1$sample_ID == mm), ]
	tmp1$IDH_status <- paste(tmp1$Hugo_Symbol, tmp1$Chromosome, tmp1$Start_position, tmp1$Variant_Classification, sep="_")
	#print(tmp1$IDH_status)
	LGG.meta[mm,"IDH_status"] <- tmp1$IDH_status
	LGG.meta[mm,"IDH_exome_VF"] <- tmp1$t_alt_count/(tmp1$t_alt_count + tmp1$t_ref_count)
}

# GBM
file1 <- "/costellolab/mazort/TCGA_data/Mutations_LGG_GBM/GBM/Somatic_Mutations/UCSC__Automated_Mutation_Calling/Level_2/ucsc.edu_GBM.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.maf"
dat1 <- read.delim(file1, as.is=TRUE, comment.char="#")
dat1$sample_ID <- sapply(sapply(dat1$Tumor_Sample_Barcode, function(x) { paste0(strsplit(x,"-")[[1]][1:4],collapse="-") } ), function(x) { substr(x,1,nchar(x)-1) })
sum(!(GBM.meta$sample_ID %in% dat1$sample_ID)) ## sanity check
tmp1 <- dat1[which(dat1$Hugo_Symbol %in% c("IDH1", "IDH2") & dat1$sample_ID %in% GBM.meta$sample_ID), ]
tmp1$IDH_status <- paste(tmp1$Hugo_Symbol, tmp1$Chrom, tmp1$Start_Position, tmp1$Variant_Classification, sep="_")
tmp1$IDH_exome_VF <- tmp1$Tumor_Alt_Count / (tmp1$Tumor_Ref_Count + tmp1$Tumor_Alt_Count)
tmp1$IDH_RNA_VF <- tmp1$RNA_Tumor_Alt_Count / (tmp1$RNA_Tumor_Ref_Count + tmp1$RNA_Tumor_Alt_Count)
GBM.meta$IDH_status <- "wt"; GBM.meta[tmp1$sample_ID, "IDH_status"] <- tmp1$IDH_status
GBM.meta$IDH_exome_VF <- NA; GBM.meta[tmp1$sample_ID, "IDH_exome_VF"] <- tmp1$IDH_exome_VF
GBM.meta$IDH_RNA_VF <- NA; GBM.meta[tmp1$sample_ID, "IDH_RNA_VF"] <- tmp1$IDH_RNA_VF


## Add copy number at IDH1/2
# LGG
file1 <- "downloaded_data/CBIO_LGG_log2.txt"
file2 <- "downloaded_data/CBIO_LGG_GISTIC.txt"
dat1 <- read.delim(file1, as.is=TRUE, skip=1); dat2 <- read.delim(file2, as.is=TRUE, skip=1)
tmp1 <- dat1[which(dat1$COMMON %in% LGG.meta$sample_ID), ]; tmp2 <- dat2[which(dat2$COMMON %in% LGG.meta$sample_ID), ]
LGG.meta$IDH1_log2 <- NA; LGG.meta[tmp1$COMMON, "IDH1_log2"] <- tmp1$IDH1
LGG.meta$IDH1_GISTIC <- NA; LGG.meta[tmp2$COMMON, "IDH1_GISTIC"] <- tmp2$IDH1
LGG.meta$IDH2_log2 <- NA; LGG.meta[tmp1$COMMON, "IDH2_log2"] <- tmp1$IDH2
LGG.meta$IDH2_GISTIC <- NA; LGG.meta[tmp2$COMMON, "IDH2_GISTIC"] <- tmp2$IDH2

# GBM
file1 <- "downloaded_data/CBIO_GBM_log2.txt"
file2 <- "downloaded_data/CBIO_GBM_GISTIC.txt"
dat1 <- read.delim(file1, as.is=TRUE, skip=1); dat2 <- read.delim(file2, as.is=TRUE, skip=1)
tmp1 <- dat1[which(dat1$COMMON %in% GBM.meta$sample_ID), ]; tmp2 <- dat2[which(dat2$COMMON %in% GBM.meta$sample_ID), ]
GBM.meta$IDH1_log2 <- NA; GBM.meta[tmp1$COMMON, "IDH1_log2"] <- tmp1$IDH1
GBM.meta$IDH1_GISTIC <- NA; GBM.meta[tmp2$COMMON, "IDH1_GISTIC"] <- tmp2$IDH1
GBM.meta$IDH2_log2 <- NA; GBM.meta[tmp1$COMMON, "IDH2_log2"] <- tmp1$IDH2
GBM.meta$IDH2_GISTIC <- NA; GBM.meta[tmp2$COMMON, "IDH2_GISTIC"] <- tmp2$IDH2


## Combine into one table
LGG.meta$source <- "LGG"
GBM.meta$source <- "GBM"
TCGA.meta <- rbind(LGG.meta, GBM.meta)
TCGA.meta$sample_ID2 <- gsub("-",".",TCGA.meta$sample_ID)
save(TCGA.meta, file="data/TCGA_metadata.Rdata")
