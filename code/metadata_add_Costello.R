
## Read in our patient metadata
costello.meta <- read.delim("data/sample.info.IDH1.txt", header=TRUE, sep="\t", as.is=TRUE)
costello.meta <- costello.meta[order(costello.meta$source, costello.meta$patient_ID, costello.meta$sample_type), ]
costello.meta$barcodes <- paste(costello.meta$beadchip, costello.meta$strip, sep="_")
costello.meta$name_uniq <- paste(costello.meta$patient_ID, costello.meta$sample_type, sep="_")
print("*** confirm that name_uniq are appropriate column names for future data frame -- output should be TRUE:")
identical(costello.meta$name_uniq, make.names(costello.meta$name_uniq, unique=TRUE))  ## should return TRUE
if(!(identical(costello.meta$name_uniq, make.names(costello.meta$name_uniq, unique=TRUE)))) {
  print("*** sample(s) that don't have matching names are: ")
  costello.meta$name_uniq[which(!(costello.meta$name_uniq == make.names(costello.meta$name_uniq, unique=TRUE)))]
}

## Save
save(costello.meta, file="data/costello.meta.Rdata")

## Rebuild metadata table to match format of TCGA metadata
cost2.meta <- data.frame(sample_ID = costello.meta$name_uniq, stringsAsFactors=FALSE); rownames(cost2.meta) <- cost2.meta$sample_ID
cost2.meta$sample_ID2 <- cost2.meta$sample_ID
cost2.meta$source <- costello.meta$source
cost2.meta$barcode <- costello.meta$barcode
cost2.meta$patient_ID <- costello.meta$patient_ID
cost2.meta$grade <- costello.meta$grade
cost2.meta$grade_of_rec1 <- costello.meta$grade_of_rec1
cost2.meta$IDH_status <- costello.meta$IDH_status
cost2.meta$IDH_exome_VF <- NA; cost2.meta$IDH_RNA_VF <- NA
cost2.meta$IDH1_log2 <- NA; cost2.meta$IDH1_GISTIC <- NA; cost2.meta$IDH2_log2 <- NA; cost2.meta$IDH2_GISTIC <- NA

## Open TCGA metdata
load("data/TCGA_metadata.Rdata")

## Combine into one table
TCGA.meta$grade <- NA; TCGA.meta$grade_of_rec1 <- NA; TCGA.meta$patient_ID <- NA
meta.all <- rbind(cost2.meta, TCGA.meta)
save(meta.all, file="data/metadata.merged.Rdata")

## Make sure there are links to all of these idat files
idatdir <- "/costellolab/mazort/TCGA_data/450k_idat_links/"
filenames <- list.files(path=idatdir)

for(bb in cost2.meta$barcode) {
  grn <- paste0(bb,"_Grn.idat");
  if(!(grn %in% filenames)) {
    cmd <- paste0("ln -s /home/mazort/LG3/infinium/methylumi/rawdata/", grn, " ", idatdir)
    print(cmd)
    system(cmd)
  }

  red <- paste0(bb,"_Red.idat")
  if(!(red %in% filenames)) {
    cmd <- paste0("ln -s /home/mazort/LG3/infinium/methylumi/rawdata/", red, " ", idatdir)
    print(cmd)
    system(cmd)
  }
}



