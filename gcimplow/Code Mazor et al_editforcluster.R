##Script
#libraries required to run script
library(methylumi)
#library(XLConnect)  ### failed installation
library(caret)
library(randomForest)
library(doMC)
library(plyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(RCurl)
##Mazor et al - IDH mut + IDH CNA at recurrence - GCIMP classification

##--------------------------------------------------
##downloaded files from Tali's box.com account.
##unzipped file and contained idats in a folder called 'data' and a key.txt file.

##Normalize Data according to TCGA metric.
Mazor_barcode <- read.delim("IDH_450k_key.txt", header = FALSE)  ##28:3
Mazor_barcode$Files <- substr(as.character(Mazor_barcode$V2),1,nchar(as.character(Mazor_barcode$V2))-9)  ##28:4
#idat <- methylumIDAT(Mazor_barcode$Files, idatPath = "data/")
idat <- methylumIDAT(Mazor_barcode$Files, idatPath= "~/LG3/infinium/methylumi/rawdata/")
proc <- stripOOB(normalizeMethyLumiSet(methylumi.bgcorr(idat)))
betaValue <- betas(proc)  ##485577:28

##--------------------------------------------------
##Random Forest IDHmutant+IDHmutantCNA recurrent project, Mazor et al. (IDH mutant/G-CIMP)
##Load TCGA LGG-GBM DNA methylation matrix
#URL <- "https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2016/LGG.GBM.meth.txt"
#destfile <- "LGG.GBM.meth.txt"
#download.file(URL, destfile, mode = "wb")
LGG.GBM <- read.table(file = "LGG.GBM.meth.txt", sep = "\t", check.names=F, header = T)
rownames(LGG.GBM) <- LGG.GBM[,1]
LGG.GBM <- LGG.GBM[,-1]

##Load TCGA LGG-GBM metadata

#URL <- "http://www.cell.com/cms/attachment/2045372863/2056783242/mmc2.xlsx"
#destfile <- "mmc2.xlsx"
#download.file(URL, destfile, mode = "wb")
#d <- loadWorkbook("mmc2.xlsx") #get from Ceccarelli et al. 2016 supplemental table 1 - http://www.cell.com/cms/attachment/2045372863/2056783242/mmc2.xlsx
#pd <- readWorksheet(d, sheet = 1, startRow=2, header = TRUE)
### having trouble installing package to do this so download manually, save first sheet as tab-delimited text (after deleting first row) and read that in
pd <- read.table(file = "mmc2.sheet1.txt", sep="\t", header=TRUE)

## Load 163 CpG probes that define G-CIMP (probes from Cell paper) 
#URL <- "https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2015/PanGlioma_MethylationSignatures.xlsx"
#destfile <- "PanGlioma_MethylationSignatures.xlsx"
#download.file(URL, destfile, mode = "wb")
#d <- loadWorkbook("PanGlioma_MethylationSignatures.xlsx") #download from https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2015/PanGlioma_MethylationSignatures.xlsx
#IDHmut_supervised_163p <- readWorksheet(d, sheet = 5, header = TRUE)
### having trouble installing package to do this so download manually, save fifth sheet as tab-delimited text and read that in
IDHmut_supervised_163p <- read.table(file = "PanGlioma_MethylationSignatures_sheet5.txt", sep="\t", header=TRUE)

##Find common probes between LGG.GBM and IDHmut_meth CpG methylation data
betaValue_163 <- betaValue[as.character(IDHmut_supervised_163p[,1]),]  ##163:28
dim(na.omit(betaValue_163))  ##163:28
betaValue_163 <- na.omit(betaValue_163) ##163:28
aux <- subset(pd, Supervised.DNA.Methylation.Cluster %in% c("G-CIMP-low","G-CIMP-high","Codel"))  ##448:51
LGG_GBM_163 <- LGG.GBM[rownames(betaValue_163), substr(colnames(LGG.GBM),1,12) %in% as.character(aux$Case)]  ##163:448
LGG_GBM_163_s <- na.omit(LGG_GBM_163)  ##159:448
colnames(LGG_GBM_163_s) <- substr(colnames(LGG_GBM_163_s),1,12)

##Train data
TrainingData <- t(LGG_GBM_163_s)  ##448:159

##Add an additional column of data that you wish to train on
##For classifying the new data add a new column called "clustM.supervised2"
TrainingData <- merge(TrainingData, aux[,c("Case","Supervised.DNA.Methylation.Cluster")], by.x = 0, by.y = "Case")  ##448:161
rownames(TrainingData) <- as.character(TrainingData$Row.names)
TrainingData <- TrainingData[,-1]  ##448 samples, 159 probes, 1 Supervised.DNA.Methylation.Cluster
TrainingData$Supervised.DNA.Methylation.Cluster <- factor(TrainingData$Supervised.DNA.Methylation.Cluster)

##Save data
save(TrainingData, file = "RF_Training_IDHmut_163.Rda")
save(LGG_GBM_163_s, betaValue_163, aux, file = "RF_obj_IDHmut_163.Rda")

##Close the R session and start from here
load("RF_Training_IDHmut_163.Rda")  ##TrainingData  ##448:160

##Register cores for doMC
registerDoMC(cores = 12)

##Set up k-fold cross validation
fitControl <- trainControl(  ##10-fold CV
  method = "repeatedcv",
  number = 10,
  repeats = 10)

##Set your seed so your work is repeatable
set.seed(42)

##Create a subset of your data to train your model on. This makes sure you have equal representation of the 'clustM.supervised2' groups in your training set
inTraining <- createDataPartition(TrainingData$Supervised.DNA.Methylation.Cluster, p = 0.8, list = FALSE, times = 1)

##Training set
myTrain <- TrainingData[inTraining, ]

##Testing set
myTest <- TrainingData[-inTraining, ]

##Confirm seed was set
set.seed(210)

##Set values for mtry
##mtry is the "number of variables randomly sampled as candidates at each split"
##Traditionally for classification we use the sqrt of the number of variables, but here we are trying a range of mtry values to find the best parameters for our model
mtryVals <- floor(c(seq(100, 2000, by = 100),
                    sqrt(ncol(TrainingData))))
mtryGrid <- data.frame(.mtry = mtryVals)

##Confirm seed again
set.seed(420)

##Set the number of cores
registerDoMC(cores = 12)

##Run training (takes up to 2 hours on a 12 core machine and 32 gig RAM)
RF_IDHmut_163p <- train(Supervised.DNA.Methylation.Cluster ~ .,  ##variable to be trained on
                        data = TrainingData,  ##data we are working on
                        method = "rf",  ##method we are using
                        trControl = fitControl,  ##how we validate (we created this object)
                        ntree = 5000,  ##number of trees which is dependent on training data size
                        importance = TRUE,  ##calculate variable importance (this step can be omitted to speed up calculation)
                        tuneGrid = mtryGrid,  ##set mtrys
                        subset = inTraining  ##define training set
)

save.image("after_train_case_163.Rda")

##--------------------------------------------------
##Random Forest prediction model (applying my model to the new cases from Costello lab)
load("after_train_case_163.Rda")

### copy code from above because it appears to be required
Mazor_barcode <- read.delim("IDH_450k_key.txt", header = FALSE)  ##28:3
Mazor_barcode$Files <- substr(as.character(Mazor_barcode$V2),1,nchar(as.character(Mazor_barcode$V2))-9)  ##28:4
idat <- methylumIDAT(Mazor_barcode$Files, idatPath= "~/LG3/infinium/methylumi/rawdata/")
proc <- stripOOB(normalizeMethyLumiSet(methylumi.bgcorr(idat)))
betaValue <- betas(proc)  ##485577:28
## Load 163 CpG probes that define G-CIMP (probes from Cell paper) 
IDHmut_supervised_163p <- read.table(file = "PanGlioma_MethylationSignatures_sheet5.txt", sep="\t", header=TRUE)
##Find common probes between LGG.GBM and IDHmut_meth CpG methylation data
betaValue_163 <- betaValue[as.character(IDHmut_supervised_163p[,1]),]  ##163:28
dim(na.omit(betaValue_163))  ##163:28
betaValue_163 <- na.omit(betaValue_163) ##163:28



Mazor_163 <- betaValue[as.character(rownames(betaValue_163)),]  ##163:28
Mazor_163_t <- t(Mazor_163)  ##28:162
Mazor_163p_pred <- predict(RF_IDHmut_163p, Mazor_163_t, type = "prob")
tmp <- predict(RF_IDHmut_163p, Mazor_163_t)
Mazor_163p_pred$GliomaSubtype <- as.character(tmp)

##--------------------------------------------------
##Add IDH mut scores to Costello's data original keys
Mazor_163p_pred$Files <- as.character(rownames(Mazor_163p_pred))
IDH_450k_key <- merge(Mazor_barcode, Mazor_163p_pred, by = "Files", sort = FALSE)  ##28:7  ##combine metadata[annotation] with predicted RF scores for this new set of cases
IDH_450k_key_scores <- IDH_450k_key[,c(2:8)]
colnames(IDH_450k_key_scores) <- c("TumorID","Red.idat","Grn.idat","IDHmut-codel","G-CIMP-high","G-CIMP-low","GliomaSubtype")
save(IDH_450k_key_scores, file = "IDH_450k_key_scores.Rda")
write.table(IDH_450k_key_scores, file = "IDH_450k_key_scores.txt", quote = FALSE, row.names = FALSE, sep = "\t")

##--------------------------------------------------
load(file = "IDH_450k_key_scores.Rda")
##Plots
d <- IDH_450k_key_scores  ##28:9
tmps <- colsplit(as.character(d$TumorID),"_",c("PatientID","SampleType","Fragment"))
#tmps <- colsplit(as.character(tmp$SampleType),"_",c("SampleType2","Fragment"));tmp <- tmp[,-2]
#d <- cbind(d, tmp, tmps)
d <- cbind(d, tmps)
#names(d)[8] <- "PatientID"
d$SampleType2 <- as.factor(d$SampleType)
levels(d$SampleType2) <- c("Primary","First Recurrence","Second Recurrence")
d$GliomaSubtype <- as.factor(d$GliomaSubtype)
d$GliomaSubtype_2 <- factor(d$GliomaSubtype, levels = c("G-CIMP-low","G-CIMP-high"))
levels(d$GliomaSubtype_2)
d$PatientID <- factor(d$PatientID)
levels(d$PatientID)

##Scatter plot G-CIMP-low score x G-CIMP-high score (colored by G-CIMP classification)
pdf(file = "Mazor_scatter_plot_G-CIMP.pdf", width = 12, height = 8)

q <- ggplot(d, aes(x = `G-CIMP-high`,y= `G-CIMP-low`, color = GliomaSubtype_2, shape = SampleType2)) + geom_point(size = 3.5) +
  scale_color_manual(values = c("#006400","#FF0000"),
  name = "G-CIMP Classification") +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + theme(axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold")) + theme(strip.text.y = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12, face = "bold")) + xlab("G-CIMP-high score") + ylab("G-CIMP-low score")

q + geom_vline(xintercept = 0.5, linetype = 2) + geom_hline(yintercept = 0.5, linetype = 2)

dev.off()

## same plot with my color scheme
library(RColorBrewer)
mut_cols <- brewer.pal(8,"Set2")
d$IDH_status <- "pre_change"
d$IDH_status[which(d$PatientID=="Patient68" & d$SampleType=="Recurrence1")] <- "gain_wt"
d$IDH_status[which(d$PatientID=="Patient17" & d$SampleType=="Recurrence1")] <- "gain_mut"
d$IDH_status[which(d$PatientID=="Patient169" & d$SampleType=="Recurrence1")] <- "del_mut"
d$IDH_status[which(d$PatientID=="Patient21" & d$SampleType=="Recurrence1")] <- "amp_mut"
d$IDH_status[which(d$PatientID=="Patient14" & d$SampleType=="Recurrence2")] <- "del_mut"
d$IDH_status <- factor(d$IDH_status, levels=c("mut", "wt", "amp_mut", "del_mut", "gain_mut", "gain_wt", "pre_change","del_wt"))

pdf("scatterplot_GCIMPscores.pdf")
plot(d[,"G-CIMP-high"], d[,"G-CIMP-low"], xaxt="n", yaxt="n", pch=c(19,17)[d$GliomaSubtype_2], xlim=c(0,1), ylim=c(0,1), xlab="GCIMP-high score", ylab="GCIMP-low score", col=mut_cols[d$IDH_status])
axis(1, at=seq(0,1,.25), las=1)
axis(2, at=seq(0,1,.25), las=1)
abline(h=0.5,v=0.5,lty=2)
legend("topright", legend=levels(d$IDH_status), col=mut_cols, pch=15)
legend("bottomleft", legend=levels(d$GliomaSubtype_2), col="black", pch=c(19,17))
dev.off()


##--------------------------------------------------
##Scatter plot G-CIMP-low score x G-CIMP-high score (colored by PatientID)
pdf(file = "Mazor_scatter_plot_G-CIMP_PatientID.pdf", width = 12, height = 8)
q <- ggplot(d, aes(x = `G-CIMP-high`,y= `G-CIMP-low`, color = PatientID, shape = GliomaSubtype_2)) + geom_point(size = 3.5) +
     scale_color_manual(values = c("#DC143C","#1C86EE","#00C957","#FFD700","#8B3626"),
     name = "Patient ID")
q + facet_grid(.~ SampleType2) + theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + theme(axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold")) + theme(strip.text.y = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12, face = "bold")) + xlab("G-CIMP-high score") + ylab("G-CIMP-low score")

dev.off()

##--------------------------------------------------
##Barplot by PatientID
dm <- melt(d[,c("PatientID", "GliomaSubtype_2")], id.vars = "PatientID")  ##28:3
patient <- c("p1","p2","p3","p1","p2","p3")
group <- c("1","2","3","1","3","2")
x <- data.frame(patient, group)
x
aux <- plyr::ddply(dm, .(PatientID), function(x) return(length(unique(x$value)) > 1))  ##TRUE means the cases that change the RF classification between P/R
aux
d.m <- merge(d, aux, by = "PatientID", all.x = TRUE, sort = FALSE)  ##28:12
d.m$V1 <- as.factor(d.m$V1)
levels(as.factor(d.m$V1))
levels(d.m$V1) <- c("No Change", "Change")
levels(d.m$V1)
pdf(file = "Mazor_barplot_G-CIMP_PatientID.pdf", width = 18)
ggplot(d.m, aes(PatientID, fill = GliomaSubtype_2)) + geom_bar() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, face = "bold"), axis.text.y = element_text(size = 12, hjust = 1, vjust = 1), strip.text.y = element_text(size = 10, face = "bold"), strip.text.x = element_text(size = 10, face = "bold")) + scale_fill_manual(values = c("#006400","#FF0000"), name = "G-CIMP Classification") + scale_y_continuous(breaks = seq(0, 10, 2)) + xlab("PatientID") + ylab("Count") + facet_grid(SampleType2~V1, scales = "free", space = "free")
dev.off()
##--------------------------------------------------