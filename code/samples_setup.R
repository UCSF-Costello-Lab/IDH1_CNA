### Prep sample information -- column numbers for selectively pulling samples
pat.IDH <- paste0("Patient", c("14","169","68","17", "21"))
idx.IDH <- which(meta.all$patient_ID %in% pat.IDH)
idx.costello <- which(meta.all$source == "Costello")
idx.notIDH <- setdiff(idx.costello, idx.IDH)

idx.luchman <- which(meta.all$source == "Luchman")
pat.luchman <- unique(meta.all$patient_ID[idx.luchman])

idx.LGG <- which(meta.all$source == "LGG")
idx.GBM <- which(meta.all$source == "GBM")
idx.TCGA.wt <- intersect(which(meta.all$IDH_status == "wt"), union(idx.LGG, idx.GBM))
idx.TCGA.mut <- setdiff(union(idx.LGG, idx.GBM), idx.TCGA.wt)

idx.adult <- which(meta.all$grade == "Adult")
idx.ped <- which(meta.all$grade == "Pediatric")
idx.fetal <- which(meta.all$grade == "Fetal")
idx.normals <- c(idx.adult, idx.ped, idx.fetal)

## IDH1 samples are split between UCSF core & USC core -- for better comparison, use pairs that were run at the same core 
## Since later analyses will use P14, 169 & 21, use USC samples for all of those
IDH.pre.all <- c(which(meta.all$sample_ID == "Patient17_Primary"),
				 which(meta.all$sample_ID == "Patient68_Primary"),
				 which(meta.all$sample_ID == "Patient169_Primary_REPEAT"),
				 which(meta.all$sample_ID == "Patient14_Recurrence1_REPEAT"),
				 which(meta.all$sample_ID == "Patient21_Primary"))
IDH.post.all <- c(which(meta.all$sample_ID == "Patient17_Recurrence1"),
				  which(meta.all$sample_ID == "Patient68_Recurrence1C"),
				  which(meta.all$sample_ID == "Patient169_Recurrence1_REPEAT"),
				  which(meta.all$sample_ID == "Patient14_Recurrence2v1_REPEAT"),
				  which(meta.all$sample_ID == "Patient21_Recurrence1v2"))
IDH.pre.sub <- c(which(meta.all$sample_ID == "Patient169_Primary_REPEAT"),
				 which(meta.all$sample_ID == "Patient14_Recurrence1_REPEAT"),
				 which(meta.all$sample_ID == "Patient21_Primary"))
IDH.post.sub <- c(which(meta.all$sample_ID == "Patient169_Recurrence1_REPEAT"),
				  which(meta.all$sample_ID == "Patient14_Recurrence2v1_REPEAT"),
				  which(meta.all$sample_ID == "Patient21_Recurrence1v2"))
				  
IDH.luchman.pre <- c(which(meta.all$sample_ID == "BT142_cell_het"), which(meta.all$sample_ID == "BT142_tumor"), which(meta.all$sample_ID == "BT142_tumor"),
				 which(meta.all$sample_ID == "BT61_cell_het"), which(meta.all$sample_ID == "BT61_tumor"), which(meta.all$sample_ID == "BT92_tumor"),
				 which(meta.all$sample_ID == "BT88_tumor"), which(meta.all$sample_ID == "BT88_tumor"),
				 which(meta.all$sample_ID == "BT54_tumor"),
				 which(meta.all$sample_ID == "BT257_tumor"), which(meta.all$sample_ID == "BT257_tumor"),	
				 which(meta.all$sample_ID == "BT257_tumor"),	 which(meta.all$sample_ID == "BT257_tumor"),
				 which(meta.all$sample_ID == "BT257_X1"), which(meta.all$sample_ID == "BT257_X1"), which(meta.all$sample_ID == "BT257_X1") )
IDH.luchman.post <- c(which(meta.all$sample_ID == "BT142_cell_hemi"), which(meta.all$sample_ID == "BT142_cell_hemi"), which(meta.all$sample_ID == "BT142_xeno"),
				  which(meta.all$sample_ID == "BT92_cell_hemi"), which(meta.all$sample_ID == "BT92_tumor"),	which(meta.all$sample_ID == "BT92_cell_hemi"),
				  which(meta.all$sample_ID == "BT88_cell_hemi"), which(meta.all$sample_ID == "BT88_xeno"),
				  which(meta.all$sample_ID == "BT54_cell_hemi"),
				  which(meta.all$sample_ID == "BT257_X1"), which(meta.all$sample_ID == "BT257_X4"),
				  which(meta.all$sample_ID == "BT257_X7"), which(meta.all$sample_ID == "BT257_X11"),
				  which(meta.all$sample_ID == "BT257_X4"), which(meta.all$sample_ID == "BT257_X7"), which(meta.all$sample_ID == "BT257_X11") )				  
IDH.luchman.names <- paste0(meta.all$sample_ID[IDH.luchman.pre], "_vs_", meta.all$sample_ID[IDH.luchman.post])

IDH.luchman.pre.sub <- c(which(meta.all$sample_ID == "BT142_cell_het"), which(meta.all$sample_ID == "BT142_tumor"),
				 which(meta.all$sample_ID == "BT61_cell_het"), which(meta.all$sample_ID == "BT92_tumor"),
				 which(meta.all$sample_ID == "BT88_tumor"), which(meta.all$sample_ID == "BT88_tumor"),
				 which(meta.all$sample_ID == "BT54_tumor"),
				 which(meta.all$sample_ID == "BT257_tumor") )
IDH.luchman.post.sub <- c(which(meta.all$sample_ID == "BT142_cell_hemi"), which(meta.all$sample_ID == "BT142_xeno"),
				  which(meta.all$sample_ID == "BT92_cell_hemi"), which(meta.all$sample_ID == "BT92_cell_hemi"),
				  which(meta.all$sample_ID == "BT88_cell_hemi"), which(meta.all$sample_ID == "BT88_xeno"),
				  which(meta.all$sample_ID == "BT54_cell_hemi"),
				  which(meta.all$sample_ID == "BT257_X11") )				  
IDH.luchman.names.sub <- paste0(meta.all$sample_ID[IDH.luchman.pre.sub], "_vs_", meta.all$sample_ID[IDH.luchman.post.sub])

IDH.luchman.tumVcell.pre <- c(which(meta.all$sample_ID == "BT142_tumor"),
                              which(meta.all$sample_ID == "BT54_tumor"),
			      which(meta.all$sample_ID == "BT88_tumor"),
			      which(meta.all$sample_ID == "BT92_tumor"))
IDH.luchman.tumVcell.post <- c(which(meta.all$sample_ID == "BT142_cell_hemi"),
                               which(meta.all$sample_ID == "BT54_cell_hemi"),
			       which(meta.all$sample_ID == "BT88_cell_hemi"),
			       which(meta.all$sample_ID == "BT92_cell_hemi"))
IDH.luchman.tumVcell.names <- paste0(meta.all$sample_ID[IDH.luchman.tumVcell.pre], "_vs_", meta.all$sample_ID[IDH.luchman.tumVcell.post])

IDH.luchman.tumVxeno.pre <- c(which(meta.all$sample_ID == "BT142_tumor"), which(meta.all$sample_ID == "BT88_tumor"), which(meta.all$sample_ID == "BT257_tumor") )
IDH.luchman.tumVxeno.post <- c(which(meta.all$sample_ID == "BT142_xeno"), which(meta.all$sample_ID == "BT88_xeno"), which(meta.all$sample_ID == "BT257_X11") )
IDH.luchman.tumVxeno.names <- paste0(meta.all$sample_ID[IDH.luchman.tumVxeno.pre],  "_vs_", meta.all$sample_ID[IDH.luchman.tumVxeno.post])

IDH.luchman.oneper.pre <- c(which(meta.all$sample_ID == "BT142_cell_het"),
                           which(meta.all$sample_ID == "BT54_tumor"),
			   which(meta.all$sample_ID == "BT88_tumor"),
			   which(meta.all$sample_ID == "BT92_tumor"),
			   which(meta.all$sample_ID == "BT257_tumor") )
IDH.luchman.oneper.post <- c(which(meta.all$sample_ID == "BT142_cell_hemi"),
                            which(meta.all$sample_ID == "BT54_cell_hemi"),
			    which(meta.all$sample_ID == "BT88_cell_hemi"),
			    which(meta.all$sample_ID == "BT92_cell_hemi"),
			    which(meta.all$sample_ID == "BT257_X11") )
IDH.luchman.oneper.names <- paste0(meta.all$sample_ID[IDH.luchman.oneper.pre], "_vs_", meta.all$sample_ID[IDH.luchman.oneper.post])


IDH.luchman.oneper61.pre <- c(which(meta.all$sample_ID == "BT142_cell_het"),
                           which(meta.all$sample_ID == "BT54_tumor"),
			   which(meta.all$sample_ID == "BT88_tumor"),
			   which(meta.all$sample_ID == "BT61_cell_het"),
			   which(meta.all$sample_ID == "BT257_tumor") )
IDH.luchman.oneper.post <- c(which(meta.all$sample_ID == "BT142_cell_hemi"),
                            which(meta.all$sample_ID == "BT54_cell_hemi"),
			    which(meta.all$sample_ID == "BT88_cell_hemi"),
			    which(meta.all$sample_ID == "BT92_cell_hemi"),
			    which(meta.all$sample_ID == "BT257_X11") )
IDH.luchman.oneper61.names <- paste0(meta.all$sample_ID[IDH.luchman.oneper61.pre], "_vs_", meta.all$sample_ID[IDH.luchman.oneper.post])


### All samples: pri & rec
pri.all <- intersect(idx.costello, grep("Primary", meta.all$sample_ID))
pri.v1.all <- pri.all[match(unique(meta.all$patient_ID[pri.all]), meta.all$patient_ID[pri.all])]
rec1.all <- intersect(idx.costello, grep("Recurrence1", meta.all$sample_ID))
rec1.v1.all <- rec1.all[match(unique(meta.all$patient_ID[rec1.all]), meta.all$patient_ID[rec1.all])]
pat.pairs <- intersect(meta.all$patient_ID[pri.v1.all], meta.all$patient_ID[rec1.v1.all])
pri.v1.pairs <- intersect(pri.v1.all, which(meta.all$patient_ID %in% pat.pairs))
rec1.v1.pairs <- intersect(rec1.v1.all, which(meta.all$patient_ID %in% pat.pairs))
pri.v1.pairs.notIDH <- intersect(pri.v1.pairs, idx.notIDH)
rec1.v1.pairs.notIDH <- intersect(rec1.v1.pairs, idx.notIDH)

pat.r1g2 <- sort(unique(meta.all$patient_ID[meta.all$source=="Costello" & meta.all$grade_of_rec1=="2"]))
pat.r1g3 <- sort(unique(meta.all$patient_ID[meta.all$source =="Costello" & meta.all$grade_of_rec1=="3"]))
pat.r1g4 <- sort(unique(meta.all$patient_ID[meta.all$source =="Costello" & meta.all$grade_of_rec1=="4"]))

pri.v1.r1g2 <- intersect(pri.v1.all, which(meta.all$patient_ID %in% pat.r1g2))
rec1.v1.r1g2 <- intersect(rec1.v1.all, which(meta.all$patient_ID %in% pat.r1g2))
pri.v1.r1g3 <- intersect(pri.v1.all, which(meta.all$patient_ID %in% pat.r1g3))
rec1.v1.r1g3 <- intersect(rec1.v1.all, which(meta.all$patient_ID %in% pat.r1g3))
pri.v1.r1g4 <- intersect(pri.v1.all, which(meta.all$patient_ID %in% pat.r1g4))
rec1.v1.r1g4 <- intersect(rec1.v1.all, which(meta.all$patient_ID %in% pat.r1g4))

pri.v1.pairs.r1g2 <- intersect(pri.v1.pairs, which(meta.all$patient_ID %in% pat.r1g2))
pri.v1.pairs.r1g3 <- intersect(pri.v1.pairs, which(meta.all$patient_ID %in% pat.r1g3))
pri.v1.pairs.r1g4 <- intersect(pri.v1.pairs, which(meta.all$patient_ID %in% pat.r1g4))
rec1.v1.pairs.r1g2 <- intersect(rec1.v1.pairs, which(meta.all$patient_ID %in% pat.r1g2))
rec1.v1.pairs.r1g3 <- intersect(rec1.v1.pairs, which(meta.all$patient_ID %in% pat.r1g3))
rec1.v1.pairs.r1g4 <- intersect(rec1.v1.pairs, which(meta.all$patient_ID %in% pat.r1g4))

pri.v1.pairs.r1g2.notIDH <- intersect(pri.v1.pairs.r1g2, idx.notIDH)
pri.v1.pairs.r1g3.notIDH <- intersect(pri.v1.pairs.r1g3, idx.notIDH)
pri.v1.pairs.r1g4.notIDH <- intersect(pri.v1.pairs.r1g4, idx.notIDH)
rec1.v1.pairs.r1g2.notIDH <- intersect(rec1.v1.pairs.r1g2, idx.notIDH)
rec1.v1.pairs.r1g3.notIDH <- intersect(rec1.v1.pairs.r1g3, idx.notIDH)
rec1.v1.pairs.r1g4.notIDH <- intersect(rec1.v1.pairs.r1g4, idx.notIDH)

### Define color scheme based on IDH status
library(RColorBrewer)
mut_cols <- brewer.pal(8,"Set2")
IDH_status_tmp <- meta.all$IDH_status
IDH_status_tmp[which(IDH_status_tmp == "")] <- "wt"
IDH_status_tmp[grep("IDH", IDH_status_tmp)] <- "mut"
IDH_status_factor <- factor(IDH_status_tmp, levels=c("mut", "wt", "amp_mut", "del_mut", "gain_mut", "gain_wt", "pre_change","del_wt"))
color_IDH <- mut_cols[IDH_status_factor]
group_cols <- brewer.pal(5, "Set1")
group_factor <- factor(meta.all$source, levels=c("Costello", "LGG", "Luchman", "GBM", "Normal"))
color_group <- group_cols[group_factor]
