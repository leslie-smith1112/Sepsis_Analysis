# Load libraries
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
library(plyr)
require(purrr)  # for map(), reduce()
library(here)


all <- readr::read_tsv(here("processed-data","Oct18_master_expression_nonnormalized.tsv"))
dim(all)
batch.ids.temp <- readr::read_tsv(here("processed-data","master_batch_ids.tsv"))
batch.ids <- batch.ids.temp$batch.ids
names(batch.ids) <- batch.ids.temp$Sample
head(batch.ids)
length(batch.ids)
metadata <- readr::read_tsv(here("processed-data", "master_sepsis_metadata_all_samples_add_run_labels.tsv"))

# make sure batch.ids, metadata, and expression matrix all have consistent samples 
batch.ids <- batch.ids[names(batch.ids) %in% colnames(all)]
length(batch.ids)
all <- column_to_rownames(all, "Gene")
all[1:5,1:5]
min(all)
# - keep metadata only from samples we have expression data for - #
dim(metadata)
meta.data <- metadata[metadata$geo_accession %in% colnames(all),]
dim(meta.data)#3335 samples 
dim(all)
all <- all[,colnames(all) %in% meta.data$geo_accession]
dim(all)
all[1:5,1:5]

#make sure all samples are common between metadata and reads
all(meta.data$geo_accession %in% colnames(all)) 
all(colnames(all) %in% meta.data$geo_accession) 

##########################################################INITIAL PCA ALL SAMPLES############################################################
all_transposed <- t(all) #trasnpose matrix for PCA function
all_transposed[1:5,1:5]

# - get rid of 0 values - #
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
library(ggplot2)
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame('Series' = meta.data$series, "dat" = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Series)) + geom_point() + ggtitle("Sepsis All Samples PCA No Batch Correction") + labs(y= "PC2", x = "PC1")
p
dtp <- data.frame('Series' = meta.data$disease_simplified, "dat" = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Series)) + geom_point() + ggtitle("Sepsis All Samples PCA No Batch Correction Disease") + labs(y= "PC2", x = "PC1")
p

##assign healthy and disease samples
# all.disease <- all[,colnames(all) %in% meta.disease$geo_accession]
# all.control <- all[,colnames(all) %in% meta.control$geo_accession]
# batch.disease <- batch.ids[colnames(all.disease)]
# batch.control <-  batch.ids[colnames(all.control)]

#write.table(meta.disease, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/ALI_ARDS_only_logged_metadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


######### METADATA FOR RUN IF NOT RUNNING CONTROL AND DISEASE SEPERATLEY#######
##################PCA ONLY CONTROL AND SEPSIS SAMPLES####################
# pca_all_meta <- rbind(not_repeat,meta.control)
# pca_all_expr <- cbind(all.disease,all.control)
# 
# pca_batch_ids <- c(batch.disease, batch.control)
# pca_batch_ids <- pca_batch_ids[colnames(pca_all_expr)]

# keep <- c("control","sepsis", "septic shock", "SIRS", "sepsis/ARDS", "septic shock/ARDS", "sepsis/ALI", "septic ARDS")
# meta.data <- meta.data[meta.data$disease_simplified %in% keep,]
# meta.data <- meta.data[meta.data$geo_accession %in% colnames(all),]
# all <- all[,colnames(all) %in% meta.data$geo_accession]
# dim(meta.data)
# dim(all)

# - update batch.ids - #
# batch.ids <- batch.ids[colnames(all)]
# length(batch.ids)
# all_transposed <- t(pca_all_expr) #trasnpose matrix for PCA function
# 
# # - get rid of 0 values - #
# all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
# all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
# dtp <- data.frame('Health' = pca_all_meta$disease_simplified, "dat" = all_pca$x[,1:2])
# r <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Health)) + geom_point() + ggtitle("Sepsis PCA No Batch Correction Only Sepsis and Control Samples") + labs(y= "PC2", x = "PC1")
# r


##################BATCH CORRECTION FOR ALL SAMPLES ####################

##
batch.ids[grep("run",batch.ids)] <- "Faheem"


# - define model matrix in metadata - #
dim(meta.data)
dim(all)
length(batch.ids)
batch.ids <- batch.ids[colnames(all)]
length(batch.ids)
all.equal(names(batch.ids), colnames(all))
all[1:5,1:5]

temp.meta <- meta.data[!is.na(meta.data$disease_simplified),] #if samples do not have a diagnosis get rid of them so we can make our model matrix
dim(temp.meta)
all <- all[,colnames(all) %in% temp.meta$geo_accession]
dim(all)
design <- model.matrix(~disease_simplified, meta.data)
batch.ids <- batch.ids[names(batch.ids) %in% temp.meta$geo_accession]
length(batch.ids)
dim(design)
# all_t <- t(all)
# A4GNT <- all_t[,1]
# ZNF211 <- all_t[,4800]
# hist(A4GNT)
# hist(ZNF211)
# 
# all_logg <- log2(all + 1)
# all[1:5,1:5]
# 
# all_t <- t(all_logg)
# all_t[1:5,1:5]
# A4GNT <- all_t[,1]
# ZNF211 <- all_t[,4800]
# hist(A4GNT)
# hist(ZNF211)
# - get rid of negative values for ComBat - #
raw_merged <- as.matrix(all - min(all))
library(sva)
dat_batch_adjusted_norm_new <- ComBat(raw_merged, batch.ids, mod = design)

## additional batch correction for batches
batch.ids.temp <- readr::read_tsv(here("processed-data","master_batch_ids.tsv"))
batch.ids <- batch.ids.temp$batch.ids
names(batch.ids) <- batch.ids.temp$Sample
head(batch.ids)
head(temp_batch)
length(batch.ids)
batch.ids <- batch.ids[colnames(all)]
length(batch.ids)
all.equal(names(batch.ids), colnames(all))

raw_merged <- as.matrix(dat_batch_adjusted_norm_new - min(dat_batch_adjusted_norm_new))
dat_batch_adjusted_norm_new <- ComBat(raw_merged, batch.ids, mod = design)


dat_batch_adjusted_norm_new <- as.data.frame(dat_batch_adjusted_norm_new)
to_write <- rownames_to_column(dat_batch_adjusted_norm_new,"Gene")
write.table(to_write,here("processed-data","Oct18_batch_corrected__all_samples.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
# dat_batch_adjusted_norm_new <- readr::read_tsv(here("processed-data","Oct9_batch_corrected_all_samples.tsv")) - # forgot to write pCA to file 
# dat_batch_adjusted_norm_new[1:5,1:5]
# dat_batch_adjusted_norm_new <- column_to_rownames(dat_batch_adjusted_norm_new, "Gene")
dim(design)# - trasnpose matrix for PCA function - #
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)

# - get rid of 0 values - 0
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame("Series"=temp.meta$series, 'Health' = temp.meta$disease_simplified, 'accession' = temp.meta$geo_accession,"dat" = dat_pca$x[,1:2])
library(ggplot2)
q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color=Series)) + geom_point() + ggtitle("Sepsis PCA Batch Corrected All Samples") + labs(y= "PC2", x = "PC1")
q

min(dat_batch_adjusted_norm_new)
#write PCA 
dat <- dat_pca$x
dat <- as.data.frame(dat_pca$x)
dat <- rownames_to_column(dat,"Sample")
write.table(dat, here("processed-data","Oct18_batch_corrected_PCA_Combat_not_logged.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
dat_batch_adjusted_norm_new <- log2(dat_batch_adjusted_norm_new + 82)
dat_batch_adjusted_norm_new <- as.data.frame(dat_batch_adjusted_norm_new)
to_write <- rownames_to_column(dat_batch_adjusted_norm_new, "Gene")
write.table(to_write, here("processed-data","Oct18_all_samples_batch_corrected_Combat_log.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

##################BATCH CORRECTION DISEASE####################

# - get rid of negative values for ComBat - #
# raw_disease <- as.matrix(all.disease) - min(all.disease)
# library(sva)
# dat_batch_adjusted_norm_new <- ComBat(raw_disease, batch.disease)
# to_write <- as.data.frame(dat_batch_adjusted_norm_new) %>% rownames_to_column("Gene") # sampels are logged wit no repeats
# write.table(to_write, here("processed-data","disease_only_no_repeats_batch_corrected.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
# 
# # - trasnpose matrix for PCA function - #
# dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)
# 
# 
# ##################BATCH CORRECTION CONTROL####################
# 
# # - get rid of negative values for ComBat - #
# raw_control <- as.matrix(all.control) - min(all.control)
# dat_batch_adjusted_norm_new_control <- ComBat(raw_control, batch.control)
# 
# to_write <- as.data.frame(dat_batch_adjusted_norm_new_control) %>% rownames_to_column("Gene") # sampels are logged wit no repeats
# write.table(to_write, here("processed-data","control_only_batch_corrected.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
# 
# # - trasnpose matrix for PCA function - #
# dat_adjusted_norm_transposed_control <- t(dat_batch_adjusted_norm_new_control)
# 
# #########COMBINE DISEASE AND CONTROL SAMPLES ########
# dat_adjusted_norm_transposed <- rbind(dat_adjusted_norm_transposed, dat_adjusted_norm_transposed_control)
# meta <- rbind(not_repeat, meta.control)
# 
# ############PCA PLOT SAMPLES ###########
# 
# # - get rid of 0 values - 0
# dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
# library(plotly)
# set.seed(417)
# dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
# dtp <- data.frame("Series"=meta$series, 'Health' = meta$disease_simplified, 'accession' = meta$geo_accession,"dat" = dat_pca$x[,1:2])
# q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color=Health)) + geom_point() + ggtitle("Sepsis PCA Batch Corrected Sepsis") + labs(y= "PC2", x = "PC1")
# q

# ###write all samples to batch##
# to_write <- t(dat_adjusted_norm_transposed)
# to_write[1:5,1:5]
# dim(to_write)
# to_write <- as.data.frame(to_write)
# to_write <- rownames_to_column(to_write, "Gene")
# write.table(to_write, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/all_samples_logged_bcorrect.tsv", sep = "\t", col.names = TRUE)
# 
# ##write sepsis samples to batch
# temp <- as.data.frame(dat_batch_adjusted_norm_new)
# temp <- rownames_to_column(temp, "Gene")
# write.table(temp, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/sepsis_only_logged_bcorrect.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
# 
# ###GET RID OF REPEATS ####
# no.repeats <- meta.disease[meta.disease$is_repeat_sample == "no",]
# temp <- column_to_rownames(temp, "Gene")
# expr.no.repeats <- temp[,colnames(temp) %in% no.repeats$geo_accession]
# expr.no.repeats <- rownames_to_column(expr.no.repeats, "Gene")
# write.table(expr.no.repeats, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/ALI_ARDS_sepsis_only_logged_bcorrect_NO_REPEATS.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
# write.table(temp, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/sepsis_only_logged_bcorrect.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)




