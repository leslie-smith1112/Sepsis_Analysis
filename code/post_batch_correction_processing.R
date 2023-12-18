
# Load libraries
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(plyr)
library(here)


## analysis ## 
all <- readr::read_tsv(here("processed-data","Oct18_batch_corrected__all_samples.tsv"))
metadata <- readr::read_tsv(here("processed-data", "master_sepsis_metadata_all_samples_add_run_labels.tsv"))
metadata$series[grep("run",metadata$series)] <- "UF_Faheem"

metadata <- metadata[!is.na(metadata$disease_simplified),] #if samples do not have a diagnosis get rid of them so we can make our model matrix
dim(metadata)
dim(all)

####get rid of meliodosis samples and repeated time points  ####
library(stringr)
#melioidosiss samples
metadata<- metadata[metadata$geo_accession %in% colnames(all),]
dim(metadata)
to_erase <- str_detect(metadata$disease, 'melioidosis')
meli_samples<- (metadata$geo_accession[to_erase])
meta.disease <- metadata[!(metadata$geo_accession %in% meli_samples),]
dim(meta.disease)

all <- all[,!(colnames(all) %in% meli_samples)]
dim(all)
batch.dat <- batch.ids[names(batch.ids) %in% colnames(all)]
length(batch.dat)

#get rid of repeat timepoints 
meta.disease$is_repeat_sample[meta.disease$is_repeat_sample == "na"] <- "no"
meta.disease$is_repeat_sample[meta.disease$is_repeat_sample == "time point: N/A"] <- "no"
to_remove <- meta.disease[meta.disease$is_repeat_sample == "yes",]
not_repeat <- meta.disease[!(meta.disease$geo_accession %in% to_remove$geo_accession),]
dim(not_repeat)


#write.table(not_repeat,here("processed-data","disease_only_metadata.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)
all<- column_to_rownames(all,"Gene")
all[1:5,1:5]
all <- all[,colnames(all) %in% not_repeat$geo_accession]
dim(all)
metadata <- not_repeat[not_repeat$disease_simplified %in% c("SIRS", "septic ARDs","sepsis","septic shock"),]
dim(metadata)
all <- all[,colnames(all) %in% metadata$geo_accession]
dim(all)
batch.disease <- batch.ids[colnames(all)]
length(batch.disease)



write.table(metadata,here("processed-data","Oct18_updated_disease_only_metadata.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)
to_write <- rownames_to_column(all,"Gene")
to_write[1:5,1:5]
dim(to_write)
write.table(to_write,here("processed-data","Oct18_updated_disease_only_expr_batch_corrected_not_logged_Faheem1.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)
batch.dat <- data.frame(batch.disease, names(batch.disease))
head(batch.dat)
colnames(batch.dat) <- c("Batch","Sample")
write.table(batch.dat,here("processed-data","Oct18_updated_disease_only_batch_ids.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)


dat_adjusted_norm_transposed <- t(all)
metadata
# - get rid of 0 values - 0
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame("Series"=metadata$series, 'Health' = metadata$disease_simplified, 'accession' = metadata$geo_accession,"dat" = dat_pca$x[,1:2])
library(ggplot2)
q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color=Health)) + geom_point() + ggtitle("Sepsis PCA Batch Corrected All Samples") + labs(y= "PC2", x = "PC1")
q






