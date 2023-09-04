# Load libraries
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(here)
##### get metadata ##### - TEMP
# new.dirs <- c("GSE32707","GSE10474","GSE66890","SRP049820")
# new.names <- paste0(df.path,'/',new.dirs,'/metadata_',new.dirs,".tsv")
# samples.meta <- lapply(new.names, readr::read_tsv)

###############################READ IN DATA###############################

##GEO & SRA DATASETS##
# - read in quantile normalized file for GEO samples- #
new.dirs <- c("EMTAB1548","EMTAB4421","GSE100159","GSE106878","GSE131411","GSE131761","GSE154918","GSE185263","GSE32707","GSE69063","GSE13015")
base_file <- "_QN_new_sepsis.tsv"
new.names <- here("datasets", new.dirs, paste0(new.dirs, "_QN_new_sepsis.tsv"))
all(file.exists(new.names))
newfiles <- lapply(new.names, readr::read_tsv)

# - get batch.ids - one is subtracted to account for the gene name column - #
new.batch.ids <- unlist(sapply(seq_along(new.dirs), function(i) {rep(new.dirs[i],ncol(newfiles[[i]])-1)} ))

# - reduce to common genes - #
new.dat <- newfiles %>% purrr::reduce(inner_join)
dim(new.dat)#6769 x 1694
#make Gene rownames so that they are used as names for the batch ids
new.dat <- column_to_rownames(new.dat, "Gene") 

##log gene values
new.dat <- log2(new.dat + 1)

# - add sample names to batch ids
names(new.batch.ids) <- colnames(new.dat)

# - we lose the below sample (commented out) because there is a type in the metadata downloaded from GEO, this should keep it - #
names(new.batch.ids)[names(new.batch.ids) == "SS_39esp2"] <- "SS_39exp2"
colnames(new.dat)[colnames(new.dat) == "SS_39esp2"] <- "SS_39exp2"

new.dat <- rownames_to_column(new.dat, "Gene")
new.dat[1:5,1:5]
dim(new.dat) #6769 x 1694

## REFINE.BIO DATASETS ##
# - read in datasets - #
df.dirs <- c("GSE10474","GSE28750","GSE33118","GSE57065","GSE66890","GSE65682","GSE74224","GSE95233","SRP132709","SRP049820")
## ALI/ARDS studies
#df.dirs <- c("GSE10474","GSE66890","SRP049820")
df.names <- here("datasets",df.dirs, paste0(df.dirs,".tsv"))
all(file.exists(df.names))
#df.names <- paste0(df.path,'/',df.dirs,'/',df.dirs,'.tsv')
myfiles <- lapply(df.names, readr::read_tsv)

# - get batch.ids - one is subtracted to account for the gene name column - #
old.batch.ids <- unlist(sapply(seq_along(df.dirs), function(i) {rep(df.dirs[i],ncol(myfiles[[i]])-1)} ))

# - reduce to common genes - #
dat <- myfiles %>% purrr::reduce(inner_join)
dat <- column_to_rownames(dat, "Gene")
dat[1:5,1:5]

dat <- log2(dat + 1)

# - get symbol id for genes from old dataset - #
library(org.Hs.eg.db)
annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = rownames(dat), column = "SYMBOL", keytype = "ENSEMBL"), dat)
annot.df$Symbols <- toupper(annot.df$Symbols)
annot.df[1:5,1:5]
# - line below is meant to replace the "TODO: don't do this", however it takes a long time to run - #
#annot.sum <- annot.df %>% group_by(Symbols) %>% summarise_all(mean)
#annot.df <- annot.sum
annot.df <- annot.df %>% as.data.frame
# annot.df <- annot.df[!duplicated(annot.df$Symbols),] # TODO: don't do this lol
# annot.df <- annot.df[complete.cases(annot.df),]
annot.df <- annot.df[!(is.na(annot.df$Symbols)),]
dat <- column_to_rownames(annot.df, "Symbols")
dat[1:5,1:5]


names(old.batch.ids) <- colnames(dat)

# # - make sure batch.ids are in the same order as expression matrix - #
# old.batch.ids <- old.batch.ids[colnames(dat)]
dat <- rownames_to_column(dat, "Gene")

# - put old and new together - #
batch.ids <- c(old.batch.ids, new.batch.ids, f.batch.ids)
all <- merge(x = dat, y = new.dat, by = "Gene")
all <- merge(all, faheem, by = "Gene")
all <- column_to_rownames(all, "Gene")

all <- all %>% dplyr::select(-SR191) # get rid of this sample becuase of low quality (from Dongyuan)
############################# METADATA MODIFICATIONS ####################################
### %% mention matching each studies attribute, we worked a lot with Faheem to read papers and is coagulation failure in one paper the same
## %% as low platelet count in another.
# - get sample metadata - #
# - we get warnings here because of bacteria explanation for one study, doesn't affect anything - #
#metadata <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/sepsismetadata_master_one_drive_2.xlsx")
metadata <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/sepsis_only_meta_NO_REPEATS_combinesALIARDS_with_Faheem_samples.xlsx")

head(new.faheem)
#must do for all faheems runs
metadata$series[metadata$geo_accession %in% colnames(run11)] <- "run11"
metadata$series[metadata$geo_accession %in% colnames(run12)] <- "run12"
metadata$series[metadata$geo_accession %in% colnames(run10)] <- "run10"
metadata$series[metadata$geo_accession %in% colnames(run9)] <- "run9"
metadata$series[metadata$geo_accession %in% colnames(run8)] <- "run8"
metadata$series[metadata$geo_accession %in% colnames(run7)] <- "run7"
metadata$series[metadata$geo_accession %in% colnames(run6)] <- "run6"
metadata$series[metadata$geo_accession %in% colnames(run5)] <- "run5"
metadata$series[metadata$geo_accession %in% colnames(run4)] <- "run4"
metadata$series[metadata$geo_accession %in% colnames(run3)] <- "run3"
metadata$series[metadata$geo_accession %in% colnames(run2)] <- "run2"
metadata$series[metadata$geo_accession %in% colnames(run1)] <- "run1"
     # - keep metadata only from samples we have expression data for - #
meta.data <- metadata[metadata$geo_accession %in% colnames(all),]

# - four samples get cut from gene expression - we can fix one noted above) but haven't found where the other 3 are in metadata - #
all <- all[,colnames(all) %in% meta.data$geo_accession]
dim(all)
batch.ids <- batch.ids[colnames(all)]
length(batch.ids)



##########################################################INITIAL PCA ALL SAMPLES############################################################
all_transposed <- t(all) #trasnpose matrix for PCA function


# - get rid of 0 values - #
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
library(ggplot2)
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame('Series' = meta.data$series, "dat" = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Series)) + geom_point() + ggtitle("Sepsis All Samples PCA No Batch Correction") + labs(y= "PC2", x = "PC1")
p

##### COMABT DISEASE AND CONTROL SEPERATLEY ######
#disease <- c("sepsis", "septic shock", "SIRS", "sepsis/ARDS", "septic shock/ARDS", "sepsis/ALI")
disease <- c("sepsis", "septic shock", "SIRS", "sepsis/ALI/ARDS", "septic shock/ARDS", "septic ARDS")

## keep
# to.keep  <- c("GSE10474", "GSE32707", "GSE66890","SRP049820")
#cut.meta <- metadata[metadata$series %in% to.keep,]
# meta.t <- meta.data
#meta.data <- cut.meta
# disease <- c("no sepsis", "sepsis", "SIRS", "sepsis/ARDS", "septic shock", "septic shock/ARDS", "sepsis/ALI")
control <- "control"

meta.disease.temp <- meta.data[meta.data$disease_simplified %in% disease,]
meta.control <- meta.data[meta.data$disease_simplified %in% control,]

####GET RID OF MELIODOSIS ANF SEPERATE DISEASE VS CONTROL ####
library(stringr)
#melioidosiss samples
to_erase <- str_detect(meta.disease.temp$disease, 'melioidosis')
meli_samples<- (meta.disease.temp$geo_accession[to_erase])
meta.disease <- meta.disease.temp[!(meta.disease.temp$geo_accession %in% meli_samples),]

##assign healthy and disease samples
all.disease <- all[,colnames(all) %in% meta.disease.temp$geo_accession]
all.control <- all[,colnames(all) %in% meta.control$geo_accession]
batch.disease <- batch.ids[colnames(all.disease)]
batch.control <-  batch.ids[colnames(all.control)]

write.table(meta.disease, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/ALI_ARDS_only_logged_metadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


######### METADATA FOR RUN IF NOT RUNNING CONTROL AND DISEASE SEPERATLEY#######

##################PCA ONLY CONTROL AND SEPSIS SAMPLES####################
keep <- c("control","sepsis", "septic shock", "SIRS", "sepsis/ARDS", "septic shock/ARDS", "sepsis/ALI", "septic ARDS")
meta.data <- meta.data[meta.data$disease_simplified %in% keep,]
meta.data <- meta.data[meta.data$geo_accession %in% colnames(all),]
all <- all[,colnames(all) %in% meta.data$geo_accession]
dim(meta.data)
dim(all)

# - update batch.ids - #
batch.ids <- batch.ids[colnames(all)]
length(batch.ids)
all_transposed <- t(all) #trasnpose matrix for PCA function

# - get rid of 0 values - #
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame('Health' = meta.data$disease_simplified, "dat" = all_pca$x[,1:2])
r <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Health)) + geom_point() + ggtitle("Sepsis PCA No Batch Correction Only Sepsis and Control Samples") + labs(y= "PC2", x = "PC1")
r



##################BATCH CORRECTION FOR ALL SAMPLES ####################
# - define model matrix in metadata - #
meta.data$matrix <- 1
meta.data$matrix[meta.data$disease_simplified == "control"] <- 0
mod = model.matrix(~as.factor(matrix), data=meta.data)

# - get rid of negative values for ComBat - #
raw_merged <- as.matrix(all) - min(all)
library(sva)
dat_batch_adjusted_norm_new <- ComBat(raw_merged, batch.ids, mod = mod)

# - trasnpose matrix for PCA function - #
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)

# - get rid of 0 values - 0
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
library(plotly)
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame("Series"=meta.data$series, 'Health' = meta.data$disease_simplified, 'accession' = meta.data$geo_accession,"dat" = dat_pca$x[,1:2])
q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color=Health)) + geom_point() + ggtitle("Sepsis PCA Batch Corrected Sepsis and Controls") + labs(y= "PC2", x = "PC1")
q

##################BATCH CORRECTION DISEASE####################

# - get rid of negative values for ComBat - #
raw_disease <- as.matrix(all.disease) - min(all.disease)
library(sva)
dat_batch_adjusted_norm_new <- ComBat(raw_disease, batch.disease)
to_write <- as.data.frame(dat_batch_adjusted_norm_new) %>% rownames_to_column("Gene") # sampels are logged wit no repeats
write.table(to_write, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/04_03_all_samples_Faheem.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

# - trasnpose matrix for PCA function - #
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)


##################BATCH CORRECTION CONTROL####################

# - get rid of negative values for ComBat - #
raw_control <- as.matrix(all.control) - min(all.control)
dat_batch_adjusted_norm_new_control <- ComBat(raw_control, batch.control)

# - trasnpose matrix for PCA function - #
dat_adjusted_norm_transposed_control <- t(dat_batch_adjusted_norm_new_control)

#########COMBINE DISEASE AND CONTROL SAMPLES ########
dat_adjusted_norm_transposed <- rbind(dat_adjusted_norm_transposed, dat_adjusted_norm_transposed_control)
meta <- rbind(meta.disease, meta.control)

############PCA PLOT SAMPLES ###########

# - get rid of 0 values - 0
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
library(plotly)
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame("Series"=meta$series, 'Health' = meta$disease_simplified, 'accession' = meta$geo_accession,"dat" = dat_pca$x[,1:2])
q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color=Health)) + geom_point() + ggtitle("Sepsis PCA Batch Corrected Sepsis") + labs(y= "PC2", x = "PC1")
q

###write all samples to batch##
to_write <- t(dat_adjusted_norm_transposed)
to_write[1:5,1:5]
dim(to_write)
to_write <- as.data.frame(to_write)
to_write <- rownames_to_column(to_write, "Gene")
write.table(to_write, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/all_samples_logged_bcorrect.tsv", sep = "\t", col.names = TRUE)

##write sepsis samples to batch
temp <- as.data.frame(dat_batch_adjusted_norm_new)
temp <- rownames_to_column(temp, "Gene")
write.table(temp, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/sepsis_only_logged_bcorrect.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

###GET RID OF REPEATS ####
no.repeats <- meta.disease[meta.disease$is_repeat_sample == "no",]
temp <- column_to_rownames(temp, "Gene")
expr.no.repeats <- temp[,colnames(temp) %in% no.repeats$geo_accession]
expr.no.repeats <- rownames_to_column(expr.no.repeats, "Gene")
write.table(expr.no.repeats, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/ALI_ARDS_sepsis_only_logged_bcorrect_NO_REPEATS.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)



write.table(temp, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/sepsis_only_logged_bcorrect.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)







