#Read in all datasets and metadata and put them into 1 common expression matrix and 1 common metadata file.

# Load libraries
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(plyr)
library(here)

###############################READ IN DATA###############################
##GEO DATASETS##
# - read in quantile normalized file for GEO samples- #
new.dirs <- c("GSE131411","GSE154918","GSE185263","EMTAB1548","EMTAB4421","GSE100159", "GSE106878","GSE131761","GSE69063","GSE13015","GSE32707")
base_file <- "_QN_new_sepsis.tsv"
new.names <- here("datasets", new.dirs, paste0(new.dirs, "_QN_new_sepsis.tsv"))
all(file.exists(new.names))
newfiles <- lapply(new.names, readr::read_tsv)

# - get batch.ids - one is subtracted to account for the gene name column - #
new.batch.ids <- unlist(sapply(seq_along(new.dirs), function(i) {rep(new.dirs[i],ncol(newfiles[[i]])-1)} ))
length(new.batch.ids)

# - reduce to common genes - #
new.dat <- newfiles %>% purrr::reduce(inner_join)
dim(new.dat)

#make Gene rownames so that they aren't used as names for the batch ids
new.dat <- column_to_rownames(new.dat, "Gene") 
new.dat[1:5,1:5]
names(new.batch.ids) <- colnames(new.dat)

# - we lose the below sample because there is a type in the metadata downloaded from GEO, this should keep it - #
names(new.batch.ids)[names(new.batch.ids) == "SS_39esp2"] <- "SS_39exp2"
colnames(new.dat)[colnames(new.dat) == "SS_39esp2"] <- "SS_39exp2"

new.dat <- rownames_to_column(new.dat, "Gene")
new.dat[1:5,1:5]
dim(new.dat) #6769 x 1694

## REFINE.BIO DATASETS ##
# - read in datasets - #
df.dirs <- c("GSE10474","GSE28750","GSE33118","GSE57065","GSE66890","GSE65682","GSE74224","GSE95233","SRP132709","SRP049820")
df.names <- here("datasets",df.dirs, paste0(df.dirs,".tsv"))
all(file.exists(df.names))
myfiles <- lapply(df.names, readr::read_tsv)

# - get batch.ids - one is subtracted to account for the gene name column - #
old.batch.ids <- unlist(sapply(seq_along(df.dirs), function(i) {rep(df.dirs[i],ncol(myfiles[[i]])-1)} ))

# - reduce to common genes - #
dat <- myfiles %>% purrr::reduce(inner_join)
dat <- column_to_rownames(dat, "Gene")
dat[1:5,1:5]

# - get symbol id for genes from old dataset - #
library(org.Hs.eg.db)
annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = rownames(dat), column = "SYMBOL", keytype = "ENSEMBL"), dat)
annot.df$Symbols <- toupper(annot.df$Symbols)
annot.df[1:5,1:5]
annot.df <- annot.df %>% as.data.frame
annot.df <- annot.df[!(is.na(annot.df$Symbols)),]
rownames(annot.df) <- NULL
dat <- column_to_rownames(annot.df, "Symbols")
dat[1:5,1:5]
# add sample names to batches
names(old.batch.ids) <- colnames(dat)

# # - make sure batch.ids are in the same order as expression matrix - #
old.batch.ids <- old.batch.ids[colnames(dat)]
dat <- rownames_to_column(dat, "Gene")
dim(dat)

# - put all datasets together - #
faheem<- readr::read_tsv(here("datasets", "Faheem","expression_QN_new_sepsis.tsv"), col_names = TRUE)
to_name <- colnames(faheem)
to_name <- to_name[-1]
all <- merge(x = dat, y = new.dat, by = "Gene")
all <- merge(all, faheem, by = "Gene")
all <- column_to_rownames(all, "Gene")
dim(all)
all <- all %>% dplyr::select(-SR191) # get rid of this sample becuase of low quality (from Dongyuan)

# - put all batch.ids together - # 
f.batch.ids <- rep("UF_Faheem",(ncol(faheem)-1)) #all faheems runs are viewed as 1 (aka not seperated by run)
names(f.batch.ids) <- to_name
batch.ids <- c(old.batch.ids, new.batch.ids,f.batch.ids)
expression_dat <- all

#########################################
## edit the ids from the datasets we reprocessed through SRA ## - they have SRA names instead of GEO names used in metadata
length(batch.ids)
batch.ids <- batch.ids[!(batch.ids %in% c("GSE131411","GSE154918","GSE185263"))]
length(batch.ids)

#for datasets we processed through SRA we used the GEO Ids in the metadata file but now our matrices have the SRR Ids - so we change that here
GSE131411_meta <- read_csv(here("datasets","GSE131411","GSE131411_metadata.txt"))
sample_map <- GSE131411_meta %>% dplyr::select(Run, `GEO_Accession (exp)`)

#several sample here map to the same GSE ID - we take the second run because there was more coverage 
dupsss <- sample_map$`GEO_Accession (exp)`[duplicated(sample_map$`GEO_Accession (exp)`)]
sample_temp <- sample_map[sample_map$`GEO_Accession (exp)` %in% dupsss,]
discardme <- seq_len(nrow(sample_temp))%%2
discardme <- sample_temp[discardme == 1, ]
expression_dat[1:5,1:5]
expression_dat <- expression_dat[,!(colnames(expression_dat) %in% discardme$Run)]
sample_map <- sample_map[!(sample_map$Run %in% discardme$Run),]
names(expression_dat) <- plyr::mapvalues(colnames(expression_dat), from = sample_map$Run, to = sample_map$`GEO_Accession (exp)`)
#expression_dat <- column_to_rownames(expression_dat,"Gene")
temp.batch  <- rep("GSE131411", length(sample_map$`GEO_Accession (exp)`))
names(temp.batch) <- sample_map$`GEO_Accession (exp)`
head(temp.batch)
batch.ids <- c(batch.ids, temp.batch)
length(batch.ids)

GSE154918_meta <- read_csv(here("datasets","GSE154918","GSE154918_metadata.txt"))
sample_map <- GSE154918_meta %>% dplyr::select(Run, `GEO_Accession (exp)`)
names(expression_dat) <- plyr::mapvalues(colnames(expression_dat), from = sample_map$Run, to = sample_map$`GEO_Accession (exp)`)
temp.batch  <- rep("GSE154918", length(sample_map$`GEO_Accession (exp)`))
names(temp.batch) <- sample_map$`GEO_Accession (exp)`
length(temp.batch)
head(temp.batch)
batch.ids <- c(batch.ids, temp.batch)
length(batch.ids)

GSE185263_meta <- read_csv(here("datasets","GSE185263","GSE185263_metadata.txt"))
sample_map <- GSE185263_meta %>% dplyr::select(Run, `Sample Name`)
names(expression_dat) <- plyr::mapvalues(colnames(expression_dat), from = sample_map$Run, to = sample_map$`Sample Name`)
temp.batch  <- rep("GSE185263", length(sample_map$`Sample Name`))
names(temp.batch) <- sample_map$`Sample Name`
length(temp.batch)
head(temp.batch)
batch.ids <- c(batch.ids, temp.batch)
length(batch.ids)
temp_batch <- batch.ids[names(batch.ids) %in% colnames(expression_dat)]
length(temp_batch)
to_write_batch <- data.frame(temp_batch, names(temp_batch))
colnames(to_write_batch) <- c("batch.ids", "Sample")
write.table(to_write_batch, here("processed-data","master_batch_ids_faheem_same_run.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
expression_dat <- rownames_to_column(expression_dat, "Gene")
write.table(expression_dat, here("processed-data","Oct18_master_expression_nonnormalized.tsv"), sep = "\t",col.names = TRUE, row.names = FALSE)

# expression_logged <- column_to_rownames(expression_dat, "Gene")
# expression_logged[1:5,1:5]
# expression_logged <- log2(expression_logged + 1)
# expression_logged <- rownames_to_column(expression_logged,"Gene")
# write.table(expression_logged, here("processed-data","Oct8_master_expression_logged.tsv"), sep = "\t",col.names = TRUE, row.names = FALSE)

############################# METADATA MODIFICATIONS ####################################
### %% mention matching each studies attribute, we worked a lot with Faheem to read papers and is coagulation failure in one paper the same
## %% as low platelet count in another.
# - get sample metadata - #
# - we get warnings on read in here because of bacteria explanation for one study, doesn't affect anything - #
#metadata <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/sepsismetadata_master_one_drive_2.xlsx")
#metadata_t <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/sepsis_only_meta_NO_REPEATS_combinesALIARDS_with_Faheem_samples.xlsx")

# metadata <- readxl::read_xlsx(here("processed-data", "master_sepsis_metadata_all_samples_Sept4_2023.xlsx"))
# head(metadata)
# dim(metadata)
# table(metadata$ARDs)
#
# #bin faheems runs as they were read in - the run dataframes come from faheem's script
# metadata$series[metadata$geo_accession %in% colnames(run11)] <- "run11"
# metadata$series[metadata$geo_accession %in% colnames(run12)] <- "run12"
# metadata$series[metadata$geo_accession %in% colnames(run10)] <- "run10"
# metadata$series[metadata$geo_accession %in% colnames(run9)] <- "run9"
# metadata$series[metadata$geo_accession %in% colnames(run8)] <- "run8"
# metadata$series[metadata$geo_accession %in% colnames(run7)] <- "run7"
# metadata$series[metadata$geo_accession %in% colnames(run6)] <- "run6"
# metadata$series[metadata$geo_accession %in% colnames(run5)] <- "run5"
# metadata$series[metadata$geo_accession %in% colnames(run4)] <- "run4"
# metadata$series[metadata$geo_accession %in% colnames(run3)] <- "run3"
# metadata$series[metadata$geo_accession %in% colnames(run2)] <- "run2"
# metadata$series[metadata$geo_accession %in% colnames(run1)] <- "run1"
# table(metadata$series)

#write.table(metadata,here("processed-data", "master_sepsis_metadata_all_samples_add_run_labels.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

### BELOW IS TO DELETE ### 
# not_included <- write_all[,!(colnames(write_all) %in% meta.data$geo_accession)]
# not_sample <- colnames(not_included)
# write.table(not_sample, here("processed-data","excluded_samples_from_merge.tsv"),sep = "\t")
##### COMABT DISEASE AND CONTROL SEPERATLEY ######
#disease <- c("sepsis", "septic shock", "SIRS", "sepsis/ARDS", "septic shock/ARDS", "sepsis/ALI")
#disease <- c("sepsis", "septic shock", "SIRS", "septic ARDs")

## keep
# to.keep  <- c("GSE10474", "GSE32707", "GSE66890","SRP049820")
#cut.meta <- metadata[metadata$series %in% to.keep,]
# meta.t <- meta.data
#meta.data <- cut.meta
# disease <- c("no sepsis", "sepsis", "SIRS", "sepsis/ARDS", "septic shock", "septic shock/ARDS", "sepsis/ALI")
#control <- "control"
# 
# meta.disease.temp <- meta.data[meta.data$disease_simplified %in% disease,]
# meta.control <- meta.data[meta.data$disease_simplified %in% control,]



