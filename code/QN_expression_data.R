## Quantile normalization script for use in pipeline consistency with refine.bio  ## 
## human reference distribution for quantile normalization downloaded from: https://api.refine.bio/v1/qn_targets/
library(preprocessCore)
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr) 

#change depending on dataset
dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/new_SRA_dowloads/snaky/DATASET_READS/"
dataset <- "GSE131411"
expression <- readr::read_tsv(paste0(dir,dataset,"/expression_summarized.tsv"))
#Faheems data located: /home/leslie.smith1/blue_kgraim/leslie.smith1/new_SRA_dowloads/

#load reference distribution
reference <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/sepsis/rawData/new/refine_bio_QN_reference.tsv", col_names = FALSE)
#get reference into a vector
reference.v <- unlist(reference)

#ensure dataset duplicate genes have been dealt with
expression <- expression %>% group_by(Symbols) %>% summarise_all(mean)
expression <- expression[!(is.na(expression$Symbols)),] #get rid of NA row
dat.matrix <- column_to_rownames(expression,"Symbols")

#dat.matrix <- expression[,-1] ## alternate 
log <- scale(dat.matrix, center = TRUE, scale = TRUE) #scaling here does not affect the output of the quantile normalization, but we do it for consisitency
qn <- normalize.quantiles.use.target(log,reference.v, copy = TRUE)
colnames(qn) <- colnames(log)
rownames(qn) <- rownames(log)
qn[1:5,1:5]
qn <- as.data.frame(qn)

#write to file 
dat <- rownames_to_column(qn,"Gene")
write.table(dat,paste0(dir,dataset,"/expression_QN.tsv"),sep = "\t", col.names = TRUE,row.names=FALSE)






