require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()

#declare path to metadata files
df.path <- 'SepsisStudy'
df.dirs <- c("GSE13015","GSE13205","GSE141864","GSE28750","GSE57065","GSE65682","GSE66890","GSE74224","GSE79962","GSE95233","SRP049820","SRP078541","SRP132709","SRP168038","GSE3037","GSE63990")
df.names <- paste0(df.path,'/',df.dirs,'/metadata_',df.dirs,'.tsvtrimmed_metadata')
#read in files
myfiles <- lapply(df.names, readr::read_tsv)
seqtype <- c("microarray","microarray","microarray","microarray","microarray","microarray","microarray","microarray","microarray","microarray","rnaseq","rnaseq","rnaseq","rnaseq","microarray","microarray")
meta.df <- data.frame(df.dirs,seqtype) 
dflist <- list()
#make each file a data frame, add on the type of sequencing
for(i in 1:length(myfiles)){
  df <- as.data.frame(myfiles[[i]]) 
  df[,c(1:ncol(df))] <- sapply(df[,c(1:ncol(df))], as.character)
  temp <- merge(x = df, y = meta.df, by.x = "experiment_accession", by.y = "df.dirs", all.x = TRUE)
  dflist[[i]] <- temp
}
combined <- bind_rows(dflist,.id=NULL)
write.table(combined,file = "combined_metdata.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

library(ggplot2)
