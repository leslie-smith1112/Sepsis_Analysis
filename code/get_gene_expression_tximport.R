### TXIMPORT FOR GENE EXPRESSION FROM SALMON OUTPUT FILES ### 
#REFERENCE: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon 
#REFERENCE FOR ENSEBLE DATABASES: https://bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/ensembldb.html
#REFERENCE: https://www.bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/ensembldb.html
library("optparse") 
library(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)

##NOTE: You must change the database you used depending on the transcriptome version you used to annotate your samples in aligner
library(EnsDb.Hsapiens.v86)
#short cut:
edb <- EnsDb.Hsapiens.v86
Tx <- transcripts(edb, return.type = "data.frame")
head(Tx)
tx2gene <- Tx[,c("tx_name","gene_id")]
head(tx2gene)
colnames(tx2gene) <- c("TX_NAME","GENE_ID")
tx2gene <- tx2gene[complete.cases(tx2gene), ]
#specify dataset, script assumes that accession file for dataset is name Sra_Acc.txt
dataset <- "GSE185263"
in_dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/new_SRA_dowloads/snaky/DATASET_READS"


# get transcript name from reference - note if you want to use ensembl database code is below. Consult the first reference here. 
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# transcripts("EnsDb.Hsapiens.v38", return.type="Dataframe")
# k <- keys(txdb, keytype = "TXID")
# tx2gene <- select(txdb, k, "GENEID", "TXID")
# tail(tx2gene)
# tx2gene <- tx2gene[complete.cases(tx2gene), ]

#read in sample files
dir <- paste0(in_dir,"/",dataset)
sample.list <- readr::read_tsv(paste0(in_dir,"/",dataset,"/Sra_Acc.txt"),col_names = FALSE)
files <- file.path(in_dir,dataset,"refine_bio_salmon_quant",sample.list$X1,"quant.sf") # dir should be to salmon wuant file 
head(files)
names(files) <- paste0(sample.list$X1)
all(file.exists(files))

#compute abundance matrix
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
names(txi)
head(txi$counts)

#prepare matrix to get symbol names
count_mat <- as.data.frame(txi$counts)
library(tidyr)
count_mat <- tibble::rownames_to_column(count_mat, "ENSEMBL")
count_mat[1:5,1:5]
# library( "biomaRt" ) #example code for mouse gene id mapping
# ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id", "mgi_symbol"), filters = "ensembl_gene_id",values = rownames(my_normalizedMatrix), mart = ensembl )

library(org.Hs.eg.db)

#annot.df$Symbols <- toupper(annot.df$Symbols)

## get gene names in symbol 
keytypes(org.Hs.eg.db)
# cols <- c("SYMBOL", "ENSEMBL") ## check with kiley about the doubles 
# map_genes <- select(org.Hs.eg.db, keys=count_mat$ENSEMBL, columns=cols, keytype="ENSEMBL")
annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = count_mat$ENSEMBL, column = "SYMBOL", keytype = "ENSEMBL"), count_mat)
# map genes our gene expression matrix
#new_final <- merge(map_genes, count_mat, by = "ENSEMBL")
annot.df[1:5,1:5]
annot.df <- annot.df[,-2]
annot.df[1:5,1:5]
dim(annot.df)
# new_final <- new_final[,-1]
# new_final[1:5,1:5]
# dim(new_final)
library(dplyr)
#Counts = aggregate(annot.df,FUN = mean,by=list(annot.df$Symbols))
temp <- annot.df[duplicated(annot.df$Symbols),]
temp[1:5,1:5]
temp <- data.frame(table(annot.df$Symbols))
temp <- temp[temp$Freq > 1,]
dim(temp)
# temp <- temp[!(is.na(temp$Symbols)),]

#mine <- annot.df[annot.df$Symbols %in% c("CLN3", "CT47A11","DEFB103B","TSPAN6"),]
annot.sum <- annot.df %>% group_by(Symbols) %>% summarise_all(mean)

#mine[mine$Symbols == "CLN3",]

#annot.sum <- annot.df %>% group_by(Symbols) %>% summarise(mean)
dim(annot.sum)
annot.sum[1:5,1:5]
write.table(annot.sum, paste0(dir,"/expression_summarized.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)


# ### FOR USE WITH ENSEMBL ###
# library(EnsDb.Hsapiens.v86) ## ussed this version because it was the one that is on hipergator 
# library(ensembldb)
# edb <- EnsDb.Hsapiens.v86
# #### 
# 
# 
