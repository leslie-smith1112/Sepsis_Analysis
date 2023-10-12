
## clustering 
library(plyr)
library(dplyr)
library(here)
library(tidyr)
library(tibble)
library(ConsensusClusterPlus)

expr.dat <- readr::read_tsv(here("processed-data","disease_only_no_repeats_batch_corrected_no_quotes.tsv"))
expr.dat[1:5,1:5]
expr.dat <- column_to_rownames(expr.dat, "Gene")
expr.dat[1:5,1:5]
mat <- as.matrix(expr.dat)
mat[1:5,1:5]

metadata <- readr::read_tsv(here("processed-data","disease_only_metadata.tsv"))

#pearson correlation distance - preprocess by median centering data 
d <- sweep(mat, 1, apply(mat, 1, median, na.rm = T))
dis.mat = as.dist(1-cor(d,method="pearson"))
length(dis.mat)
class(dis.mat)

library(ConsensusClusterPlus)
title = "Sept6Clustering"

results = ConsensusClusterPlus(dis.mat,maxK=6,reps=100,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388,plot="png")
icl = calcICL(results,title=title,plot="png")
cluster_scores <- icl[["clusterConsensus"]]
consensus <- icl[["itemConsensus"]][1:5,]
cluster_id <- results[[6]][["consensusClass"]]
head(cluster_id)
class(results)
table(cluster_id)
clust.dat <- as.data.frame(cluster_id)
cluster.dat <- rownames_to_column(clust.dat, "sample_id")
colnames(cluster.dat) <- c("sample_id", "6_clusters")
head(cluster.dat)

metadata <- merge(metadata,cluster.dat,by.x = "geo_accession", "sample_id")
metadata$`2_clusters` <- paste0("cluster",metadata$`2_cluster`)

write.table(metadata, here("processed-data", "disease_metadata_with_cluster.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

table(cluster.dat$cluster_id)
library(cluster)
library(stats)
clus_vect <- as.vector(cluster_id)
names(clus_vect) <- names(cluster_id)
sil <- silhouette(clus_vect, dis.mat)
summary(sil)

##did not use anything below this
expr.dat <- as.data.frame(expr.dat)

### seperate clusters #### 
clus1 <- cluster.dat[cluster.dat$cluster_id == 1,]
clus2 <- cluster.dat[cluster.dat$cluster_id == 2,]
clus3 <- cluster.dat[cluster.dat$cluster_id == 3,]
#clus4 <- cluster.dat[cluster.dat$cluster_id == 4,]
expr1 <- expr.dat[,colnames(expr.dat) %in% clus1$sample_id]
expr2 <- expr.dat[,colnames(expr.dat) %in% clus2$sample_id]
expr3 <- expr.dat[,colnames(expr.dat) %in% clus3$sample_id]
#expr4 <- expr.dat[,colnames(expr.dat) %in% clus4$sample_id]

meta1 <- meta.disease.temp[meta.disease.temp$geo_accession %in% clus1$sample_id,]
meta2 <- meta.disease.temp[meta.disease.temp$geo_accession %in% clus2$sample_id,]
meta3 <- meta.disease.temp[meta.disease.temp$geo_accession %in% clus3$sample_id,]
#meta4 <- meta.disease.temp[meta.disease.temp$geo_accession %in% clus4$sample_id,]

cluster.dat$cluster_id[cluster.dat$cluster_id == 1] <- "cluster1"
cluster.dat$cluster_id[cluster.dat$cluster_id == 2] <- "cluster2"
cluster.dat$cluster_id[cluster.dat$cluster_id == 3] <- "cluster3"
#cluster.dat$cluster_id[cluster.dat$cluster_id == 4] <- "cluster4"
write.table(cluster.dat, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/05_22_23_all_genes_clusterID19batches.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
#######
t.met <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/05_22_23_all_genes_clusterID19batches.tsv")
# for chi squared test 
mets <- metadata[metadata$geo_accession %in% t.met$sample_id,]
merged <- merge(mets,t.met, by.x = "geo_accession", by.y = "sample_id")
dim(merged)
table(merged$endotype_class, merged$cluster_id)
chisq.test(merged$endotype_class == "Hypo", merged$cluster_id == 'cluster3')

chisq.test(merged$endotype_class == "Normo", merged$cluster_id == 'cluster3')


chisq.test( dat$hypo == 'HYPO', dat$cluster =='1') 



# for tumor map with ARDS. and HYpo normo
metadata$ARDs[metadata$ARDs == "UF ARDS"] <- "UF ARDs"
hypnorm <- metadata[metadata$endotype_class == "Hypo" | metadata$endotype_class == "Normo",]
hypnorm <- hypnorm[!is.na(hypnorm$endotype_class),]
hypnormacc <- hypnorm %>% select(geo_accession, endotype_class)

new <- data.frame(ards.samples$geo_accession)
new$cluster_id <- rep("ARDs",116)
colnames(new) <- c("sample_id", "cluster_id")
t.met <- t.met[!(t.met$sample_id %in% hypnormacc$sample_id),]
colnames(hypnormacc) <- c("sample_id", "cluster_id")
new.t <- rbind(t.met, hypnormacc)

ards.samples <- metadata[!is.na(metadata$ARDs),]
dim(ards.samples)

write.table(new.t, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/05_22_23_clusterID_withHypoNormo3.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
question <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/05_22_23_clusterID_withHypoNormo3.tsv")
#######################################################################
expression.dat <- cbind(expr1, expr2)
meta.dat <- rbind(clus1, clus2)
colnames(meta.dat) <- c("sample_id", "subtype")

all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- as.factor(meta.dat$subtype) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object1 <- compute_DE(expression.dat, meta.dat)

deseq_object1[deseq_object1$Gene == "CYP51A1",]#FALSE
deseq_object1[deseq_object$Gene == "CYP4F2",] #TRUE

deseq_object1[deseq_object1$Gene == "LDLR",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR24",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR7",] #FALSE
deseq_object1[deseq_object1$Gene == "MSMO1",] #FALSE
deseq_object1[deseq_object1$Gene == "HMGCR",] #TRUE
deseq_object1[deseq_object1$Gene == "PCSK9",] #TRUE
dim(expression.dat)
#######################################################################

expression.dat <- cbind(expr1, expr3)
meta.dat <- rbind(clus1, clus3)
colnames(meta.dat) <- c("sample_id", "subtype")

all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- factor(meta.dat$subtype, levels = c(3,1)) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
levels(meta.dat$subtype)

## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object2 <- compute_DE(expression.dat, meta.dat)

deseq_object2[deseq_object$Gene == "CYP51A1",]
deseq_object2[deseq_object$Gene == "CYP4F2",] 

deseq_object2[deseq_object2$Gene == "LDLR",] 
deseq_object2[deseq_object2$Gene == "DHCR24",] 
deseq_object2[deseq_object2$Gene == "DHCR7",] #TRUE
deseq_object2[deseq_object2$Gene == "MSMO1",] #
deseq_object2[deseq_object2$Gene == "HMGCR",] 
deseq_object2[deseq_object2$Gene == "PCSK9",] #

#######################################################################
expression.dat <- cbind(expr1, expr4)
meta.dat <- rbind(clus1, clus4)
colnames(meta.dat) <- c("sample_id", "subtype")

all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- as.factor(meta.dat$subtype) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object1 <- compute_DE(expression.dat, meta.dat)

deseq_object1[deseq_object1$Gene == "CYP51A1",]#FALSE
deseq_object1[deseq_object$Gene == "CYP4F2",] #TRUE

deseq_object1[deseq_object1$Gene == "LDLR",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR24",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR7",] #FALSE
deseq_object1[deseq_object1$Gene == "MSMO1",] #FALSE
deseq_object1[deseq_object1$Gene == "HMGCR",] #TRUE
deseq_object1[deseq_object1$Gene == "PCSK9",] #TRUE

#######################################################################

expression.dat <- cbind(expr2, expr3)
meta.dat <- rbind(clus2, clus3)
colnames(meta.dat) <- c("sample_id", "subtype")

all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- as.factor(meta.dat$subtype) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
levels(meta.dat$subtype)
## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object3 <- compute_DE(expression.dat, meta.dat)

deseq_object3[deseq_object$Gene == "CYP51A1",]#FALSE
deseq_object3[deseq_object$Gene == "CYP4F2",] #TRUE
deseq_object3[deseq_object3$Gene == "LDLR",] #TRUE
deseq_object3[deseq_object3$Gene == "DHCR24",] #TRUE
deseq_object3[deseq_object3$Gene == "DHCR7",] #TRUE
deseq_object3[deseq_object3$Gene == "MSMO1",] #TRUE
deseq_object3[deseq_object3$Gene == "HMGCR",] #TRUE
deseq_object3[deseq_object3$Gene == "PCSK9",] #TRUE
#######################################################################

expression.dat <- cbind(expr2, expr4)
meta.dat <- rbind(clus2, clus4)
colnames(meta.dat) <- c("sample_id", "subtype")

all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- as.factor(meta.dat$subtype) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object1 <- compute_DE(expression.dat, meta.dat)

deseq_object1[deseq_object1$Gene == "CYP51A1",]#FALSE
deseq_object1[deseq_object$Gene == "CYP4F2",] #TRUE

deseq_object1[deseq_object1$Gene == "LDLR",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR24",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR7",] #FALSE
deseq_object1[deseq_object1$Gene == "MSMO1",] #FALSE
deseq_object1[deseq_object1$Gene == "HMGCR",] #TRUE
deseq_object1[deseq_object1$Gene == "PCSK9",] #TRUE

#######################################################################

expression.dat <- cbind(expr3, expr4)
meta.dat <- rbind(clus3, clus4)
colnames(meta.dat) <- c("sample_id", "subtype")

all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- as.factor(meta.dat$subtype) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object1 <- compute_DE(expression.dat, meta.dat)

deseq_object1[deseq_object1$Gene == "CYP51A1",]#FALSE
deseq_object1[deseq_object$Gene == "CYP4F2",] #TRUE

deseq_object1[deseq_object1$Gene == "LDLR",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR24",] #FALSE
deseq_object1[deseq_object1$Gene == "DHCR7",] #FALSE
deseq_object1[deseq_object1$Gene == "MSMO1",] #FALSE
deseq_object1[deseq_object1$Gene == "HMGCR",] #TRUE
deseq_object1[deseq_object1$Gene == "PCSK9",] #TRUE














