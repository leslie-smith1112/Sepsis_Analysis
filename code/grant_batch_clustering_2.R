### batch correction and clustering for grant ###### 
## lipid genes 
genes <- c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9")
########################################################
## newly added geo sets, quantile normalized ##
df.path <- '/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS'
base_file <- "_QN_new_sepsis.tsv"

#new.dirs <- c("EMTAB1548","EMTAB4421","GSE100159","GSE106878","GSE131411","GSE131761","GSE154918","GSE185263","GSE32707","GSE69063","GSE13015")

#new.dirs <- c("EMTAB1548","EMTAB4421","GSE100159","GSE106878","GSE131411","GSE131761","GSE154918","GSE185263","GSE32707","GSE69063","GSE13015")

# only datasets that have all genes
new.dirs <- c("GSE131411","GSE131761","GSE154918","GSE185263","GSE69063")

## only datasets that have all genes or only missing first 
# new.dirs <- c("GSE131411","GSE131761","GSE154918","GSE185263","GSE69063")
# 
# ## only datasets that have all genes or only missing first and/or 6th 
# new.dirs <- c("EMTAB4421","GSE100159","GSE131411","GSE131761","GSE154918","GSE185263","GSE32707","GSE69063","GSE13015")

new.names <- paste0(df.path,'/',new.dirs,'/',new.dirs,base_file)
newfiles <- lapply(new.names, readr::read_tsv)
length(newfiles)

new.batch.ids <- unlist(sapply(seq_along(new.dirs), function(i) {rep(new.dirs[i],ncol(newfiles[[i]])-1)} ))
new.dat <- newfiles %>% purrr::reduce(inner_join)
dim(new.dat)#6743 x 1685
#make Gene rownames so that they are used as names for the batch ids
new.dat <- column_to_rownames(new.dat, "Gene") 

##log gene values 
new.dat <- log2(new.dat + 1)

# - put batch.ids and gene expression matrix in same order - #
names(new.batch.ids) <- colnames(new.dat)
new.batch.ids <- new.batch.ids[colnames(new.dat)]
new.dat <- rownames_to_column(new.dat, "Gene")
new.dat[1:5,1:5]
dim(new.dat)

########################################################
## old datasets from refine.bio ##
df.path <- '/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS'
#df.dirs <- c("GSE10474","GSE28750","GSE33118","GSE57065","GSE66890","GSE65682","GSE74224","GSE95233","SRP132709","SRP049820")

#df.dirs <- c("GSE10474","GSE28750","GSE33118","GSE57065","GSE66890","GSE65682","GSE74224","GSE95233","SRP132709","SRP049820")

##only datatsets that keep all genes ## 
df.dirs <- c("SRP132709","SRP049820")

## only datatsets that keep all genes  or missing first  ## 
#df.dirs <- c("GSE28750","GSE33118","GSE57065","GSE66890","GSE65682","GSE74224","GSE95233","SRP132709","SRP049820")

## ALI/ARDS studies 
#df.dirs <- c("GSE10474","GSE66890","SRP049820")
df.names <- paste0(df.path,'/',df.dirs,'/',df.dirs,'.tsv')
myfiles <- lapply(df.names, readr::read_tsv)

old.batch.ids <- unlist(sapply(seq_along(df.dirs), function(i) {rep(df.dirs[i],ncol(myfiles[[i]])-1)} ))
dat <- myfiles %>% purrr::reduce(inner_join)
dat <- column_to_rownames(dat, "Gene")
dat[1:5,1:5]
dat <- log2(dat + 1)

annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = rownames(dat), column = "SYMBOL", keytype = "ENSEMBL"), dat)
annot.df$Symbols <- toupper(annot.df$Symbols)
annot.df[1:5,1:5]

annot.df <- annot.df %>% as.data.frame
annot.df <- annot.df[!duplicated(annot.df$Symbols),] # TODO: don't do this lol
annot.df <- annot.df[complete.cases(annot.df),]
rownames(annot.df) <- annot.df$Symbols
annot.df <- annot.df[,-1]
annot.df[1:5,1:5]
dat <- annot.df
names(old.batch.ids) <- colnames(dat)

old.batch.ids <- old.batch.ids[colnames(dat)]
dat <- rownames_to_column(dat, "Gene")
dim(dat)
f.samples <- faheem()
faheem <- f.samples$expression.mat
f.batch.ids <- f.samples$batch.ids.dat


batch.ids <- c(old.batch.ids, new.batch.ids, f.batch.ids)
all <- merge(x = dat, y = new.dat, by = "Gene")
all <- merge(all, faheem, by = "Gene")
all <- column_to_rownames(all, "Gene")
# get rid of this sample becuase of low quality (from Dongyuan) 
all <- all %>% dplyr::select(-SR191) 
all[1:5,1:5]
dim(all)
############################# METADATA MODIFICATIONS ####################################

## metadata are added together and edited in excel before beling loaded in here ## 
metadata <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/rewritten_Faheem_sepsis_only_meta_NO_REPEATS_combinesALIARDS.tsv")
dim(metadata)
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
## getting rid of these runs for now becayse they mess up the batch correction ## 
# erase <- c("run10", "run11", "run12")
# meta.data <- meta.data[!(meta.data$series %in% erase),]
# - four samples get cut from gene expression (i.e. they aren't in the metadata)- we can fix one noted above) but haven't found where the other 3 are in metadata - #
all <- all[,colnames(all) %in% meta.data$geo_accession]
dim(all)
batch.ids <- batch.ids[colnames(all)]
length(batch.ids)
to.return <- list(batch.id.dat = batch.ids, expr.dat = all)


# read.val <- read.files()
# 
# all <- read.val$expr.dat
# batch.ids <- read.val$batch.id.dat
meta.data$disease_simplified[meta.data$disease_simplified == "septic  ARDs"] <- "septic ARDs"##TODO can't find why this sample is like this 
# old.all <- all
# all <- old.all[,colnames(old.all) %in% meta.norun$geo_accession]
##########################################################INITIAL PCA ALL SAMPLES############################################################
#trasnpose matrix for PCA function
all_transposed <- t(all) 

# - get rid of 0 values - #
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
library(ggplot2)
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame('Series' = meta.data$disease_simplified, "dat" = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Series)) + geom_point() + ggtitle("Sepsis All Samples PCA No Batch Correction") + labs(y= "PC2", x = "PC1")
p

##### COMBAT DISEASE  ######
disease <- c("sepsis", "septic shock", "SIRS", "sepsis/ALI/ARDs", "septic shock/ARDs", "septic ARDs")
meta.data$disease_simplified[meta.data$geo_accession == "K011"] <- "sepsis"
###TODO 
meta.disease.temp <- meta.data[meta.data$disease_simplified %in% disease,]
dim(meta.disease.temp)
table(meta.disease.temp$series)
####GET RID OF MELIODOSIS ANF SEPERATE DISEASE VS CONTROL ####
library(stringr)
##assign healthy and disease samples 
all.disease <- all[,colnames(all) %in% meta.disease.temp$geo_accession]
dim(all.disease)
batch.disease <- batch.ids[colnames(all.disease)]
length(batch.disease)

#write.table(meta.disease, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/pca_checks_all_studies/ALI_ARDS_only_metadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


################BATCH CORRECTION DISEASE####################
# - get rid of negative values for ComBat - #
# add disease info to combat 
#meta.disease.temp <- meta.disease.temp[!(meta.disease.temp$series == "run10"),]

#mod.mat <- model.matrix(~as.factor(disease_simplified), data = meta.disease.temp)
mod.mat <- model.matrix(~as.factor(disease_simplified), data = meta.disease.temp)
raw_disease <- as.matrix(all.disease) - min(all.disease)

library(sva)
#dat_batch_adjusted_norm_new <- cbind(dat_batch_adjusted_norm_new.sepsis,dat_batch_adjusted_norm_new.sepsis.ARDs,dat_batch_adjusted_norm_new.septic.shock,dat_batch_adjusted_norm_new.sepsis.SIRS)
dim(raw_disease)

dat_batch_adjusted_norm_new <- ComBat(raw_disease, batch.disease)
to_write <- as.data.frame(dat_batch_adjusted_norm_new) %>% rownames_to_column("Gene") # sampels are logged wit no repeats

# - trasnpose matrix for PCA function - #
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new) 
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
library(plotly)
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
#meta.temp <- met.temp[!is.na(met.temp$disease_simplified),]
dtp <- data.frame("Series"=meta.disease.temp$series, 'Health' = meta.disease.temp$disease_simplified, 'accession' = meta.disease.temp$geo_accession,"dat" = dat_pca$x[,1:2])
q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color=Health)) + geom_point() + ggtitle("Sepsis PCA Batch Corrected Sepsis with No Mod.Mat") + labs(y= "PC2", x = "PC1")
q
write.table(to_write, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/SEPSIS/05_21_23_disease_samples_Faheem_withNOmodmat_30batch2_missingCYP51A1andMSMO1.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
#####

######################################################
## clustering ## 

expr.dat <- dat_batch_adjusted_norm_new
mat <- as.matrix(expr.dat)
mat[1:5,1:5]

#pearson correlation distance - preprocess by median centering data 
d <- sweep(mat, 1, apply(mat, 1, median, na.rm = T))
dis.mat = as.dist(1-cor(d,method="pearson"))
length(dis.mat)
class(dis.mat)

library(ConsensusClusterPlus)
title = "5_22_23_Cluster_all_NOmodmat"

results = ConsensusClusterPlus(dis.mat,maxK=6,reps=100,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388,plot="png")
icl = calcICL(results,title=title,plot="png")
cluster_scores <- icl[["clusterConsensus"]]
consensus <- icl[["itemConsensus"]][1:5,]
cluster_id <- results[[3]][["consensusClass"]]
head(cluster_id)
class(results)

clust.dat <- as.data.frame(cluster_id)
cluster.dat <- rownames_to_column(clust.dat, "sample_id")
table(cluster.dat$cluster_id)
library(cluster)
library(stats)
clus_vect <- as.vector(cluster_id)
names(clus_vect) <- names(cluster_id)
sil <- silhouette(clus_vect, dis.mat)
summary(sil)

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







































