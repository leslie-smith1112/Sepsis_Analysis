## CLUSTER PROFILER RESULTS ## 
## FOR NETWORKING FIGURES ## 

source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breast.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/network_grid_search/expression_elasticnet.R")

########## read in network files  ########
# temp.net <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/figures/base_eur_afr_networking_fig.tsv")
# eur.t <- temp.net[temp.net$ancestry == "eur",]
# afr.t <- temp.net[temp.net$ancestry == "afr",]
# e.genes <- unique(c(eur.t$Gene1, eur.t$Gene2))
# e.genes <- e.genes[order(e.genes),]
# a.genes <- unique(c(afr.t$Gene1, afr.t$Gene2))
# a.genes <- a.genes[order(a.genes)]
# 
# ####### get list of common genes, unique afr genes, and eur genes #######
# common <- a.genes[a.genes %in% e.genes]
# u.a.genes <- a.genes[!(a.genes %in% e.genes)]
# u.e.genes <- e.genes[!(e.genes %in% a.genes)]
# length(unique(common))

## load libraries for DE and GSEA ## 
library(clusterProfiler)
library("org.Hs.eg.db")
library(DESeq2)
library(msigdbr)
library(magrittr)
library(EnhancedVolcano)
library(apeglm)
library(ggplot2)
set.seed(12345)

## read in sample expression and ancestry data ## 
# expression <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_data_mrna_seq_v2_rsem.txt", col_names = TRUE)
# clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_tcga_pan_can_atlas_2018_clinical_data.tsv", col_names = TRUE)
# estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
# samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "BRCA",]
# samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
# samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
# samples.ancestry  <- samples.ancestry[samples.ancestry$patient %in% colnames(expression),]
# 
# ############## organize metadata #############
# ####### pull out sample and subtype from clinical data #######
# temp <- clinical %>% dplyr::select(`Sample ID`, Subtype)
# colnames(temp) <- c("sample_id", "subtype")
# 
# #######keep only luminal and basal samples #######
# tcga.binomial <- temp[temp$subtype == "BRCA_Basal",] #160 patients
# tcga.binomial1 <- temp[temp$subtype == "BRCA_LumA",] #20 patients
# binomial <- rbind(tcga.binomial,tcga.binomial1)
# tcga.binomial2 <- temp[temp$subtype == "BRCA_LumB",]
# binomial2 <- rbind(binomial,tcga.binomial2)
# binomial2 <- na.omit(binomial2)
# pattern1 <- "BRCA_LumA"
# pattern2 <- "BRCA_LumB"
# pattern3 <- "BRCA_Basal"
# 
# #### replace subtypes for logistic regression model####
# binomial2$subtype <- stri_replace_all_fixed(binomial2$subtype, pattern2,"Luminal")
# binomial2$subtype <- stri_replace_all_fixed(binomial2$subtype, pattern1,"Luminal")
# binomial2$subtype <- stri_replace_all_fixed(binomial2$subtype, pattern3,"Basal")
# binomial2$subtype <- as.factor(binomial2$subtype)
# 
# ####### keep only eur and afr samples #######
# cut.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr"| samples.ancestry$consensus_ancestry == "eur" ,] # CHANGE SAMPLES USED IN ANALYSIS HERE
# common.meta <- binomial2[binomial2$sample_id %in% cut.samples$patient,]
# # common 
# # u.a.genes 
# # u.e.genes
# ####### trim expression matrix to only genes we want for given signature ####### 
# expression <- trim.expr.matrix(expression, NULL, NULL) # get rid of null values for genes 
# common.expr <- expression[,colnames(expression) %in% common.meta$sample_id | colnames(expression) == "Hugo_Symbol"]
# common.meta$subtype <- as.factor(common.meta$subtype)
# common.expr <- common.expr[common.expr$Hugo_Symbol %in% common,]
# common.expr[1:5,1:5]
# dim(common.expr)
# com.expr <- column_to_rownames(common.expr, "Hugo_Symbol")
# common.df <- com.expr %>%
#   dplyr::select(common.meta$sample_id)
# 
# ####### Check if this is in the same order #######
# all.equal(colnames(common.df), common.meta$sample_id)
# common.meta$subtype <- as.factor(common.meta$subtype)


################ START OF DIFFERENTIAL EXPRESSION ################
####### differential expression code copied from: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html 
####### the code in put into function compute_DE defined at bottom of script 
# 
# genes.in.both <-  # this is a self defined function defined below, return deseq object 
# # afr.deseq_df <- genes.in.both
# # eur.deseq_df <- genes.in.both
# # common.deseq_df <- genes.in.both
# write.table(afr.deseq_df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/PAPER/afr_sig_differential_expression_new.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(eur.deseq_df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/PAPER/eur_sig_differential_expression_new.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(common.deseq_df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/PAPER/common_sig_differential_expression_new.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)
# 
# volcano_plot <- EnhancedVolcano::EnhancedVolcano(
#   deseq_df,
#   lab = deseq_df$Gene,
#   x = "log2FoldChange",
#   y = "padj",
#   title = "Differentially Expressed Eur Genes",
#   pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
# )
# 
# volcano_plot
# 
# msigdbr_species()
# 
# mm_hallmark_sets <- msigdbr(
#   species = "Homo sapiens", # Replace with species name relevant to your data
#   category = "C7"
# )
# head(mm_hallmark_sets)
# 
# 
# keytypes(org.Hs.eg.db)
# 
# 
# # eur.deseq_df
# # afr.deseq_df - for afr genes 
# # common.deseq_df 
# 
# any(duplicated(eur.deseq_df$Gene))
# any(duplicated(afr.deseq_df$Gene))
# any(duplicated(genes.in.both$Gene))
# 
# # Let's create a named vector ranked based on the log2 fold change values
# lfc_vector <- genes.in.both$log2FoldChange
# names(lfc_vector) <- genes.in.both$Gene
# head(lfc_vector)
# 
# # We need to sort the log2 fold change values in descending order here
# lfc_vector <- sort(lfc_vector, decreasing = TRUE)
# 
# 
# set.seed(2020)
# gsea_results <- GSEA(
#   geneList = lfc_vector, # Ordered ranked gene list
#   minGSSize = 0, # Minimum gene set size
#   maxGSSize = 500, # Maximum gene set set
#   pvalueCutoff = 0.05, # p-value cutoff
#   eps = 0, # Boundary for calculating the p value
#   seed = TRUE, # Set seed to make results reproducible
#   pAdjustMethod = "BH", # Benjamini-Hochberg correction
#   TERM2GENE = dplyr::select(
#     mm_hallmark_sets,
#     gs_name,
#     gene_symbol
#   )
# )
# 
# head(gsea_results@result)
# commongsea_result_df <- data.frame(gsea_results@result)
# 
# write.table(afrgsea_result_df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/PAPER/afr_sig_GSEA.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(eurgsea_result_df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/PAPER/eur_sig_GSEA.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(commongsea_result_df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/PAPER/common_sig_GSEA.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)

##



#################################################################################
## FUNCTION FOR DIFFERENTIAL EXPRESSION ANALYSIS - 
#################################################################################
library(clusterProfiler)
library("org.Hs.eg.db") #to get genes
library(DESeq2) #computeDE
library(msigdbr)
library(magrittr)
library(EnhancedVolcano)
library(apeglm)
library(ggplot2)
set.seed(12345)

#################################################################################
## READ IN ALL EXPRESSION DATA AND METADATA INFORMATION 
meta.dat<- readr::read_tsv(here("processed-data", "Oct91_updated_disease_only_metadata_clusters.tsv"))
expression.dat <- readr::read_tsv(here("processed-data","Oct91_updated_disease_only_expr_batch_corrected_logged.tsv"))
expression.dat[1:5,1:5]
expression.dat <- column_to_rownames(expression.dat, "Gene")
expression.dat[1:5,1:5]
# expression.dat <- as.matrix(expression.dat)
# colnames(meta.dat)
# expression.dat<- t(expression.dat)
# expression.dat[1:5,1:5]

meta_hold <- meta.dat
meta.dat <- meta_hold %>% dplyr::select("geo_accession","2_clusters")
head(meta.dat)
colnames(meta.dat) <- c("sample_id","subtype")

reorder_idx <- match(colnames(exression.dat),meta.dat$sample_id) 
meta.dat <- meta.dat[reorder_idx,]

## MAKE MODIFICATION TO DATA IF YOU NEED TO, WE NEED:
## gene expression matrix with samples as column names and gene as rownames 
## metadata is n x 2 (rows x columns) containing sample ID in the first column and subtype (EX: disease vs not)


###check that the column names of your expression matrix are same order as the samples in the metadata matrix
all.equal(colnames(expression.dat), meta.dat$sample_id) ## CHANGE HERE
# make sure subtype is a factor
meta.dat$subtype <- as.factor(meta.dat$subtype) ## CHANGE HERE - note if you don't call you subtype subtype, make sure you change "design = ~subtype" in  the function below,
## replacing subtye with whatever name you use

## do differentual expression analysis 
deseq_object <- compute_DE(expression.dat, meta.dat)



#################################################################################
## INPUT TO FUNCTION:
## dat: gene expression matrix with samples as column names and gene as rownames 
## dat.meta should be a matrix that is n x 2 (rows x columns) containing sample ID in the first column and subtype (EX: disease vs not)
#################################################################################
compute_DE <- function(dat, dat.meta){
  filtered_expression_df <- dat %>%
    dplyr::filter(rowSums(.) >= 10)
  gene_matrix <- round(filtered_expression_df)
  if(any(gene_matrix < 0) == TRUE){
    gene_matrix <- gene_matrix + abs(gene_matrix) # get rid of negative values 
  }
  ddset <- DESeqDataSetFromMatrix(
    # Here we supply non-normalized count data
    countData = gene_matrix,
    # Supply the `colData` with our metadata data frame
    colData = dat.meta,
    # Supply our experimental variable to `design`
    design = ~subtype ## POSSIBLE CHANGE HERE 
  )
  deseq_object <- DESeq(ddset)
  deseq_results <- results(deseq_object)
  deseq_results <- lfcShrink(
    deseq_object, # The original DESeq2 object after running DESeq()
    coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
    res = deseq_results # The original DESeq2 results table
  )
  
  deseq_df <- deseq_results %>%
    # make into data.frame
    as.data.frame() %>%
    # the gene names are row names -- make them a column for easy display
    tibble::rownames_to_column("Gene") %>%
    # add a column for significance threshold results
    dplyr::mutate(threshold = padj < 0.05) %>%
    # sort by statistic -- the highest values will be genes with
    # higher expression in RPL10 mutated samples
    dplyr::arrange(dplyr::desc(log2FoldChange))
  return(deseq_df)
}

write.table(deseq_object, here("processed-data","disease_deseq_diseasesimplified_results.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

deseq_object1 <- deseq_object
# see if faheems lipid genes are differentially expressed

faheem_genes <- deseq_object[deseq_object$Gene %in% c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9"),]
faheem_genes
#DHCR7 TRUE
#DHCR24 FALSE 

mm_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "H"
)
any(duplicated(deseq_object$Gene))




pf.og.gene.list <- deseq_object$log2FoldChange 
names(pf.og.gene.list) <- deseq_object$Gene
gene.list <- na.omit(pf.og.gene.list)
gene.list <- gene.list[order(gene.list, decreasing = TRUE)]
head(gene.list)
set.seed(2020)
gsea_results <- GSEA(
  geneList = gene.list, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  nPermSimple = 100000,
  TERM2GENE = dplyr::select(
    mm_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

head(gsea_results@result)

Hresults <- data.frame(gsea_results@result)
head(Hresults)
dim(Hresults)
write.table(Hresults, here("processed-data","H_GSEA_results_clusters.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
##STOP 


organism <- "org.Hs.eg.db"
keytypes(org.Hs.eg.db)
gse <- gseGO(geneList=gene.list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

gseaGO2 <- GSEA(gene.list, 
                TERM2GENE=c5, 
                minGSSize = 10,
                nPerm = 1000, 
                pvalueCutoff = 0.05,
                verbose=FALSE)
