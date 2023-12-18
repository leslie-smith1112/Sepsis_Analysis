
library(tidymodels)
library(workflows) #tidy models package for bundling model specs 
library(parsnip) #for modeling 
library(workflows) #put model in workflow 
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi) #string mutations 
library(rpart.plot)
library(vip)
set.seed(4831)
library(here)

######################################################################################################
# ## PHYLOFRAME PREPROCESSING BEFORE ELASTICNET RUN + ELASTICNET RUN##
######################################################################################################
## -- set up output directory -- ##
set.directory <- function(in.dir){
  if(!(file.exists(in.dir))){
    dir.create(in.dir)
  }
}

meta.dat<- readr::read_tsv(here("processed-data", "Oct18_updated_disease_only_metadata_clusters.tsv"))
expression.dat <- readr::read_tsv(here("processed-data","Oct18_updated_disease_only_expr_batch_corrected_not_logged_Faheem1.tsv"))
expression.dat[1:5,1:5]

expression.dat <- column_to_rownames(expression.dat, "Gene")
expression.dat <- t(expression.dat)
expression.dat[1:5,1:5]

#get train samples
temp <- c("Hypo", "Normo")
meta <- meta.dat[meta.dat$endotype_class %in% temp,]
dim(meta)

train_meta <- meta %>% dplyr::select(geo_accession, endotype_class)
dim(train_meta)

expression_train <- expression.dat[rownames(expression.dat) %in% train_meta$geo_accession,]
dim(expression_train)
expression_train[1:5,1:5]
expression_train <- rownames_to_column(as.data.frame(expression_train), "geo_accession")
expression_train[1:5,1:5]
training_dat <- merge(expression_train, train_meta, by="geo_accession")
training_dat[1:5,1:5]
head(training_dat$endotype_class)
training_dat <- column_to_rownames(training_dat, "geo_accession")

##test set 
test_meta <- meta.dat[!(meta.dat$geo_accession %in% train_meta$geo_accession),]
dim(test_meta)
expression_test <- expression.dat[rownames(expression.dat) %in% test_meta$geo_accession,]
expression_test[1:5,1:5]


## expects matrix with samples as row names 
## -- trim expression matrix to the selected samples and genes for current run --  ##
## result is genes x sample matrix

## expects matrix with samples as row names 
elasticnet.run <- function(in.matrix, directory, out_file, en.mix, seed=4831, endotype_class){
  print(en.mix)
  set.seed(seed)
  
  split <- initial_split(in.matrix, strata = endotype_class)
  train <- training(split)
  test <- testing(split)
  
  ## -- check that there is no intersection between train and test set -- ##
  length(intersect(rownames(train), rownames(test)))
  
  #validation set 
  val_set <- validation_split(train,
                              strata = endotype_class,
                              prop = 0.80)
  ### logistic regression ###
  ### set model - we will tune the penalty mixture of elasticnet is passed in ###
  ###ELASTIC NET ####
  lr_mod <-
    logistic_reg(penalty = tune(), mixture = en.mix) %>%
    set_engine("glmnet") #provides variable importance scores
  
  ### create recipe, remove indicator values that only contain 0 and noramlize the predictors ###
  rec <-
    recipe(endotype_class ~ ., data = train) 
  
  #### create workflow ###
  lr_wf <- workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(rec)
  
  #grid for tuning 
  lr_reg_grid <- tibble(penalty = 10^seq(-4,-1,length.out = 30))
  
  lr_res <-
    lr_wf %>% 
    tune_grid(val_set,
              grid = lr_reg_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(yardstick::roc_auc)) #first use roc metric to get tune9) parameter 
  
  ### view top models predicted with tuning ###
  lr_res %>%
    show_best("roc_auc", n = 15) %>%
    arrange(penalty)
  
  ### printing out all tested penalty to get correct slice number ###
  print(lr_res %>% 
          collect_metrics(), n = 30)
  
  ### get the slice with the best roc_auc ###
  best <- lr_res %>% 
    collect_metrics()
  
  ### selecting the penalty with best ROC, if there is a tie one is randomly selected ###
  w.score <- max(best$mean)
  chosen <- best[best$mean == w.score,]
  row <- sample(1:nrow(chosen), 1)
  the.one <- chosen[row,]
  penal <- the.one$penalty
  
  lr_best <-
    lr_res %>% 
    collect_metrics() %>%
    arrange(penalty) %>% 
    dplyr::slice(row)
  
  
  #========================LAST FIT =======================================
  #last fit 
  last_lr_mod <- 
    logistic_reg(penalty = penal, mixture = en.mix) %>% 
    set_engine("glmnet", importance = "impurity") #provides variable importance scores 
  
  #last workflow
  last_wf <- 
    lr_wf %>%
    update_model(last_lr_mod)
  
  last_fitt <- last_fit(last_wf, split) #split is our train and test data from beginning - it trains on both train set and validation set this time and test on test set  
  pred <- collect_metrics(last_fitt)
  write.table(pred, paste0(directory, "/",out_file,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  desc.scores <- last_fitt %>%
    extract_fit_engine() %>%
    vi()
  desc.scores <- desc.scores[desc.scores$Importance > 0,]
  print(desc.scores)
  ## -- added run information to check genes in ancestry signatures -- ## 
  #out_scores <- paste0(directory, "/", out_file,"_general_info.txt")
  print(paste0("Writing to ",out_file))
  write.table(desc.scores, file = paste0(directory, "/", out_file,"_all_sig.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  #look at model
  my_model <- extract_workflow(last_fitt)
  dat <- tidy(my_model)
  write.table(dat, paste0(directory, "/",out_file,"_model_coefficients.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  write_rds(my_model, paste0(directory,"/",out_file, "_EN_model.rds"))
  return(my_model)
}
in.matrix <- training_dat
en.mix <- 1
directory <- here("model_results")
out_file <- "HypoNormo_prediction_1_penalty"

my_model <- elasticnet.run(training_dat, directory, out_file = out_file, en.mix, seed=4831,"endotype_class")
result <- predict(my_model, expression_test)

predicted_pheno <- data.frame("geo_accession"  = rownames(expression_test), "endotype_class" = result$.pred_class)
write.table(predicted_pheno, here("model_results",paste0(out_file,"phenotype_predictions.tsv")), sep= "\t", col.names = TRUE, row.names = FALSE)
sig <- readr::read_tsv(here("model_results",paste0(out_file,"_all_sig.txt")))

genes <- c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9","FDFT1","SQLE","ALOX5", "LSS","LCAT","LBR","CYP4Z1",
           "KCNH7", "CYP4X1","CYP46A1","APOA1","EBP","CYP4V2", "CYP39A1","SCARB1","HSD17B7","CYP4F22","ALOX15","ABCG1","NSDHL", "CYP4F12","ALOX12",
           "ABCA1","CYP4F11", "LOX","LIPE","TM7SF2","CYP4F8","PTGS1","LPA","LBR","CYP4F3","APOB","LBP","CETP","CYP4B1","CYP4A11",
           "FDFT1","LCAT","CYP4A22","PLTP","PTGS2")

common <- sig[sig$Variable %in% genes,]
dim(common)


predicted <- readr::read_tsv(here("model_results","HypoNormo_prediction_05_penaltyphenotype_predictions.tsv"))
head(train_meta)
head(predicted)

all <- rbind(predicted, train_meta)
write.table(all, here("model_results","all_phenotypes.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)





