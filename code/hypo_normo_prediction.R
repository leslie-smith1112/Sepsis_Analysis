
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

expression_train <- expression.dat[rownames(expression.dat) %in% train_meta$geo_accession,]
dim(expression_train)
expression_train[1:5,1:5]
expression_train <- rownames_to_column(as.data.frame(expression_train), "geo_accession")
expression_train[1:5,1:5]
training_dat <- merge(expression_train, train_meta, by="geo_accession")
training_dat[1:5,1:5]
head(training_dat$)

test_meta <- meta.dat[!(meta.dat$geo_accession %in% train$geo_accession),]
dim(test_meta)
table(test_meta$ARDs)
test_meta$ARDs[test_meta$ARDs %in% temp] <- "ARDs"
test_meta <- test_meta  %>% select(geo_accession, ARDs) ## test metadata here 
test_meta$ARDs[is.na(test_meta$ARDs)] <- "noARDS"
table(test_meta$ARDs)
expression.test <- expression.test[rownames(expression.test) %in% test_meta$geo_accession,]
dim(expression.test)
expression.test[1:2, 4875:4878]
expression.test <- as.data.frame(expression.test)
expression.test[1:5,1:5]
expression.test <- rownames_to_column(expression.test, "Gene")
expression_matrix <- merge(expression.test, test_meta, by.x = "Gene", by.y = "geo_accession")
dim(expression_matrix)
expression_matrix <- column_to_rownames(expression_matrix, "Gene")
expression_matrix[1:5,1:5]
dim(expression_matrix)

expression.dat <- as.data.frame(expression.train)
expression.dat[1:5,1:5]
expression.dat <- rownames_to_column(expression.dat, "Gene")
expression_matrix <- merge(expression.dat, train, by.x = "Gene", by.y = "geo_accession")
dim(expression_matrix)
expression_matrix <- column_to_rownames(expression_matrix, "Gene")
expression_matrix[1:5,1:5]

dim(expression_matrix)
expression_matrix[1:2, 4875:4879]
## expects matrix with samples as column names 
## -- trim expression matrix to the selected samples and genes for current run --  ##
## result is genes x sample matrix
in.matrix <- expression_matrix
en.mix <- 1
directory <- here("model_results")
out_file <- "ARDs_prediction"

elasticnet.run(expression_matrix, directory, out_file = out_file, en.mix, seed=4831,"ARDs")
## expects matrix with samples as row names 
elasticnet.run <- function(in.matrix, directory, out_file, en.mix, seed=4831, ARDs){
  print(en.mix)
  set.seed(seed)
  
  split <- initial_split(in.matrix, strata = ARDs)
  train <- training(split)
  test <- testing(split)
  
  ## -- check that there is no intersection between train and test set -- ##
  length(intersect(rownames(train), rownames(test)))
  
  #validation set 
  val_set <- validation_split(train,
                              strata = ARDs,
                              prop = 0.80)
  ### logistic regression ###
  ### set model - we will tune the penalty mixture of elasticnet is passed in ###
  ###ELASTIC NET ####
  lr_mod <-
    logistic_reg(penalty = tune(), mixture = en.mix) %>%
    set_engine("glmnet") #provides variable importance scores
  
  ### create recipe, remove indicator values that only contain 0 and noramlize the predictors ###
  rec <-
    recipe(ARDs ~ ., data = train) 
  #%>%
  # step_zv(all_predictors()) %>%
  # step_normalize(all_predictors())
  
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
  
  #------commented out for library issues, only code for confusion matrix------------- After tuning to get confusion matrix ---------------#
  ####TEMPORARILY COMMENTED OUT BECAUSE OF SPEC MASKING ISSUES #####
  # lr_res <-
  #   lr_wf %>%
  #   tune_grid(val_set,
  #             grid = lr_reg_grid,
  #             control = control_grid(save_pred = TRUE), #save model predictions
  #             metrics = yardstick::metric_set(spec)) #run this second round to get specificity
  # print("best")
  # #for spec
  # lr_auc <-
  #   lr_res %>%
  #   collect_predictions(parameters = lr_best)
  # 
  # predicted <- lr_auc$.pred_class
  # actual <- lr_auc$subtype
  # confusion_matrix <- table(predicted,actual) #want this printed out
  # out_matrix <- paste0(directory, "/", out_file,"_confusion_matrix.txt")
  # write.table(confusion_matrix, file = out_matrix, col.names = T)
  
  
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

model.metrics(my_model = my_model, dat = expression.test, directory = here("model_results"), out_file = "ARDsvalidation", "ARDs","noARDS")
model.metrics <- function(my_model, dat, directory, out_file, subtype1, subtype2)
{
  pred_prob <- predict(my_model, dat, type = "prob")
  pred_class <- predict(my_model, dat, type = "class")
  results <- dat %>% dplyr::select(ARDs) %>% bind_cols(pred_class, pred_prob)
  results <- rownames_to_column(results, "sample_id")
  write.table(results, paste0(directory,"/", out_file,"_results.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  # new.res1 <- bind_cols(dat$subtype, results$.pred_Basal, results$.pred_Luminal, results$.pred_class)
  # colnames(new.res) <- c("subtype", ".pred_Basal", ".pred_Luminal", ".pred_class")
  new.res <- bind_cols(dat$ARDs, results[,3], results[,4], results$.pred_class)
  colnames(new.res) <- c("ARDs", paste0(".pred_ARDs"), paste0(".pred_noARDs"), ".pred_class")
  confusion <- conf_mat(results, truth = ARDs,estimate = .pred_class)
  results$ARDs <- as.factor(results$ARDs)
  confusion
  conf.df <- as.data.frame(confusion$table)
  write_delim(conf.df, paste0(directory,"/", out_file,"_confusion_matrix.tsv"), delim  = "\t")
  #auc <- roc_auc(new.res, truth = subtype, estimate = .pred_Basal)
  auc <- roc_auc(new.res, ARDs, paste0(".pred_ARDs")) 
  new.res$ARDs <- as.factor(new.res$ARDs)
  #- all can be put in 1 matrix - #
  senss <- yardstick::sens(results,  ARDs, .pred_class)
  specc <- yardstick::spec(results,  ARDs,  .pred_class)
  acc <- accuracy(results,  ARDs, .pred_class)
  prec <- yardstick::precision(results,  ARDs, .pred_class)
  re <- yardstick::recall(results, ARDs, .pred_class)
  f <- yardstick::f_meas(results, ARDs, .pred_class)
  kapp <- kap(results,  ARDs, .pred_class)
  mccc <- mcc(results, ARDs, .pred_class)
  metrics <- rbind(auc, acc, senss, specc, prec, re, f, kapp, mccc)
  print(metrics)
  write.table(metrics,paste0(directory,"/", out_file,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  #return(metrics)
}
## MODEL METRICS FOR SUBSAHARAN AFRICAN VALIDATION SET:
# model.metrics <- function(my_model, dat, directory, out_file, subtype1, subtype2)
# {
#   pred_prob <- predict(my_model, dat, type = "prob")
#   pred_class <- predict(my_model, dat, type = "class")
#   results <- dat %>% dplyr::select(subtype) %>% bind_cols(pred_class, pred_prob)
#   results <- rownames_to_column(results, "sample_id")
#   write.table(results, paste0(directory,"/", out_file,"_results.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
#   false_negative <- nrow(results[results$.pred_class == "Luminal",])
#   true_positive <- nrow(results[results$.pred_class == "Basal",])
#   rec <- true_positive/(true_positive + false_negative)
# 
#   prec <- true_positive/(true_positive + 0)
#   to.write <- data.frame(c("Recall", "Precision") ,c(rec, prec))
#   write.table(to.write, paste0(directory,"/", out_file,"_results.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
# }



