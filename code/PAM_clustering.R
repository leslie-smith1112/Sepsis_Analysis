## PAM clustering 

library()

meta.dat<- readr::read_tsv(here("processed-data", "Oct18_updated_disease_only_metadata_clusters.tsv"))
expression.dat <- readr::read_tsv(here("processed-data","Oct18_updated_disease_only_expr_batch_corrected_not_logged_Faheem1.tsv"))
expression.dat[1:5,1:5]