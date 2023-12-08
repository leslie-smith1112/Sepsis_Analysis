## Look at datasets to see where genes get cut out 
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(plyr)
library(here)

#basing all 
faheem<- readr::read_tsv(here("datasets", "Faheem","expression_QN_new_sepsis.tsv"), col_names = TRUE)
to_name <- colnames(faheem)

new.dirs <- c("GSE131411","GSE154918","GSE185263","EMTAB1548","EMTAB4421","GSE100159", "GSE106878","GSE131761","GSE69063","GSE13015","GSE32707")
base_file <- "_QN_new_sepsis.tsv"
new.names <- here("datasets", new.dirs, paste0(new.dirs, "_QN_new_sepsis.tsv"))
all(file.exists(new.names))
newfiles <- lapply(new.names, readr::read_tsv)
lapply(newfiles,dim)


df.dirs <- c("GSE10474","GSE28750","GSE33118","GSE57065","GSE66890","GSE65682","GSE74224","GSE95233","SRP132709","SRP049820")
df.names <- here("datasets",df.dirs, paste0(df.dirs,".tsv"))
all(file.exists(df.names))
myfiles <- lapply(df.names, readr::read_tsv)

#check gene intersection
for(i in 1:length(newfiles)){
  print(length(intersect(faheem$Gene,newfiles[i][[1]]$Gene)))
}
#add datasets to see if there are any major drops
for(i in 1:length(newfiles)){
  print(ncol(newfiles[i][[1]])-1)
}
#check number of patients 
#"GSE131411","GSE154918","GSE185263","EMTAB1548","EMTAB4421","GSE100159", "GSE106878","GSE131761","GSE69063","GSE13015","GSE32707"
all <- faheem
for(i in 1:length(newfiles)){
  if(i != 4){
    all<- merge(all,newfiles[i][[1]], by = "Gene")
    print(dim(all))
  }
  
}

#get refine.bio genes in the correct format
library(org.Hs.eg.db)
annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = dat$Gene, column = "SYMBOL", keytype = "ENSEMBL"), dat)
annot.df$Symbols <- toupper(annot.df$Symbols)
for(i in 1:length(myfiles)){
  annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = myfiles[i][[1]]$Gene, column = "SYMBOL", keytype = "ENSEMBL"), myfiles[i][[1]])
  myfiles[i][[1]]$Symbols <- toupper(annot.df$Symbols)
}
#check intersection with faheems data
for(i in 1:length(myfiles)){
  print(length(intersect(faheem$Gene,myfiles[i][[1]]$Symbols)))
}
#check number of patients 
for(i in 1:length(myfiles)){
  print(ncol(myfiles[i][[1]])-1)
}

# make refine.bio genes symbols
new <- myfiles
for(i in 1:length(myfiles)){
  new[i][[1]]$Gene <- new[i][[1]]$Symbols
}

##all genes 
genes <- c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9","FDFT1","SQLE","ALOX5", "LSS","LCAT","LBR","CYP4Z1",
           "KCNH7", "CYP4X1","CYP46A1","APOA1","EBP","CYP4V2", "CYP39A1","SCARB1","HSD17B7","CYP4F22","ALOX15","ABCG1","NSDHL", "CYP4F12","ALOX12",
           "ABCA1","CYP4F11", "LOX","LIPE","TM7SF2","CYP4F8","PTGS1","LPA","LBR","CYP4F3","APOB","LBP","CETP","CYP4B1","CYP4A11",
           "FDFT1","LCAT","CYP4A22","PLTP","PTGS2")
dat <- c(newfiles, new)
all <- unique(faheem$Gene)
for(i in 1:length(dat)){
  if(i != 5){
    all<- intersect(all, unique(dat[i][[1]]$Gene))
    print(length(all))
  }
}
print(length(intersect(genes,all))) 


