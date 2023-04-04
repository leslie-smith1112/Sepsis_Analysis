#### ADDING FAHEEMS RNA SEQ DATA FOR TUMORMAP ### 

df.path <- '/blue/fguirgis/dongyuanwu/RNASeqExpMat'
#new.dirs <- c("GSE32707","GSE10474","GSE66890","SRP049820")
f.new.dirs <- list(1:12)
f.new.dirs <- unlist(f.new.dirs)
# f.new.dirs2 <- list(7:9)
# f.new.dirs2 <- unlist(f.new.dirs2)
first.new.names <- paste0(df.path,'/run',f.new.dirs,'_matrix_single.txt')
#second.new.names <- paste0(df.path,'/run',f.new.dirs2,'1_matrix.txt')

## rename columns 
run1 <- readr::read_table(first.new.names[1],skip= 1)
run1 <- run1[,-2:-6]
run1[1:5,1:5]
colnames(run1) ## 17:20
temp <- substr(colnames(run1), 17, 20)
temp <- temp[-1]
temp <- c("Gene",temp)
colnames(run1) <- temp


run2 <- readr::read_table(first.new.names[2],skip= 1)
run2 <- run2[,-2:-6]
colnames(run2)
temp <- substr(colnames(run2), 17, 20)
temp <- temp[-1]
temp <- c("Gene",temp)
colnames(run2) <- temp
run2[1:5,1:5]

run3 <- readr::read_table(first.new.names[3],skip= 1)
run3 <- run3[,-2:-6]
run3[1:5,1:5]
colnames(run3)
temp <- colnames(run3)
temp <- temp[-1]
temp[1] <- substr(temp[1], 17, 21)
temp[2:3] <- substr(temp[2:3], 17, 20)
temp[4] <- substr(temp[4], 17, 21)
temp[5:16] <- substr(temp[5:16], 17, 20)

temp <- c("Gene",temp)
colnames(run3) <- temp
run3[1:5,1:5]

run4 <- readr::read_table(first.new.names[4],skip= 1)
run4 <- run4[,-2:-6]
colnames(run4)
temp <- substr(colnames(run4), 17, 20)
temp <- temp[-1]
temp <- c("Gene",temp)
colnames(run4) <- temp
run4[1:5,1:5]

run5 <- readr::read_table(first.new.names[5],skip= 1)
run5 <- run5[,-2:-6]
colnames(run5) # 1:10= 17:20, 11:12 = 17:21, 13 = 17:22, 14:15 = 17:21, 16 = 17:20 

temp <- colnames(run5)
temp <- temp[-1]
temp[1:10] <- substr(temp[1:10], 17, 20)
temp[11:12] <- substr(temp[11:12], 17, 21)
temp[13] <- substr(temp[13], 17, 22)
temp[14:15] <- substr(temp[14:15], 17, 21)
temp[16] <- substr(temp[16], 17, 20)
temp <- c("Gene",temp)
colnames(run5) <- temp
run5[1:5,1:5]

run6 <- readr::read_table(first.new.names[6],skip= 1)
run6 <- run6[,-2:-6]
colnames(run6) # 1:12 = 17:20, 13:14 = 17:21, 15:16 = 17:20
temp <- colnames(run6)
temp <- temp[-1]
temp[1:12] <- substr(temp[1:12], 17, 20)
temp[13:14] <- substr(temp[13:14], 17, 21)
temp[15:16] <- substr(temp[15:16], 17, 20)
temp <- c("Gene",temp)
colnames(run6) <- temp
run6[1:5,1:5]

run7 <- readr::read_table(first.new.names[7],skip= 1)
run7 <- run7[,-2:-6]
colnames(run7) # 18:22
temp <- substr(colnames(run7), 18, 22)
temp <- temp[-1]
temp <- c("Gene",temp)
colnames(run7) <- temp
run7[1:5,1:5]


run8 <- readr::read_table(first.new.names[8],skip= 1)
run8 <- run8[,-2:-6]
colnames(run8)# 18:22
temp <- substr(colnames(run8), 18, 22)
temp <- temp[-1]
temp <- c("Gene",temp)
colnames(run8) <- temp
run8[1:5,1:5]

run9 <- readr::read_table(first.new.names[9],skip= 1)
run9 <- run9[,-2:-6]
colnames(run9) #1:2 = 18:22, 3:8 = 18:23, 9:20 = 18:22
temp <- colnames(run9)
temp <- temp[-1]
temp[1:2] <- substr(temp[1:2], 18, 22)
temp[3:8] <- substr(temp[3:8], 18, 23)
temp[9:20] <- substr(temp[9:20], 18, 22)
temp <- c("Gene",temp)
colnames(run9) <- temp
run9[1:5,1:5]

run10 <- readr::read_table(first.new.names[10],skip= 1)
run10[1:5,1:8]
run10 <- run10[,-2:-6] 
colnames(run10)
temp <- colnames(run10)
temp <- temp[-1]
temp <- substr(temp, 19, 23)
temp <- c("Gene",temp)
colnames(run10) <- temp
run10[1:5,1:5]

run11 <- readr::read_table(first.new.names[11],skip= 1)
run11[1:5,1:8]
run11 <- run11[,-2:-6] 
colnames(run11)
run11[1:5,1:8]
temp <- colnames(run11)
temp <- temp[-1]
temp[1:4] <- substr(temp[1:4], 19, 24)
temp[5:22] <- substr(temp[5:22], 19, 23)
temp <- c("Gene", temp)
colnames(run11) <- temp

run12 <- readr::read_table(first.new.names[12],skip= 1)
run12[1:5,1:8]
run12 <- run12[,-2:-6] 
colnames(run12)
run12[1:5,1:8]
temp <- colnames(run12)
temp <- temp[-1]
temp[1] <- substr(temp[1], 19, 23)
temp[2:24] <- substr(temp[2:24], 19, 24)
temp <- c("Gene", temp)
colnames(run12) <- temp



## before merging all files resolve duplicate samples
## run8 and run 6 duplicate sample : SR104
all.equal(run8$Gene, run6$Gene)
SR104 <- data.frame(run8$Gene, run8$SR104, run6$SR104)
SR104 <- column_to_rownames(SR104, "run8.Gene")
SR104$mean <- rowMeans(SR104)
## add mean to row 6 and remove from row 8
run6$SR104 <- SR104$mean
run8 <- run8 %>% dplyr::select(-SR104)

## run 7 and 5 
all.equal(run7$Gene, run5$Gene)
t.temp <- data.frame(run7$Gene, run7$SR130, run5$SR130 , run7$SR127, run5$SR127)
t.temp <- column_to_rownames(t.temp, "run7.Gene")
t.temp$mean <- rowMeans(t.temp[1:2])
t.temp$mean127 <- rowMeans(t.temp[3:4])
## add mean to row 5 and remove from row 7
run5$SR130 <- t.temp$mean
run5$SR127 <- t.temp$mean127
run7 <- run7 %>% dplyr::select(-SR130, -SR127)

f.newfiles <- list(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10, run11, run12)

#7,8
dat <- f.newfiles %>% purrr::reduce(inner_join)


f.batch.ids <- c(rep("run1", ncol(run1)-1), rep("run2", ncol(run2)-1), rep("run3", ncol(run3)-1),
                 rep("run4", ncol(run4)-1), rep("run5", ncol(run5)-1), rep("run6", ncol(run6)-1),
                 rep("run7", ncol(run7)-1), rep("run8", ncol(run8)-1), rep("run9", ncol(run9)-1),rep("run10",ncol(run10)-1),
                  rep("run11",ncol(run11)-1),rep("run12",ncol(run12)-1))

new.dat <- column_to_rownames(dat, "Gene") #make Gene rownames so that they are used as names for the batch ids
new.dat<- log2(new.dat + 1)

# - put batch.ids and gene expression matrix in same order - #
names(f.batch.ids) <- colnames(new.dat)
f.batch.ids <- f.batch.ids[colnames(new.dat)]

new.dat <- rownames_to_column(new.dat, "Gene")
new.dat[1:5,1:5]
dim(new.dat)

faheem <- new.dat #ready to be merged with others. (from batch_correction_pca.R)

## merge with the other dats from ensemble study 








DF1$Activity[DF2$NAME == DF1$NAME] <- DF2$Activity[DF2$NAME == DF1$NAME]





