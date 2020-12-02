## Packages ##
library(dplyr)
library(tidyr)
library(glmnet)
library(DT)
## Step 1: Set working environment ##
# dir <- The/Path/to/TLHLASSO/
# setwd(dir)


## Step 2: Load datasets and divide samples into 10 groups ##
source("./R/load_data_and_divide_group.R")


## Step 3: predict mRNAs and microRNAs based only on methylatomic markers ##

source("./R/cross_lasso.R")
source("./R/R_square.R")
###### The prediction micro RNA ##########################
predicted_mir <- c()
for (i in 1:ncol(mir)){
  sub <- cross_lasso(as.matrix(Methyl), as.matrix(as.matrix(mir)[,i]), Group)
  predicted_mir <- cbind(predicted_mir,sub)
}
colnames(predicted_mir) <- colnames(mir)
###### View cross_lasso predicted micro RNA ########################
datatable(predicted_mir)

###### The prediction mRNA ##########################
predicted_mrna <- c()
for (j in 1:ncol(mrna)){
  sub <- cross_lasso(as.matrix(Methyl), as.matrix(as.matrix(mrna)[,j]), Group)
  predicted_mrna <- cbind(predicted_mrna,sub)
}
colnames(predicted_mrna) <- colnames(mrna)
###### View cross_lasso predicted  mRNA ########################
datatable(predicted_mrna)

## Layer 1 selection ###
## Set threshold of predicted_mrna to do selection and do second regression based on selected predicted mrna to gleason score
#### Get the R square value of each predicted mRNA
R_square_vector <- c()
for(i in 1:ncol(mrna)){
  R_square_vector[i] <- R_square(as.matrix(predicted_mrna[,i]), as.matrix(as.matrix(mrna)[,i]))
}
R_square_vector_matrix <- as.matrix(R_square_vector)
rownames(R_square_vector_matrix) <- colnames(mrna)

## Select R square and predict gleason

L1_P_value = c()
L1_R_square = c()
input_pred_mRNA_number <- c()

for (i in seq(0, 0.95, by=0.05)) {
  submRNA <- predicted_mrna[, rownames(R_square_vector_matrix)[which(R_square_vector_matrix>i)]]
  if(sum(R_square_vector_matrix>i)!=0){
    pred <- cross_lasso(submRNA, gleason_score, Group)
    R2 <- R_square(pred, gleason_score)
    L1_P_value <- c(L1_P_value,i)
    L1_R_square <- c(L1_R_square,R2)
    input_pred_mRNA_number <- c(input_pred_mRNA_number,ncol(submRNA))
  }
}

Layer_1_result <- data.frame(L1_P_value, L1_R_square, input_pred_mRNA_number)
print(paste0("The best R square is ",max(Layer_1_result$L1_R_square), " when P = ", Layer_1_result$L1_P_value[which.max(Layer_1_result$L1_R_square)]))
L1_mrna_names <- rownames(R_square_vector_matrix)[which(R_square_vector_matrix>max(Layer_1_result$L1_R_square))]
UL1_mrna_names <- rownames(R_square_vector_matrix)[which(R_square_vector_matrix<max(Layer_1_result$L1_R_square))]


## Step 4: combine methylation and estimated mirco RNA profile to improve predictability of mRNAs in UL1##
combined_predictor <- merge(Methyl, predicted_mir, by = 'row.names')
rownames(combined_predictor) <- combined_predictor$Row.names
combined_predictor <- combined_predictor[,-1]
predicted_UL1 <- c()
for (j in UL1_mrna_names){
  sub <- cross_lasso(as.matrix(combined_predictor), as.matrix(mrna[,j]), Group)
  predicted_UL1 <- cbind(predicted_UL1,sub)
}
colnames(predicted_UL1) <- UL1_mrna_names
  
R_square_vector_UL1 <- c()
for(i in UL1_mrna_names){
  R_square_vector_UL1[i] <- R_square(as.matrix(predicted_UL1[,i]), as.matrix(mrna[,i]))
}

R_square_vector_UL1 <- as.matrix(R_square_vector_UL1)
R_square_vector_UL1_original <- as.matrix(R_square_vector_matrix[UL1_mrna_names,])
UL1_compare <- cbind(R_square_vector_UL1_original, R_square_vector_UL1)
UL1_compare <- cbind(UL1_compare, UL1_compare[,2] - UL1_compare[,1])

L2_candadiate_names <-rownames(UL1_compare)[UL1_compare[,3]>0]
R_square_vector_UL1_L2_candadiate <- as.matrix(R_square_vector_UL1[L2_candadiate_names, ])
## Tuning parameter Q ##
L2_Q_value <- c()
L2_R_square <- c()
input_pred_mRNA_number <- c()
for (i in seq(0, 0.95, by=0.05)) {
  submRNAL2 <- predicted_UL1[, rownames(R_square_vector_UL1_L2_candadiate)[which(R_square_vector_UL1_L2_candadiate>i)]]
  submRNAL2_n_L1 <- merge(submRNAL2, predicted_mrna[,L1_mrna_names], by = 'row.names') 
  rownames(submRNAL2_n_L1) <- submRNAL2_n_L1[,1]
  submRNAL2_n_L1 <-  submRNAL2_n_L1[, -1]
  if(sum(R_square_vector_UL1_L2_candadiate>i)!=0){
    pred <- cross_lasso(as.matrix(submRNAL2_n_L1), gleason_score, Group)
    R2 <- R_square(pred, gleason_score)
    L2_Q_value <- c(L2_Q_value,i)
    L2_R_square <- c(L2_R_square,R2)
    input_pred_mRNA_number <- c(input_pred_mRNA_number,ncol(submRNAL2_n_L1))
  }
}
Layer_2_result <- data.frame(L2_Q_value, L2_R_square, input_pred_mRNA_number)
print(paste0("The best R square is ",max(Layer_2_result$L2_R_square), " when Q = ", Layer_2_result$L2_Q_value[which.max(Layer_2_result$L2_R_square)]))
L2_mrna_names <- rownames(R_square_vector_UL1_L2_candadiate)[which(R_square_vector_UL1_L2_candadiate>max(Layer_1_result$L1_R_square))]

saveRDS(L1_mrna_names,'./output/L1_mrna_ids.rds')
saveRDS(L2_mrna_names,'./output/L2_mrna_ids.rds')
