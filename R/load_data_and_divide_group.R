library(dplyr)
library(tidyr)
print("This code is wrote by Lei Yu to do lasso regression on micro RNAs vs Phenotype")

# Read the micro RNA dataset into program
print("Reading the micro RNA data")
mir <- readRDS("./dataset/mir_demo.rds")
# Read the  MRNA dataset into program
print("Reading the mRNA data")
mrna <- readRDS("./dataset/mrna_demo.rds")
# Read the Methyl dataset into program
print("Reading the Methyl data")
Methyl <- readRDS("./dataset/methyl_demo.rds")
print("Reading the Gleason score")
clinicalDa <- readRDS("./dataset/clinicalDa.rds")
clinicalDa <- as.matrix(clinicalDa)

# Because we only use the cancer cells transportome to make the prediction, therefore, we should add 
# "-01" to clinicalDa's row name

row.names(clinicalDa) <- paste(row.names(clinicalDa), '-01', sep="")
gleason_score <- clinicalDa[,"gleason_score"]
gleason_score <- as.matrix(gleason_score)
gleason_score <- apply(gleason_score, 2, as.numeric)
rownames(gleason_score) <- rownames(clinicalDa)

print("Data input finished")
print("Patients number in clinical Data, mir RNA, mRNA and Methyl are:")
print(paste(nrow(gleason_score), nrow(mir), nrow(mrna), nrow(Methyl),sep = ', '))
l <- list(rownames(gleason_score), rownames(mir), rownames(mrna), rownames(Methyl))
Patients <- Reduce(intersect,l)

gleason_score <- gleason_score[Patients,]
gleason_score <- as.matrix(gleason_score)

Methyl <- Methyl[Patients,]
mir <- mir[Patients,]
mrna <- mrna[Patients,]
number_of_groups <- 10
if(!file.exists('./output/Group.rds')){

number_of_groups <- 10
print("Split Gleason score into groups")
name_of_six <- rownames(gleason_score)[(which(gleason_score == 6))]
name_of_seven <- rownames(gleason_score)[(which(gleason_score == 7))]
name_of_eight<- rownames(gleason_score)[(which(gleason_score == 8))]
name_of_nine <- rownames(gleason_score)[(which(gleason_score == 9))]
name_of_ten <- rownames(gleason_score)[(which(gleason_score == 10))]

print("Constructing groups")
Group <- list()
Group_six <- split(name_of_six, 1:number_of_groups)
Group_seven <- split(name_of_seven, 1:number_of_groups)
Group_eight <- split(name_of_eight, 1:number_of_groups)
Group_nine <- split(name_of_nine, 1:number_of_groups)
Group_ten <- split(name_of_ten, 1:number_of_groups)

for (i in 1:number_of_groups){
    Group[[i]] <- c(Group_six[[i]], Group_seven[[i]], Group_eight[[i]], Group_nine[[i]], Group_ten[[i]])
  }
saveRDS(Group, './output/Group.rds')
}
Group <- readRDS('./output/Group.rds')
