# Part 3 R Code: Read counting with Feature Counts, combining counts into one matrix  

library("dplyr")
library("ggplot2")
library("DESeq2")
library("reshape2")



#load in featureCount samples from SCC to create csv file // featureCounts directory
SRR1177998 <- read.table("featureCounts_outputSRR1177998.txt")
SRR1178001 <- read.table("featureCounts_outputSRR1178001.txt")
SRR1178003 <- read.table("featureCounts_outputSRR1178003.txt")
SRR1177993 <- read.table("featureCounts_outputSRR1177993.txt")
SRR1177994 <- read.table("featureCounts_outputSRR1177994.txt")
SRR1177995 <- read.table("featureCounts_outputSRR1177995.txt")
SRR1177966 <- read.table("featureCounts_outputSRR1177966.txt")
SRR1177969 <- read.table("featureCounts_outputSRR1177969.txt")
SRR1177970 <- read.table("featureCounts_outputSRR1177970.txt")

#making the first column the gene ID
genes <- SRR1177970[, 1]

#combining data from all featureCounts into one data frame
count_matrix <- data.frame(Geneid = genes, SRR1177998[, 7], SRR1178001 [, 7], SRR1178003[, 7],
                       SRR1177993[, 7], SRR1177994[, 7], SRR1177995[, 7], SRR1177966[, 7], SRR1177969[, 7], SRR1177970[, 7])
colnames(count_matrix) <- c("Geneid", "SRR1177998", "SRR1178001","SRR1178003",
                        "SRR1177993", "SRR1177994", "SRR1177995", "SRR1177966", "SRR1177969","SRR1177970")
                        
count_matrix = count_matrix[-c(1),]
write.csv(count_matrix, "count_matrix.csv")



#create box plot

graphmatrix = count_matrix
graphmatrix[graphmatrix==0] = 1
boxplot(as.integer(paste(graphmatrix$SRR1177998)), 
        as.integer(paste(graphmatrix$SRR1178001)), 
        as.integer(paste(graphmatrix$SRR1178003)), 
        as.integer(paste(graphmatrix$SRR1177993)),
        as.integer(paste(graphmatrix$SRR1177994)),
        as.integer(paste(graphmatrix$SRR1177995)),
        as.integer(paste(graphmatrix$SRR1177966)),
        as.integer(paste(graphmatrix$SRR1177969)),
        as.integer(paste(graphmatrix$SRR1177970)),
        main = "feature counts for SRR1177998", names = names(count_matrix)[-1], log = "y", xlab = "Samples", ylab = "Feature Counts")






