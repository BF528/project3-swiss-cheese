
#Programmer 
library("DESeq2")

BiocManager::install(c("apeglm"))
BiocManager::available()
BiocManager::valid()

#Part 4: RNA-Seq Differential Expression with DESeq2


#combined matrix from part3 
write.csv(count_matrix, "count_matrixnew.csv")

count_matrixnew = read.csv("count_matrixnew.csv", header = TRUE, row.names=1)

#edit control matrix to only have relevant samples 
control = read.csv("control_counts.csv", header = TRUE)
control = control[,-1]
control = read.csv('control_counts.csv',)[ ,c('Geneid', 'SRR1178030', 'SRR1178040', 'SRR1178056', 'SRR1178024', 'SRR1178035', 'SRR1178045', 'SRR1178004', 'SRR1178006', 'SRR1178013'),]



#Three different groups + controls 
AhR_samples = cbind(count_matrixnew[,c("SRR1177998", "SRR1178001", "SRR1178003")], control[,c("SRR1178030", "SRR1178040", "SRR1178056")])
CARPXR_samples = cbind(count_matrixnew[,c("SRR1177993", "SRR1177994", "SRR1177995")], control[,c("SRR1178024", "SRR1178035", "SRR1178045")])
Cytotoxic_samples = cbind(count_matrixnew[,c("SRR1177966", "SRR1177969", "SRR1177970")], control[,c("SRR1178004", "SRR1178006", "SRR1178013")])




#4.3 

# load counts for AhR
AhR_counts = AhR_samples

# filter out rows that have any zeros for funzies
AhR_counts = subset(AhR_counts,rowSums(AhR_counts==0)==0)

# sample information
info = read.csv('group_2_rna_info.csv', header = TRUE, row.names=1)
AhR_info = info[c("SRR1177998", "SRR1178001", "SRR1178003", "SRR1178030", "SRR1178040", "SRR1178056"), ]
# create the DESeq object
dds1 = DESeqDataSetFromMatrix(
  countData = AhR_counts,
  colData = AhR_info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds1$mode_of_action = relevel(dds1$mode_of_action, ref='Control')

# run DESeq
dds1 = DESeq(dds1)
res1 = results(dds1, contrast=c('mode_of_action','AhR','Control'))
res1 = lfcShrink(dds1, coef=2)
# write out DE results
write.csv(res1,'AhR_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds1,normalized=TRUE),'AhR_deseq_norm_counts.csv')



#For CARPXR

#load counts for CARPXR
CARPXR_counts = CARPXR_samples

# filter out rows that have any zeros
CARPXR_counts = subset(CARPXR_counts,rowSums(CARPXR_counts==0)==0)
CARPXR_info = info[c("SRR1177993", "SRR1177994", "SRR1177995", "SRR1178024", "SRR1178035", "SRR1178045"), ]
dds2 = DESeqDataSetFromMatrix(
  countData = CARPXR_counts,
  colData = CARPXR_info,
  design= ~ mode_of_action
)
dds2 = DESeq(dds2)
res2 = results(dds2, contrast=c('mode_of_action','CAR/PXR','Control'))
res2 = lfcShrink(dds2, coef=2)
write.csv(res2,'CARPXR_deseq_results.csv')
write.csv(counts(dds2,normalized=TRUE),'CARPXR_deseq_norm_counts.csv')


#For Cytotoxic
Cytotoxic_counts = Cytotoxic_samples
Cytotoxic_counts = subset(Cytotoxic_counts,rowSums(Cytotoxic_counts==0)==0)
Cytotoxic_info = info[c("SRR1177966", "SRR1177969", "SRR1177970", "SRR1178004", "SRR1178006", "SRR1178013"), ]
dds3 = DESeqDataSetFromMatrix(
  countData = Cytotoxic_counts,
  colData = Cytotoxic_info,
  design= ~ mode_of_action
)
dds3 = DESeq(dds3)
res3 = results(dds3, contrast=c('mode_of_action','Cytotoxic','Control'))
res3 = lfcShrink(dds3, coef=2)
write.csv(res3,'Cytotoxic_deseq_results.csv')
write.csv(counts(dds3,normalized=TRUE),'Cytotoxic_deseq_norm_counts.csv')



#4.4 
res1_sorted = res1[order(res1$padj),]
res1_sorted = res1_sorted[complete.cases(res1_sorted), ]
res2_sorted = res2[order(res2$padj),]
res2_sorted = res2_sorted[complete.cases(res2_sorted), ]
res3_sorted = res3[order(res3$padj),]
res3_sorted = res3_sorted[complete.cases(res3_sorted), ]

AhR_sig_genecounts = dim(res1_sorted[res1_sorted$padj<0.05,])[1]
AhR_top10 = row.names(res1_sorted)[1:10]
CARPXR_sig_genecounts = dim(res2_sorted[res2_sorted$padj<0.05,])[1]
CARPXR_top10 = row.names(res2_sorted)[1:10]
Cytotoxic_sig_genecounts = dim(res3_sorted[res3_sorted$padj<0.05,])[1]
Cytotoxic_top10 = row.names(res3_sorted)[1:10]



#4.5 
Result = data.frame(col1=c(AhR_sig_genecounts, AhR_top10), col2=c(CARPXR_sig_genecounts, CARPXR_top10), col3=c(Cytotoxic_sig_genecounts, Cytotoxic_top10))
names(Result) = c("AhR", "CARPXR", "Cytotoxic")
rownames(Result) = c("Significant genes count", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th")
write.csv(Result, "Toxgroup significant gene count and top 10 most significant genes.csv")


#4.6 Histogram and Scatter Plots: (group2)
