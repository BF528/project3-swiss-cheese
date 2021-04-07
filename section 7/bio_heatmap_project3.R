ECONAZOLE<-read.csv('/projectnb/bf528/users/swiss_cheese2/project_3/programmer/DESeq_results/ECONAZOLE_deseq_norm_counts.csv')


Beta_Nap<-read.csv('/projectnb/bf528/users/swiss_cheese2/project_3/programmer/DESeq_results/Beta_Nap_deseq_norm_counts.csv')


THIOACETAMIDE<-read.csv('/projectnb/bf528/users/swiss_cheese2/project_3/programmer/DESeq_results/THIOACETAMIDE_deseq_norm_counts.csv')

ECON<-as.matrix(ECONAZOLE[1:11329,2:7])
BETA<-as.matrix(Beta_Nap[1:11329,2:7])
THIO<-as.matrix(THIOACETAMIDE[1:11329,2:7])
JOIN<-cbind(ECON,BETA,THIO)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(pheatmap)


colors = brewer.pal(n = 11, name = "GnBu")
colors = colorRampPalette(colors)(50)
colors = rev(colors)

heatmap <- pheatmap(JOIN, scale = "row", color = colors,fontsize_row = 4,border_color = NA, clustering_distance_rows="euclidean",
                    clustering_distance_cols="euclidean", main = "Clustered MOA ")


e<-c()
ei<-c()
b<-c()
bi<-c()
t<-c()
ti<-c()

ee<-ECONAZOLE

for (i in 2:nrow(ee)){
  if(sd(ee[i,2:7])/(sum(ee[i,2:7])/6)>0.186){
    ee[-c(i),]
    e<-c(e,sd(ee[i,2:7])/(sum(ee[i,2:7])/6))
    ei<-c(ei,i)
  }
}


bb<-Beta_Nap

for (i in 2:nrow(bb)){
  if(sd(bb[i,2:7])/(sum(bb[i,2:7])/6)>0.186){
    bb[-c(i),]
    b<-c(b,sd(ee[i,2:7])/(sum(ee[i,2:7])/6))
  }
}


tt<-THIOACETAMIDE
for (i in 2:nrow(tt)){
  if(sd(tt[i,2:7])/(sum(tt[i,2:7])/6)>0.186){
    tt[-c(i),]
    t<-c(t,sd(ee[i,2:7])/(sum(ee[i,2:7])/6))
  }
}