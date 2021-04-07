#Volcano Plots 

# add a column of NAs
carpxr_deseq$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
carpxr_deseq$diffexpressed[carpxr_deseq$log2FoldChange > 0.6 & carpxr_deseq$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
carpxr_deseq$diffexpressed[carpxr_deseq$log2FoldChange < -0.6 & carpxr_deseq$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=carpxr_deseq, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


carpxr_deseq$carpxr_deseqlabel <- NA
carpxr_deseq$carpxr_deseqlabel[carpxr_deseq$diffexpressed != "NO"] <- carpxr_deseq$gene_symbol[carpxr_deseq$diffexpressed != "NO"]



library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=carpxr_deseq, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=carpxr_deseqlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")







ggplot(data=Cytotoxic_deseq, aes(x=log2FoldChange, y=pvalue)) + geom_point()
p <- ggplot(data=Cytotoxic_deseq, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()+ theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# add a column of NAs
Cytotoxic_deseq$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Cytotoxic_deseq$diffexpressed[Cytotoxic_deseq$log2FoldChange > 0.6 & Cytotoxic_deseq$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Cytotoxic_deseq$diffexpressed[Cytotoxic_deseq$log2FoldChange < -0.6 & Cytotoxic_deseq$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=Cytotoxic_deseq, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


Cytotoxic_deseq$Cytotoxic_deseqlabel <- NA
Cytotoxic_deseq$Cytotoxic_deseqlabel[Cytotoxic_deseq$diffexpressed != "NO"] <- Cytotoxic_deseq$gene_symbol[Cytotoxic_deseq$diffexpressed != "NO"]




library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=Cytotoxic_deseq, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=Cytotoxic_deseqlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

