'''
This script is modified by the Analyst Yichi Zhang to run limma analysis.
'''
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages('plotly')
install.packages('ggrepel')
library(glue)
library(limma)
library(plotly)
library(ggrepel)

setwd("/home/eikthedragonslayer/Desktop/BF528/project03")
# sample info dataframe with array_id and chemical columns
samples <- read.csv('group_2_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

# create a function for repeat calls
visualization <- function(chemical, samples){ # chemical is the name of the treatment
  # subset the full expression matrix to just those in this comparison
  rma.subset <- rma[paste0('X',samples[samples$chemical == 'Control' | samples$chemical == chemical,]$array_id)]
  
  # construct a design matrix modeling treatment vs control for use by limma
  design <- model.matrix(
    ~factor(
      samples$chemical,
      levels=c('Control',chemical)
    )
  )
  colnames(design) <- c('Intercept',chemical)
  
  # run limma
  fit <- lmFit(rma.subset, design)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH', sort.by = 'p') # sort results by adj-p-values and select at p-adjust < 0.05
  t <- t[t[,5] < 0.05,]
  top10 <- t[1:10,]
  
  # write out the results to file
  write.csv(t,glue('{chemical}_limma_results.csv'))
  write.csv(top10,glue('top10_{chemical}_limma_results.csv'))
  
  # plot a scatter plot of fold change values and scatter plots of fold change vs nominal-pvalue from the significant DE genes and save
  indexes = row.names(t)
  plot <- t %>% 
    ggplot(aes(size = AveExpr, x = adj.P.Val, y = logFC, color = B)) + 
    geom_point(aes(text = indexes)) +
    geom_text_repel(
      mapping = aes(x = adj.P.Val, y = logFC, label = indexes),
      size=4, box.padding = unit(0.5, "lines")
    ) +
    theme_bw(base_size=20) + 
    geom_hline(yintercept = 0, linetype = "dotted", colour= 'red') +
    scale_color_continuous(low='skyblue', high='midnightblue') +
    labs(x=paste0("Nominal P-values"),
         y=glue("Log2 Fold Change (Control vs. {chemical})")) +
    theme(legend.position="none") 
  plot+theme(legend.position = "right") 
  #ggplotly(plot) %>% layout(showlegend=T)
  # save plot in tiff for high qualities
  ggsave(filename = glue('{chemical}_scatter.tiff'), width = 16, height = 9, device='tiff', dpi=300)
  
  # plot a histogram of FC values
  histo <- t %>%
    ggplot( aes(x = logFC)) +
    geom_histogram( binwidth=0.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("Bin size = 0.1") +
    geom_vline(aes(xintercept = 0),
               color="blue", linetype="dashed", size=1) +
    geom_density(alpha=.2, fill="#FF6666") +
    labs(x=glue("Log2 Fold Change (Control vs. {chemical})"),
         y=paste0("Frequency")) +
    theme(
      plot.title = element_text(size=20)
    )
  #histo
  ggsave(filename = glue('{chemical}_hist.tiff'), width = 16, height = 9, device='tiff', dpi=300)
}
for (chemical in c('ECONAZOLE','BETA-NAPHTHOFLAVONE','THIOACETAMIDE')) {
  visualization(chemical, samples)
}
