import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
from scipy import stats

# 1. create two dict to map refseq id and probe id or vice versa
def get_gene_maps(infile):
    '''
    returns gene maps corresponding to different ids
    '''
    id_map = pd.read_csv(infile)
    id_map.dropna(inplace=True)
    refseq_to_gene = {}
    affy_to_gene = {}
    for index, row in id_map.iterrows():
        if row['REFSEQ'] not in refseq_to_gene:
            refseq_to_gene[row['REFSEQ']] = row['SYMBOL']
        if row['PROBEID'] not in affy_to_gene:
            affy_to_gene[row['PROBEID']] = row['SYMBOL']
    return refseq_to_gene, affy_to_gene
    
# 2. Compute concordance for the significantly DE genes
def concordance(refseq_genes, affy_genes, x=0):
    '''
    The concordance is adjusted by random chance prior to the computation.
    concordance = 2*intersect(ref_genes,affy_genes)/(refseq_genes + affy_genes)
    where intersect = ref_genes * affy_genes / (ref_gens + affy_genes)
    input:
    two lists of significant genes in each method
    x: background-corrected intersection. In two independent sets, it is 0
    output: a float 
    '''
    return 2*intersect(refseq_genes, affy_genes, x)/(len(refseq_genes)+len(affy_genes))

def intersect(refseq_genes, affy_genes, x):
    '''
    A helper function that calculates background-corrected intersection between 2 sets of DEgenes
    '''
    n0 = len(list(set(refseq_genes)&set(affy_genes))) # observed intersections
    n1 = len(refseq_genes)
    n2 = len(affy_genes)
    N = n1+n2-n0
    return x + (n1-x)*(n2-x)/(N-x)
    
def get_concordance(df_deseq, df_limma, refseq_to_gene, affy_to_gene):
    '''
    inputs: dataframes, gene mapping dicts
    return: concordance of a single treatment across two platforms, number of DEGs from both platforms
    '''
    # we are only interested in first column names
    refseq_genes = df_deseq.iloc[:, 0].tolist()
    affy_genes = df_limma.iloc[:, 0].tolist()
    # replace refseq acc/probe id to gene symbols
    refseq_genes = list(set([refseq_to_gene[x] for x in refseq_genes if x in refseq_to_gene]))
    affy_genes = list(set([affy_to_gene[x] for x in affy_genes if x in affy_to_gene]))
    # compute concordance
    return [concordance(refseq_genes, affy_genes), len(refseq_genes), len(affy_genes)]

def subdivide_concordance(df_deseq, df_limma, refseq_to_gene, affy_to_gene):
    '''
    for each treatment, get the concordnaces of above-median genes and below-median genes
    '''
    # subset each dataframe into above-median and blow-median dataframes
    df_deseq_above = df_deseq.loc[df_deseq['baseMean'] >= df_deseq['baseMean'].median()]
    df_deseq_below = df_deseq.loc[df_deseq['baseMean'] < df_deseq['baseMean'].median()]
    df_limma_above = df_limma.loc[df_limma['AveExpr'] >= df_limma['AveExpr'].median()]
    df_limma_below = df_limma.loc[df_limma['AveExpr'] < df_limma['AveExpr'].median()]
    
    c_above, deseq_genes_above, limma_genes_above = get_concordance(df_deseq_above, df_limma_above, refseq_to_gene, affy_to_gene)
    c_below, deseq_genes_below, limma_genes_below = get_concordance(df_deseq_below, df_limma_below, refseq_to_gene, affy_to_gene)
    return [c_above, deseq_genes_above, limma_genes_above],[c_below, deseq_genes_below, limma_genes_below]

def label_point(x, y, val, ax):
    '''
    This function labels points on a plot of matplotlib
    '''
    for a, b, label in zip(x,y,val):
        ax.text(a-120, b-5, label[:3])

def plot_concordance(overall_concordance, num_deseq, num_limma, ylabel, xlabel1, xlabel2, fname):
    fig, axes = plt.subplots(1,2,sharey=True, figsize=(25, 16))
    sns.regplot(x=num_deseq, y=overall_concordance, ax=axes[0], marker = 'x', scatter_kws={"s": 80},line_kws={'label':"Linear Reg"})
    axes[0].lines[0].set_linestyle("--")
    sns.regplot(x=num_limma, y=overall_concordance, ax=axes[1], marker = 'x', scatter_kws={"s": 80},line_kws={'label':"Linear Reg"})
    axes[1].lines[0].set_linestyle("--")
    axes[0].set_ylim(0,100)
    axes[0].yaxis.set_major_formatter(mtick.PercentFormatter())
    axes[0].set_xlim(min(num_deseq)-150,max(num_deseq)+150)
    axes[1].set_xlim(min(num_limma)-150,max(num_limma)+150)
    axes[0].legend()
    leg = axes[0].get_legend()
    L_labels = leg.get_texts()
    slope, intercept, r_value, p_value, std_err = stats.linregress(overall_concordance,num_deseq)
    label_line = r'$R^2:{0:.2f}$'.format(r_value*r_value)
    L_labels[0].set_text(label_line)
    axes[1].legend()
    leg = axes[1].get_legend()
    L_labels = leg.get_texts()
    slope, intercept, r_value, p_value, std_err = stats.linregress(overall_concordance,num_limma)
    label_line = r'$R^2:{0:.2f}$'.format(r_value*r_value)
    L_labels[0].set_text(label_line)
    label_point(num_deseq, overall_concordance, chemicals, axes[0])
    label_point(num_limma, overall_concordance, chemicals, axes[1])
    axes[0].set_ylabel(ylabel)
    axes[0].set_xlabel(xlabel1)
    axes[1].set_xlabel(xlabel2)
    #plt.show()
    sns.set(font_scale=3)
    plt.savefig(fname, dpi=100)
    
if __name__ == '__main__':
    #get the concordance and plot for each treatment
    refseq_to_gene, affy_to_gene = get_gene_maps('refseq_affy_map.csv')
    overall_concordance = []
    num_deseq = []
    num_limma = []
    data = []
    chemicals = ['BETA-NAPHTHOFLAVONE','ECONAZOLE','THIOACETAMIDE']
    for chemical in chemicals:
        df_deseq = pd.read_csv(f'{chemical}_deseq_results.csv')
        df_limma = pd.read_csv(f'{chemical}_limma_results.csv')
        c = get_concordance(df_deseq, df_limma, refseq_to_gene, affy_to_gene)
        print(f"for treatment {chemical}:")
        print(f"overall concordance: {c[0]}")
        print(f"# of DEGs in RNA-seq: {c[1]}")
        print(f"# of DEGs in microarray: {c[2]}")
        overall_concordance.append(c[0]*100)
        num_deseq.append(c[1])
        num_limma.append(c[2])
        data.append([chemical, 'overall', c[0]*100])
        
        c1,c2 = subdivide_concordance(df_deseq, df_limma, refseq_to_gene, affy_to_gene)
        print(f"above-median concordance: {c1[0]}")
        print(f"# of DEGs in RNA-seq: {c1[1]}")
        print(f"# of DEGs in microarray: {c1[2]}")
        data.append([chemical, 'above-median', c1[0]*100])
        
        print(f"below-median concordance: {c2[0]}")
        print(f"# of DEGs in RNA-seq: {c2[1]}")
        print(f"# of DEGs in microarray: {c2[2]}")
        data.append([chemical, 'below-median', c2[0]*100])
        print("\n")

    plot_concordance(overall_concordance, 
                     num_deseq, num_limma, 
                     "Concordance of DEG", 
                     "Treatment Effect\n(number of DEGs from RNA-seq)", 
                     "Treatment Effect\n(number of DEGs from Microarray)",
                     '6a.png')

    df = pd.DataFrame(data, columns = ['Treatment', 'Expression Level', 'Concordance']) 
    fig, ax = plt.subplots(figsize=(25, 16))
    splot = sns.barplot(x='Treatment', y='Concordance', data=df, hue='Expression Level', ax=ax)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    for p in splot.patches:
        splot.annotate(format(p.get_height(), '.1f'), 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'center', 
                       xytext = (0, 9), 
                       textcoords = 'offset points')
    plt.savefig('6b.png', dpi=100)
    
    # get top 10 expression genes from limma results
    affy_genes = df_limma.iloc[:, 0].tolist()
    affy_genes = list(set([affy_to_gene[x] for x in affy_genes if x in affy_to_gene]))
    outfile = open(f'top10_{chemical}_limma.txt','w')
    for gene in affy_genes[:10]:
        outfile.write(gene+'\n')
    outfile.close()
