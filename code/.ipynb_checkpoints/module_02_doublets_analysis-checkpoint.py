import pandas as pd
import numpy as np
import scipy
import seaborn as sns
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import warnings
import scanpy as sc
import sklearn
from adjustText import adjust_text
import anndata as ad

warnings.filterwarnings("ignore")

def get_CTprops_barplot_df(CTprobs, clusts):
    ## barplot of cell type proportions per cluster
    tmp=CTprobs.copy()
    tmp['clusters']=clusts
    tmp=tmp[tmp.iloc[:,0:len(tmp.columns)-1].sum(axis=1)>0]
    tmp=tmp.groupby(['clusters']).sum().T/tmp.clusters.value_counts()
    tmp=tmp.T

    biggestVals_clust=pd.DataFrame([np.sort(tmp.iloc[i,:])[-1] for i in range(len(tmp.index))], columns=['Value'])
    #biggestVals['order']=1
    secondVals_clust=pd.DataFrame([np.sort(tmp.iloc[i,:])[-2] for i in range(len(tmp.index))], columns=['Value'])
    #secondVals['order']=2
    biggestVals_clust=pd.concat([biggestVals_clust, secondVals_clust], axis=1)
    biggestVals_clust.columns=['firstVal', 'secondVal']
    biggestVals_clust.index=tmp.index
    tmp=tmp.T

    for i in tmp.columns:
        tmp[i][(tmp[i]!=biggestVals_clust.loc[i][0])&(tmp[i]!=biggestVals_clust.loc[i][1])]=0

    other=1-tmp.sum()
    tmp=pd.concat([tmp.T, other], axis=1).T
    tmp.rename(index={0:'other'},inplace=True)
    return tmp


def doubletsDEA_niche(adata, niche_col, niche_name, condition_col, disease_cond_name, fdr=0.05, logFC=1.5, plot=True, title=None):
    ## differential expression per cluster for volcano plots

    tmp=adata[adata.obs[niche_col]==niche_name]
    
    tmp=tmp[:,(pd.DataFrame(tmp.to_df(), index=tmp.obs_names)[tmp.obs[condition_col]!=disease_cond_name].mean()>0)&(pd.DataFrame(tmp.to_df(), index=tmp.obs_names)[tmp.obs[condition_col]==disease_cond_name].mean()>0)]
    sc.pp.filter_genes(tmp, min_cells = 10)
    
    sc.tl.rank_genes_groups(tmp, condition_col, method='wilcoxon')
    sc.pl.rank_genes_groups(tmp, n_genes=25, sharey=False)
    t=tmp.uns['rank_genes_groups']
    DE_DF=pd.concat([pd.DataFrame(t['names'])[disease_cond_name], pd.DataFrame(t['pvals'])[disease_cond_name], pd.DataFrame(t['logfoldchanges'])[disease_cond_name]], axis=1)
    DE_DF.columns=['gene_name', 'p-value', 'logFC']
    ###just positive pvals
    DE_DF=DE_DF[DE_DF['p-value']>0]
    
    
    if plot == True:
        DE_DF["-logQ"] = -np.log10(DE_DF["p-value"].astype("float"))
        lowqval_de = DE_DF.loc[abs(DE_DF["logFC"]) > logFC]
        other_de = DE_DF.loc[abs(DE_DF["logFC"]) <= logFC]
        
        fig, ax = plt.subplots()
        sns.regplot(
            x=other_de["logFC"],
            y=other_de["-logQ"],
            fit_reg=False,
            #scatter_kws={"s": 6},
        )
        sns.regplot(
            x=lowqval_de["logFC"],
            y=lowqval_de["-logQ"],
            fit_reg=False,
            #scatter_kws={"s": 6},
        )
        ax.set_xlabel("log2 FC")
        ax.set_ylabel("-log Q-value")
            
        txt=[ax.text(DE_DF.sort_values('-logQ', ascending=False, ignore_index=True)['logFC'][i], 
                    DE_DF.sort_values('-logQ', ascending=False, ignore_index=True)['-logQ'][i], 
                    DE_DF.sort_values('-logQ', ascending=False, ignore_index=True)['gene_name'][i]) for i in range(20)]
        adjust_text(txt)
        plt.title(title)
    return DE_DF


def getViolinData(adata, adata_imp, g, annot_col, EXP_stage, WT_stage, savename='test'):
    ## Get data for violin plots
    gene_imp=pd.DataFrame(adata_imp.X, index=adata_imp.obs_names, columns=adata_imp.var_names)[g]
    gene=pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)[g]
    
    
    d=pd.DataFrame({'gene':gene, 'gene_imp':gene_imp, 'stage':adata_imp.obs.stage, 'annotation': adata_imp.obs[annot_col].astype(str)})
    
    p=pd.Series([scipy.stats.ranksums(d.gene[(d.stage==EXP_stage)&(d.annotation==cell)], 
                                      d.gene[(d.stage==WT_stage)&(d.annotation==cell)]).pvalue for cell in d.annotation.unique()], index=d.annotation.unique())

    FC=pd.Series([(0.000001+d.gene_imp[(d.annotation==c) & (d.stage==EXP_stage)].mean())/(0.000001+d.gene_imp[(d.annotation==c) & (d.stage==WT_stage)].mean()) for c in p.index], index=p.index)
    st=pd.concat([p, FC], axis=1)
    st.columns=['p-value', 'FC']
    st.sort_values('p-value')

    d.to_csv(savename+'.csv')
    return d