### Coloc DF

import sys
import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad
import nichesphere
from sklearn.neighbors import NearestNeighbors

data=sc.read('../data/CODEX/final_data/CODEX_run45_anndata_final.h5ad')

for rad in [600]:
    ## In for loop
    CTcolocalizationP_=pd.DataFrame()
    
    for smpl in list(data.obs['sample'].unique()):
        print(smpl)
        anndata=data[data.obs['sample']==smpl].copy()
        colocHM=nichesphere.coloc.compute_sc_spatial_radius_coloc_matrix(adata=anndata, cluster_col='clusters_named', radius=rad)
        testDistribution_df=nichesphere.tl.get_spatial_radius_BG_OEratios_DF(anndata, cluster_col='clusters_named', radius=rad, n_permutations=1000)
        testDistribution_df.to_csv('../res_data/revision/CODEX/'+smpl+'_'+str(rad)+'_testDistribution_df_final.csv')
        colocHM['sample']=smpl
        CTcolocalizationP_=pd.concat([CTcolocalizationP_, colocHM])
    
    CTcolocalizationP_.to_csv('../res_data/revision/CODEX/CODEX_run45_CTcolocalizationP_radius'+str(rad)+'_final.csv')



















