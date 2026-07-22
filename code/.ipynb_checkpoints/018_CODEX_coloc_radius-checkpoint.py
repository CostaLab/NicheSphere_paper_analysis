### Coloc DF

import sys
import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.neighbors import NearestNeighbors

def compute_colocalization_matrix(adata, cluster_col, k=5):
    """
    Compute a colocalization matrix showing the number of kNN interactions between clusters.
    
    Parameters:
        adata (anndata): A spatial anndata obj with cell centroids as obsm['spatial'] and a cluster column in the obs.
        cluster_col (str): The name of the column containing cluster labels. (must be categorical)
        k (int): Number of nearest neighbors to consider.

    Returns:
        pd.DataFrame: A colocalization matrix where entry (i, j) represents 
                      the number of times cluster j appears in the kNN of cluster i.
    """
    # Extract centroids of polygons
    coords=adata.obsm['spatial']

    # Fit kNN model
    ### Try radius??
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(coords)
    _, indices = nbrs.kneighbors(coords)

    # Get cluster labels
    cluster_labels = adata.obs[cluster_col].values
    #unique_clusters = np.unique(cluster_labels)
    unique_clusters = adata.obs[cluster_col].cat.categories

    # Initialize colocalization matrix
    cluster_index = {c: i for i, c in enumerate(unique_clusters)}
    colocalization_matrix = np.zeros((len(unique_clusters), len(unique_clusters)), dtype=int)

    # Populate the matrix
    for i, neighbors in enumerate(indices):
        cluster_i = cluster_labels[i]
        for neighbor in neighbors[1:]:  # Exclude self (first neighbor)
            cluster_j = cluster_labels[neighbor]
            colocalization_matrix[cluster_index[cluster_i], cluster_index[cluster_j]] += 1
    W = colocalization_matrix
    W = W/W.sum() ## Prob to find an boundary
    # Convert to a Pandas DataFrame for readability
    colocalization_df = pd.DataFrame(W, index=unique_clusters, columns=unique_clusters)
    
    return colocalization_df


def compute_radius_colocalization_matrix(adata, cluster_col, radius=40.0):
    coords = adata.obsm['spatial']
    nbrs = NearestNeighbors(radius=radius, algorithm='kd_tree').fit(coords)
    indices = nbrs.radius_neighbors(coords, return_distance=False)

    cluster_labels = adata.obs[cluster_col].values
    categories = adata.obs[cluster_col].cat.categories
    cluster_to_idx = {cat: i for i, cat in enumerate(categories)}
    
    n = len(categories)
    mat = np.zeros((n, n))

    for i, neighbors in enumerate(indices):
        # map current cell label to index
        row = cluster_to_idx[cluster_labels[i]]
        for nb in neighbors:
            if i == nb: continue # ignore self
            col = cluster_to_idx[cluster_labels[nb]]
            mat[row, col] += 1

    # Symmetrize and Normalize
    mat = (mat + mat.T) / 2
    total = mat.sum()
    if total > 0:
        mat = mat / total

    return pd.DataFrame(mat, index=categories, columns=categories)


def generate_spatial_background(adata, cluster_col, radius=40.0, n_permutations=1000):
    """
    Generates a distribution of O/E ratios by shuffling labels 
    while keeping the spatial graph fixed.
    """
    coords = adata.obsm['spatial']
    # 1. Pre-calculate the neighbors to save time during shuffling
    nbrs = NearestNeighbors(radius=radius, algorithm='kd_tree').fit(coords)
    adj_list = nbrs.radius_neighbors(coords, return_distance=False)

    categories = adata.obs[cluster_col].cat.categories
    n_clusters = len(categories)
    cluster_to_idx = {cat: i for i, cat in enumerate(categories)}
    
    # Store results for each permutation
    null_distributions = []

    print(f"Generating {n_permutations} permutations...")
    for _ in range(n_permutations):
        # 2. Shuffle the labels
        shuffled_labels = adata.obs[cluster_col].sample(frac=1).values
        
        # 3. Compute interactions on the shuffled graph
        mat = np.zeros((n_clusters, n_clusters))
        for i, neighbors in enumerate(adj_list):
            row = cluster_to_idx[shuffled_labels[i]]
            for nb in neighbors:
                if i == nb: continue
                col = cluster_to_idx[shuffled_labels[nb]]
                mat[row, col] += 1
        
        # Symmetrize and normalize as in the observed calculation
        mat = (mat + mat.T) / 2
        if mat.sum() > 0:
            mat = mat / mat.sum()
            
        # Flatten to a Series where index is 'CT1_CT2' to match your OvsE function
        flat_probs = []
        for i, ct1 in enumerate(categories):
            for j, ct2 in enumerate(categories):
                flat_probs.append(mat[i, j])
        
        null_distributions.append(flat_probs)

    # Create the testDistribution_df (rows = permutations, columns = cell type pairs)
    colnames = [f"{c1}-{c2}" for c1 in categories for c2 in categories]
    return pd.DataFrame(null_distributions, columns=colnames)

data=sc.read('../data/CODEX/final_data/CODEX_run45_anndata_final.h5ad')

for rad in [600]:
    ## In for loop
    CTcolocalizationP_=pd.DataFrame()
    
    for smpl in list(data.obs['sample'].unique()):
        print(smpl)
        anndata=data[data.obs['sample']==smpl].copy()
        colocHM=compute_radius_colocalization_matrix(adata=anndata, cluster_col='clusters_named', radius=rad)
        testDistribution_df=generate_spatial_background(anndata, cluster_col='clusters_named', radius=rad, n_permutations=1000)
        testDistribution_df.to_csv('../res_data/revision/CODEX/'+smpl+'_'+str(rad)+'_testDistribution_df_final.csv')
        colocHM['sample']=smpl
        CTcolocalizationP_=pd.concat([CTcolocalizationP_, colocHM])
    
    CTcolocalizationP_.to_csv('../res_data/revision/CODEX/CODEX_run45_CTcolocalizationP_radius'+str(rad)+'_final.csv')



















