library(Seurat)
library(sctransform)
library(harmony)
library(tidyverse)
library(hdf5r)
library(patchwork)

library(cowplot)
library(RColorBrewer)
library(viridis)
library(Matrix)

options(future.globals.maxSize = 5000 * 1024^2)

## Folders
fig_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run46_Bella/visual_output/"
out_data_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run46_Bella/"

## resolution to use
chosen_res='SCT_snn_res.0.3'

## cluster color palette
cluster_palette="Paired"

#######################

cur_seurat<-readRDS(paste0(out_data_dir, 'seurat_object_harmony_clusts_filt_bg.rds'))
cur_seurat@meta.data$seurat_clusters<-cur_seurat@meta.data[,chosen_res]	
Idents(object = cur_seurat) <- "seurat_clusters"
# Classic palette BuPu, with 4 colors
newpal <- brewer.pal(12, cluster_palette) 
# Add more colors to this palette :
newpal <- colorRampPalette(newpal)(20)

pdf(paste0(fig_dir, "umap_clusters_sct_harmony_filt_bg.pdf"), width=7, height=7)
p1 <- DimPlot(cur_seurat, group.by='seurat_clusters', reduction='umap', cols = newpal, raster = FALSE) +
  ggtitle('annotation')
print(p1)
dev.off()

pdf(paste0(fig_dir, "umap_sample_sct_harmony_filt_bg.pdf"), width=7, height=7)
p3 <- DimPlot(cur_seurat, group.by='sample', reduction='umap', raster = FALSE) +
  ggtitle('Sample')
print(p3)
dev.off()

pdf(paste0(fig_dir, "umap_run_sct_harmony_filt_bg.pdf"), width=7, height=7)
p4 <- DimPlot(cur_seurat, group.by='run', reduction='umap', raster = FALSE) +
  ggtitle('Run')
print(p4)
dev.off()

#cluster markers with Wilcoxon test
clustMarkers <- FindAllMarkers(cur_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(clustMarkers, paste0(out_data_dir, 'clusterMarkers_filt_bg.csv'))

pdf(paste0(fig_dir, 'posClustMarkers_filt_bg.pdf'), height=4, width=9)
#png(paste0(fig_dir, 'pngs/posClustMarkers.png'), height=4, width=9, res=250, units='in')
DotPlot(cur_seurat, features = unique(as.character(unlist(tapply(clustMarkers$gene, clustMarkers$cluster, function(x){head(x, 2)}))))) + RotatedAxis()
dev.off()

pdf(paste0(fig_dir, 'allMarkers_bg.pdf'), height=4, width=9)
#png(paste0(fig_dir, 'pngs/posClustMarkers.png'), height=4, width=9, res=250, units='in')
DotPlot(cur_seurat, features = rownames(cur_seurat)) + RotatedAxis()
dev.off()

write.csv(cur_seurat@meta.data, paste0(out_data_dir,'seurat_object_harmony_clusts_filt_meta_bg.csv'))
saveRDS(cur_seurat, paste0(out_data_dir,'seurat_object_harmony_clusts_filt_bg_id.rds'))
counts_matrix <- GetAssayData(object = cur_seurat, assay='RNA')
writeMM(counts_matrix, paste0(out_data_dir,'seurat_object_harmony_clusts_filt_bg_id_counts.mtx'))
write.csv(Features(cur_seurat), paste0(out_data_dir,'seurat_object_harmony_clusts_filt_bg_id_mk_names.csv'))

print('END!!!')
