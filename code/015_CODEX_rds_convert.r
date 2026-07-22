### Data convert
## conda activate sceasy_env

library(reticulate)
# Replace 'sceasy_env' with the exact name you used
use_condaenv("sceasy_env", required = TRUE) 
py_config()
reticulate::py_module_available(module = 'anndata')
## TRUE

library(scCustomize)
library(Seurat)

obj=readRDS('../data/CODEX/intermediate_steps/seurat_object_harmony_clusts_filt_bg_id.rds')
# This function natively supports Seurat v5 objects
as.anndata(x = obj, file_path = "../data/CODEX/intermediate_steps/", file_name = "seurat_object_harmony_clusts_filt_bg_id.h5ad")