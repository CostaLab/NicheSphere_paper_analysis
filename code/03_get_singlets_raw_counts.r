## conda activate sceasy_env
library(Seurat)
library(Matrix)

dat=readRDS('../data/PICseq/intermediate_steps/rds/scrna_phase_comparing_singlet_final.Rds')
counts_matrix <- GetAssayData(object = dat, assay='RNA')
writeMM(counts_matrix, '../data/PICseq/intermediate_steps/rds/singlets_raw/singlets_raw.mtx')
write.csv(Features(dat), '../data/PICseq/intermediate_steps/rds/singlets_raw/singlets_raw_genes.csv')
write.csv(dat@meta.data, '../data/PICseq/intermediate_steps/rds/singlets_raw/singlets_raw_meta.csv')