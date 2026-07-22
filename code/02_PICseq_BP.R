#!/usr/bin/env Rscript

#### Run BayesPrism on PIC-seq dataset

library(CIMseq)
library(BayesPrism)

PICseq_singlets<-readRDS("../data/PICseq/intermediate_steps/rds/scrna_phase_comparing_singlet_final_noSMCs_Xist.Rds")
PICseq_doublets<-readRDS("../data/PICseq/intermediate_steps/rds/scrna_phase_clustering_doublets.Rds")

PICseqDoublets_counts<-as.matrix(PICseq_doublets@assays$RNA@counts)
PICseqSinglets_counts<-as.matrix(PICseq_singlets@assays$RNA@counts)
PICseqDoublets_metadata<-PICseq_doublets@meta.data
PICseqSinglets_metadata<-PICseq_singlets@meta.data

rm(PICseq_singlets)
rm(PICseq_doublets)

sc.stat <- plot.scRNA.outlier(
    input=t(PICseqSinglets_counts), 
    cell.type.labels=PICseqSinglets_metadata$annotation,
    species="mm",
    return.raw=TRUE, #return the data used for plotting. 
    pdf.prefix="gbm.sc.stat" #specify pdf.prefix if need to output to pdf
)

bk.stat <- plot.bulk.outlier(
    bulk.input=t(PICseqDoublets_counts),
    sc.input=t(PICseqSinglets_counts),
    cell.type.labels=PICseqSinglets_metadata$annotation,  
    species="mm", #currently only human(hs) and mouse(mm) annotations are supported
    return.raw=TRUE,
    pdf.prefix="gbm.bk.stat" #specify pdf.prefix if need to output to pdf
)

sc.dat.filtered <- cleanup.genes (input=t(PICseqSinglets_counts),
    input.type="count.matrix",
    species="mm",
    gene.group=c("Rb","Mrp","other_Rb","chrM","chrX","chrY"),
    exp.cells=5)

commonGenes<-intersect(rownames(PICseqDoublets_counts), colnames(sc.dat.filtered))
PICseqDoublets_counts<-t(PICseqDoublets_counts)
PICseqSinglets_counts<-t(PICseqSinglets_counts)

diff.exp.stat <- get.exp.stat(sc.dat=PICseqSinglets_counts[,colSums(PICseqSinglets_counts>0)>3],
    cell.type.labels=PICseqSinglets_metadata$annotation,
    cell.state.labels=PICseqSinglets_metadata$annotation,
    psuedo.count=0.1, 
    cell.count.cutoff=20, 
    n.cores=10) #number of threads

sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered,
    stat=diff.exp.stat,
    pval.max=0.01,
    lfc.min=0.1)

myPrism <- new.prism(
    reference=sc.dat.filtered.pc.sig, 
    mixture=PICseqDoublets_counts[,intersect(colnames(PICseqDoublets_counts),   colnames(sc.dat.filtered.pc.sig))],
    input.type="count.matrix", 
    cell.type.labels = PICseqSinglets_metadata$annotation, 
    cell.state.labels = PICseqSinglets_metadata$annotation,
    key=NULL,
    outlier.cut=0.01,
    outlier.fraction=0.1)

print('checkpoint!')

bp.res <- run.prism(prism = myPrism, n.cores=25)

theta <- get.fraction (bp=bp.res,
    which.theta="final",
    state.or.type="type")

write.csv(theta, '../res_data/revision/deconvolution/PICseq_doublets_CTprops_BP_annotation_noMSCs_Xist.csv')
print('done!')

