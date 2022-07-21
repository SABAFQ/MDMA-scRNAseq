args <- commandArgs(trailingOnly = TRUE)

name <- args[1]

dir.create(paste0('diffexp/',name))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))
# suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SingleCellExperiment))

sce_prep <- function (counts,meta){
    ### Load data as `SingleCellExperiment` object 
    sce <- SingleCellExperiment(list(counts=counts),colData=meta)
    
    # ## low count filter - at least 10 cells with count of 5 or more
    # keep <- rowSums(counts(sce) >= 1) >= 10 
    # sce <- sce[keep,] 
    
    assay(sce) <- as.matrix(assay(sce))

    return (sce)
}

plotDisp <- function(dds){
    p1 = plotDispEsts(dds)
    keep <- rowSums(counts(dds) >= 10) >= 25
    dds2 <- estimateDispersionsFit(dds[keepForDispTrend,])
    p2 = plotDispEsts(dds2, ylim=c(1e-3,1))
}


correct_batch <- function (dds,gr,out){
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    p0 <- plot_PCA(vsd, gr, 'Before removeBatchEffect')
    mat <- assay(vsd)
    
    mat <- limma::removeBatchEffect(mat, vsd$batch) 
    assay(vsd) <- mat
    p1 <- plot_PCA(vsd, gr, 'After removeBatchEffect',legend="right")
    counts_batch_corrected <- assay(vsd)
    
    if (out == 'plot'){return (p0 + p1)}
    if (out == 'vsd') {return (vsd)}
    if (out == 'cbc') {return (counts_batch_corrected)}
} 

get_res<- function(fit,coef){
    fc <- fit[,coef]$coef [,coef] %>% data.frame 
    fc <- fc %>% dplyr::rename(log2FoldChange=colnames(fc))
    pv <- fit[,coef]$p.value [,coef] %>% data.frame 
    pv <- pv %>% dplyr::rename(pvalue=colnames(pv))
    res <- cbind(fc, pv)
    return(res)
}


counts <- fread(paste0('diffexp/',name,'_counts.txt'),sep='\t') %>% 
    data.frame %>% column_to_rownames('V1')

meta <- fread(paste0('diffexp/',name,'_meta.txt'),sep='\t') %>% 
    data.frame %>% column_to_rownames('V1') %>% dplyr::rename('cond' = 'group')

rownames(meta) <- colnames(counts) 

meta$cond <- factor (meta$cond)
meta$cond <- relevel(meta$cond, ref="MDMA")
conds <- meta$cond

sc <- sce_prep(counts,meta)
de <- convertTo(sc, type="edgeR")

design <- model.matrix(~0 + conds)
colnames(design) <- gsub("conds", "", colnames(design))

contr.matrix <- makeContrasts(
    MDMAvsMeth= MDMA - Meth, 
    MDMAvsSaline= MDMA - Saline,
    levels = colnames(design)
)

# Log tranform & Quantile normalize
de2 <- normalizeBetweenArrays(log(de$counts + 0.001))
vfit2 <- lmFit(de2, design)
vfit2 <- contrasts.fit(vfit2, contrasts=contr.matrix)
efit2 <- eBayes(vfit2)

res_MDMAvsMeth <- get_res(efit2,'MDMAvsMeth')
write.table(
    res_MDMAvsMeth,
    paste0('diffexp/',name,'_delta_exp_MDMAvsMeth.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_MDMAvsSaline <- get_res(efit2,'MDMAvsSaline')
write.table(
    res_MDMAvsSaline,
    paste0('diffexp/',name,'delta_exp_MDMAvsSaline.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)
