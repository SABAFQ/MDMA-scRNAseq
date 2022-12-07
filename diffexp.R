args <- commandArgs(trailingOnly = TRUE)

name <- args[1]

suppressWarnings(dir.create(paste0('diffexp/',name)))

message('---\tload libraries\t---')
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


prep_de <- function(de,design,contr.matrix){
    # Log tranform & Quantile normalize
    de2 <- normalizeBetweenArrays(log(de$counts + 0.001))
    vfit2 <- lmFit(de2, design)
    vfit2 <- contrasts.fit(vfit2, contrasts=contr.matrix)
    efit2 <- eBayes(vfit2)
    return(efit2)
}


get_res<- function(fit,coef){
    # https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_March_UCSF_mRNAseq_Workshop/master/differential_expression/DE_Analysis.Rmd
    
    res <- topTable(fit = fit, coef = coef, sort.by = "none", n = Inf)
    # fc <- fit[,coef]$coef [,coef] %>% data.frame 
    # fc <- fc %>% dplyr::rename(log2FoldChange=colnames(fc))
    # pv <- fit[,coef]$p.value [,coef] %>% data.frame 
    # pv <- pv %>% dplyr::rename(pvalue=colnames(pv))
    # res <- cbind(fc, pv)
    return(res)
}


message('---\tload counts\t---')
counts <- fread(paste0('diffexp/',name,'_counts.txt'),sep='\t') %>% 
    data.frame %>% column_to_rownames('V1')

message('---\tprep meta\t---')

meta <- fread(paste0('diffexp/',name,'_meta.txt'),sep='\t') %>% 
    data.frame %>% column_to_rownames('V1') %>% dplyr::rename('cond' = 'group')

rownames(meta) <- colnames(counts) 

meta$cond <- factor (meta$cond)
meta$cond <- relevel(meta$cond, ref="MDMA")
conds <- meta$cond

meta$treat <- factor (meta$treat)
meta$treat <- relevel(meta$treat, ref="notAcute")
treat <- meta$treat

message('---\trun DE analysis\t---')
sc <- sce_prep(counts,meta)
de <- convertTo(sc, type="edgeR")

### design a
design_a <- model.matrix(~0 + conds)
colnames(design_a) <- gsub("conds", "", colnames(design_a))

contr.matrix_a <- makeContrasts(
    MDMAvsMeth= MDMA - Meth,
    MethvsSaline= Meth - Saline,
    MDMAvsSaline= MDMA - Saline,
    levels = colnames(design_a)
)

efit_a <- prep_de(de,design_a,contr.matrix_a)

res_MDMAvsMeth <- get_res(efit_a,'MDMAvsMeth')

write.table(
    res_MDMAvsMeth,
    paste0('diffexp/',name,'/delta_exp_MDMAvsMeth.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_MDMAvsSaline <- get_res(efit_a,'MDMAvsSaline')
write.table(
    res_MDMAvsSaline,
    paste0('diffexp/',name,'/delta_exp_MDMAvsSaline.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_MethvsSaline <- get_res(efit_a,'MethvsSaline')
write.table(
    res_MethvsSaline,
    paste0('diffexp/',name,'/delta_exp_MethvsSaline.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

### design b
design_b <- model.matrix(~0 + treat)
colnames(design_b) <- gsub("treat", "", colnames(design_b))

contr.matrix_b <- makeContrasts(
    MDMA_Acute_vs_notAcute = Acute - notAcute,
    levels = colnames(design_b)
)

efit_b <- prep_de(de,design_b,contr.matrix_b)

res_MDMA_Acute_vs_notAcute <- get_res(efit_b,'MDMA_Acute_vs_notAcute')

write.table(
    res_MDMA_Acute_vs_notAcute,
    paste0('diffexp/',name,'/delta_exp_MDMA_Acute_vs_notAcute.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)
