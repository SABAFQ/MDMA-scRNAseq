suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))

plot_Volcano <- function(res, lfc.cutoff  = 1,pval.cutoff = 0.05, title='', x_min=-20,x_max=20, y_max=FALSE){
    res$sig <- as.factor(res$P.Value < pval.cutoff & abs(res$logFC) > lfc.cutoff)
    relevel(res$sig, ref=TRUE)
    res <- res %>% mutate(name=rownames(res))
    
    data = subset(res[order(res$P.Value),], sig == TRUE)[1:7,]
    
    if (y_max == FALSE){-log10(min(res$P.Value))}
    
    vol = res %>% ggplot(
        aes(x=logFC, y=-log10(P.Value), colour=sig, fill=sig)) +
        geom_point(aes(color = sig),alpha = 1/10) +
            geom_hline(yintercept=-log10(pval.cutoff), linetype="dashed", alpha = 4/10) +
            geom_vline(xintercept=lfc.cutoff, linetype="dashed", alpha = 4/10) +
            geom_vline(xintercept=(-1)*lfc.cutoff, linetype="dashed", alpha = 4/10) +
            scale_color_manual(values = c("grey", "red")) +
            theme_bw() + 
            theme(
                legend.position="none",plot.title = element_text(hjust = 0.5)
            ) +
            ggtitle (title) + 
            ggrepel::geom_text_repel(
                data = data,
                aes(label = rownames(data)),
                size = 3,
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            ) + 
            xlim(c(x_min,x_max)) +
            if (y_max != FALSE){ylim(c(0,y_max))}

     return (vol)
} 


plot_PCA <- function(vsd, gr, title,legend="none"){
    colData <- colData(vsd)
    z <- plotPCA(vsd,intgroup=gr, returnData=TRUE)
    percentVar <- round(100 * attr(z, "percentVar"))
    pca <- ggplot(z, aes(PC1, PC2)) + 
            geom_point(aes(size = 0.01, colour=group,fill=group), alpha = 0.1,shape=18) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) +
            guides (size = FALSE) +
            ggtitle (title)+ 
            theme(legend.position=legend)
    return (pca)
}
