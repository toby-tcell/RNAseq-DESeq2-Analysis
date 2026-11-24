library(SummarizedExperiment)
library(tidyverse)
library(ComplexHeatmap)
library(DESeq2) #cite ApeGLM
library(vsn)
library(ggrepel)
library(org.Hs.eg.db) #cite this
library(msigdbr) #cite this
library(clusterProfiler) #cite this
library(enrichplot) #cite this


dds <- readRDS('Data/DESeqDataset.rds')
plotpath <- 'Plots/DESeq'
dir.create(plotpath, showWarnings = F)
msig <- msigdbr(species = 'human')

levels(dds$Line)
levels(dds$Line) <- gsub('\\+','', levels(dds$Line))
levels(dds$Line) <- gsub('PD1-','', levels(dds$Line))
levels(dds$Line)

levels(dds$Condition) <- gsub('\\+','_', levels(dds$Condition))
levels(dds$Condition)

dds$Group <- factor(paste0(dds$Line, '_', dds$Condition))
design(dds) <- formula(~ Group)
levels(dds$Group)



# Run DESeq for all refs required

dir.create('Data/DESeqRDS', showWarnings = F)

dds$Group <- relevel(dds$Group, ref = 'Truncated_0_hours')
dds_Trunc_0 <- DESeq(dds)
saveRDS(dds_Trunc_0, 'Data/DESeqRDS/Trunc_0.rds')

dds$Group <- relevel(dds$Group, ref = 'TIM3_0_hours')
dds_TIM3_0 <- DESeq(dds)
saveRDS(dds_TIM3_0, 'Data/DESeqRDS/TIM3_0.rds')

dds$Group <- relevel(dds$Group, ref = 'TIM3_265272F_0_hours')
dds_TIM3_265272_0 <- DESeq(dds)
saveRDS(dds_TIM3_265272_0, 'Data/DESeqRDS/265272_0.rds')

dds$Group <- relevel(dds$Group, ref = 'TIM3_277trunc_0_hours')
dds_TIM3_277trunc_0 <- DESeq(dds)
saveRDS(dds_TIM3_277trunc_0, 'Data/DESeqRDS/277trunc_0.rds')

dds$Group <- relevel(dds$Group, ref = 'Truncated_24_hours_block')
dds_Trunc_24b <- DESeq(dds)
saveRDS(dds_Trunc_24b, 'Data/DESeqRDS/Trunc_24b.rds')

dds$Group <- relevel(dds$Group, ref = 'TIM3_24_hours_block')
dds_TIM3_24b <- DESeq(dds)
saveRDS(dds_TIM3_24b, 'Data/DESeqRDS/TIM3_24b.rds')

dds$Group <- relevel(dds$Group, ref = 'TIM3_24_hours')
dds_TIM3_24 <- DESeq(dds)
saveRDS(dds_TIM3_24, 'Data/DESeqRDS/TIM3_24.rds')

dds$Group <- relevel(dds$Group, ref = 'Truncated_24_hours')
dds_Trunc_24 <- DESeq(dds)
saveRDS(dds_Trunc_24, 'Data/DESeqRDS/Trunc_24.rds')

# Function to run standard workflow

DESeq2_analysis <- function(filepath, coef){
    
    results <- lfcShrink(dds, coef = coef)
    results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    summary(results)
    
    plotpath <- paste0('Plots/DESeq/', filepath)
    dir.create(plotpath, showWarnings = F)
    
    gseacolors <- c(
        "blue4",  
        "red3",  
        "green4"
    )
    
    
    # Plot volcano plot
    volcano_df <- as.data.frame(results)
    volcano_df$group <- "Not significant"
    volcano_df$group[
        volcano_df$padj < 0.05 &
            volcano_df$log2FoldChange > 0.25
    ] <- "Positive"
    volcano_df$group[
        volcano_df$padj < 0.05 &
            volcano_df$log2FoldChange < -0.25
    ] <- "Negative"
    
    positive <- sum(volcano_df$log2FoldChange > 0.25)
    negative <- sum(volcano_df$log2FoldChange < -0.25)
    label_df <- volcano_df %>%
        filter(
            abs(log2FoldChange) > 0.25,
            padj < 0.01
        )
    
    # Remove excess labels
    if(nrow(label_df>25)){
        label_df <- rbind(
            head(label_df %>%
                     arrange(desc(-log10(padj)*abs(log2FoldChange))),
                 n=20),
            head(label_df %>%
                     arrange(desc(abs(log2FoldChange))), n=5),
            head(label_df %>%
                     arrange(padj), n=5)
        ) %>%
        distinct()
            
            
            
            label_df %>%
            arrange(desc(-log10(padj)*abs(log2FoldChange))) %>%
            head(n=20)
    }
    
    ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(data = volcano_df, colour = "grey90") +
        geom_point(data = subset(volcano_df, group == "Negative"), colour = "blue4", , alpha = 0.5) + 
        geom_point(data = subset(volcano_df, group == "Positive"), colour = "red3", alpha = 0.5) + 
        geom_vline(xintercept = c(0.25,-0.25), linetype = "dashed", color = 'grey40') +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = 'grey40') +
        geom_text_repel(
            data = label_df,
            aes(label = gene_symbol),
            color = 'black', segment.size = 0.2,
            max.overlaps = Inf, 
            point.padding = 0.2, box.padding = 0.2, min.segment.length = 0.2
        ) +
        theme_classic(base_size = 14) 
    ggsave(file.path(plotpath, 'Volcanoplot.pdf'), height = 6, width = 7)
    
    
    # Hallmark gene enrichment
    genelist <- results %>% as.data.frame() %>%
        filter(!is.na(log2FoldChange), !is.na(gene_symbol)) %>%
        distinct(gene_symbol, .keep_all = T) %>%
        arrange(desc(log2FoldChange)) %>%
        with(setNames(log2FoldChange, gene_symbol))
    
    hallmark <- msigdbr(
        species = 'Homo sapiens',
        collection = 'H'
    ) %>%
        select(gs_name, gene_symbol)
    
    hallmark$gs_name <- gsub('HALLMARK_','',hallmark$gs_name)
    
    gsea <- GSEA(
        genelist,
        TERM2GENE = hallmark,
        pvalueCutoff = 0.05,
        minGSSize = 20,
        maxGSSize = 500,
        verbose = FALSE
    )
    
    gsea_res <- as.data.frame(gsea) %>%
        arrange(abs(NES))
    
    # Dotplot top 10
    
    restop10 <- tail(gsea_res, n=10)

    restop10$ID <- factor(restop10$ID, levels = restop10$ID)
    
    ggplot(restop10, aes(x = ID, y = NES, colour = NES)) +
        geom_point(aes(size = -log10(p.adjust))) +
        scale_colour_gradient2(
            low = "blue4",
            mid = "white",
            high = "red3",
            midpoint = 0
        ) +
        theme_bw(base_size = 14) + 
        coord_flip() +
        labs(y='Normalized Enrichment Score (NES)',x = '', colour = 'NES', size = '-log10(padj)') +
        theme(axis.text.y = element_text(size = 10, colour = 'black')) +
        ylim(min(restop10$NES)-0.5, max(restop10$NES)+0.5)
    
    
    ggsave(file.path(plotpath, 'pathwayenrichmentscores_hallmark.pdf'), height = 5, width = 6)
    
    # GSEA plot for top 3
    p<- NULL
    p<- gseaplot2(gsea, geneSetID = 1:3, pvalue_table = F, color = gseacolors, base_size = 14)
    ggsave(plot = p, filename = file.path(plotpath, 'top3enrichment_hallmark.pdf'), height = 8, width = 9)
    p <- NULL
    
    # GSEA plot for NFKB, IL2 and IL6 signalling
    pathways <- c('TNFA_SIGNALING_VIA_NFKB','INFLAMMATORY_RESPONSE','IL2_STAT5_SIGNALING')
    if (any(pathways %in% gsea_res$ID)){
        present <- pathways[pathways %in% gsea@result$ID]
        q <- NULL
        q<- gseaplot2(gsea, geneSetID = present, pvalue_table = F, color = gseacolors, base_size = 14)
        ggsave(plot = q, filename = file.path(plotpath, 'TCR_enrichment_hallmark.pdf'), height = 8, width = 9)
        q <- NULL
    }
    print('Completed analysis')
    
    
}

# Run Trunc_0 v Trunc_24+block

dds <- readRDS('Data/DESeqRDS/Trunc_0.rds')
resultsNames(dds)

coef = 'Group_Truncated_24_hours_block_vs_Truncated_0_hours'
filepath = 'Trunc_activation'


DESeq2_analysis(coef = coef, filepath = filepath)


# Run TIM3_0 v TIM3_24+block

dds <- readRDS('Data/DESeqRDS/TIM3_0.rds')
resultsNames(dds)

coef = 'Group_TIM3_24_hours_block_vs_TIM3_0_hours'
filepath = 'TIM3_activation'
DESeq2_analysis(coef = coef, filepath = filepath)

# Run TIM3_277trunc_0 v TIM3_277trunc_24+block

dds <- readRDS('Data/DESeqRDS/277trunc_0.rds')
resultsNames(dds)

coef = 'Group_TIM3_277trunc_24_hours_block_vs_TIM3_277trunc_0_hours'
filepath = 'TIM3_277trunc_activation'
DESeq2_analysis(coef = coef, filepath = filepath)

# Run TIM3_265+272_0 v TIM3_265+272_24+block

dds <- readRDS('Data/DESeqRDS/265272_0.rds')
resultsNames(dds)

coef = 'Group_TIM3_265272F_24_hours_block_vs_TIM3_265272F_0_hours'
filepath = 'TIM3_265+272F_activation'
DESeq2_analysis(coef = coef, filepath = filepath)

# Run Trunc_24block v Trunc_24

dds <- readRDS('Data/DESeqRDS/Trunc_24b.rds')
resultsNames(dds)

coef = 'Group_Truncated_24_hours_vs_Truncated_24_hours_block'
filepath = 'Trunc_engagedvnonengaged'
DESeq2_analysis(coef = coef, filepath = filepath)

# Run TIM3_24block v TIM3_24

dds <- readRDS('Data/DESeqRDS/TIM3_24b.rds')
resultsNames(dds)

coef = 'Group_TIM3_24_hours_vs_TIM3_24_hours_block'
filepath = 'TIM3_engagedvnonengaged'
DESeq2_analysis(coef = coef, filepath = filepath)

# Run Trunc_24 v TIM3_24

dds <- readRDS('Data/DESeqRDS/Trunc_24.rds')
resultsNames(dds)

coef = 'Group_TIM3_24_hours_vs_Truncated_24_hours'
filepath = 'TIM3_vTrunc_engaged'
DESeq2_analysis(coef = coef, filepath = filepath)

# Run TIM3_277trunc_24 AND TIM3_265+272F v TIM3_24 

dds <- readRDS('Data/DESeqRDS/TIM3_24.rds')
resultsNames(dds)

coef = 'Group_TIM3_277trunc_24_hours_vs_TIM3_24_hours'
filepath = 'TIM3277trunc_vTIM3_engaged'
DESeq2_analysis(coef = coef, filepath = filepath)

coef = 'Group_TIM3_265272F_24_hours_vs_TIM3_24_hours'
filepath = 'TIM3265+272F_vTIM3_engaged'
DESeq2_analysis(coef = coef, filepath = filepath)



