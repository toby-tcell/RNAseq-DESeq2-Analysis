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


activationgenes <- c('IL2RA','CD69','IL2','NFKBIA','MYC','SELL','MAL','KLF2')

#  6.1 Trunc Activation ----

dds <- readRDS('Data/DESeqRDS/Trunc_0.rds')
resultsNames(dds)

coef = 'Group_Truncated_24_hours_block_vs_Truncated_0_hours'
filepath = 'Trunc_activation'

results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)

# Activation genes
lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)




activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'logFC Activated v Control', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-1, max(activation_df$log2FoldChange) +1) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5)

#  6.1 TIM3 Activation ----

dds <- readRDS('Data/DESeqRDS/TIM3_0.rds')
resultsNames(dds)

coef = 'Group_TIM3_24_hours_block_vs_TIM3_0_hours'
filepath = 'TIM3_activation'

results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)

# Activation genes
lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)
   



activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'logFC Activated v Control', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-1, max(activation_df$log2FoldChange) +1) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5)


# 6.2 Trunc ± block - activation genes

dds <- readRDS('Data/DESeqRDS/Trunc_24b.rds')
resultsNames(dds)

coef = 'Group_Truncated_24_hours_vs_Truncated_24_hours_block'
filepath = 'Trunc_engagedvnonengaged'

results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)

# Activation genes
lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)


activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'LogFC Engaged v Non-engaged', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-0.5, max(activation_df$log2FoldChange) +0.5) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5)

truncactivation_df <- activation_df






# 6.2 TIM3 ± block - epithelial transition ----

dds <- readRDS('Data/DESeqRDS/TIM3_24b.rds')
resultsNames(dds)

coef = 'Group_TIM3_24_hours_vs_TIM3_24_hours_block'
filepath = 'TIM3_engagedvnonengaged'

results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)


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




epithelial_genes <- unlist(strsplit(gsea_res["EPITHELIAL_MESENCHYMAL_TRANSITION","core_enrichment"], "/"))
coagulation_genes <- unlist(strsplit(gsea_res["COAGULATION","core_enrichment"], "/"))

# Epithelial
plot_df <- data.frame(gene = epithelial_genes,
               hallmark = "Epithelial")
    
lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

plot_df <- plot_df %>%
    left_join(lfc_df, by = "gene")

plot_df <- plot_df %>%
    arrange(abs(log2FoldChange)) %>%
    mutate(rank = row_number())

label_df <- plot_df %>%
    slice_max(order_by = abs(log2FoldChange), n = 8)

ggplot(plot_df,
       aes(x = rank,
           y = log2FoldChange)) +
    geom_point(size = 2) +
    scale_color_manual(values = c('red3','blue3')) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "Leading edge gene rank",
        y = "LogFC Engaged v Non-engaged",
    ) +
    geom_text_repel(data = label_df, aes(label = gene), colour = 'black',
                    max.overlaps = 9, direction = 'both', point.padding = 0.3, nudge_y = 0.01, box.padding = 0.6, force = 3, segment.alpha = 0.5) +
    theme_bw(base_size = 14)
ggsave(file.path(plotpath, 'Epithelial_genes.pdf'), height = 5, width = 5)

# Coagulation
plot_df <- data.frame(gene = coagulation_genes,
                      hallmark = "Coagulation")

lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

plot_df <- plot_df %>%
    left_join(lfc_df, by = "gene")

plot_df <- plot_df %>%
    arrange(abs(log2FoldChange)) %>%
    mutate(rank = row_number())

label_df <- plot_df %>%
    group_by(hallmark) %>%
    slice_max(order_by = abs(log2FoldChange), n = 8)

ggplot(plot_df,
       aes(x = rank,
           y = log2FoldChange)) +
    geom_point(size = 2) +
    scale_color_manual(values = c('red3','blue3')) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "Leading edge gene rank",
        y = "LogFC Engaged v Non-engaged",
    ) +
    geom_text_repel(data = label_df, aes(label = gene), colour = 'black',
                    max.overlaps = 9, direction = 'both', point.padding = 0.3, nudge_y = 0.01, box.padding = 0.6, force = 3, segment.alpha = 0.5) +
    theme_bw(base_size = 14)
ggsave(file.path(plotpath, 'Coagulation_genes.pdf'), height = 5, width = 5)


# Activation genes

activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'LogFC Engaged v Non-engaged', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-0.3, max(activation_df$log2FoldChange) +0.3) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5)

# Comparison with trunc

truncactivation_df$line <- 'Truncated'

activation_df$line <- 'TIM-3'

comparison <- rbind.data.frame(activation_df, truncactivation_df) %>%
    select(!padj) %>%
    pivot_wider(names_from = 'line', values_from = 'log2FoldChange', ) %>%
    mutate(delta = comparison$`TIM-3`-comparison$Truncated)

ggplot(comparison, aes(x = delta, y = gene, fill = delta)) +
    geom_col() +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = '\u0394 LogFC TIM-3 v Truncated', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes_comparison.pdf'), height = 5, width = 5, device = cairo_pdf)



# 6.3 TIM3 v Trunc engaged Hypoxia and myogenesis ----

dds <- readRDS('Data/DESeqRDS/Trunc_24.rds')
resultsNames(dds)

coef = 'Group_TIM3_24_hours_vs_Truncated_24_hours'
filepath = 'TIM3_vTrunc_engaged'

results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)



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



hypoxia_genes <- unlist(strsplit(gsea_res["HYPOXIA","core_enrichment"], "/"))
cholesterol_genes <- unlist(strsplit(gsea_res["CHOLESTEROL_HOMEOSTASIS","core_enrichment"], "/"))

# Hypoxia
plot_df <- data.frame(gene = hypoxia_genes,
                      hallmark = "Hypoxia")

lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

plot_df <- plot_df %>%
    left_join(lfc_df, by = "gene")

plot_df <- plot_df %>%
    arrange(abs(log2FoldChange)) %>%
    mutate(rank = row_number())

label_df <- plot_df %>%
    group_by(hallmark) %>%
    slice_max(order_by = abs(log2FoldChange), n = 8)

ggplot(plot_df,
       aes(x = rank,
           y = log2FoldChange)) +
    geom_point(size = 2) +
    scale_color_manual(values = c('red3','blue3')) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "Leading edge genes",
        y = "LogFC PD-1-TIM-3 v PD-1-Truncated",
    ) +
    geom_text_repel(data = label_df, aes(label = gene), colour = 'black',
                    max.overlaps = 9, direction = 'both', point.padding = 0.3, nudge_y = 0.01, box.padding = 0.6, force = 3, segment.alpha = 0.5) +
    theme_bw(base_size = 14)
ggsave(file.path(plotpath, 'Hypoxia_genes.pdf'), height = 5, width = 5)

# Cholesterol
plot_df <- data.frame(gene = cholesterol_genes,
                      hallmark = "Cholesterol")

lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

plot_df <- plot_df %>%
    left_join(lfc_df, by = "gene")

plot_df <- plot_df %>%
    arrange(abs(log2FoldChange)) %>%
    mutate(rank = row_number())

label_df <- plot_df %>%
    group_by(hallmark) %>%
    slice_max(order_by = abs(log2FoldChange), n = 8)

ggplot(plot_df,
       aes(x = rank,
           y = log2FoldChange)) +
    geom_point(size = 2) +
    scale_color_manual(values = c('red3','blue3')) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "Leading edge genes",
        y = "LogFC PD-1-TIM-3 v PD-1-Truncated",
    ) +
    geom_text_repel(data = label_df, aes(label = gene), colour = 'black',
                    max.overlaps = 9, direction = 'both', point.padding = 0.3, nudge_y = 0.01, box.padding = 0.6, force = 3, segment.alpha = 0.5) +
    theme_bw(base_size = 14)
ggsave(file.path(plotpath, 'Cholesterol_genes.pdf'), height = 5, width = 5)


activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'LogFC PD-1-TIM-3 v PD-1-Truncated', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-0.3, max(activation_df$log2FoldChange) +0.3) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5)




# 6.4 277trunc v TIM3 cholesterol and E2F ----

dds <- readRDS('Data/DESeqRDS/TIM3_24.rds')
resultsNames(dds)

coef = 'Group_TIM3_277trunc_24_hours_vs_TIM3_24_hours'
filepath = 'TIM3277trunc_vTIM3_engaged'


results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)


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




cholesterol_genes <- unlist(strsplit(gsea_res["CHOLESTEROL_HOMEOSTASIS","core_enrichment"], "/"))
E2F_genes <- unlist(strsplit(gsea_res["E2F_TARGETS","core_enrichment"], "/"))


# Cholesterol
plot_df <- data.frame(gene = cholesterol_genes,
                      hallmark = "Cholesterol")

lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

plot_df <- plot_df %>%
    left_join(lfc_df, by = "gene")

plot_df <- plot_df %>%
    arrange(abs(log2FoldChange)) %>%
    mutate(rank = row_number())

label_df <- plot_df %>%
    group_by(hallmark) %>%
    slice_max(order_by = abs(log2FoldChange), n = 8)

ggplot(plot_df,
       aes(x = rank,
           y = log2FoldChange)) +
    geom_point(size = 2) +
    scale_color_manual(values = c('red3','blue3')) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "Leading edge genes",
        y = "LogFC PD-1-TIM-3 \u0394 277  v PD-1-TIM-3",
    ) +
    geom_text_repel(data = label_df, aes(label = gene), colour = 'black',
                    max.overlaps = 10, direction = 'both', point.padding = 0.3, nudge_y = 0.01, box.padding = 0.6, force = 3, segment.alpha = 0.5) +
    theme_bw(base_size = 14)
ggsave(file.path(plotpath, 'Cholesterol_genes.pdf'), height = 5, width = 5, device = cairo_pdf)

# NFKB
plot_df <- data.frame(gene = NFKB_genes,
                      hallmark = "NFKB")

lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

plot_df <- plot_df %>%
    left_join(lfc_df, by = "gene")

plot_df <- plot_df %>%
    arrange(abs(log2FoldChange)) %>%
    mutate(rank = row_number())

label_df <- plot_df %>%
    group_by(hallmark) %>%
    slice_max(order_by = abs(log2FoldChange), n = 8)

ggplot(plot_df,
       aes(x = rank,
           y = log2FoldChange)) +
    geom_point(size = 2) +
    scale_color_manual(values = c('red3','blue3')) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "Leading edge genes",
        y = "LogFC PD-1-TIM-3 \u0394 277  v PD-1-TIM-3",
    ) +
    geom_text_repel(data = label_df, aes(label = gene), colour = 'black',
                    max.overlaps = 9, direction = 'both', point.padding = 0.3, nudge_y = 0.01, box.padding = 0.6, force = 3, segment.alpha = 0.5) +
    theme_bw(base_size = 14)
ggsave(file.path(plotpath, 'NFKB_genes.pdf'), height = 5, width = 5, device = cairo_pdf)


activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'LogFC PD-1-TIM-3 \u0394 277  v PD-1-TIM-3', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-0.3, max(activation_df$log2FoldChange) +0.3) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5, device = cairo_pdf)



# 6.5 265+272F v TIM3 signalling ----

dds <- readRDS('Data/DESeqRDS/TIM3_24.rds')
resultsNames(dds)

coef = 'Group_TIM3_265272F_24_hours_vs_TIM3_24_hours'
filepath = 'TIM3265+272F_vTIM3_engaged'

results <- lfcShrink(dds, coef = coef)
results$gene_symbol <- mapIds(org.Hs.eg.db,keys=row.names(results),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
summary(results)

plotpath <- paste0('Plots/DESeq/', filepath)
dir.create(plotpath, showWarnings = F)


# Add LFC from DESeq2 results
lfc_df <- as.data.frame(results) %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    rename(gene = gene_symbol)

activation_df <- lfc_df %>%
    filter(gene %in% activationgenes) %>%
    mutate(significance = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
activation_df$gene <- factor(activation_df$gene, levels = activationgenes)
ggplot(activation_df, aes(y = gene, x = log2FoldChange, fill = log2FoldChange, label = significance)) +
    geom_col() +
    geom_text(hjust = ifelse(activation_df$log2FoldChange>0, -0.5, 1.5)) +
    scale_fill_gradient2(low = 'blue4', mid = 'grey',high= 'red3') +
    theme_classic(base_size = 12) +
    labs(x = 'LogFC PD-1-TIM-3 265+272F  v PD-1-TIM-3', y = NULL) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    xlim(min(activation_df$log2FoldChange)-0.01, max(activation_df$log2FoldChange) +0.01) +
    theme(legend.position = 'none')
ggsave(file.path(plotpath, 'Activation_genes.pdf'), height = 5, width = 5)


