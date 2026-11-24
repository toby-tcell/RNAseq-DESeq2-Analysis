library(SummarizedExperiment)
library(tidyverse)
library(ComplexHeatmap)
library(DESeq2)
library(vsn)
library(ggrepel)
library(viridis)
library(org.Hs.eg.db)

plotpath <- 'Plots/QC/'
dir.create(plotpath, showWarnings = F)


se <- readRDS('Data/summarizedexp.rds')
dds <- readRDS('Data/DESeqDataset.rds')
vsd <- readRDS('Data/Normalized_se.rds')
mycolors <- c(
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # green
    "#CC79A7",  # purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000"   # black
)


# Size factors/lib sizes ----




colData(dds) |>
    as.data.frame() |>
    ggplot(aes(x = Sample_name, y = libSize / 1e6, fill = Line)) + 
    geom_bar(stat = "identity") + theme_bw() + 
    scale_fill_manual(values = mycolors) +
    labs(x = "Sample", y = "Total count (millions)") + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 1)) +
    coord_flip()
ggsave(file.path(plotpath, 'libsizes.pdf'), height = 6, width = 6)


colData(dds) %>%
    as.data.frame() %>%
    ggplot(aes(x = Sample_name, y = sizeFactor, fill = Line)) + 
    geom_bar(stat = "identity") + theme_bw() + 
    scale_fill_manual(values = mycolors) +
    labs(x = "Sample", y = "Size Factor") + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 1)) +
    coord_flip()
ggsave(file.path(plotpath, 'Sizefactors.pdf'), height = 6, width = 6)

colData(dds) %>%
    as.data.frame() %>%
    ggplot(aes(x = libSize/1e6, y = sizeFactor)) + 
    scale_color_manual(values = mycolors) +
    geom_smooth(se = F, linewidth = 0.5, color = 'grey50') +
    geom_point(aes(size = libSize/1e6, color = Line)) + theme_bw() + 
    labs(x = "Total count (millions)", y = "Size Factor", size = 'Total count (millions)') + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 1)) +
    coord_flip()
ggsave(file.path(plotpath, 'libsizevSizefactors.pdf'), height = 6, width = 7)




# PCA analysis ----

vsd <- readRDS('Data/Normalized_se.rds')

pcaData <- DESeq2::plotPCA(vsd, intgroup = c("Line","Condition"),
                           returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


centroids <- pcaData %>%
    group_by(Sample_name) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2))

ggplot(pcaData, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Line, shape = Condition), size = 4, show.legend = T) +
    scale_color_manual(values = mycolors) +
    theme_bw(base_size = 14) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) 
ggsave(file.path(plotpath, 'PCA.pdf'), height = 5, width = 7)



# Sample density vsd counts ----

assay(vsd) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene,
                 names_to = "sample",
                 values_to = "expression") %>%
    as.data.frame() %>%
    ggplot(aes(x = expression, colour = sample)) +
    scale_colour_viridis_d(option = 'D') +
    geom_density() + 
    labs(x='Normalized Expression',y='Density') +
    theme_bw(base_size = 14) +
    theme(legend.position = 'none') 
ggsave(file.path(plotpath, 'Sample_density.pdf'), height = 6, width = 6)


# Sample distance heatmap ----

mat <- assay(vsd)
sampleDists <- dist(t(mat)) %>%
    as.matrix()
samplenames<- colData(vsd) %>%
    as.data.frame() %>%
    pull(Sample_Type)
rownames(sampleDists) <- samplenames
colnames(sampleDists) <- samplenames

pdf(file.path(plotpath, 'Sample_distance.pdf'), height = 6, width = 9)
Heatmap(
    sampleDists,
    name = 'Sample \nDistances',
    show_column_dend = F, show_column_names = F
)    
dev.off()


# Sample correlation heatmap ----

mat <- assay(vsd)

samplecorr <- cor(mat, method = 'spearman')


samplenames<- colData(vsd) %>%
    as.data.frame() %>%
    pull(Sample_Type)
rownames(samplecorr) <- samplenames
colnames(samplecorr) <- samplenames

pdf(file.path(plotpath, 'Sample_corr.pdf'), height = 6, width = 9)
Heatmap(
    samplecorr,
    name = 'Sample \nCorrelation',
    show_column_dend = F, show_column_names = F
)    
dev.off()


