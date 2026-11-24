library(SummarizedExperiment)
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(DESeq2)
library(vsn)
library(ggrepel)
library(SingleCellExperiment)
library(org.Hs.eg.db)

setwd("~/Desktop/260604_JurkatRNASeq/Sequencing_analysis")
plotpath <- 'Plots/'
dir.create(plotpath, showWarnings = F)

# To-Do list



# Assemble SE ----

# Read data
data <- read.delim('Data/MXMVMV-expression-matrix.tsv', row.names = 1)


# Create count matrix
colnamestoget <- grep('count',colnames(data), value = T)

counts <- data[,colnamestoget]
counts <- round(counts)
colnames(counts) <- gsub('_count','',colnames(counts))

# Create RowData
rowdata <- DataFrame(
    row.names = rownames(data),
    gene_symbol = data$gene_name,
    genetype = data$gene_biotype
)


# Create ColData
coldata <- read.csv('Data/ColData.csv', row.names = 1)
coldata$Line <- gsub(' ','_', coldata$Line)
coldata$Sample_Type <- str_c(str_sub(coldata$Line, 5), gsub('_hours','',coldata$Condition), sep = '_')
coldata$Sample_name <- str_c(str_sub(coldata$Line, 5), gsub('_hours','',coldata$Condition), coldata$Replicate, sep = '_')
coldata <- coldata[colnames(counts),]


all.equal(colnames(counts), rownames(coldata))
all.equal(rownames(counts), rownames(rowdata))


# Create SummarizedExperiment
se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    colData = coldata,
    rowData = rowdata
)

head(rowData(se))
head(colData(se))
head(assay(se))

# Remove low expression
'HAVCR2' %in% rowData(se_filtered)$gene_symbol

se_filtered <- se[rowMeans(assay(se, "counts")) > 5, ]
dim(se)
dim(se_filtered)
head(colData(se_filtered))

saveRDS(se_filtered, 'Data/summarizedexp.rds')

# Switch to DESeqDataSet
se$libSize <- colSums(assay(se))
se <- readRDS('Data/summarizedexp.rds')
dds <- DESeq2::DESeqDataSet(se, design = ~ Condition + Line)
saveRDS(dds, 'Data/DESeqDataset.rds')

# Transform data to account for increased variance with increased count, giving normalized assay
vsd <- DESeq2::vst(dds, blind = TRUE)
saveRDS(vsd, 'Data/Normalized_se.rds')

