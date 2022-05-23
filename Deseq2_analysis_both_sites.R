library(DESeq2)
library(splines)
library(dplyr)
library(gplots)
library(ggplot2)
library(pheatmap)


# Read the data with raw counts from the table.
countData <- read.table("./Pele_raw_counts_2_genes.all.tsv", header=TRUE, sep="\t", row.name=1)
countData[1:24] <- lapply(countData[1:24], as.integer)

# design table
design <- data.frame(Point = colnames(countData),
                     Site = c(rep('head', 12), rep('tail', 12)),
                     TPA = rep(c(0, 0, 4, 4, 12, 12, 24, 24, 48, 48, 96, 96), 2),
                     Replicates = rep(c(1:12), each = 2)
                     )

# Time Spline creation
sites_matTimeSplineBasis <- ns(design$TPA, df=4)
colnames(sites_matTimeSplineBasis) <- paste0("spline", seq(1, dim(sites_matTimeSplineBasis)[2]))
design <- cbind(design, sites_matTimeSplineBasis)

# factor variables
design$Site <- factor(design$Site)
design$Replicates <- factor(design$Replicates)
design$TPA <- factor(design$TPA)

# Deseq dataset creation
dds <- DESeqDataSetFromMatrix(countData, 
                              colData = design, 
                              design = ~Site + Site:spline1 + Site:spline2 + Site:spline3 + Site:spline4)
# filter genes with low expression level
dds <- dds[rowSums(counts(dds)) > 10, ]

# manually run deseq
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds)
dds <- estimateDispersionsMAP(dds, outlierSD = 10)

dds <- nbinomLRT(dds, full=~Site + Site:spline1 + Site:spline2 + Site:spline3 + Site:spline4,
                                    reduced = ~spline1 + spline2 + spline3 + spline4)

# Exploratory data analysis
## PCA
vsd <- vst(dds, blind = F)
# rld <- rlog(dds, blind = F) # recommended for small data
pcaData <- plotPCA(vsd, intgroup=c("Site", "TPA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=TPA, shape=Site)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

## Correlation plot
rld_cor <- cor(assay(vsd))
pheatmap(rld_cor)

# results
res <- results(dds[which(mcols(dds)$fullBetaConv),])

# Results transformation for visualization as heatmap
sorted <- res[with(res, order(padj, -log2FoldChange)), ]
sorted.df <- data.frame("id"=rownames(sorted),sorted)
selected_genes <- sorted.df %>% 
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) >= 0.58) %>%
  #slice_max(order_by = abs(log2FoldChange) + padj, n = 50) %>% 
  select(id)
selected_genes <- selected_genes$id

norm_cnt <- counts(dds, normalized=TRUE)
dt <- data.frame("id"=rownames(norm_cnt), norm_cnt)
short_dt <- dt[match(selected_genes, dt$id), ]


# heatmap
gene <- short_dt[,1]
vals <- as.matrix(short_dt[, 2:25])
vals <- jitter(vals, factor = 1, amount=0.00001)
score <- NULL
for (i in 1:nrow(vals)) {
  row <- vals[i,]
  zscore <- (row-mean(row))/sd(row)
  score <- rbind(score,zscore)
}
row.names(score) <- gene
zscore <- score
mat <- as.matrix(zscore)
colors <- colorRampPalette(c("blue4","cornsilk2","red"),space="rgb")(256)
heatmap.2(mat, col=colors, density.info="none", trace="none", 
          margins=c(8,8), lhei=c(1,5), Colv=FALSE, cexRow=0.3)


# # heatmap for log2(foldchange against time 0)
# ntd <- normTransform(dds)
# ntd_dt <- data.frame("id"=rownames(assay(ntd)), assay(ntd))
# short_ntd <- ntd_dt[match(selected_genes, ntd_dt$id),]
# gene = short_ntd[,1]
# vals = as.matrix(short_ntd[,14:25])
# vals = jitter(vals, factor = 1, amount=0.00001)
# score = NULL
# for (i in 1:nrow(vals)) {
#   row=vals[i,]
#   point0 = mean(c(row[1], row[2]))
#   zscore=row-point0
#   score =rbind(score,zscore)
# }
# row.names(score) = gene
# zscore=score
# mat = as.matrix(zscore)
# colors = colorRampPalette(c("blue4","cornsilk2","red"),space="rgb")(256)
# heatmap.2(mat,col=colors,density.info="none",trace="none", margins=c(8,8),lhei=c(1,5), Colv=FALSE, cexRow=0.3)
# dev.off()
