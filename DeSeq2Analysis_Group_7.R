##########################################################################################
# This script performs differential gene expression analysis with DESeq2 #
##########################################################################################

### LOAD REQUIRED LIBRARIES
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("vsn")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("genefilter")
library("biomaRt")
library("IHW")
library("ggplot2")
library("apeglm")
library("ashr")

### SET WORKING DIRECTORY
# note: this directory should be populated with the raw counts file
setwd("C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs")

### Import count table and details on experimental design
# NB: Make sure column names in the sample(table) file and counts file are exactly the same and in the same order
CountTable <- read.table("counts.txt", header=TRUE, row.names=1) 
CountTable <- subset(CountTable, select = -c(Chr, Start, End, Strand, Length))
names(CountTable) <- substring(names(CountTable), 8, 17)

samples <- read.table("samples.txt", header=TRUE)
Dataset <- DESeqDataSetFromMatrix(countData = CountTable, colData=samples, design=~Condition) # incorporates batch information (if available)

### PRELIMINARY ANALYSES ###
# The first steps in your analysis should focus on better understanding the relationship of the datasets being studied. 
# This can be simply achieved by generating a PCA plot showing the relationship of your samples.
# First we transform our raw count data using a variance stabilizing transformation (VST) that roughly mirrors how DeSeq2 models the data.
vsd1 <- varianceStabilizingTransformation(Dataset, blind=FALSE)

# Then we plot a PCA, grouping and coloring our datasets according to batch
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_PCA.png")
plotPCA(vsd1, "Condition")
dev.off()

### note that you can attach additional information based on the column headers in your sample table
# plotPCA(vsd1, c("Condition","Batch")) 
    # looks like we don't need because we don't have batch info

# we can also attempt to replicate the batch effect correction performed by DeSeq2 using the limma::removeBatchEffect function
vsd2 <- varianceStabilizingTransformation(Dataset, blind=FALSE)
assay(vsd2) <- limma::removeBatchEffect(assay(vsd2), vsd2$batch)

png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_PCA_2.png")
plotPCA(vsd2, "Condition")
dev.off()

### We can also calculate and plot sample distances using either the batch corrected (vsd2) or uncorrected (vsd1) data. 

# uncorrected
sampleDists <- dist( t( assay(vsd1) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_DistMat.png")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main = "Feature Counts Sample Distance Matrix")
dev.off()

# corrected
sampleDistsCorr <- dist( t( assay(vsd2) ) )
sampleDistsCorr
sampleDistCorrMatrix <- as.matrix( sampleDistsCorr )
colnames(sampleDistCorrMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_DistMat_corrected.png")
pheatmap(sampleDistCorrMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main = "Feature Counts Sample Distance Matrix (Corrected)")
dev.off()

# At this stage, you should have a good sense of how your samples cluster and the effect of batch correction (if used)
# In simple RNA-Seq situations (control vs treatment, 3-5 bioreps each), this is all that should be required. 
# For more complex situations, you will need to dive deep into working of DeSeq2
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


### BASIC DGE ANALYSIS USING DESEQ2 ###

Dataset <- DESeqDataSetFromMatrix(countData = CountTable, colData=samples, design=~Condition) # incorporates batch information (if available)

# Run DESEQ and generate a simple plot showing the distribution of regulated and unregulated genes
DatasetProcessed <- DESeq(Dataset) # runs DESEQ
par(mfrow=c(1,1))
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_MeanCounts.png")
plotMA(DatasetProcessed, main="Feature Counts DESeq2", ylim=c(-5,5))
dev.off()

# Next we perform a contrast analysis to produce a list of differentially regulated genes between our two conditions
# First we set CTRL dataset as baseline
Dataset$Condition <- relevel(Dataset$Condition, "Ctrl")

# Next we create our results object while performing shrinkage of effect size 
# (this reduces the impact of apparent gross changes in low expressed genes)
resultsNames(DatasetProcessed)
res1 <- lfcShrink(DatasetProcessed, contrast=c("Condition","dsDNA","Ctrl"), type = "ashr")

# Here we modify our output data to include two additional columns that contain the baseMeans (a proxy for counts)
# This is useful for downstream filtering of lowly expressed genes
baseMeanCtrl = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$Condition == "Ctrl"])
baseMeandsDNA = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$Condition == "dsDNA"])
res1 = cbind(as.data.frame(res1), baseMeanCtrl, baseMeandsDNA)

# Here we add two further columns, the gene symbol (common name) and entrez ID - both of which may be useful downstream
res1$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="SYMBOL", keytype="SYMBOL", multiVals="first") # MAPS GENE IDs
res1$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="ENTREZID", keytype="SYMBOL", multiVals="first")

# Finally we write the complete results object to an outfile
write.csv(res1, "DGEanalysis.csv", row.names=TRUE)

############## PLOTTING ############## 

### VOLCANO PLOTS
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_volcano.png")
Dataset$Condition <- relevel(Dataset$Condition, "Ctrl")
res1 <- lfcShrink(DatasetProcessed, contrast=c("Condition","dsDNA","Ctrl"), type = "ashr")
with(res1, plot(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, main="Feature Counts Volcano Plot", xlim=c(-5,5)))
with(subset(res1, padj>0.01), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="gray"))
with(subset(res1, padj<0.01 & log2FoldChange> 2), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="red"))
with(subset(res1, padj<0.01 & log2FoldChange< -2), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="blue"))
### ADD BELLS AND WHISTLES
pval = 0.01
abline(h = -log10(pval), col = "black", lty = 2, lwd=4)
mtext(paste("pval = 0.01", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
dev.off()

### Plot heatmap - sorting by padj and log2fc
# PICK ALL GENES WITH pADJ < 0.01 AND THEN SUBSET FOR THOSE WITH Log2FC > 1 THEN PICK TOP 25 HITS
subset <- (subset(res1, res1$padj < 0.01 & abs(res1$log2FoldChange) > 1))
subset <- head(subset[order(-subset$log2FoldChange),],25)
sigGenes <- rownames(subset)
rows <- match(sigGenes, row.names(vsd1))
mat <- assay(vsd1)[rows,]
mat<-mat-rowMeans(mat) ##plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
df <- as.data.frame(colData(DatasetProcessed))

png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/03_Crash_Course_in_Gene_Expression/Assignments/bp_outputs/Figs/FeatCounts_heatmap.png")
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, main = "Feature Counts Heatmap")
dev.off()
