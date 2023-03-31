library( "DESeq2" )
library(ggplot2)
setwd("/t1-data/user/aangeles/fastaq_files/DESEQ")
countsTable<-read.delim("name_output_expression_matrix_full.out", header=T, row.names=1)
countData<-matrix(countsTable)
sampleInfo<-read.csv("DESEQ_samples_flipped.csv", header=T)
attach(sampleInfo)
meta<-data.frame(row.names=colnames(countsTable), type=factor(condition))
dds<-DESeqDataSetFromMatrix(countData=countsTable,meta,formula(~type))
#removal of all entries with 0 counts
idx <- which(rowSums(counts(dds)) > 0)
dds <- dds[idx,]
dds <- estimateSizeFactors(dds)
temp2<-counts(dds, normalized=TRUE)
write.table(temp2, file="normalized_counts_reduced.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

###Principal Component Analysis
###The following code was used:
###GGFORTIFY Method (actual using) (taken from http://rstudio-pubs-static.s3.amazonaws.com/53162_cd16ee63c24747459ccd180f69f07810.html; https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html)

vsd <- vst(dds, blind=FALSE)
pcaResult<-prcomp(t(assay(vsd)))
summary(pcaResult)
write.table(pcaResult$x, file="pca_results_SW.txt", sep="\t", quote=F, col.names = TRUE, row.names = TRUE)
head(pcaResult$x)
pdf("PCA_DESeq_ggplot2.pdf")
library(ggfortify)
autoplot(pcaResult, data = meta, colour = 'type')
dev.off()

###Heatmap Analysis (taken from https://informatics.fas.harvard.edu/differential-expression-with-deseq2-mouse-immune-cells.html)
###The following code was used:

rld <- rlogTransformation(dds, blind=TRUE)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor) 
pheatmap(rld_cor)
###Generate MA Plot
dds2 <- DESeq(dds)
result <- results(dds2, contrast=c('type', 'iPSC', 'iPSC-neuron'))
head(result)
plotMA(result, alpha = .05, main=paste0('Condition: iPSC vs. iPSC-neurons'), ylim=c(-5,5))
List of top 250 up-regulated and top 250 down-regulated genes, sorted by p-value (taken from https://informatics.fas.harvard.edu/differential-expression-with-deseq2-mouse-immune-cells.html)
n = 250 
resOrdered <- result[order(result$padj),]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,],
                     +                      resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]  #print results for top and bottom 5 genes
write.table(topResults, file="topResults_neuron.txt", sep="\t", quote=F, col.names = TRUE, row.names = TRUE)

###Export DESEQ2 results to table:

write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
(from https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exporting-results-to-csv-files)

###Plot Analysis of Single Genes
###The following code was used:

d<-plotCounts(dds, gene = "ENSG00000141837.19", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,250,400,1000,1300)) +
  ggtitle("CACNA1A") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000148408.12", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1000,3000)) +
  ggtitle("CACNA1B") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000151067", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1000)) +
  ggtitle("CACNA1C") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000157388.16", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNA1D") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000198216.11", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1200)) +
  ggtitle("CACNA1E") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000102001.12", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNA1F") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000006283.17", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNA1G") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000196557.12", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1000,2000)) +
  ggtitle("CACNA1H") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000100346.17", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNA1I") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000081248.10", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(1,5,25,100,400)) +
  ggtitle("CACNA1S") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000153956.15", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,800,1600,3200,4000,6000)) +
  ggtitle("CACNA2D1") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000007402.11", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(5,25,100,400,800,1200,1600,2000,5000)) +
  ggtitle("CACNA2D2") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000157445.14", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1000)) +
  ggtitle("CACNA2D3") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000151062.14", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNA2D4") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000067191.15", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1000,2000)) +
  ggtitle("CACNB1") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000165995.19", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNB2") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000167535.7", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1000,2000)) +
  ggtitle("CACNB3") + 
  theme(plot.title = element_text(hjust = 0.5))

d<-plotCounts(dds, gene = "ENSG00000182389.19", intgroup='type', returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("CACNB4") + 
  theme(plot.title = element_text(hjust = 0.5))