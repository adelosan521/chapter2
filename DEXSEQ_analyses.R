##This code uses the DEXSeq library to generate normalized counts (for analysis of differential exon usage (DEU)) and plotting of DEU for individual genes.

## Load the DEXSeq library
library(DEXSeq)

## Set working directory
setwd("/t1-data/user/aangeles/fastaq_files/revisedcounts2")

## List all count files with specified pattern and store names in "countFiles"
countFiles = list.files("/t1-data/user/aangeles/fastaq_files/revisedcounts2", pattern = "CACNA1C.counts.txt$", full.names = TRUE)

## Print the base names of the countFiles
basename(countFiles)

## Set the path for the flattened annotation file
flattenedFile = "/t1-data/user/aangeles/fastaq_files/new_annotation.gff"

## Read the CSV file containing sample information
sampleTable=read.csv("samples_revised3.csv", row.names=1)

## Display first few rows of the sampleTable
head(sampleTable)

## Create a DEXSeqDataSet object using the count files and sample information
dxd = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design=~sample + exon + condition:exon, flattenedfile=flattenedFile)

## Show dimensions of DEXSeqDataSet
dim(dxd)

## Estimate size factors for the DEXSeqDataSet object
dxd=estimateSizeFactors(dxd)

##Generate normalized counts for DEU file

test<-counts(dxd, normalized=T)
write.table(test, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/normalized_counts_DEXSeq.txt", sep="\t",  row.names = TRUE, col.names = TRUE)

##Make table for VGCC subunit DEU and plotting of individual genes. First code section is for CACNA1C (ENSG00000151067) but can be adapted accordingly for the other gens.

##Grep only the gene of interest (example ENSG00000151067) from the DEXSeqDataSet object "dxd")
dxd2<-dxd[grep("ENSG00000151067", rownames(dxd)), ]

## Display dimensions of DEXSeqDataSet object "dxd2"
dim(dxd2)

## Estimate dispersion for dxd2
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))

## Test for differential exon usage in dxd2
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))

##Estimate exon fold changes in 'dxd2' with respect to the 'condition' variable
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))

##Extract DEXseqResults from "dxd2" and store in "drx1"
drx1=DEXSeqResults(dxd2)

## Display the first few rows of drx1
head(drx1)

## Write the drx1 results to a text file
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000151067.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

## Create a PDF file to store the DEXSeq plot
pdf("CACNA1C_DEXSeq.pdf")

## Generate the DEXSeq plot for the gene of interest
plotDEXSeq(drx1, "ENSG00000151067", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

## Save drx1 as an Rdata file
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

## Repeat of lines 40 - 72 for other genes indicated below (for example line 76 - 87 is for CACNA1D (ENSG00000157388))

dxd2<-dxd[grep("ENSG00000157388", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000157388.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1D_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000157388.16", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000102001", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000102001.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1F_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000102001.12", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000007402", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000007402.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA2D2_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000007402.11", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000198216", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000198216.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1E_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000198216.11", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000157445", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000157445.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA2D3_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000157445.14", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000153956", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000153956.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA2D1_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000153956.15", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000148408", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000148408.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1B_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000148408.12", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000151062", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000151062.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA2D4_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000151062.14", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000167535", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000167535.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNB3_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000167535.7", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000196557", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000196557.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1H_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000196557.12", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000067191", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000067191.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNB1_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000067191.15", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000006283", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000006283.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1G_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000006283.17", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000141837", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000141837.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1A_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000141837.19", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()

dxd2<-dxd[grep("ENSG00000100346", rownames(dxd)), ]
dim(dxd2)
dxd2=estimateDispersions(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=testForDEU(dxd2, BPPARAM=MulticoreParam(workers=16))
dxd2=estimateExonFoldChanges(dxd2, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))
drx1=DEXSeqResults(dxd2)
head(drx1)
write.table(drx1, file="/t1-data/user/aangeles/fastaq_files/revisedcounts2/test_DEU_ENSG00000100346.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
pdf("CACNA1I_DEXSeq.pdf")
plotDEXSeq(drx1, "ENSG00000100346.17", splicing = TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
save(drx1, file = "/t1-data/user/aangeles/fastaq_files/revisedcounts2/iPSC_Neuron.Rdata")
dev.off()
