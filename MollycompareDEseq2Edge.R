##  Trying to run DESeq and EdgeR.
##  Running it on 2 different tables generated of the same data.
#   Then do a Volcano plo, PCA and p-value table using all the data!!!!

rm(list = ls())

#source("http://www.bioconductor.org/biocLite.R") # This is old and full of terrors
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

install.packages('NMF')
install.packages("htmltools")
BiocManager::install("tximport")
BiocManager::install('EnhancedVolcano')
library("tximport")
library(htmltools)
library(DESeq2)
library(ggplot2)
library(stringr)
library(EnhancedVolcano)

setwd("/Users/mollydavis333/Desktop/Internship/Tomato_Meta")
#dir="/Users/robreid/Dropbox (UNC Charlotte)/r/tomato/seedling-pollen/quants"
filemeta="/Users/mollydavis333/Desktop/Internship/Tomato_Meta/metadata.csv"
filecounts= "/Users/mollydavis333/Desktop/Internship/Tomato_Meta/counts.txt"
filesalmon= "/Users/mollydavis333/Desktop/Internship/Tomato_Meta/salmon.merged.gene_counts.tsv"
filesorel= "/Users/mollydavis333/Desktop/Internship/Tomato_Meta/BAM_Files/loop1-counts.csv"


dat <- read.table(file.path(filecounts), header=TRUE)
dat

dat <- read.table(file.path(filesalmon), header=TRUE)
dat

dat <- read.csv(file.path(filesorel), header=TRUE)
dat

## Remove the gene desc column but save it for later
#save the gene desc
genedesc <- subset(dat,select=c(2))
head(genedesc)

dat <- subset(dat,select=-c(2))
head(dat)


#Data for salmon Heinz
geneID<- dat[,1]
geneID
dat<- dat[,2:7]
dat<-round(dat)
dat
dat<- cbind(geneID, dat)
dat

##Data for Heiz Hot pollen counts
Heinzdat<- dat[,-(2:5)]
Heinzdat<- Heinzdat[,-(8:25)]
Heinzdat

##Data for Malintka
Malintkadt<- dat[,-(2:5)]
Malintkadt
Malintkadt<- Malintkadt[,-(2:7)]
Malintkadt
Malintkadt<-Malintkadt[,0:7]
Malintkadt

##Data for Nagcarlang
NagcarlangDat<- dat[,-(2:5)]
NagcarlangDat
NagcarlangDat<- NagcarlangDat[,-(2:13)]
NagcarlangDat
NagcarlangDat<-NagcarlangDat[,0:7]
NagcarlangDat

##Data for Tamualipas
Tamualipas<- dat[,-(2:5)]
Tamualipas
Tamualipas<- Tamualipas[,-(2:19)]
Tamualipas


## Metadata for Heinz
metaData <- read.csv(file = filemeta, header = TRUE,strip.white=TRUE,stringsAsFactors=FALSE)
metaData<- data.frame(metaData)
HeinzMeta<-metaData[1:6,]
HeinzMeta

##DESeq for Heinz unfiltered
dds <- DESeqDataSetFromMatrix(countData=dat, colData=HeinzMeta, design=~Condition, tidy = TRUE)
dds
#ddsZero <- dds[ rowSums(counts(dds)) > 0, ] #removes rows with zero counts
#dds <- DESeq(ddsZero)
dds <- DESeq(dds)
#Cloudy equal bad 
plotDispEsts(dds)
#needs to run results function to get mean p-values from DESeq analysis
resNoFilt <- results(dds, independentFiltering=FALSE) #turn off independent filtering from results() function
summary(resNoFilt)



##DESeq for Heinz filtered
library("genefilter")
dds <- DESeqDataSetFromMatrix(countData=dat, colData=HeinzMeta, design=~Condition, tidy = TRUE)
dds
ddsZero <- dds[ rowSums(counts(dds)) > 0, ] #removes rows with zero counts
dds <- estimateSizeFactors(ddsZero)
tmp = varFilter(counts(dds), var.func=IQR, var.cutoff=0.80, filterByQuantile=TRUE)
#idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
ddsvar <- DESeqDataSetFromMatrix(tmp,colData = HeinzMeta, design =~ Condition)
ddsZero <- ddsvar[rowSums(counts(ddsvar)) > 0, ]  ## Remove rows that = 0, aka no reads
dds <- DESeq(ddsZero)

#Cloudy equal bad 
plotDispEsts(dds)
res <- results(dds,cooksCutoff=FALSE, independentFiltering=FALSE) #results() function independently filters the counts < 0.1 (alpha) 
summary(res)

##Significant Genes with pvalue & padj <= 0.05 No filter
ResSigpval<- resNoFilt$pvalue <= 0.05
length(ResSigpval[ResSigpval== TRUE])

ResSigpadj<- resNoFilt$padj <= 0.05 
length(ResSigpadj[ResSigpadj== TRUE])


##Significant Genes with padj & pvalue <= 0.05 with filter
ResSigpval<- res$pvalue <= 0.05
length(ResSigpval[ResSigpval== TRUE])

ResSigpadj<- res$padj <= 0.05 
length(ResSigpadj[ResSigpadj== TRUE])



##DESeq top 200 Varying Genes
library("DESeq2")
library("RColorBrewer")
library( "genefilter" )
#Regularized log transformation
rld <- rlog( dds, fitType='mean', blind=TRUE)
#Get 200 top varying genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 200)
#make a subset of the log transformed counts for just the top 200 varying genes
top200Counts<-assay(rld)[topVarGenes,]
write.csv(top200Counts, file="top200counts_SalmonIQR.rld.csv", quote=FALSE)

##DESeq Heat Map
#install.packages("pheatmap")
library("pheatmap")
#Read in the regularized log transformed counts for the top 100 varying genes
top100Counts<-read.table("top200counts_hotpollen.rld.csv", sep=",",header=TRUE, row.names=1)

#PLOT HEATMAP USING PHEATMAP
#INCLUDE NEXT LINE IF YOU WANT TO SAVE THE FIGURE IN A FILE
pdf(file="gene_hotpollen_200.heatmap.pdf")
#heatmap using top 25 varying genes
pheatmap(top100Counts,scale = "row")
#INCLUDE NEXT LINE IF YOU WANT TO SAVE THE FIGURE IN A FILE
dev.off()







##EdgeR for Heinz
# load EdgeR library
#BiocManager::install("edgeR")
library(edgeR)

dat <- read.csv(file.path(filesorel), header = TRUE)
dat

## Remove the gene desc column but save it for later
#save the gene desc
genedesc <- subset(dat,select=c(2))
head(genedesc)

dat <- subset(dat,select=-c(2))
head(dat)
##Data for Heiz
Heinzdat<- dat[,-(2:5)]
Heinzdat<- Heinzdat[,-(8:25)]
Heinzdat

#Data for salmon Heinz
geneID<- dat[,1]
geneID
dat<- dat[,2:7]
dat<-round(dat)
dat
dat<- cbind(geneID, dat)
dat


dds <- DESeqDataSetFromMatrix(countData=Heinzdat, colData=HeinzMeta, design=~Condition, tidy = TRUE)
dds

#Into the Edge object.
rownames(Heinzdat) <- Heinzdat[, 1]
Heinzdat[, 1] <- NULL
Heinzdat
d<-DGEList(counts=Heinzdat,group=factor(dds$Condition))
dim(d)

#Backup
d.full <- d # keep the old one in case we mess up
head(d$counts)

apply(d$counts, 2, sum) # total gene counts per sample before filter


##Use if you want to filter 
keep <- rowSums(cpm(d)>1) >= 2  #Keep only rows with data #cpm() only done once. This keeps those genes which have minimum cpm of 1 in at least 2 samples.
d <- d[keep,]
dim(d)
apply(d$counts, 2, sum) # total gene counts per sample after filter

#This reduces the dataset. For the filtered tags,
#there is very little power to detect differential expression, so little information is lost by filtering.
#After filtering, it is a good idea to reset the library sizes:
d$samples$lib.size <- colSums(d$counts)
d$samples

####   edgeR is concerned with differential expression analysis rather than with the quantification of expression levels.
####   It is concerned with relative changes in expression levels between conditions,
####   but not directly with estimating absolute expression levels.

#
#   The calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors
#   for the library sizes that minimize the log-fold changes between the samples for most genes.
#    The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
#    We call the product of the original library size and the scaling factor the effective library size.
#   The effective library size replaces the original library size in all downsteam analyses.
#
d <- calcNormFactors(d)


#Quick plot to compare covarying  ((  multidimensional scaling ))
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
#estimate the dispersion parameter for each tag
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)

## plotBCV() plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
plotBCV(d1)

## GLM estimates of dispersion
# First, you must fit the common dispersion. Then you need to fit a trended model
#  we can also estimate a generalized linear model (glm) fit using edgeR. In the same way that we've been doing above,
#  we will just add these as additional data to the object we've been working with.
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="bin.loess")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

#Differential Expression
#Once the dispersions are estimated, we can proceed with testing procedures for determining differential expression.
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
SigEdgePVal<-et12$table$PValue <= 0.05
length(SigEdgePVal[SigEdgePVal== TRUE])

topTags<- topTags(et12, n=200)
topTags


##Print Summary of EdgeR Results
#BVC plots
pdf("BCV_Heinz.pdf")
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
plotBCV(d1)
plotBCV(d2)
dev.off()

#Differential Expression with P-values
head(topTags)
outname=paste("heinz", "EdgeR_top200_sorel.csv", sep = "_")
head(outname)
write.table(topTags, file=outname, sep='\t')



###Venn Diagram of data from all counts and top 100 genes each 
#Installs
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

install.packages("VennDiagram")
install.packages("venn") 

##DESeq
library(VennDiagram)
DESeq<- read.csv(file.path('/Users/mollydavis333/Desktop/Internship/Tomato_Meta/Top_200_genes_IQR_DESeq\ .csv'), header = TRUE)
DESeq
venn.diagram(DESeq, main = "Top 200 Genes from IQR DESeq Analysis", sub= "Author: Molly Davis (mdavi258@uncc.edu)" ,fill = c("red", "blue", "yellow"), filename = "venn200-DESeq_ IQR-dimensions.png")

#Show the results for each intersection 
library(gplots)
salmon<- DESeq[,1]
Sorel<- DESeq[,2]
HotPollen<- DESeq[,3]
lst<- list(salmon,Sorel,HotPollen)
lst
names(lst)<- c('Salmon', 'Sorel', 'HotPollen')
ItemsList <- venn(lst, show.plot = FALSE)
Deseq_results<-attributes(ItemsList)$intersections
Deseq_results
#Save the results in file
capture.output(Deseq_results, file = "DESeq_top200_venn__IQR_results.txt")


##EdgeR
library(VennDiagram)
EdgeR<- read.csv(file.path('/Users/mollydavis333/Desktop/Internship/Tomato_Meta/Top_200_genes_EdgeR.csv'), header = TRUE)
EdgeR
venn.diagram(EdgeR,main = "Top 200 Genes from EdgeR Analysis", sub= "Author: Molly Davis (mdavi258@uncc.edu)", fill = c("#999999", "#E69F00", "#56B4E9"), filename = "venn200-EdgeR-dimensions.png")

#Show the results for each intersection 
salmonR<- EdgeR[,1]
SorelR<- EdgeR[,2]
HotPollenR<- EdgeR[,3]
lst<- list(salmonR,SorelR,HotPollenR)
lst
names(lst)<- c('Salmon', 'Sorel', 'HotPollen')
ItemsList <- venn(lst, show.plot = FALSE)
EdgeR_results<-attributes(ItemsList)$intersections
EdgeR_results
#Save the results in file
capture.output(EdgeR_results, file = "EdgeR_top200_venn_results.txt")



###Compare both DESeq and EdgeR in Venn Diagram






###################################################



#Get just the data we want from the larger set.
exp1=grep("H(1|2|3)-C",samples[,1],value=TRUE)
head(exp1)
exp2=grep("H(1|2|3)-S",samples[,1],value=TRUE)
exp2
#changed to define experiments verbatim (did not use this,#)
#exp1<-c("N1-C", "N2-C", "N3-C")
#exp2<-c("N1-S", "N2-S", "N3-S")
#### Get sub samples
sub1 <- str_subset(samples$id,exp1)
sub2 <- str_subset(samples$id,exp2)
#samp <- append(sub2,sub1)   #(this will reverse the direction of a volcano plot, LF directions)
samp <- append(exp1,exp2) # by doing them in this order a + fc means exp2 is greater than exp1
head(samp)


#changing this line to samples$id instead of samp brings in all samples; not just those identified by samp (above)
#files <- file.path(dir,paste(samp,"_quant.sf",sep = ""))
#files <- file.path(dir,paste(samples$id,"_quant.sf",sep = ""))
#head(files)
#tx_salmon = tximport(files, type = "salmon", ignoreTxVersion=T,txOut=TRUE)
#head(tx_salmon)
#data are now in place in R as an object

### Need a metafile
metaall <- read.csv(file = metafile, header = TRUE,strip.white=TRUE,stringsAsFactors=FALSE)
meta1 <- metaall[grep("H(1|2|3)-C", metaall$id),]
meta2 <- metaall[grep("H(1|2|3)-S", metaall$id),]
metboth <- rbind(meta1,meta2)
head(metboth)

## Make the Deseq Object
##   design = ~Strain + Time
# for this first run, we're simply comparing condition within variety
# Rob adds pairing - "And then this line gets changed in the code":
tx_sdds = DESeqDataSetFromTximport(tx_salmon, colData = metboth, design = ~Variety + Condition)
dds <- DESeqDataSetFromMatrix(countData=dat, colData=metaData, design=~Variety, tidy = TRUE)
dds
ddsZero <- dds[ rowSums(counts(dds)) > 0, ]
#tx_sdds = DESeqDataSetFromTximport(tx_salmon, colData = metboth, design = ~Condition)
data4edge = tx_sdds
#this is the original, not paired version:
#tx_sdds = DESeqDataSetFromTximport(tx_salmon, colData = metboth, design =~ Condition)
#tx_sdds = DESeqDataSetFromTximport(tx_salmon, colData = metaall, design =~Variety + Condition)
#this line can be used to run with all samples

txddsZero <- tx_sdds[rowSums(counts(tx_sdds)) > 0, ]  ## Remove rows that = 0, aka no reads
#lhc <- 40 (could use light harvesting complex as a lower minimum)
#txddsLHC <- tx_sdds[rowSums(counts(tx_sdds)) > lhc, ]  ## Remove rows that < light harvest complex genes

#the following is to have a look at the data in a different way - need to install Glimpse package
#glimpse(txddsZero)
#str(tx_sdds) # a way to peak at # of rows
#str(txddsZero@rowRanges)


#### Save the NA values - by seeing the minreplicatesforreplace to Inf, the NA will stay!!!
####  Now we are keeping every P value we can of P-adjust and pvalue.
dds <- DESeq(txddsZero, minReplicatesForReplace=Inf)
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)

#Will only look at data that have some sort of reads
dds <- DESeq(txddsZero)
res <- results(dds)
summary(res)

library("genefilter")
# varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
#tmp = varFilter(counts(tx_sdds), var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

tx_sdds = DESeqDataSetFromTximport(tx_salmon, colData = metboth, design =~ Variety + Condition)
txddsZero <- tx_sdds[rowSums(counts(tx_sdds)) > 0, ]
dds <- estimateSizeFactors(txddsZero)
tmp = varFilter(counts(dds), var.func=IQR, var.cutoff=0.80, filterByQuantile=TRUE)
#idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
ddsvar <- DESeqDataSetFromMatrix(tmp,colData = metboth, design =~ Variety + Condition)
txddsZero <- ddsvar[rowSums(counts(ddsvar)) > 0, ]  ## Remove rows that = 0, aka no reads
dds <- DESeq(txddsZero)
res <- results(dds,cooksCutoff=FALSE, independentFiltering=FALSE)
summary(res)
dim(tmp)



#Print Summary of comparison
res <- res[order(res$padj),]  ## Reordering the table by adjusted Pvalue
head(res)
outname=paste("heinz", "padjust-iqr05_rob.txt", sep = "_")
head(outname)
write.table(res, file=outname, sep='\t')

### PCA
#Need rld
rld <- rlog( dds )

#Principal component plot of the samples
png(file="basic-PCA-nag.png",width=600, height=350)
plotPCA(rld, intgroup=c("Condition"))
dev.off()

#PCA with shapes
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Variety"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
outpng=paste("heinz", "PCA-CookOff.png", sep = "_")
png(file=outpng,width=600, height=550)
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Variety)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + scale_colour_brewer(palette="Set2") +
  coord_fixed()
dev.off()
#End of PCA shapes





### VOLCANO PLOT
#Enhanced Volcano PLots

## Get the min and max values for LOg fold change, get the max for Pvalues!
lfcmin <- min(res$log2FoldChange)
lfcmax <- max(res$log2FoldChange)
padjmax <- -log10(min(res$padj))

outname=paste("heinz","Volcano", sep = "-")

#R hex color
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
EnhancedVolcano(res,
                lab = rownames(res),
                #selectLab = c('Solyc06g036290.3'),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 5e-2,
                col = colorBlindGrey8,
                xlim = c((lfcmin-1), lfcmax+2),
                ylim = c(0,padjmax),
                ylab = "neg. LOG of Adjusted Pvalue",
                title = outname)
outname=paste(exp1,exp2,"log2FoldChange-cookoff", sep = "_")
outpng=paste(exp1,exp2,"enhancevolcano-cookoff-padj.png", sep = "_")
png(file=outpng,width=600, height=550)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 5e-2,
                selectLab = c('TMEM176B','ADH1A'),
                col = colorBlindGrey8,
                xlim = c((lfcmin-1), lfcmax+2),
                ylim = c(0,padjmax),
                ylab = "neg. LOG of Adjusted Pvalue",
                title = outname)
dev.off()





#######  EDGER   Run on same data
## Based on this site:   https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
## Comparing the models in DESeq and edgeR
## DESeq always only uses a gamma glm as its model. Since edgeR does not have gamma glm as an option,
## we cannot produce the same glm results in edgeR as we can in DESeq and vice versa.

# load EdgeR library
library(edgeR)

#INto the Edge object.
d <- DGEList(counts=tx_salmon$counts,group=factor(data4edge$Condition))
dim(d)

#backup
d.full <- d # keep the old one in case we mess up
head(d$counts)

##Removing low counts via CPM (counts per million)
head(cpm(d))

apply(d$counts, 2, sum) # total gene counts per sample before filter

keep <- rowSums(cpm(d)>100) >= 2  #Keep only rows with data!
d <- d[keep,]
dim(d)
apply(d$counts, 2, sum) # total gene counts per sample after filter

#This reduces the dataset. For the filtered tags,
#there is very little power to detect differential expression, so little information is lost by filtering.
#After filtering, it is a good idea to reset the library sizes:
d$samples$lib.size <- colSums(d$counts)
d$samples

####   edgeR is concerned with differential expression analysis rather than with the quantification of expression levels.
####   It is concerned with relative changes in expression levels between conditions,
####   but not directly with estimating absolute expression levels.

#
#   The calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors
#   for the library sizes that minimize the log-fold changes between the samples for most genes.
#    The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
#    We call the product of the original library size and the scaling factor the effective library size.
#   The effective library size replaces the original library size in all downsteam analyses.
#
d
d <- calcNormFactors(d)
d

#Quick plot to compare covarying  ((  multidimensional scaling ))
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

#estimate the dispersion parameter for each tag
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)

## plotBCV() plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
plotBCV(d1)

## GLM estimates of dispersion
# First, you must fit the common dispersion. Then you need to fit a trended model
#  we can also estimate a generalized linear model (glm) fit using edgeR. In the same way that we've been doing above,
#  we will just add these as additional data to the object we've been working with.
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="bin.loess")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

#Differential Expression
#Once the dispersions are estimated, we can proceed with testing procedures for determining differential expression.
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
topTags(et12, n=25)



