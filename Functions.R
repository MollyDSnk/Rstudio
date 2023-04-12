library(stringr)


makeExperimentMetaDataFrame = function(fname="../results/muday-144-SL5_counts-salmon.txt") {
  experiment_data = read.csv(fname,
                             sep = "\t",
                             stringsAsFactors = F, 
                             header=T, 
                             check.names = FALSE, 
                             row.names = "gene_name")
  sample_names = colnames(experiment_data)
  boolean_vector = str_detect(sample_names,'[VFA]\\.\\d\\d\\.\\d\\d\\.\\d')
  sample_names = sample_names[boolean_vector]
  genotype = sapply(strsplit(sample_names,"\\."),function(x){x[[1]]})
  time = sapply(strsplit(sample_names,"\\."),function(x){x[[3]]})
  temperature = sapply(strsplit(sample_names,"\\."),function(x){x[[2]]})
  to_return = data.frame(genotype,time,temperature)
  rownames(to_return) = sample_names
  return(to_return)
}



# Creating and running a Function
Gentoype_DE_Analysis <- function(data, sampleData, CSV, title) {
  cts <- as.matrix(data) # Count matrix data: the function DESeqDataSetFromMatrix can be used if you already have a matrix of read counts prepared from another source. 
  coldata <- sampleData
  coldata <- coldata[,c("genotype", "time", "temperature")]
  coldata$genotype <- factor(coldata$genotype, levels = c("A", "F", "V"))
  coldata$time <- factor(coldata$time, levels = c("15", "30", "45", "75"))
  coldata$temperature <- factor(coldata$temperature, levels = c("28", "34"))
  cts <- cts[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts)) # Columns of the count matrix and the rows of the column data (information about samples) are in the same order
  
  # DESeq dataset creation : No pre-filtering or cooks cutoff was used
  dds <- DESeqDataSetFromMatrix(countData = round(cts), # Needed to round to get rid of decimals for DESeq to work
                                colData = coldata,
                                design = ~ temperature)
  head(dds)
  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  mcols(dds)
  
  # Running DESeq
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  
  # group the factors creates interactions to design the formula
  #dds$group <- factor(paste0(dds$genotype,dds$time, dds$temperature))
  #design(dds) <- ~ group
  #dds <- DESeq(dds)
  
  #Pre-Filtering
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  resultsNames(dds)
  # By default the argument alpha is set to 0.1. If the adjusted p value cutoff will be a value other than 0.1,alpha should be set to that value:
  res05 <- results(dds, alpha=0.05)
  print(res05)
  summary(res05)
  
  numSignGenes<- sum(res05$pvalue < 0.05, na.rm=TRUE) # Not adjusted
  
  print("Number of Significant Genes with pvalue < 0.05:")
  print(numSignGenes)
  
  # A p-value adjustment is necessary when one performs multiple comparisons or multiple testing in a more general sense: performing multiple tests of significance where only one significant result will lead to the rejection of an overall hypothesis.
  numSignGenes<- sum(res05$padj < 0.05, na.rm=TRUE) # Adjusted
  
  print("Number of Significant Genes with pvalue-Adjusted < 0.05:")
  print(numSignGenes) 
  
  resOrdered <- res05[order(res05$pvalue),]
  write.csv(as.data.frame(resOrdered), 
            file=CSV)
  
  #VSD Transformed values 
  vsd <- vst(dds, blind=FALSE)
  #rld <- rlog(dds, blind=FALSE) # Used for heat maps 
  
  #par(mar=c(8,5,2,2))
  #boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
  
  #plotMA(res05)
  #plotDispEsts(dds)
  
  
  #R hex color
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  lfcmin <- min(res05$log2FoldChange,na.rm=TRUE)
  lfcmax <- max(res05$log2FoldChange,na.rm=TRUE)
  #padjmax <- max(res50$pvalue,na.rm=TRUE)
  padjmax <- -log10(min(res05$pvalue,na.rm=TRUE))
  # create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
  # this can be achieved with nested ifelse statements
  keyvals <- ifelse(
    res05$log2FoldChange < -1, 'royalblue',
    ifelse(res05$log2FoldChange > 1, 'orange',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'orange'] <- 'Up-regulated'
  names(keyvals)[keyvals == 'black'] <- 'Mid-regulated'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'
  EnhancedVolcano(res05,
                  lab = rownames(res05),
                  #selectLab = c('Solyc06g036290.3'),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  selectLab = rownames(res05)[which(names(keyvals) %in% c('Up-regulated', 'Down-regulated'))],
                  colCustom = keyvals,
                  colAlpha = 1,
                  shape = c(4, 1, 6, 3),
                  pCutoff = 0.05, # will show up as 1.4 on y-axis but it 0.05 significance. 
                  FCcutoff = 1.0,
                  col = colorBlindGrey8,
                  xlim = c((lfcmin-1), lfcmax+1),
                  ylim = c(0,padjmax),
                  ylab = "Neg. LOG of Pvalue",
                  title = title, 
                  legendPosition = "right")
}






# Creating and running a Function
PCA_plots <- function(data, sampleData, title) {
  cts <- as.matrix(data) # Count matrix data: the function DESeqDataSetFromMatrix can be used if you already have a matrix of read counts prepared from another source. 
  coldata <- sampleData
  coldata <- coldata[,c("genotype", "time", "temperature")]
  coldata$genotype <- factor(coldata$genotype, levels = c("A", "F", "V"))
  coldata$time <- factor(coldata$time, levels = c("15", "30", "45", "75"))
  coldata$temperature <- factor(coldata$temperature, levels = c("28", "34"))
  cts <- cts[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts)) # Columns of the count matrix and the rows of the column data (information about samples) are in the same order
  
  # DESeq dataset creation : No pre-filtering or cooks cutoff was used
  dds <- DESeqDataSetFromMatrix(countData = round(cts), # Needed to round to get rid of decimals for DESeq to work
                                colData = coldata,
                                design = ~ time + temperature)
  head(dds)
  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  mcols(dds)
  
  # Running DESeq
  dds <- DESeq(dds)
  res <- results(dds)
  summary(res)
  
  # group the factors creates interactions to design the formula
  #dds$group <- factor(paste0(dds$genotype,dds$time, dds$temperature))
  #design(dds) <- ~ group
  #dds <- DESeq(dds)
  #resultsNames(dds)
  
  
  # How many adjusted p-values were less than 0.1
  sum(res$padj < 0.1, na.rm=TRUE)
  # By default the argument alpha is set to 0.1. If the adjusted p value cutoff will be a value other than 0.1,alpha should be set to that value:
  res05 <- results(dds, alpha=0.05)
  summary(res05)
  # How many adjusted p-values were less than 0.05
  sum(res05$padj < 0.05, na.rm=TRUE)
  
  #VSD Transformed values 
  vsd <- vst(dds, blind=FALSE)
  head(vsd)
  
  #rld <- rlog(dds, blind=FALSE) # Used for heat maps 
  
  #PCA Plots for all 3 genotypes at control and stress temperatures in each time point (so, this would be 4 PCA plots).
  pcaData <- plotPCA(vsd, intgroup=c("temperature", "time"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=temperature, shape=time)) +
    geom_point() +
    scale_color_manual(values =  c("28" = "blue",
                                   "34" = "orange")) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(title) +
    coord_fixed()
}