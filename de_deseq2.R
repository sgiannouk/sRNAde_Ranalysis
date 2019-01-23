### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 5) p-value 6) path/to/packrat_directory
### Example Command: Rscript DESeq2.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj 0.05 /projects/packrat


# Input arguments and error control
args <- commandArgs(TRUE)
if (length(args) == 4) {
  if (!file.exists(args[1])) {
    cat("ERROR - The input matrix does NOT exist...\nEXITING!\n")
    quit()
  }
  matfile <- read.delim(args[1], header=TRUE, row.names=1)   # Input the input delimited text file containing the count matrix
  groups <- unlist(strsplit(args[2], ","))  # Sample description
  sampletypevalues <- rev(unique(groups))  # Getting the group levels
  if (!dir.exists(args[3])) {
    cat("ERROR - The output directory does NOT exist...\nEXITING!\n")
    quit()
  }
  outdir <- args[3]  # Output directory
  basename <- args[4]  # Base name
  pvalue <- 0.05  # Default differential expression cutoff
  # if (dir.exists("/opt/sRNAtoolboxDB/packrat")) {
  #   packrat_path <- "/opt/sRNAtoolboxDB/packrat"  # Alu/Epigenoma
  # } else {
  #   cat("ERROR - The packrat directory could not be found...\nEXITING!\n")
  #   quit()
  # }
} else if (length(args) == 5) {
  if (!file.exists(args[1])) {
    cat("ERROR - The input matrix does NOT exist...\nEXITING!\n")
    quit()
  }
  matfile <- read.delim(args[1], header=TRUE, row.names=1)   # Input the input delimited text file containing the count matrix
  groups <- unlist(strsplit(args[2], ","))  # Sample description
  sampletypevalues <- rev(unique(groups))  # Getting the group levels
  if (!dir.exists(args[3])) {
    cat("ERROR - The output directory does NOT exist...\nEXITING!\n")
    quit()
  }
  outdir <- args[3]  # Output directory
  basename <- args[4]  # Base name
  pvalue <- as.numeric(args[5])  # Differential expression cutoff
  # packrat_path <- args[6]  # Path of packrat environment
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- rev(unique(groups))
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# pvalue <- 0.05

# Initiating packrat environment and 
# loading all necessary libraries
# library("packrat")
# packrat::init(packrat_path)
library("vsn")
library("DESeq2")
library("ggplot2")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/deseq2_", basename, ".log", sep=""), append = TRUE)

print("DESeq2 Differential Expression analysis of RNA-Seq expression profiles is now running...")

# Designing the data's factors which indicate the experimental group for each sample
samplefactors <- data.frame(row.names=colnames(matfile), condition = factor(groups, levels=sampletypevalues))
# Constructing the data set object that DESeq2 needs to perform the analysis
dds <- DESeqDataSetFromMatrix(countData = matfile, colData = samplefactors, design = ~condition)
# Calling DESeq method 
dds <- DESeq(dds, fitType="local")
# Obtaining the 'normalised' matrix
normalised_counts <- counts(dds, normalized=TRUE)

# Exporting the 'normalised' table
print(paste("Exporting DESeq2 normalised table containing all genes to: ", outdir,"/",basename,"_deseq2_normTable.csv", sep=""))
write.table(data.frame("name"=rownames(normalised_counts), normalised_counts), file=paste(outdir,"/",basename,"_deseq2_RLEnormTable.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Plotting the per-gene dispersion estimates together with the fitted mean-dispersion relationship.
print("Plotting the estimated dispersion...") 
png(paste(outdir,"/",basename,"_deseq2_dispersionPlot.png",sep=""), units='px', height=900, width=1600, res=90)
plotDispEsts(dds, xlab = "Mean of normalized counts", ylab = "Dispersion")
title(main = "Dispersion Estimates")
dev.off()

# Standard deviation of the raw data, across samples, against the mean
print("Plotting the mean stdv of the raw data...")
png(paste(outdir,"/",basename,"_deseq2_meanSD.png",sep=""), units='px', height=900, width=900, res=90)
meanSdPlot <- meanSdPlot(log2(counts(dds, normalized=TRUE) + 1))
meanSdPlot$gg + ggtitle("Standard deviation against the mean") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


for(i in 1:(length(sampletypevalues)-1)) {
    for(j in (i+1):length(sampletypevalues)){ 

        # Extracting the table containing the test results for each gene in the original count table
        results <- results(dds, contrast=c("condition", sampletypevalues[i], sampletypevalues[j]))
        
        if(sum(is.na(results))>0) {
          print("WARNING: missing values (NA) appear in DESeq2 result, probably due to 0 counts in both samples")
          print("Rows with NA values will be omitted...")}    
        
        # Selecting the samples
        selected_samples <- (which(groups==sampletypevalues[i] | groups==sampletypevalues[j]))
        # Obtaining the list of genes with absolute log2FoldChange greater than 1 and p adjusted value lower than the user-input value
        selected <- row.names(results)[(abs(results$log2FoldChange)>=1) & (results$padj<=pvalue)]
        # Removing rows with missing values on columns
        selected <- na.omit(selected)
        # Combing the normalised data along with statistical analysis results ("log2FoldChange", "pval", "padj")
        ncounts_selected <- as.data.frame(cbind(normalised_counts[ , selected_samples], results$log2FoldChange, results$pvalue, results$padj))
        # Naming the new columns
        colnames(ncounts_selected) <- c(head(colnames(ncounts_selected), n=-3), "log2FoldChange", "pvalue", "padj")
        
        # Generating a matrix containing the mean normalised results per group per (ALL) gene
        mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(ncounts_selected[ ,which(groups==sampletypevalues[j])]), 
                                                     mean_Bgroup=rowMeans(ncounts_selected[ ,which(groups==sampletypevalues[i])]),
                                                     ncounts_selected[,c("log2FoldChange", "pvalue", "padj")]))
        colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[j],sep=""), paste("mean_",sampletypevalues[i],sep=""))
        # Insert FoldChange calculations
        mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[,1:2],
                                                     transform(mean_ncounts_selected[0], FoldChange = (mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                                     mean_ncounts_selected[,3:5]))
        
        ncounts_selected <- as.data.frame(cbind(normalised_counts[ , selected_samples], mean_ncounts_selected$FoldChange, results$log2FoldChange, results$pvalue, results$padj))
        # Naming the new columns
        colnames(ncounts_selected) <- c(head(colnames(ncounts_selected), n=-4), "FoldChange", "log2FoldChange", "pvalue", "padj")
        # Order dataframe based on the padj
        mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj),]
        
        # Obtaining the final matrix of selected genes
        result <- ncounts_selected[selected, ]
        
        # Exporting the normalised results table containing the selected genes with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
        print(paste("Exporting DESeq2 normalised table containing ONLY selected genes (adjusted P value < ", pvalue, ") to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_topGenesBelow", gsub("[.]", "", pvalue), ".csv", sep=""))
        write.table(data.frame("name"=rownames(result), result), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_topGenesBelow", gsub("[.]", "", pvalue), ".csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        # Exporting the normalised results table containing ALL genes
        print(paste("Exporting DESeq2 normalised table containing all genes to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_DEseq2_allGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(ncounts_selected), ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        # Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
        print(paste("Exporting DESeq2 mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_meanGroupsAllGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        
        # Expression plot
        expression_plot <- function (res, thresh=pvalue, labelsig=TRUE, ...) {
          plot(log2(res[ ,1]+1),log2(res[ ,2]+1), pch = 20, cex = 0.5, col = 1,
               xlab = colnames(res)[1],
               ylab = colnames(res)[2],
               main="Average expression values of each condition\n(differentially expressed features are highlighted)")
          gn.selected <- abs(res$log2FoldChange)>=1 & res$padj<=pvalue
          points(log2(res[ ,1]+1)[gn.selected], log2(res[ ,2]+1)[gn.selected], col="red", pch=20, cex=.8)
        }
        print("Average expression values of each condition and highlight the features declared as differentially expressed...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_expressionPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        expression_plot(mean_ncounts_selected)
        dev.off()
        
        # MA-plot
        # Plotting the log2 fold changes against the mean normalised counts, colouring in red those 
        # genes that are significant
        print("Plotting the log2 fold changes against the mean normalised counts...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[i],"VS",sampletypevalues[j],"_deseq2_maPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        plotMA(results, alpha=pvalue, cex=.5, main="MAplot")
        abline(h=c(-1,1), col="dodgerblue", lwd=2)
        dev.off()
        
        # Volcano plot
        # Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
        resultDESeq2 <- ncounts_selected[order(ncounts_selected$padj), ]
        resultDESeq2$sig <- (-log10(resultDESeq2$padj))
        resultDESeq2[is.infinite(resultDESeq2$sig),"sig"] <- 350
        # Select genes with a defined p-value (DESeq2 assigns NA to some genes)
        genes.to.plot <- !is.na(resultDESeq2$pvalue)
        ## Volcano plot of adjusted p-values
        cols <- densCols(resultDESeq2$log2FoldChange, resultDESeq2$sig)
        cols[resultDESeq2$pvalue==0] <- "purple"
        resultDESeq2$pch <- 19
        resultDESeq2$pch[resultDESeq2$padj==0] <- 6
        print("Generating the volcano plot...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_deseq2_volcanoPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        plot(resultDESeq2$log2FoldChange,
             resultDESeq2$sig,
             col=cols, panel.first=grid(),
             main="Volcano plot",
             xlab="Effect size: log2(fold-change)",
             ylab="-log10(adjusted p-value)",
             pch=resultDESeq2$pch, cex=0.4)
        abline(v=0)
        abline(v=c(-1,1), col="brown")
        abline(h=-log10(pvalue), col="brown")
        ## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
        gn.selected <- (abs(resultDESeq2$log2FoldChange)>=1 & resultDESeq2$padj<=pvalue)
        text(resultDESeq2$log2FoldChange[gn.selected],
             resultDESeq2$sig[gn.selected],
             lab=rownames(resultDESeq2)[gn.selected ], cex=0.6)
        dev.off()
        
    }
}

## DESeq2 DE analysis finished ##
print("DESeq2 Differential Expression analysis has finished successfully!")

### DESeq2 analysis is now being terminated...
print(" ---------------------------------------------- DESeq2 ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()