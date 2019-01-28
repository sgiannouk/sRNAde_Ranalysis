### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 5) p-value 6) path/to/packrat_directory
### Example Command: Rscript edgeR.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj 0.05 /projects/packrat

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
# names(matfile) <- gsub(x = names(matfile), pattern = "\\.", replacement = " ")
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- rev(unique(groups))
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# pvalue <- 0.05

# Initiating packrat environment and 
# loading all necessary libraries
# library("packrat")
# packrat::init(packrat_path)
library("edgeR")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/edger_", basename, ".log", sep=""), append = TRUE)

print("EdgeR Differential Expression analysis of RNA-Seq expression profiles is now running...")
# Storing the raw read counts table in a simple list-based data object called a DGEList.
edgeR_table <- DGEList(counts=matfile, group=factor(groups))  # From now on, all the necessary info will be stored in this variable
# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most genes. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
edgeR_table <- calcNormFactors(edgeR_table)
# The first major step in the analysis of DGE  data using the NB model is to estimate the dispersion
# parameter for each tag, a measure of the degree of inter-library variation for that tag. 
# Estimating the common dispersion gives an idea of overall variability across the genome for this dataset.
edgeR_table <- estimateCommonDisp(edgeR_table)
# For routine differential expression analysis, we use empirical Bayes tagwise dispersions. 
# Note that common dispersion needs to be estimated before estimating tagwise dispersions.
edgeR_table <- estimateTagwiseDisp(edgeR_table)

# Exporting the TMM normalised table
print(paste("Exporting edgeR normalised table containing all genes to: ", outdir,"/",basename,"_edger_TMMnormTable.csv", sep=""))
write.table(data.frame("name"=rownames(edgeR_table$pseudo.counts), edgeR_table$pseudo.counts), file=paste(outdir,"/",basename,"_edger_TMMnormTable.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


for(i in 1:(length(sampletypevalues)-1)) {
    for(j in (i+1):length(sampletypevalues)){
        
        print(paste("Performing an exact test between:", sampletypevalues[i], "and", sampletypevalues[j]))
        # Performing an exact test for the difference in expression between each pair of conditions. 
        # More explicitly, computing gene-wise exact tests for differences in the means between two 
        # groups of negative-binomially distributed counts.
        dgeTest <- exactTest(edgeR_table, pair=c(sampletypevalues[i], sampletypevalues[j]))
        # The test results for the most significant tags are conveniently in the following function.
        topTags <- topTags(dgeTest, n=Inf)
        # Selecting the samples
        selected_samples <- (which(groups==sampletypevalues[i] | groups==sampletypevalues[j]))
        # Obtaining the pseudo-counts for each sample
        z <- as.data.frame(edgeR_table$pseudo.counts)[ ,selected_samples]
        # Obtaining the statistical results from the edge test
        t <- as.data.frame(topTags)
        # Merging the above info into one table
        data <- as.data.frame(merge(z, t, by="row.names"))
        row.names(data) <- data$Row.names  # row names manipulation
        data$Row.names <- NULL  # row names manipulation
        # Combing the normalised data along with statistical analysis results ("log2FoldChange", "pval", "padj")
        data <- as.data.frame(cbind(data[ , selected_samples], data$logFC, data$PValue, data$FDR))
        # Naming the new columns
        colnames(data) <- c(head(colnames(data), n=-3), "log2FoldChange", "pvalue", "padj")
        # Selecting only genes with padj lower than the input p-value
        selected <- which(data$padj<=pvalue) 
        
        # Generating a matrix containing the mean normalised results per group per (ALL) gene
        mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(data[ ,which(groups==sampletypevalues[j])]), 
                                                     mean_Bgroup=rowMeans(data[ ,which(groups==sampletypevalues[i])]),
                                                     data[,c("log2FoldChange", "pvalue", "padj")]))
        colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[j],sep=""), paste("mean_",sampletypevalues[i],sep=""))
        # Inserting FoldChange calculations
        mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[,1:2],
                                                     transform(mean_ncounts_selected[0], FoldChange = (mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                                     mean_ncounts_selected[,3:5]))
        # Inserting FoldChange calculations
        data <- as.data.frame(cbind(data[ ,selected_samples], mean_ncounts_selected$FoldChange, data$log2FoldChange, data$pvalue, data$padj))
        # Naming the new columns
        colnames(data) <- c(head(colnames(data), n=-4), "FoldChange", "log2FoldChange", "pvalue", "padj")
        # Order dataframe based on the padj
        mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj), ]
        
        # Obtaining the final matrix of selected genes
        result <- data[selected, ] 
        
        # Exporting the normalised results table containing the selected genes below the chosen threshold
        print(paste("Exporting edgeR normalised table containing ONLY selected genes (adjusted P value <", pvalue, ") to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_topGenesBelow", gsub("[.]", "", pvalue), ".csv", sep=""))
        write.table(data.frame("name"=rownames(result), result), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_topGenesBelow", gsub("[.]", "", pvalue), ".csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        # Exporting the normalised results table containing ALL genes
        print(paste("Exporting edgeR normalised table containing all genes to: ", outdir,"/",basename,"_edger_allGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(data), data), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        # Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
        print(paste("Exporting edgeR mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_meanGroupsAllGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        
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
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_expressionPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        expression_plot(mean_ncounts_selected)
        dev.off()
        
        # MA-plot
        # Plotting the log2 fold changes against the mean normalised counts, colouring in red those
        # genes that are significant
        mean_ncounts_selected <- as.data.frame(cbind(mean=rowMeans(mean_ncounts_selected[ ,c(1,2)]), mean_ncounts_selected))
        maplot <- function (res, thresh=pvalue, labelsig=TRUE, ...) {
          with(res, plot(mean, log2FoldChange, pch=20, cex=.5, log="x", ...))
          with(subset(res, padj<=thresh), points(mean, log2FoldChange, col="red", pch=20, cex=.8))
        }
        print("Plotting the log2 fold changes against the mean normalised counts...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_maPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        maplot(mean_ncounts_selected, 
               xlab = "mean of normalized counts", 
               ylab = "log2 fold change", 
               main="MAplot")
        abline(h=c(-1,1), col="dodgerblue", lwd=2)
        abline(h=0, col="red3", lwd=3)
        dev.off()
        
        # Volcano plot
        # Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
        resultedgeR <- data[order(data$padj), ]
        resultedgeR$sig <- (-log10(resultedgeR$padj))
        resultedgeR[is.infinite(resultedgeR$sig),"sig"] <- 350
        # Select genes with a defined padj (edgeR assigns NA to some genes)
        genes.to.plot <- !is.na(resultedgeR$pvalue)
        ## Volcano plot of adjusted padj
        cols <- densCols(resultedgeR$log2FoldChange, resultedgeR$sig)
        cols[resultedgeR$pvalue == 0] <- "purple"
        resultedgeR$pch <- 19
        resultedgeR$pch[resultedgeR$padj == 0] <- 6
        print("Generating the volcano plot...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_edger_volcanoPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        plot(resultedgeR$log2FoldChange,
             resultedgeR$sig,
             col=cols, panel.first=grid(),
             main="Volcano plot",
             xlab="Effect size: log2(fold-change)",
             ylab="-log10(adjusted p-value)",
             pch=resultedgeR$pch, cex=0.4)
        abline(v=0)
        abline(v=c(-1,1), col="brown")
        abline(h=-log10(pvalue), col="brown")
        ## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
        gn.selected <- abs(resultedgeR$log2FoldChange)>=1 & resultedgeR$padj<=pvalue
        text(resultedgeR$log2FoldChange[gn.selected],
             resultedgeR$sig[gn.selected],
             lab=rownames(resultedgeR)[gn.selected ], cex=0.6)
        dev.off()
    }
}

## edgeR DE analysis finished ##
print("edgeR Differential Expression analysis has finished successfully!")

### edgeR analysis is now being terminated...
print(" ---------------------------------------------- edgeR ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()