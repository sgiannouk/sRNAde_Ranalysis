### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 5) NOISeq probab. value 6) path/to/packrat_directory
### Example Command: Rscript NOISeq.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj 0.8 /projects/packrat


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
  noiseq_prob <- 0.8  # Default differential expression cutoff
  # if (dir.exists("/opt/sRNAtoolboxDB/packrat")) {
  #   packrat_path <- "/opt/sRNAtoolboxDB/packrat"  # Alu/Epigenoma
  # } else {
  #   cat("ERROR - The packrat directory could not be found...\nEXITING!\n")
  #   quit()
  # }
} else if (length(args) == 6) {
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
  noiseq_prob <- as.numeric(args[5])  # Differential expression cutoff
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
# noiseq_prob <- 0.8

# Initiating packrat environment and
# loading all necessary libraries
# library("packrat")
# packrat::init(packrat_path)
library("NOISeq")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/noiseq_", basename, ".log", sep=""), append = TRUE)

print("NOISeq Differential Expression analysis of RNA-Seq expression profiles is now running...")
  
# Designing the data's factors which indicate the experimental group for each sample
samplefactors <- data.frame(row.names=colnames(matfile), condition = factor(groups, levels=sampletypevalues))
# Importing all necessary information into a NOISeq object
countdata <- readData(data = matfile, factors = samplefactors)
# Trimmed Mean of M values (TMM) normalisation to correct the sequencing depth bias
# No length is taken into account
TMMvalues = tmm(assayData(countdata)$exprs, long = 1000, lc = 0)

# Exporting the TMM normalised table
print(paste("Exporting NOISeq TMM normalised table containing all genes to: ", outdir,"/",basename,"_noiseq_TMMnormTable.csv", sep=""))
write.table(data.frame("name"=rownames(TMMvalues), TMMvalues), file=paste(outdir,"/",basename,"_noiseq_TMMnormTable.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


for(i in 1:(length(sampletypevalues)-1)) {
    for(j in (i+1):length(sampletypevalues)){ 
      
        # Computing the differential expression between experimental conditions from the filtered read count data
        NOISeq <- noiseq(countdata, k=0.5, lc=0, norm="tmm", factor="condition", conditions=c(sampletypevalues[i], sampletypevalues[j]))

        # Extract the table containing the test results for each gene of the original count table
        curresult <- NOISeq@results[[1]]
        
        if(sum(is.na(curresult))>0){
          print("WARNING: NA values appear in NOISeq result, probably due to 0 counts in both samples")
          print("Rows with NA values will be omitted...")}
        
        # Selecting the samples
        selected_samples <- (which(groups==sampletypevalues[i] | groups==sampletypevalues[j]))
        # Obtaining the list of genes with probability higher than the user-input value (default 0.8)
        selected <- rownames(degenes(NOISeq, q = noiseq_prob, M = NULL))
          
        selected <- na.omit(selected)
        # Combing the normalised data along with statistical analysis results ("M", "prob", "padj")
        TMMvalues_selected <- as.data.frame(cbind(TMMvalues[, selected_samples], curresult$M, curresult$prob, (1-curresult$prob)))
        # Naming the new columns
        colnames(TMMvalues_selected) <- c(head(colnames(TMMvalues_selected),n=-3), "log2FoldChange", "prob", "padj")
        
        # Generating a matrix containing the mean normalised results per group per (ALL) gene
        mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(TMMvalues_selected[ ,which(groups==sampletypevalues[j])]), 
                                                     mean_Bgroup=rowMeans(TMMvalues_selected[ ,which(groups==sampletypevalues[i])]),
                                                     TMMvalues_selected[,c("log2FoldChange", "prob", "padj")]))
        colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[j],sep=""), paste("mean_",sampletypevalues[i],sep=""))
        # Insert FoldChange calculations
        mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[,1:2],
                                                     transform(mean_ncounts_selected[0], FoldChange = (mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                                     mean_ncounts_selected[,3:5]))
        
        # Combing the normalised data along with statistical analysis results ("log2FoldChange", "prob", "padj")
        TMMvalues_selected <- as.data.frame(cbind(TMMvalues[, selected_samples], mean_ncounts_selected$FoldChange, curresult$M, curresult$prob, (1-curresult$prob)))
        # Naming the new columns
        colnames(TMMvalues_selected) <- c(head(colnames(TMMvalues_selected),n=-4), "FoldChange", "log2FoldChange", "prob", "padj")
        # Order dataframe based on the FDR
        mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj), ]
        
        # Obtaining the final matrix of selected genes
        result <- TMMvalues_selected[selected, ]
        
        # Exporting the normalised results table containing the selected genes below the chosen threshold
        print(paste("Exporting NOISeq normalised table containing ONLY selected genes (NOISeq probability >", noiseq_prob, ") to: ",outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_topGenesBelow", gsub("[.]", "", noiseq_prob), ".csv", sep=""))
        write.table(data.frame("name"=rownames(result), result), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_topGenesBelow", gsub("[.]", "", noiseq_prob), ".csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        # Exporting the normalised results table containing ALL genes
        print(paste("Exporting NOISeq normalised table containing all genes to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_allGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(TMMvalues_selected), TMMvalues_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        # Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
        print(paste("Exporting NOISeq mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_meanGroupsAllGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        
        # Expression plot
        print("Average expression values of each condition and highlight the features declared as differentially expressed...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_expressionPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        DE.plot(NOISeq, q = noiseq_prob, graphic = "expr", log.scale = TRUE, main="Average expression values of each condition\n(differentially expressed features are highlighted)")
        dev.off()
        
        # MA-plot
        # Plotting the log2 fold changes against the mean normalised counts, colouring in red those 
        # genes that are significant
        mean_ncounts_selected <- as.data.frame(cbind(mean=rowMeans(mean_ncounts_selected[ ,c(1,2)]), mean_ncounts_selected))
        maplot <- function (res, thresh=noiseq_prob, labelsig=TRUE, ...) {
          with(res, plot(mean, log2FoldChange, pch=20, cex=.5, log="x", ...))
          with(subset(res, prob>=thresh), points(mean, log2FoldChange, col="red", pch=20, cex=.8))
        }
        print("Plotting the log2 fold changes against the mean normalised counts...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_MAPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        maplot(mean_ncounts_selected, 
               xlab = "mean of normalized counts", 
               ylab = "log2 fold change", 
               main="MAplot")
        abline(h=c(-1,1), col="dodgerblue", lwd=2)
        abline(h=0, col="red3", lwd=3)
        dev.off()

        # Volcano plot
        # Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
        resultNOISeq <- TMMvalues_selected[order(TMMvalues_selected$prob), ]
        resultNOISeq$sig <- (-log10(resultNOISeq$prob))
        resultNOISeq[is.infinite(resultNOISeq$sig),"sig"] <- 350
        # Select genes with a defined p-value (DESeq2 assigns NA to some genes)
        genes.to.plot <- !is.na(resultNOISeq$prob)
        ## Volcano plot of adjusted p-values
        cols <- densCols(resultNOISeq$log2FoldChange, resultNOISeq$sig)
        cols[resultNOISeq$prob ==0] <- "purple"
        resultNOISeq$pch <- 19
        resultNOISeq$pch[resultNOISeq$prob == 0] <- 6
        print("Generating the volcano plot...")
        png(paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_volcanoPlot.png",sep=""), units='px', height=900, width=1600, res=100)
        plot(resultNOISeq$log2FoldChange,
             resultNOISeq$sig,
             col=cols, panel.first=grid(),
             main="Volcano plot",
             xlab="Effect size: log2(fold-change)",
             ylab="-log10(probability)",
             pch=resultNOISeq$pch, cex=0.4)
        abline(v=0)
        abline(v=c(-1,1), col="brown")
        abline(h=-log10(noiseq_prob), col="brown")
        ## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
        gn.selected <- resultNOISeq$prob>=noiseq_prob
        text(resultNOISeq$log2FoldChange[gn.selected],
             resultNOISeq$sig[gn.selected],
             lab=rownames(resultNOISeq)[gn.selected ], cex=0.6)
        dev.off()

    }
} 

## NOISeq DE analysis finished ##
print("NOISeq Differential Expression analysis has finished successfully!")

### NOISeq analysis is now being terminated...
print(" ---------------------------------------------- NOISeq ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()
