### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 5) percentage/top_n_genes 6) path/to/packrat_directory
### Example Command: Rscript visualisation.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj 0.1 /projects/packrat

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
  n_top <- 0.1  # # Cutoff for plotting the 10% most variable genes
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
  n_top <- as.numeric(args[5])  # Cutoff for plotting the top_n most variable genes
  # packrat_path <- args[6]  # Path of packrat environment
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# packrat_path <- "/Users/stavris/R/packrat"
# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- rev(unique(groups))
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# n_top <- 0.125
# n_top <- 25

# Initiating packrat environment and 
# loading all necessary libraries
# library("packrat")
# packrat::init(packrat_path)
library("edgeR")
library("heatmaply")
library("RColorBrewer")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/visualisation_", basename, ".log", sep=""), append = TRUE)

print("Visualisation of the RNA-seq experiment is now running...")
data <- DGEList(counts=matfile, group=factor(groups))  # Summarise the input data

check.integer <- function(N){
  !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}
# Getting log2 of cmp data, 2 is being added to raw data
log2counts <- cpm(data$counts, prior.count=2, log=TRUE)

# Colouring the different conditions
if (length(sampletypevalues) == 2) {
  col_condition <- c("#7FC97F", "#BEAED4")[data$samples$group]
} else {
  col_condition <- c(brewer.pal(n = length(sampletypevalues), name = "Accent"))[data$samples$group]
}

# Examine the distributions of the raw counts by plotting the log2CPM of the counts
print("Checking the distribution of the read counts on the log2 scale...")
png(paste(outdir,"/",basename,"_visualisation_Log2DistPlot.png",sep=""), units='px', height=900, width=1600, res=90)
# Check distributions of samples using boxplots
par(mar=c(8.1, 4.1, 4.1, 2.1))
boxplot(cpm(data$counts, prior.count=2, log=TRUE),col=col_condition, xlab="", ylab="Log2 counts per million", las=2)
# Adding a blue horizontal line that corresponds to the median log2CPM
abline(h=median(log2counts), col="slategrey", lwd=2)
title("Boxplots of log2CPMs (unnormalised)")
dev.off()

# An MDSplot is a visualisation of a principle components analysis, which determines the greatest sources of variation in the data. 
# If the experiment is well controlled and has worked well, what we hope to see is that the greatest sources of variation are the 
# treatments/groups we are interested in. It is also an incredibly useful tool for quality control and checking for outliers.
print("Creating the MultiDimensional Scaling plot...") 
png(paste(outdir,"/",basename,"_visualisation_plotMDS.png",sep=""), units='px', height=900, width=1600, res=90)
plotMDS(data, col=col_condition)
title(main = "MultiDimensional Scaling plot\n(distances approximate the log2 fold changes between the samples)")
dev.off()

# Estimate the variance for each row in the log2counts matrix
var_genes <- apply(log2counts, 1, var)


if (check.integer(n_top)) {
  # Obtaining the percentage of top most variable genes
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:n_top]
} else {
  # Obtaining the n top most variable genes
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:(nrow(log2counts)*n_top)]
}

# Subset log2counts matrix
highly_variable <- log2counts[select_var, ]

# Creating the heatmaps. If samples are more than 100, then sample labels are being omitted
if (length(groups) < 100) {
	# Heatmap of top selected normalised genes
	heatmaply(highly_variable, file = paste(outdir,"/",basename,"_heatmap_normalLog2CPM.html",sep =""),
          limits = NULL, xlab = "Samples", ylab = "Genes", colors = brewer.pal(11,"Spectral"), scale = "row",
	        main = paste("Top",nrow(highly_variable),"most variable genes across samples",sep=" "), key.title=NULL, 
	        col_side_colors = data.frame(groups), hide_colorbar = FALSE, column_text_angle=60, fontsize_col = 9,
	        fontsize_row = 8,  showticklabels=c(TRUE,TRUE))
  h <- heatmaply(highly_variable, #file = paste(outdir,"/",basename,"_heatmap_normalLog2CPM.png",sep =""),
          limits = NULL, xlab = "Samples", ylab = "Genes", colors = brewer.pal(11,"Spectral"), scale = "row",
          main = paste("Top",nrow(highly_variable),"most variable genes across samples",sep=" "), key.title=NULL, 
          col_side_colors = data.frame(groups), hide_colorbar = FALSE, column_text_angle=60, fontsize_col = 9,
          fontsize_row = 8,  showticklabels=c(TRUE,TRUE))
  h$width <- 1200
  h$height <- 800
  export(h, file = paste(outdir,"/",basename,"_heatmap_normalLog2CPM.png",sep =""))
  
  } else {
	# Heatmap of top selected genes, No sample-names 
	heatmaply(highly_variable, # file = paste(outdir,"/",basename,"_heatmap_normalLog2CPM.html",sep =""),
	          limits = NULL, xlab = "Samples", ylab = "Genes", colors = brewer.pal(11, "Spectral"), scale = "row",
	          main = paste("Top",nrow(highly_variable),"most variable genes across samples",sep=" "), key.title=NULL, 
	          col_side_colors = data.frame(groups), hide_colorbar = FALSE, showticklabels = c(FALSE, TRUE), fontsize_row = 8)
  h <- heatmaply(highly_variable, file = paste(outdir,"/",basename,"_heatmap_normalLog2CPM.png",sep =""),
            limits = NULL, xlab = "Samples", ylab = "Genes", colors = brewer.pal(11, "Spectral"), scale = "row",
            main = paste("Top",nrow(highly_variable),"most variable genes across samples",sep=" "), key.title=NULL, 
            col_side_colors = data.frame(groups), hide_colorbar = FALSE, showticklabels = c(FALSE, TRUE), fontsize_row = 8)  
  h$width <- 1200
  h$height <- 800
  export(h, file = paste(outdir,"/",basename,"_heatmap_normalLog2CPM.png",sep =""))
  
}



### Heatmap generation is now being terminated...
print(" ---------------------------------------------- HEATMAP HAS BEEN GENERATED ---------------------------------------------- ")


# Deactivating sink(). Revert output back to the console
sink()