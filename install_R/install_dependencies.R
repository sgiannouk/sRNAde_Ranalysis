#!/usr/bin/env Rscript

dependencies <- c("ggplot2", "RColorBrewer", "viridis", "heatmaply", "plotly", "vsn", "DESeq", "DESeq2", "edgeR", "NOISeq", "UpSetR")
install.packages("BiocManager", repos="http://cran.us.r-project.org")
BiocManager::install(dependencies)
