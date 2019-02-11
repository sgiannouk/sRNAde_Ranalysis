#Input mappingStat.txt

if (!require("ggplot2")){ 
  install.packages("ggplot2", repos='http://cran.us.r-project.org') 
} 

if (!require("scales")){ 
  install.packages("scales", repos='http://cran.us.r-project.org') 
} 

args <- commandArgs(TRUE)
STATFILE <- args[1]
GRAPHFILE <- args[2]

th<-theme_bw()+theme(axis.title.x = element_text(size=18), axis.text.x = element_text(size=12, angle=0), axis.title.y=element_text(size=18), axis.text.y = element_text(size=12), panel.border=element_rect(linetype="dashed"))+theme(plot.title=element_text( size=24, vjust = 2.5))


table_data <- read.table(STATFILE, header = TRUE, sep = "\t", na.strings="---", colClasses=c("character","numeric","numeric","numeric","numeric","numeric"))
top<-table_data

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




png(
  filename = GRAPHFILE, 
  width=1200, 
  height=1200,
  res=90
)



p1<-ggplot(top, aes(x=name, y=wMean)) + geom_bar(stat="identity",color="black", fill="#009E73")+ylab("weighted isomiR fraction")+xlab("")+scale_x_discrete(limits=(as.vector(top$name)),label=comma)+th+ggtitle("Fraction of 3’ and 5’ Length Variants")
p2<-ggplot(top, aes(x=name, y=mean)) + geom_bar(stat="identity",color="black", fill="#045FB4")+ylab("mean isomiR fraction per microRNA")+xlab("")+scale_x_discrete(limits=(as.vector(top$name)),label=comma)+th

multiplot(p1, p2)


dev.off()
