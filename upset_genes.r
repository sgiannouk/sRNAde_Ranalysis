# input_path = "/Users/ernesto/Desktop/toolbox/cellVSexounder_FC.venn"
# title = "this is a test"

library("UpSetR")

args = commandArgs(trailingOnly=TRUE)
input_path = args[1]
# output_path = gsub(".venn", "_venn.jpg", input_path)
output_path = args[2]
  
input_table = read.table(input_path, sep="\t", 
                         row.names = 1, stringsAsFactors = FALSE)
l = list()
v = c()
i <- 1
for (val in row.names(input_table)) {
  current_string<-input_table[val,"V2"]
  current_vector <-unlist(strsplit(current_string, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE))
  # print(length(current_vector))
  print(current_vector)
  fake_list<- list()
  fake_list[[val]]<- current_vector
  l<- append(l,fake_list)
  # l[[ val]] <- current_vector
  v <- c(v,toString(val))
  # print(l[[val]])
  i <- i + 1
}

jpeg(output_path,width = 8, height = 5, units = 'in', res = 300)
upset(fromList(l), nsets = 4,sets.bar.color = "#56B4E9", order.by = "freq")
dev.off()

# 
# library("grid")
# grid.text("My_Title",x = 0.65, y=0.95, gp=gpar(fontsize=20))

# upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
#       order.by = "freq", empty.intersections = "on")
# 
# upset(mutations,  sets.bar.color = "#56B4E9",
#       , empty.intersections = "on")
# 
# 
# input_table[,1] = apply(input_table["V2"],1,toString)
# 
# 
# example_string = toString(input_table[1,1]) 
# 
# strsplit(example_string, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)
