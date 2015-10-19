## 43140 OTUs and 2898 samples


rm(list=ls())

require(reshape2)
require(ggplot2)
require(grid)

setwd('/Users/alifaruqi/Desktop/Projects/Development_Tools/Github_Scripts/biological_convergence/hmp_data')


otutable <- "otu_table_hmp_qiime_mapped.txt"

MYdata <- read.table(otutable,header = T, sep = "\t", check.names = F, row.names =1, comment.char= "", skip =1,quote="")


mapfile <- "mapfile_hmp_samples_selected.txt"
MYmeta <- read.table(mapfile,header = T, sep = "\t", check.names = F, comment.char= "")
