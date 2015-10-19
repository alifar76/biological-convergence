rm(list=ls())

require(reshape2)
require(ggplot2)
require(grid)

setwd('/Users/alifaruqi/Desktop/Projects/Kei_Project/URECA_7_year_outcome/asthma_7yr_022315/non_atopic_exclusion_081815/asthma_vs_nonasthma/elastic_net_090815/plots_elastic_net_101215')


infile <- "en_result.txt"

MYmeta <- read.table(infile,header = T, sep = "\t", check.names = F, comment.char= "")



for (val in seq(0.1, 1, by =0.1)){
	coln <- paste("en",val,sep="")
	MYmeta[,coln] <- ifelse(is.na(MYmeta[,coln]),0,val)
}


MYmeta[,"name_val"] <- paste(MYmeta[,"OTU_ID"],MYmeta[,"Genus"],MYmeta[,"Species"],sep="_")
allen <- MYmeta[,c("name_val",paste("en",seq(0.1,1,by=0.1),sep=""))]
allen <- melt(allen,id.vars=c("name_val","threemod"))
allen <- allen[order(allen[,"name_val"]), ]

allen <- subset(allen, value > 0)


ggplot(allen, aes_string(x = "name_val",y="value")) + geom_point(aes(colour = threemod),size=2.5) +
xlab("Taxa") + ylab("Alpha") +
scale_y_continuous(breaks=seq(0.1,1,by=0.1)) + coord_flip() + 
ggtitle("Taxa in final LOOCV model (of elastic nets) at different levels of alpha that are also significantly enriched based on three-model approach") +
theme(plot.title = element_text(lineheight=.8, face="bold",size=20),axis.text.x=element_text(size=20),
        axis.title=element_text(size=40,face="bold"),legend.key.size = unit(1, "cm"),legend.text=element_text(size=30))














dat <- subset(melt(MYmeta,id.vars="OTU_ID"), variable %in% c(paste("en",seq(0.1,1,by=0.1),sep="")))




ggplot(dat, aes_string(x = "OTU_ID",y="value")) + geom_dotplot()


test <- matrix(nrow=20,ncol=2)
test[,1] <-rep(c("OTU_A","OTU_B"),1,each=10)
test[,2] <- c(0,0.2,0.3,0.4,0,0.6,0,0,0.9,1.0,0.1,0,0,0,0.5,0,0.7,0,0,0)
test <- as.data.frame(test)
colnames(test) <- c("otus","alpha")






ggplot(test, aes(x = otus, y = alpha)) +
  geom_point()