#Rscript nb_regression_outlier_filtering.R high_vs_low_otu_table.txt high_low_mapfile.txt High Low Treatment ZINB_NB_Output_result.txt 2

start.time <- Sys.time()

library(pscl)
library(MASS)
library(foreach)
library(doMC)

data_load <- function(MYdata,mapfile,treatment,trt1,trt2){    
  colnames(MYdata)
  MYmeta <- read.table(mapfile,header = T, sep = "\t", check.names = F, comment.char= "") #change Group header to Treatment
  colnames(MYmeta)
  allcols <- length(colnames(MYdata))
  MYdata <- MYdata[,order(names(MYdata))]
  MYdata2 <- MYdata[,1:(allcols-1)]
  MYdata2 <- MYdata2[,order(names(MYdata2))]
  MYmeta <- MYmeta[order(MYmeta[,"#SampleID"]), ]
  matrix1 <- MYdata2[c(as.factor(MYmeta[,"#SampleID"][MYmeta[,treatment]==trt1]))] # change whichever group you are testing
  matrix2 <- MYdata2[c(as.factor(MYmeta[,"#SampleID"][MYmeta[,treatment]==trt2]))] # change whichever group you are testing
  matrix3 <- cbind(matrix1,matrix2)
  mat.array <- t(matrix3)
  both <- merge(mat.array,MYmeta,by.x=0,by.y="#SampleID")
}
 
zinb_nb_test <- function(both,MYdata,trt,categ1,categ2){
  all_data <- foreach(i=2:(length(rownames(MYdata))+1), .combine = rbind) %dopar% {
  final_vec <- c()
  formula1 <- as.formula(paste("both[,i] ~ ",trt," | 1",sep=""))
  formula2 <- as.formula(paste("both[,i] ~ ",trt,sep=""))
  result.pois <- tryCatch(glm(formula2, family="poisson", data = both),error=function(e) NA)
  result.zinb <- tryCatch(zeroinfl(formula1, data = both, dist = "negbin"),error=function(e) NA)			
  result.nb <- tryCatch(glm.nb(formula2, data = both),error=function(e) NA)
  pois.coeff <- tryCatch(exp(summary(result.pois)$coefficients[2,1]),error=function(e) NA)			# Column 1
  pois.pval <- tryCatch(summary(result.pois)$coefficients[2,4],error=function(e) NA)					# Column 2
  nb.coeff <- tryCatch(exp(summary(result.nb)$coefficients[2,1]),error=function(e) NA)				# Column 3
  nb.pval <- tryCatch(summary(result.nb)$coefficients[2,4],error=function(e) NA)						# Column 4
  zinb.coeff <- tryCatch(exp(summary(result.zinb)$coefficients$count[2,1]),error=function(e) NA)		# Column 5
  zinb.pval <- tryCatch(summary(result.zinb)$coefficients$count[2,4],error=function(e) NA)			# Column 6
  aic.pois <- tryCatch(AIC(result.pois),error=function(e) NA)											# Column 7
  aic.nb <- tryCatch(AIC(result.nb),error=function(e) NA)											# Column 8
  aic.zinb <- tryCatch(AIC(result.zinb),error=function(e) NA)											# Column 9
  # BIC calculated by passing k=log(n) as argument to AIC because BIC function causing problems
  bic.pois <- tryCatch(AIC(result.pois,k=log(length(both[,i]))),error=function(e) NA)											# Column 10	
  bic.nb <- tryCatch(AIC(result.nb,k=log(length(both[,i]))),error=function(e) NA)											# Column 11
  bic.zinb <- tryCatch(AIC(result.zinb,k=log(length(both[,i]))),error=function(e) NA)											# Column 12
  final_vec <- c(pois.coeff,pois.pval,nb.coeff, nb.pval, zinb.coeff, zinb.pval, aic.pois, aic.nb, aic.zinb, bic.pois, bic.nb, bic.zinb)						# Appended data from Columns 1-12
    
  shap_wilk_pval <- tryCatch(shapiro.test(both[,i])$p.value,error=function(e) NA)       # Significant p-value indicates data is not normally distributed. (Column 15)
  pval_ttest <- tryCatch(t.test(formula2, data=both)$p.value,error=function(e) NA)		# Column 14
  estimate_tab <- tryCatch(t.test(formula2, data=both)$estimate,error=function(e) NA)		# Column 17-20 drawn from this object
  heading <- paste(gsub(" ","",strsplit(names(estimate_tab)[1],"mean in group")[[1]][2],fixed=TRUE),"_minus_",gsub(" ","",strsplit(names(estimate_tab)[2],"mean in group")[[1]][2],fixed=TRUE),"_mean",sep="")		# Column 16
  mean_diff <- tryCatch((estimate_tab[1][[1]] - estimate_tab[2][[1]]),error=function(e) NA)			# Column 13
  kwtest <- tryCatch(kruskal.test(formula2,data=both)$p.value,error=function(e) NA)					# Column 21
  warn.nb <- tryCatch(glm.nb(formula2, data = both),error=function(e) NA,warning=function(w) w)		
  valwarn.nb <- ifelse(class(warn.nb)[1] == "simpleWarning", 'yes', 'no')								# Column 22
  trt1vals <- both[,i][which(as.vector(both[,trt]) == categ1)]
  trt2vals <- both[,i][which(as.vector(both[,trt]) == categ2)]
  zerotrt1 <- sum(trt1vals == 0)																		# Column 23
  zerotrt2 <- sum(trt2vals == 0)																		# Column 24
  nonzerotrt1 <- length(trt1vals) - zerotrt1															# Column 25
  nonzerotrt2 <- length(trt2vals) - zerotrt2															# Column 26
  totaltrt1 <- sum(trt1vals)																			# Column 27
  totaltrt2 <- sum(trt2vals)																			# Column 28
  mean.otu <- mean(both[,i])																			# Column 29
  var.otu <- var(both[,i])																			# Column 30
  var.mean.ratio <- var.otu/mean.otu																	# Column 31
  newd <- both[,i][-(which(both[,i] > 5*IQR(both[,i])))]       # Select indices of values that are not greater than 5 times the IQR (i.e., values > 5*IQR will be removed)
  treatment <- as.vector(both[,trt])
  newmeta <- treatment[-(which(both[,i] > 5*IQR(both[,i])))]    # Select values greater than 5 times the IQR
  newpval_pois <- tryCatch(summary(glm(newd ~ newmeta, family="poisson"))$coefficients[2,4],error=function(e) NA)		# Column 32
  newpval_nb <- tryCatch(summary(glm.nb(newd ~ newmeta))$coefficients[2,4],error=function(e) NA)						# Column 33
  newpval_zinb <-  tryCatch(summary(zeroinfl(newd ~ newmeta | 1, dist = "negbin"))$coefficients$count[2,4],error=function(e) NA)	# Column 34
  aic.filt.pois <- tryCatch(AIC(glm(newd ~ newmeta, family="poisson")),error=function(e) NA)										# Column 35
  aic.filt.nb <- tryCatch(AIC(glm.nb(newd ~ newmeta)),error=function(e) NA)														# Column 36
  aic.filt.zinb <- tryCatch(AIC(zeroinfl(newd ~ newmeta | 1, dist = "negbin")),error=function(e) NA)								# Column 37
  bic.filt.pois <- tryCatch(AIC(glm(newd ~ newmeta, family="poisson"),k=log(length(both[,i]))),error=function(e) NA)										# Column 38
  bic.filt.nb <- tryCatch(AIC(glm.nb(newd ~ newmeta),k=log(length(both[,i]))),error=function(e) NA)														# Column 39
  bic.filt.zinb <- tryCatch(AIC(zeroinfl(newd ~ newmeta | 1, dist = "negbin"),k=log(length(both[,i]))),error=function(e) NA)								# Column 40

  bestmod <- c("Poisson","NB","ZINB")
  all.aic.nonfilt <- c(aic.pois,aic.nb,aic.zinb)
  all.bic.nonfilt <- c(bic.pois,bic.nb,bic.zinb)
  all.aic.filt <- c(aic.filt.pois,aic.filt.nb,aic.filt.zinb)
  all.bic.filt <- c(bic.filt.pois,bic.filt.nb,bic.filt.zinb)
	
  aic.nonfilt.best <- bestmod[which(all.aic.nonfilt == min(all.aic.nonfilt, na.rm=TRUE))][1]											# Column 41
  bic.nonfilt.best <- bestmod[which(all.bic.nonfilt == min(all.bic.nonfilt, na.rm=TRUE))][1]											# Column 42
  aic.filt.best <- bestmod[which(all.aic.filt == min(all.aic.filt, na.rm=TRUE))][1]													# Column 43
  bic.filt.best <- bestmod[which(all.bic.filt == min(all.bic.filt, na.rm=TRUE))][1]													# Column 44
  bestmodel <- c(aic.nonfilt.best,bic.nonfilt.best,aic.filt.best,bic.filt.best)

  final_vec <- c(final_vec,mean_diff,pval_ttest,shap_wilk_pval,heading,estimate_tab[1][[1]],names(estimate_tab)[1],estimate_tab[2][[1]],names(estimate_tab)[2],kwtest,valwarn.nb,zerotrt1,zerotrt2,nonzerotrt1,nonzerotrt2,totaltrt1,totaltrt2,mean.otu,var.otu,var.mean.ratio,newpval_pois,newpval_nb,newpval_zinb,aic.filt.pois,aic.filt.nb,aic.filt.zinb,bic.filt.pois,bic.filt.nb,bic.filt.zinb,bestmodel)			# 
  final_vec
  }
  return (all_data)
}
 

final_steps <- function(otutable,mapfile,categ1,categ2,trt,outputname){
  MYdata <- read.table(otutable,header = T, sep = "\t", check.names = F, row.names =1, comment.char= "", skip =1,quote="")
  both <- data_load(MYdata,mapfile,trt,categ1,categ2)
  all_data <- zinb_nb_test(both,MYdata,trt,categ1,categ2)
  allcols <- length(colnames(MYdata))
  pois.qval <- p.adjust(all_data[,2], method = "fdr")					# Column 2 of all_data 
  nb.qval <-  p.adjust(all_data[,4], method = "fdr")					# Column 4 of all_data
  zinb.qval <- p.adjust(all_data[,6], method = "fdr")					# Column 6 of all_data
  ttest.qval <- p.adjust(all_data[,14], method = "fdr")					# Column 14 of all_data
  kw.qval <- p.adjust(all_data[,21], method = "fdr")					# Column 21 of all_data
  pois.filt.qval <- p.adjust(all_data[,32], method = "fdr")				# Column 32 from all_data
  nb.filt.qval <- p.adjust(all_data[,33], method = "fdr")				# Column 33 from all_data
  zinb.filt.qval <- p.adjust(all_data[,34], method = "fdr")				# Column 34 from all_data
  taxonomy <- MYdata[allcols]
  otuids <- colnames(both)[2:(length(colnames(both))-1)]
  taxlabels <- as.vector(taxonomy[,1])
  difflabel <- unique(all_data[,16])						# It's the unique entry of 'heading' column in all_data
  mean1_head <- unique(all_data[,18])					# It's the unique 'names(estimate_tab)[1]' column in all_data
  mean2_head <- unique(all_data[,20])					# It's the unique entry of 'names(estimate_tab)[2]' column in all_data
  zrtrt1 <- paste("# of 0's in ",categ1,sep="")
  zrtrt2 <- paste("# of 0's in ",categ2,sep="")
  nzrtrt1 <- paste("# of non-zeroes in ",categ1,sep="")
  nzrtrt2 <- paste("# of non-zeroes in ",categ2,sep="")
  tottrt1 <- paste("Total count in ",categ1,sep="")
  totttrt2 <- paste("Total count in ",categ2,sep="")
  all_data <- cbind(otuids,all_data[,1:2],pois.qval,all_data[,3:4],nb.qval,all_data[,5:6],zinb.qval,all_data[,17],all_data[,19],all_data[,13],all_data[,14],ttest.qval,all_data[,21],kw.qval,all_data[,22:31],all_data[,15],taxlabels,all_data[,32],pois.filt.qval,all_data[,33],nb.filt.qval,all_data[,34],zinb.filt.qval,all_data[,7:12],all_data[,35:44])
  colnames(all_data) <- c("OTU_IDs","Poiss_Coeff","Poiss_pval","Poiss_qval","NB_Coeff","NB_pval","NB_qval","ZINB_Coeff","ZINB_pval","ZINB_qval",mean1_head,mean2_head,difflabel,"ttest_pval","ttest_qval","KW_pval","KW_qval","NB_Coeff_Estimate_Error",zrtrt1,zrtrt2,nzrtrt1,nzrtrt2,tottrt1,totttrt2,"mean_otu","variance_otu","var/mean ratio","Shapiro_Wilk_Normality_pvalue","taxonomy","pois_filt_pval","pois_filt_qval","nb_filt_pval","nb_filt_qval","zinb_filt_pval","zinb_filt_qval","aic.pois", "aic.nb", "aic.zinb", "bic.pois", "bic.nb", "bic.zinb","aic.filt.pois","aic.filt.nb","aic.filt.zinb","bic.filt.pois","bic.filt.nb","bic.filt.zinb","aic.nonfilt.best","bic.nonfilt.best","aic.filt.best","bic.filt.best") #change Difference to groups being tested
  suppressWarnings(write.table(as.matrix(all_data),file=outputname,sep="\t",append = TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE))
}

# Run functions using CLI
argv <- commandArgs(TRUE)
registerDoMC(as.numeric(argv[7]))   # number of CPU cores to use
final_steps(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6])
print (Sys.time() - start.time)
