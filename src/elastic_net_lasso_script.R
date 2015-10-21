# Elastic Net and Lasso Wrapper Script
# Currently set to perform logistic regression (binary response variable)
# Rscript elastic_net_lasso_script.r high_vs_low_otu_table.txt high_low_mapfile.txt en_lasso_output.csv auc binomial 4

# http://stackoverflow.com/questions/17032264/big-matrix-to-run-glmnet
# Cannot use model.matrix for this dataset as we have over 40,000 OTUs to test
# family=binomial for two-level factor. type.measure can be 'auc' on Mark Segal's suggestion

#######
require(Matrix)
require(glmnet)
require(doMC)
require(methods)

lasso_enet <- function(x,y,num,methodtype,compid,outputname,typemeasure,familydist,taxnames,one_vec,two_vec){
	#fitted = glmnet(x, y,standardize=FALSE,alpha=num)
	#plot(fitted, xvar = "lambda", label = TRUE, main=paste("fitted_",compid,"_",methodtype,sep=""))
	cat(paste(compid,"_cv_started_",methodtype,sep=""),"\n",sep="")
	cv.dat = cv.glmnet(x,y,grouped=FALSE,nfolds=length(y),alpha = num,parallel=TRUE,type.measure=typemeasure,family=familydist)
	# standardize=TRUE if explanatory variables not on same scale. lambda=seq(500,0,length.out=100), [for cv.glmnet]
	cat(paste(compid,"_cv_completed_",methodtype,sep=""),"\n",sep="")
	plot(cv.dat,main=paste("cv_",compid,"_",methodtype,sep=""))
	coefval <- coef(cv.dat)
	for (i in dimnames(coefval)[[1]]){
		all_dat <- c()
		if (coefval[i,] != 0){
			if (coefval[i,] < 0){
				if (as.vector(i) != "(Intercept)"){
					all_dat <- c(compid,methodtype,coefval[i,],i,as.vector(taxnames[i,]),one_vec)
				} else {
					all_dat <- c(compid,methodtype,coefval[i,],i,"(Intercept)",one_vec)
				}
			}
			if (coefval[i,] > 0){
				if (as.vector(i) != "(Intercept)"){
					all_dat <- c(compid,methodtype,coefval[i,],i,as.vector(taxnames[i,]),two_vec)
				} else {
					all_dat <- c(compid,methodtype,coefval[i,],i,"(Intercept)",two_vec)
				}	
			}
			write.table(as.matrix(t(all_dat)),file=outputname,sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		}
	}
}


data_load_run <- function(otuinput,mapfile,outputname,typemeasure,familydist){
	otutab <- read.table(otuinput,header = T, sep = "\t", check.names = F, row.names =1, comment.char= "",skip=1)
	MYmetaEF <- read.table(mapfile,header = T, sep = "\t", check.names = T, row.names =1, comment.char= "")
	taxnames <- matrix(nrow=nrow(otutab),ncol=1)
	taxnames[,1] <- as.vector(otutab$taxonomy)
	rownames(taxnames) <- rownames(otutab)
	otutab$taxonomy <- NULL
	otutab <- subset(as.matrix(otutab),select=rownames(MYmetaEF))
	otutab <- t(otutab)
	# All cases
	for (i in 1:length(colnames(MYmetaEF))){
		respname <- colnames(MYmetaEF)[i]
		respvariable <- subset(MYmetaEF,select=respname)
		# http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
		# Positive coeff increased in two_vec
		# Negative coeff increased in one_vec
		y <- as.numeric(respvariable[,1])
		one_vec <- c()
		two_vec <- c()
		for (kn in unique(as.vector(respvariable[,1]))){
			if (length(which(as.vector(y==1))) == length(which(as.vector(respvariable[,1])==kn))){
				one_vec <- c(kn)
			}
			if (length(which(as.vector(y==2))) == length(which(as.vector(respvariable[,1])==kn))){
				two_vec <- c(kn)
			}
		}
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),1,"Lasso",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.1,"EN_0.1",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.2,"EN_0.2",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.3,"EN_0.3",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.4,"EN_0.4",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.5,"EN_0.5",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.6,"EN_0.6",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.7,"EN_0.7",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.8,"EN_0.8",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		lasso_enet(as.matrix(otutab),as.numeric(respvariable[,1]),0.9,"EN_0.9",respname,outputname,typemeasure,familydist,taxnames,one_vec,two_vec)
		}
}


argv <- commandArgs(TRUE)
otutable <- argv[1]
mapfile <- argv[2]
outputfile <- argv[3]
measuretype <- argv[4]
pdf <- argv[5]
ncor <- argv[6]
registerDoMC(cores = as.numeric(ncor))

start.time <- Sys.time()
data_load_run(otutable,mapfile,outputfile,measuretype,pdf)
print (Sys.time() - start.time)
