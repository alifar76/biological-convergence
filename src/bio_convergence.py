import os
from decimal import *
import argparse
from datetime import datetime


def zero_check(checkv,dictn):
	
	val = [Decimal(checkv)]
	result = [x/abs(x) for x in val if x !=0]
	if len(result) != 0:
		return dictn[str(result[0])]


def filter_results(iname,outname,qval,filtinds,qvals):
	""" Filters file from output of 3 model test """
	infile = open(iname,'rU')
	outfile = open(outname,"w")
	grps = {}
	for line in infile:
		spline = line.strip().split("\t")
		if line.startswith('OTU_IDs'):
			# Model name, OTU ID, Mean Difference, Significant Group, Lowest BIC value, q-value for selected model, taxonomy
			outfile.write("Selected Model"+"\t"+spline[0]+"\t"+spline[12]+"\t"+"Significant Group"+"\t"+"BIC value of Selected Model (Lowest)"+"\t"+
				"q-value for selected model"+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")		
			diffg = spline[12].split("_minus_")
			grps['1'] = diffg[0]																	# Treatment group 1 name
			grps['-1'] = diffg[1].replace("_mean","")												# Treatment group 2 name
		else:
			if spline[filtinds[0]:filtinds[1]] != ['NA','NA','NA']:									# Indices of BIC values of Poisson, NB and ZINB with outliers filtered/not filtered
				newc = filter(lambda a: a != 'NA', spline[filtinds[0]:filtinds[1]])					# Remove NAs from BIC values using lambda expression
				minbic = str(min([Decimal(e) for e in newc]))										# Convert the values to decimals and obtain the minimum value
				indexic = spline[filtinds[0]:filtinds[1]].index(minbic)								# Index of lowest value of BIC of filtered values
				if spline[qvals[2]] != 'NaN':
					if (indexic == 0) and (Decimal(spline[qvals[0]]) < qval):									# Filter OTUs with Poisson q-value less than 0.05
						outfile.write("Poisson"+"\t"+spline[0]+"\t"+spline[12]+"\t"+zero_check(spline[12],grps)+"\t"+spline[filtinds[0]]+
							"\t"+spline[qvals[0]]+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")
				if spline[qvals[2]] != 'NaN':
					if (indexic == 1) and (Decimal(spline[qvals[1]]) < qval):									# Filter OTUs with NB q-value less than 0.05
						outfile.write("NB"+"\t"+spline[0]+"\t"+spline[12]+"\t"+zero_check(spline[12],grps)+"\t"+spline[filtinds[0]+1]+
							"\t"+spline[qvals[1]]+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")
				if spline[qvals[2]] != 'NaN':
					if (indexic == 2) and (Decimal(spline[qvals[2]]) < qval):									# Filter OTUs with ZINB q-value less than 0.05
						outfile.write("ZINB"+"\t"+spline[0]+"\t"+spline[12]+"\t"+zero_check(spline[12],grps)+"\t"+spline[filtinds[0]+2]+
							"\t"+spline[qvals[2]]+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")
	outfile.close()
	return outfile


def compare_en_negbin(negbin,enout,outputn):
	infile = open(negbin,'rU')
	negbindict = {}
	meandiff = {}
	for line in infile:
		if not line.startswith("Selected Model"):
			spline = line.strip().split("\t")
			negbindict[spline[1]] = spline[3]
			meandiff[spline[1]] = spline[2]
	infile = open(enout,'rU')
	outfile = open(outputn,"w")
	outfile.write("OTU_ID"+"\t"+"Mean Difference"+"\t"+"NegBin Group"+"\t"+"EN Group"+"\t"+"Taxonomy"+"\n")
	for line in infile:
		spline = line.strip().split(",")
		if spline[3] in negbindict.keys():
			outfile.write(spline[3]+"\t"+meandiff[spline[3]]+"\t"+"In NegBin: "+negbindict[spline[3]]+
				"\t"+"In "+spline[1]+": "+spline[5]+"\t"+spline[4]+"\n")
	outfile.close()
	return outfile


def final_comparing_output(iname,oname):
	otuid = []
	infile = open(iname,'rU')
	for line in infile:
		if not line.startswith("OTU_ID"):
			spline = line.strip().split("\t")
			otuid.append(spline[0])
	otuids = list(set(otuid))
	header = ["EN %s (enriched group)" %str(float(x)/10) for x in range(1,11)]
	outfile = open(oname,"w")
	outfile.write('\t'.join(["OTU_ID"]+["Mean Difference"]+
		["Three-model significant group"]+header+["taxonomy"])+"\n")
	for x in otuids:
		headlist = []
		group = []
		relabun = []
		taxon = []
		infile = open(iname,'rU')
		for line in infile:
			if not line.startswith("OTU_ID"):
				spline = line.strip().split("\t")
				if spline[0] == x:
					group.append(spline[2].split("In NegBin:")[1].strip())
					relabun.append(spline[1])
					taxon.append(spline[-1])
					if "EN" in spline[3]:
						enlevel = spline[3].split("In EN_")[1].split(":")
						headlist.append([header.index(i) for i in header if enlevel[0] in i][0])
					elif "Lasso" in spline[3]:
						headlist.append(len(header)-1)
		datlist = ['NA']*10
		for u in headlist:
			datlist[u] = list(set(group))[0]
		tax = list(set(taxon))[0]
		mean = list(set(relabun))[0]
		threemod = list(set(group))[0]
		outfile.write(x+"\t"+mean+"\t"+threemod+"\t"+'\t'.join(datlist)+"\t"+tax+"\n")
	outfile.close()
	return outfile
	


def main_function(threemod,signif_threemod,enlasso,compare_en_threemod):
	# Infile name, output file name, q-value, indices of BIC values, indices of q-values
	notfilt = [38,41]			# Indices of BIC values with outliers not filtered
	qvalnonfilt = [3,6,9]		# Indices of q-values of 3 models with outliers not filtered
	# Get only significant OTUs from dataset
	filter_results(threemod,signif_threemod,0.05,notfilt,qvalnonfilt)
	## Output from Python filter script for Negative Binomial analysis, output from MicrobeNet, 
	## output name of comparing file
	compare_en_negbin(signif_threemod,enlasso,"en_vs_negbin_compare.txt")
	# Output from compare_en_negbin to give final output showing convergent taxa
	final_comparing_output("en_vs_negbin_compare.txt",compare_en_threemod)
	os.system("rm en_vs_negbin_compare.txt")
	return


if __name__ == '__main__':
	startTime = datetime.now()
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
	description='Biological Convergence: A confluence of statistical methods',
	epilog='''
An example to run the pipeline:
python bio_convergence.py -otu lean_obese/219_otu_table.txt -map lean_obese/lean_obese_mapfile.txt -trt1 Lean -trt2 Not_Lean -metavar obesitycat -proc 4
	''')
	parser.add_argument('-otu', metavar='OTU table', nargs=1, help='Path to OTU table',
						required=True)
	parser.add_argument('-map', metavar='Mapping file', nargs=1, help='Path to mapping file',
						required=True)
	parser.add_argument('-trt1', metavar='Treatment 1', nargs=1, help='Value of 1st treatment group as it appears in \
						mapping file (such as Lean)', required=True)
	parser.add_argument('-trt2', metavar='Treatment 2', nargs=1, help='Value of 2nd treatment group as it appears in \
						mapping file (such as Not_Lean)', required=True)
	parser.add_argument('-metavar', metavar='Heading of meta-variable', nargs=1, help='Header of the meta-varaible as it \
						appears in the mapping file (of the variable to be tested)', required=True)
	parser.add_argument('-proc', metavar='Processors', type=int ,nargs=1, help='Number of \
							processors to use',required=True)
	args = parser.parse_args()

	otutab = args.otu[0]
	mapfile = args.map[0]
	trt1 = args.trt1[0]
	trt2 = args.trt2[0]
	metavar = args.metavar[0]
	coren = args.proc[0]
	cv = 'deviance'
	dist = 'binomial'
	threemod = "three_model_ouput_%s_%s.txt"	% (trt1, trt2)		# Name of output files from 3 model approach
	signif_threemod = "neg_bin_sig_output_%s_%s.txt" % (trt1, trt2)	# Name of significant OTU file from 3 model approach filtering
	enlasso = "en_lasso_result_%s_%s.csv" % (trt1, trt2)
	compare_en_threemod = "convergence_%s_%s.txt" % (trt1, trt2)
	os.system("Rscript nb_regression_outlier_filtering.R %s %s %s %s %s %s %s" % (otutab,mapfile,trt1,trt2,metavar,threemod,str(coren)))
	os.system("Rscript elastic_net_lasso_script.R %s %s %s %s %s %s" % (otutab,mapfile,enlasso,cv,dist,str(coren)))
	main_function(threemod,signif_threemod,enlasso,compare_en_threemod)
	os.system("mkdir output_bio_convergence")
	os.system("mv %s output_bio_convergence/" % threemod)
	os.system("mv %s output_bio_convergence/" % signif_threemod)
	os.system("mv %s output_bio_convergence/" % enlasso)
	os.system("mv %s output_bio_convergence/" % compare_en_threemod)
	os.system("mv Rplots.pdf output_bio_convergence/")
	print "\n"+"Task Completed! Completion time: "+ str(datetime.now()-startTime)
