allmeta = []

infile = open("mapfile_hmp_test.txt",'rU')
for line in infile:
	if not line.startswith("#SampleID"):
		allmeta.append(line.strip().split("\t")[0])
		
		
#print len(allmeta)



### Creates an OTU table by filtering those samples not present in original mapping file.
outfile = open("otu_table_hmp_qiime_mapped.txt","w")
infile = open("v13_psn_otu.genus.fixed.txt",'rU')
sampind = []
samplenames = []
for line in infile:
	if line.startswith("# Constructed from biom"):
		outfile.write("# Constructed from biom"+"\n")
	elif line.startswith("#OTU"):
		spline = line.strip().split("\t")
		for x in spline:
			if x in allmeta:
				sampind.append(spline.index(x))
				samplenames.append(x)
		outfile.write('\t'.join(["#OTU ID"]+[spline[j] for j in sampind]+["taxonomy"])+'\n')
	else:
		spline = line.strip().split("\t")
		outfile.write('\t'.join([spline[0]]+[spline[j] for j in sampind]+[spline[-1]])+'\n')
outfile.close()



### Creates a new mapping file by removing those samples that are present in the filtered OTU table. 

outfile = open("mapfile_hmp_samples_selected.txt",'w')
infile = open("v13_map_uniquebyPSN.txt",'rU')
for line in infile:
	if line.startswith("#SampleID"):
		spline = line.strip().split("\t")
		outfile.write(spline[0]+"\t"+spline[3]+"\n")
	else:
		spline = line.strip().split("\t")
		if spline[0] in samplenames:
			outfile.write(spline[0]+"\t"+spline[3]+"\n")
outfile.close()