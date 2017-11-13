import sys

# argument list
# 1) disease genes
# 2) differentially expressed genes
# 3) output file: genes in disease
# 4) output file: genes in both files

def readFile(fi):
	data = []
	'''
	with open(filename, 'r') as fi:
		for line in fi:
			data.append(line.strip())
	'''
	for line in fi:
		data.append(line.strip())
	data = list(set(data))
	return data

def writeToFile(data, filename):
	with open(filename, 'w') as fo:
		fo.write('\n'.join(data))


def cleanupGenesets(diseasegenefile, diffgenefile, outputdir):
	disease_genes = readFile(diseasegenefile)
	diff_exp_genes = readFile(diffgenefile)

	diseasegeneonlyfile = outputdir + "/pc_disease_genes.txt"
	diffgeneonlyfile = outputdir + "/pc_diff_exp_genes.txt"
	genesinbothfile = outputdir + "/pc_overlapping_genes.log"

	only_disease = []
	both = []

	for g in disease_genes:
		if g in diff_exp_genes:
			both.append(g)
		else:
			only_disease.append(g)

	writeToFile(only_disease, diseasegeneonlyfile)
	writeToFile(both, genesinbothfile)
	writeToFile(diff_exp_genes, diffgeneonlyfile) 
	
