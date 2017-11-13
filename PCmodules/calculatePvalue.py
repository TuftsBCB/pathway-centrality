import sys
import random
import math
import calculateBetweenness as cb
import networkx as nx
import csv

from collections import defaultdict
from time import gmtime, strftime

# this version does permutation test by randomly sampling a gene set of which size is same as input pathway
# the population where the gene set is sampled from is: intersection of k-core of network and superset of pathway members

# argument list
# 1) network file
# 2) shortest paths file
# 3) pathway protein file
# 4) group centrality file
# 5) output file
# 6) number of trials

# group centrality
# = sum of centrality score (varies depending on measurement) of all members of a set with a certain property

# command line example
# python ../script/8permutationTestFractionalExclusivePathwayKcore.py shortest.paths.rop.A protein.index.rop.A pathway.rop.A.kegg.v5.0.txt group.centrality.frac.inc.rop.A.kegg.v5.0.txt pvalue.frac.inc.rop.A.kegg.v5.0.txt



# load pathway information
def loadPathwayInformation(pathwayproteinfile):
	pathwayproteins = []
	with open(pathwayproteinfile, 'r') as fi:
		for line in fi:
			tokens = line.strip().split('\t')
			pathwayproteins += tokens[2:]

	return list(set(pathwayproteins))



# get group centrality (simple sum up)
def getGroupCentrality(proteins, proteinCentrality):
	sumCentrality = 0
	for p in proteins:
		sumCentrality += proteinCentrality[p]
	if len(proteins) == 0:
		return 0
	else:
		return sumCentrality / float(len(proteins))



# run permutation test
def permutationTest(population, numSamples, score, numTrials, proteinCentrality):
	goodones = 0
	for i in xrange(numTrials):
		# select random sample
		sampleProteins = random.sample(population, numSamples)
		sampleScore = getGroupCentrality(sampleProteins, proteinCentrality)
		if sampleScore >= score:
			goodones += 1
	return float(goodones) / float(numTrials)


# get 2-core of input network
def getKcore(nwfile):
	layer = 2
	
	# load network file
        G = nx.Graph()
        
        with open(nwfile, 'r') as fi:
                for tokens in csv.reader(fi, delimiter='\t'):
                        G.add_node(tokens[0])
                        G.add_node(tokens[1])
                        G.add_edge(tokens[0], tokens[1])
        
	kcore = nx.k_core(G, k=layer)
	
	return kcore.nodes()


# run permutation test to calculate p_cent
def significanceAssessment(nwfile, shortestpathfile, pathwayproteinfile, groupcentralityfile, outputdir, num_trials):

	kcoregenes = getKcore(nwfile)	
	outputfile = outputdir + "/pc_p_cent.txt"
	logfile = outputdir + "/pc_p_cent.log"

	flog = open(logfile, 'w')

	# load shortest paths between disease genes and differentially expressed genes of a disease
	(numShortestPaths, shortestPaths, sources, destinations) = cb.loadShortestPaths(shortestpathfile)

	# find out population where we are randomly selecting
	population = loadPathwayInformation(pathwayproteinfile)

	# logging
	flog.write("total number of proteins in pathways: %d\n" % len(population))

	# remove genes that are not in k-core from population
	population = list(set(population).intersection(kcoregenes))
	flog.write("population genes: %s\n" % ', '.join(population))
	
	# logging
	flog.write("total number of pathway proteins in 2-core: %d\n" % len(population))
	
	# calculate betweenness score for each protein
	allProteinsCentrality = cb.calculateBetweenness(population, numShortestPaths, shortestPaths)	

	# while loading group centrality score file, run permutation test
	with open(groupcentralityfile, 'r') as fi, open(outputfile, 'w') as fo:
		# skip the header line
		fi.readline()

		# add header to outputfile
		fo.write("#pathway\t#p_value\t#num_genes\t#num_pairs\t#sum_fractional_group_centrality\t#averaged_fractional_group_centrality\n")

		for line in fi:
			tokens = line.strip().split('\t')
			pathwayname = tokens[0]
			numSamples = int(tokens[1])
			score = float(tokens[4])

			# run permutation test
			pvalue = permutationTest(population, numSamples, score, num_trials, allProteinsCentrality)		

			# logging 
			flog.write("permutation test for %s started: %s\n" % (pathwayname, strftime("%Y-%m-%d %H:%M:%S", gmtime()))) 

			# print out results
                        digit = int(math.log10(num_trials))
                        formatstr = "%." + str(digit) + "f"
                        pvaluestr = formatstr % pvalue
                        fo.write("%s\t%s\t%d\t%s\t%s\t%.4f\n" % (pathwayname, pvaluestr, numSamples, tokens[2], tokens[3], score))       

	flog.close()

