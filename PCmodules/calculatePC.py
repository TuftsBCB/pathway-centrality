import sys
import calculateBetweenness as cb
from collections import defaultdict

# argument list
# 1) shorteat path file
# 2) pathway file: .gmt file format
# 3) output file
# 4) pathway gene output file

# command line example:
# python /r/bcb/jpark/GroupCentrality/script/7calculateDiseaseCentralityFractionalInclusive.py asthma.shortest.paths asthma.genes.index /r/bcb/jpark/GroupCentrality/data/pathways/c5.bp.v5.0.symbols.gmt group.centrality.frac.inc.asthma.bp.v5.0.txt pathway.asthma.bp.v5.0.txt

# group centrality
# = sum of centrality score (varies depending on measurement) of all members of a set with a certain property
# this version uses betweenness score for the centrality measurement
# pathway geneset excludes disease genes and differentially expressed genes


def readPathwayInformation(pathwayfile, sources, destinations, outputpathwayfile):
        # load pathway gene sets
        inputdata = defaultdict(list)
        allpathwayproteins = []

        with open(outputpathwayfile, 'w') as fo:

		# add header to the output file
		fo.write("#pathway_name\t#(num_genes_original->num_genes_filtered)\t#genes\n")

                for line in pathwayfile:
                        tokens = line.strip().split('\t')
                        pathwayname = tokens[0]
                        genes = tokens[2:]

			# pathway gene set should not include disease genes or differentially expressed genes
			inputdata[pathwayname] = [x for x in genes if x not in sources and x not in destinations]

			# print out to a file
			fo.write("%s\t(%d->%d)\t%s\n" % (pathwayname, len(genes), len(inputdata[pathwayname]), '\t'.join(inputdata[pathwayname])))
			# collect all proteins in pathways
			allpathwayproteins += inputdata[pathwayname]

        allpathwayproteins = list(set(allpathwayproteins))
        return (inputdata, allpathwayproteins)



# calculate group centrality for each pathway and print out the results to a file
def calculateGroupCentrality(inputdata, centrality, numpairs, outputfile):
	with open(outputfile, 'w') as fo:
		fo.write("#pathway\t#num_genes\t#num_pairs\t#sum_fractional_group_centrality\t#averaged_fractional_group_centrality\n")
	
		for pathway in inputdata.keys():
			# set of genes in a pathway
			pathwayGenes = inputdata[pathway]
			numGenes = len(pathwayGenes)
			
			# sum up all fractional centrality score of each gene
			sumFractionalCentrality = 0
			for gene in pathwayGenes:
				sumFractionalCentrality += centrality[gene]

			if numGenes == 0:
				averagedNumShortestPaths = 0
			else:
				averagedNumShortestPaths = sumFractionalCentrality / float(numGenes)
			
			fo.write("%s\t%d\t%d\t%f\t%f\n" % (pathway, numGenes, numpairs, sumFractionalCentrality, averagedNumShortestPaths))


# main function
def calculatePathwayCentrality(shortestpathfile, pathwayfile, outputdir):

	outputfile = outputdir + "/pc_scores.txt"
	pathwayoutputfile = outputdir + "/pc_pathway_genes.txt"

	# load shortest path file
	(numShortestPaths, shortestPaths, sources, destinations) = cb.loadShortestPaths(shortestpathfile)

	# load pathway members and clean up the input data
	# remove proteins from pathway if they don't have string id or are not in network
	allPathwayProteins = []
	(inputdata, allPathwayProteins) = readPathwayInformation(pathwayfile, sources, destinations, pathwayoutputfile)

	# calculate betweenness score for all proteins in experiment
	allProteinsCentrality = cb.calculateBetweenness(allPathwayProteins, numShortestPaths, shortestPaths)
	
	# summarize centrality score for proteins in pathway and print out to a file
	numPairs = len(numShortestPaths.keys())
	calculateGroupCentrality(inputdata, allProteinsCentrality, numPairs, outputfile)


if __name__ == "__main__":
    main()
