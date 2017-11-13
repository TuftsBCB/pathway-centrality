import sys
import networkx as nx
import multiprocessing
import re

# argument list
# 1) network file
# 2) disease gene file
# 3) differentially expressed gene file
# 4) shortest path output file
# standard output will have genes not in the network

# networkx library should be availabe on the working machine


def getAllShortestPaths(argument):
	protein1 = argument[0]
	protein2 = argument[1]
	G = argument[2]

	resultStr = "%s\t%s" % (protein1, protein2)

	try: 
		paths =  nx.all_shortest_paths(G, source=protein1, target=protein2)
		for p in paths:
			p = map(str, p)
			resultStr += ("\t%s" % '-'.join(p))
		resultStr += "\n"

	except nx.NetworkXNoPath:
		resultStr += "\tNone\n"
	
	return resultStr


# load gene files
def loadGeneFile(filename):
	genes = []
	with open(filename, 'r') as fi:
		for line in fi:
			genes.append(line.strip())
	return genes


def calculateShortestPaths(networkfile, diseasegenefile, diffexpressedgenefile, outputdir):
	sp_outputfile = outputdir + "/pc_shortest_paths.txt"
	sp_logfile_disease = outputdir + "/pc_disease_genes_not_in_lcc.log"
	sp_logfile_diff = outputdir + "/pc_diff_exp_genes_not_in_lcc.log"


	# load ppi networks and build a graph
	G = nx.Graph()

	with open(networkfile, 'r') as fi:
		for line in fi:
			tokens = line.strip().split()
			protein1 = tokens[0]
			protein2 = tokens[1]

			# add edge between two proteins
			G.add_edge(protein1, protein2)

	numProteins = nx.number_of_nodes(G)


	# load genes
	diseasegeneset = loadGeneFile(diseasegenefile)
	diffgeneset = loadGeneFile(diffexpressedgenefile)


	# there should not any overlap between two gene lists
	assert len(set(diseasegeneset).intersection(diffgeneset)) == 0


	# collect disease and differentially expressed genes that are in the input network
	diseaseGenes = []
	diffGenes = []
	for n in G.nodes():
		if n in diseasegeneset:
			diseaseGenes.append(n)
		elif n in diffgeneset:
			diffGenes.append(n)
	

        # identify disease genes and differentially expressed genes that are not in the input network and record it out a file
	with open(sp_logfile_disease, 'w') as fo:
		fo.write("\n".join([x for x in diseasegeneset if x not in diseaseGenes]))

	with open(sp_logfile_diff, 'w') as fo:
        	fo.write("\n".join([x for x in diffgeneset if x not in diffGenes]))
	

	# calculate all shortest paths between all possible pairs of nodes and print out
	numCPU = multiprocessing.cpu_count()
	p = multiprocessing.Pool(numCPU/8)

	arguments = []
	for p1 in diseaseGenes:
		for p2 in diffGenes:
			arguments.append([p1, p2, G])

	with open(sp_outputfile, 'w') as fo:
		fo.write("#source\ttarget\tshortest_paths\n")
		
	for result in p.imap(getAllShortestPaths, arguments):
		with open(sp_outputfile, 'a') as fo:
			fo.write(result)


if __name__ == "__main__":
    main()
