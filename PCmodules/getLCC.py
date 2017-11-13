import sys
import networkx as nx

# argument lists
# 1) string-db network file
# 2) output directory

delimiter = '\t'

def getLCC(networkfile, outputdirectory):
	'''
	networkfile = sys.argv[1]
	outputdirectory = sys.argv[2]
	'''

	# initiate a networkx graph
	G = nx.Graph()
	for line in networkfile:
        	tokens = line.strip().split(delimiter)
                protein1 = tokens[0]
                protein2 = tokens[1]

		# exclude self-loops
		if protein1 != protein2:
                        # add edge between two proteins
                        G.add_edge(protein1, protein2)

	cc = sorted(nx.connected_components(G), key = len, reverse=True)		

	outputfile = outputdirectory + "/pc_network_lcc.txt"

	# only print out the largest connected component in the given ppi networks
	connectedNetwork = G.edges(cc[0])
	with open(outputfile, 'w') as fo:
		for cedge in connectedNetwork:
			fo.write("%s\t%s\n" % (cedge[0], cedge[1]))



