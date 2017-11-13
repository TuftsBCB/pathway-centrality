'''
PCmain.py -- main function to calculate pathway centrality for pre-defined functional gene sets between a set of disease genes and a set of differentially expressed genes in given protein-protein interaction networks.

usage: python PCmain.py 

'''

import sys
import argparse
import os

import PCmodules.cleanUpInput as ci
import PCmodules.getLCC as lc
import PCmodules.getSPDisease as sp
import PCmodules.calculatePC as pc
import PCmodules.calculatePvalue as pv

def main():
	parser = argparse.ArgumentParser(description="Pathway Centrality calculation")

	parser.add_argument("-d", "--disease_gene_file", type=argparse.FileType('r'),
				help="one column file, gene identification should match with other input files.")
	parser.add_argument("-e", "--diff_exp_gene_file", type=argparse.FileType('r'),
				help="one column file, gene identification should match with other input files.")

	parser.add_argument("-o", "--output_directory", help="A directory where all output files are stored: default will be named output and placed in the working directory.",)

	parser.add_argument("-p", "--ppi_network_file", type=argparse.FileType('r'), 
				help="two column file delimited by tab, gene identification should match with other input files.")
	parser.add_argument("-g", "--pathway_gmt_file", type=argparse.FileType('r'),
				help=".gmt file format, gene identification should match with other input files.")

	option = parser.parse_args()

	
	# set output directory
	outputdir = "./output"
	if option.output_directory is not None:
		outputdir = option.output_directory + "/"

	if not os.path.exists(outputdir):
    		os.makedirs(outputdir)

	'''
        # read input files 1: disease genes and differentially expressed genes
        # create three new files in output directory: pc_disease_genes.txt, pc_diff_exp_genes.txt, pc_overlapping_genes.txt 
	# pc_disease_genes.txt: filtered disease gene set --- differentially expressed genes removed
	# pc_diff_exp_genes.txt: same as initial input
	# pc_overlapping_genes.txt: for record purposes --- genes removed from disease gene set
	ci.cleanupGenesets(option.disease_gene_file, option.diff_exp_gene_file, outputdir)	


	# read input files 2: ppi network file
	# create one new file in output directory: pc_network_lcc.txt
	# pc_network_lcc.txt: the largest connected component of the given ppi network, self-loop removed
	lc.getLCC(option.ppi_network_file, outputdir)


	# calculate shortest paths between disease genes and differentially expressed genes
	# create one new file in output directory: pc_shortest_paths.txt
	# pc_shortest_paths.txt: all possible shortest paths between disease genes and differentially expressed genes in the claulcated largest connected component
	'''
	nwfile = outputdir + "/pc_network_lcc.txt"
	'''
	diseasegenefile = outputdir + "/pc_disease_genes.txt"
	diffexpgenefile = outputdir + "/pc_diff_exp_genes.txt"
	sp.calculateShortestPaths(nwfile, diseasegenefile, diffexpgenefile, outputdir)
	'''

	# calculate pathway centrality score for input pathway gene sets
	# calculatePC.py
	shortestpathfile = outputdir + "/pc_shortest_paths.txt"
	pc.calculatePathwayCentrality(shortestpathfile, option.pathway_gmt_file, outputdir)


	# calculate p-value from permutation test using 2-core genes of the input network
	num_trial = 10000
	pathwayproteinfile = outputdir + "/pc_pathway_genes.txt"
	groupcentralityfile = outputdir + "/pc_scores.txt"
	pv.significanceAssessment(nwfile, shortestpathfile, pathwayproteinfile, groupcentralityfile, outputdir, num_trial) 	



if __name__ == "__main__":
	main()
