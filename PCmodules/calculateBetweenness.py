import sys
import math
from collections import defaultdict


# load shortest path file
def loadShortestPaths(shortestpathfile):
        # load shortest paths between disease genes and differentially expressed genes of a disease
        
        # total number of shortest paths between a disease gene and a differentially expressed gene
        # the length of the -numShortestPaths- will be the number of pairs between disease genes and differentially expressed genes
        numShortestPaths = {}
        shortestPaths = defaultdict(list)
        
        # collect input genes: sources and destinations
        sources = []
        destinations = []
        
        with open(shortestpathfile, 'r') as fi:
                # skip the header line
                fi.readline()

                for line in fi:
                        tokens = line.strip().split('\t')
                        
                        sources.append(tokens[0])
                        destinations.append(tokens[1])
                        
                        pair = (tokens[0], tokens[1])
                        
                        # first two tokens are source and destination
                        # from third token, it is shortest path
                        if tokens[2] == "None":
                                numShortestPaths[pair] = 0
                        else:
                                numShortestPaths[pair] = len(tokens) - 2
        
                                for i in xrange(2, len(tokens)):
                                        path = tokens[i]
                                        nodes = path.split('-')
                                
                                        # counting the number of the shortest paths between a disease gene and a differentially expressed gene for each node
                                        for j in xrange(1, len(nodes)-1):
                                                shortestPaths[nodes[j]].append(path)
        
        sources = list(set(sources))
        destinations = list(set(destinations))
        
        return (numShortestPaths, shortestPaths, sources, destinations)


# calculate fractional betweenness score for each protein
def calculateBetweenness(allPathwayProteins, numShortestPaths, shortestPaths):
        betweenness = {}  # size of betweenness is number of proteins in allPathwayProteins

        for gene in allPathwayProteins:
                # count number shortest paths going through "gene"
                sumNumShortestPaths = {}
                # collect all shortest paths between sources and destinations going through "gene"
                allShortestPathsThrough = shortestPaths[gene]

                for p in allShortestPathsThrough:
                        # p = 'x-u-v-y'
                        tokens = p.split('-')
                        pair = (tokens[0], tokens[-1])  # pair is (source, destination)
                        
                        try:
                                sumNumShortestPaths[pair] += 1
                        except KeyError:
                                sumNumShortestPaths[pair] = 1
                
                # for each gene, calculate fractional betweenness score by summing up score for all possible (s, d) paris
                sumFractionalCentrality = 0
                for pair in sumNumShortestPaths.keys():
                        totalNumShortestPaths = float(numShortestPaths[pair])

                        if totalNumShortestPaths != 0:
                                sumFractionalCentrality += float(sumNumShortestPaths[pair]) / float(numShortestPaths[pair])

                # At least one shortest path should exist to calculate Pathway Centrality. ################################
                assert len(numShortestPaths.keys()) > 0, "No shortest path found between two input gene sets: check if the input genes are in input PPI networks."  
                ###########################################################################################################
                
                betweenness[gene] = sumFractionalCentrality / float(len(numShortestPaths.keys()))

        return betweenness



