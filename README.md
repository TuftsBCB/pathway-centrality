# pathway-centrality

This program implements a module that calculates pathway-centrality scores for pre-defined pathway gene sets. Pathway centrality measures the amount of disease-specific communication passing through each pathway gene set, by counting the number of shortest paths between disease genes and differentially expressed genes. Significance of observed pathway-centrality scores for pathways are assessed via permutation tests using 10,000 pathway genes randomly selected from 2-core of the input network.

The program requires 5 arguments:
1) a file containing genes with known mutation associated with disease of interest (-d): i.e., sample_data/bpd.disease.genes.txt
2) a file containing genes differentially expressed within disease of interest (-e): i.e., sample_data/bpd.diff.exp.genes.txt
3) a file containing protein-protein interaction pairs (-p): i.e., sample_data/hippie_high_ppi.txt
4) a file containing pathway gene sets in .gmt file format (-g): i.e., sample_data/c2.cp.kegg.v6.0.entrez.gmt
5) output directory where all output files will be placed (-o): i.e., sample_data/output/

Genes should use exactly same identifications across all the input files. In our sample_data, genes are identified using Entrez Gene IDs. 

The example command_line to run the program is:
python PCmain.py -d sample_data/bpd.disease.genes.txt -e sample_data/bpd.diff.exp.genes.txt -p sample_data/hippie_high_ppi.txt -g sample_data/c2.cp.kegg.v6.0.entrez.gmt -o sample_data/output/

The program will create 11 files:
1) pc_disease_genes.txt: input disease genes, except those that also exist in differentially expressed gene set
2) pc_diff_exp_genes.txt: duplicated copy of input differentially expressed genes
3) pc_overlapping_genes.log: genes that exist in both disease gene set and differentially expressed gene set - these genes are removed from the disease gene set
4) pc_network_lcc.txt: protein-protein interaction pairs in the largest connected component of the given ppi networks
5) pc_disease_genes_not_in_lcc.log: diseaes genes that are not in 4), excluded from the experiment
6) pc_diff_exp_genes_not_in_lcc.log: differentially expressed genes that are not in 4), excluded from the experiment
7) pc_shortest_paths.txt: all possible shortest paths from input disease genes to differentially expressed genes in the largest connected component.
8) pc_pathway_genes.txt: input pathway gene sets, excluding disease genes and differentially expressed genes
9) pc_scores.txt: pathway centrality score calculated for all pathway gene sets
10) pc_p_cent.txt: p-value calculated for observed pathway centrality score for each pathway gene set using permutation tests
11) pc_p_cent.log: log file for permutation test, contains genes in the pool for random sampling and time records for progress
  
