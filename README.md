# clusterperiod_GEE
The current folder contains files for implementing the cluster-period GEE approach introduced in  "Marginal Modeling of Cluster-Period Means and Intraclass Correlations in Stepped Wedge Designs with Binary Outcomes" by F. Li et al. (Biostatistics, 2021)

The binary responses are generated from a marginal mean model with specified correlation structures to mimic a cross sectional stepped wedge 
cluster randomized trial. The treatment is assiged at the cluster level. 

List of Files:
1) binGEN_ed.R = R file for simulating correlated binary outcome data with exponential decay correlation structure
2) binGEN_nex.R = R file for simulating correlated binary outcome data with nested exchangeable correlation structure
3) cpgee_ed.R = cluster-period GEE program for cluster-period means of binary outcomes with exponential decay correlation structure
4) cpgee_nex.R = cluster-period GEE program for cluster-period means of binary outcomes with nested exchangeable correlation structure
3) example_analysis.R = R file that illustrates the use of cpgee_nex/cpgee_ed main program based on a simulated dataset using binGEN_nex/binGEN_ed

NOTES:  1) The example program demonstrates the computation for a trial with 12 clusters, 20~30 individuals per cluster at each period and 5 periods (1 baseline). 
	    2) You will need to change path/directory names to import the functions before running the example program. 
