# LEANR
Implements the method described in [Gwinner et al., Network-based analysis of omics data: The LEAN method, Bioinformatics 2016].
Given a protein interaction network and a list of p-values describing a measure of interest 
(as e.g. differential gene expression) this method  computes an enrichment p-value for the 
protein neighborhood of each gene and compares it to a background distribution of randomly drawn p-values.
The resulting scores are corrected for multiple testing and significant hits are returned in tabular format.



See help page of <b>run.lean</b> for a more detailed description of how to use this package 
(type "?run.lean" in R prompt to do so)</br>
Type vignette("CCM-data") for an example showing the application of LEAN to the CCM knockout data set discussed in the paper.</br>
Type vignette("subnet-sim") for an example showing the application of LEAN to simulated subnetwork data discussed in the paper.</br>
