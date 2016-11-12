\name{LEANR-package}
\alias{LEAN-package}
\alias{LEANR}
\alias{LEAN}
\alias{Local enrichement analysis}
\docType{package}
\title{
Finds "local subnetworks" within an interaction network which show enrichment for differentially expressed genes
}
\description{
Implements the method described in "Network-based analysis of omics data: The LEAN method".
Given a protein interaction network and a list of p-values describing a measure of interest (as e.g. differential gene expression) 
this method computes an enrichment p-value for the protein neighborhood of each gene and compares it to a background distribution 
of randomly drawn p-values. The resulting scores are corrected for multiple testing and significant hits are returned in tabular format.
}
\details{
\tabular{ll}{
Package: \tab LEANR\cr
Type: \tab Package\cr
Version: \tab 1.4.8\cr
Date: \tab 2016-11-11\cr
License: \tab GPL-3\cr
}
 See help page of run.lean for a more detailed description of how to use this package. 
 Type vignette("CCM-data") for an example showing the application of LEAN to the CCM knockout data set discussed in the paper.
 Type vignette("subnet-sim") for an example showing the application of LEAN to simulated subnetwork data discussed in the paper.
}
\author{
Frederik Gwinner

Maintainer: Frederik Gwinner <frederik.gwinner@gmail.com>
}
\references{
Gwinner et al., Network-based analysis of omics data: The LEAN method, Bioinformatics 2016
}
\keyword{ package }
\seealso{
 \code{\link[LEANR:run.lean]{run.lean}}
 vignette("CCM-data")
 vignette("subnet-sim")
}
