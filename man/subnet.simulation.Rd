\name{subnet.simulation}
\alias{subnet.simulation}
\title{
 Simulate subnetworks
}
\description{
 Simulate subnetworks (also called modules) and gene p-values to be then used in a ROC performance evaluation study.
}
\usage{
subnet.simulation(g, nmods=10, mod_lims=c(10,50), pval_scaling=0.1, mod_enrich_perc=0.5, 
         spec='', prob_function=function(degs){degs/sum(degs)}, create.files=T)
}
\arguments{
  \item{g}{
    igraph graph representing the network in which subnetworks are supposed to be simulated
  }
  \item{nmods}{
    number of subnetworks/modules to simulate
  }
  \item{mod_lims}{
    minimum and maximum size (number of genes) of each module
  }
  \item{pval_scaling}{
    parameter value for <p_scale>
  }
  \item{mod_enrich_perc}{
    parameter value for <p_enrich>
  }
  \item{spec}{
    string, specifier appended to the created pvalue files (if create.files=T)
  }
  \item{prob_function}{
    probability function used for picking attachment point in iterative construction of subnetworks. defaults to preferential attachment based on node degree. To disable preferential attachment, use prob_function=function(degs){rep(1/length(degs),length(degs))}
  }
  \item{create.files}{
    whether to write subnetwork simulation results to file so external approaches can be run and evaluated on them
  }
}
\value{
A list object containing the simulated subnetworks. The list encompasses the following elements:
\item{mods }{List of simulated modules/subnetworks. Each module is given by the igraph indices of the nodes contained in it.}
\item{pvals }{Result table containing for each gene in the graph its simulated pvalue (column P.Value) and its association to subnetworks or background (column NodeType)}
\item{pvalfile }{String containing the name of the file containing the equivalent information to <pvals> created in this run if create.files=T}
}
\references{
Gwinner et al., Network-based analysis of omics data: The LEAN method, MS submitted to Bioinformatics
}
\author{
Frederik Gwinner
}
\examples{
### See vignette("subnet-sim") for a use case.
}
