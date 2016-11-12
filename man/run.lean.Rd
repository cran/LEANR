\name{run.lean}
\alias{lean}
\alias{run.lean}
\title{
 Run the LEAN approach
}
\description{
 Apply the LEAN approach to a given network and a list of pvalues
}
\usage{
run.lean(ranking, network, ranked = F, 
    add.scored.genes = F, keep.nodes.without.scores = F, 
    verbose = F, n_reps = 10000, bootstrap = F, ncores = NULL)
}
\arguments{
  \item{ranking}{
  Either a file containing gene p-values or a named numerical vector of p-values with names matching node names used in the network
}
  \item{network}{
  Either a file containing the network in sif format or an igraph graph object representing the network
}
  \item{ranked}{
  whether to transform input p-values into a uniformly distributed list of p-values based on the genes' rank before p* calculation
}
  \item{add.scored.genes}{
  whether to create one singleton node for each gene with a score but not occurring in the graph
}
  \item{keep.nodes.without.scores}{
  whether to keep nodes of the graph that have no recorded score. For those nodes it is still possible to compute enrichment scores if at least one of their network neighbors has a recorded score.
}
  \item{verbose}{
  whether to print additional status messages
}
  \item{n_reps}{
  the number of samples each background distribution should consist of. Largely influences the run-time, but higher values needed for meaningful empirical pvalues !
}
  \item{bootstrap}{
  whether to draw the pvalues of the background distributions with or without replacement
}
  \item{ncores}{
  number of cores to be used in parallel computation. Default (NULL) leads to automatic guessing of max number of cores to be used (depending on operating system).
}
}
\value{
A list object containing the results of the LEAN run. The list encompasses the following elements:
\item{restab }{Result table of applying LEAN to the real data}
\item{randtab }{Result table of applying LEAN to a permuted p-value list}
\item{indGraph }{igraph graph representing the input network after adapting it according to parameters <add.scored.genes>, <keep.nodes.without.scores> and the presence of gene scores in the input scores}
\item{nhs }{The extracted local subnetworks. Encoded as a named (by protein/gene ids) list of igraph node indices detailing each evaluated local subnetwork}
\item{gene.scores }{The gene p-values extracted from the input scores. Encoded as a numeric vector named with protein/gene ids}
}
\references{
Gwinner et al., Network-based analysis of omics data: The LEAN method, MS submitted to Bioinformatics
}
\author{
Frederik Gwinner
}
\seealso{
    \code{\link[LEANR:LEANR-package]{LEANR-package}}
}
\examples{
## Simple use case starting from a test network and p-value list
\dontrun{ 
# compute LEAN p-values starting from a p-value file and a network file
rank_file<-system.file('extdata/pvals_red.txt.gz', package='LEANR')
net_file<-system.file('extdata/g_red.sif.gz', package='LEANR')
system.time(res<-run.lean(ranking=rank_file, network=net_file, 
  add.scored.genes=T, verbose=T, n_reps=1000, ncores=3))

# compute LEAN p-values starting from a list of gene scores and a graph
data(pvals_red)
data(g_red)
system.time(res2<-run.lean(ranking=pvals_red, network=g_red, verbose=T, 
   n_reps=1000, ncores=3))

} 
}
