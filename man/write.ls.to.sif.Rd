\name{write.ls.to.sif}
\alias{write.ls.to.sif}
\alias{local.subnetwork.export}
\title{
 Extract the "local subnetwork" around a given protein
}
\description{
 Extract the "local subnetwork" around a given protein and write it to a Cytoscape-readable .sif file
}
\usage{
write.ls.to.sif(prot_id, LEANres, outfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prot_id}{
    protein id compatible with node names used in graph \code{g}
}
  \item{LEANres}{
    LEAN result object (list) returned by <run.lean> or <run.lean.fromdata>
}
  \item{outfile}{
  character string describing the location of an output file. Should end in .sif to be able to load it in Cytoscape.
}
}
\author{
Frederik Gwinner
}

\seealso{
 \code{\link[LEANR:run.lean]{run.lean}}
}