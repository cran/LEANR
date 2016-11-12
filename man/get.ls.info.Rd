\name{get.ls.info}
\alias{get.ls.info}
\alias{local.subnetwork.info}
\title{
 Extract the genes of a "local subnetwork"" around a given protein 
}
\description{
 Extract the genes of a "local subnetwork"" around a given protein and present in tabular format
}
\usage{
 get.ls.info(prot_id, LEANres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prot_id}{
    Protein id compatible with node names used in graph.
}
  \item{LEANres}{
    LEAN result object (list) returned by <run.lean> or <run.lean.fromdata>
}
}
\author{
Frederik Gwinner
}

\seealso{
 \code{\link[LEANR:run.lean]{run.lean}}
}