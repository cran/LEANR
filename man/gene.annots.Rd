\name{gene.annots}
\alias{gene.annots}
\docType{data}
\title{
  Annotation for STRING protein Ids
}
\description{
 Annotation table giving gene names and descriptions for each protein contained in the STRING network
}
\usage{data("gene.annots")}
\format{
  A data frame with 7342 observations on the following 4 variables.
  \describe{
    \item{\code{ensembl_gene_id}}{a character vector}
    \item{\code{mgi_symbol}}{a character vector}
    \item{\code{entrezgene}}{a character vector}
    \item{\code{description}}{a character vector}
  }
}
\details{
  Row.names of the data.frame are STRING protein Ids
}

\examples{
data(gene.annots)
str(gene.annots)
}
\keyword{datasets}
