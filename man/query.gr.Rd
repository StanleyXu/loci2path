% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loci2path.demo-data.R
\docType{data}
\name{query.gr}
\alias{query.gr}
\title{eQTL geneset enrichment query demo data}
\format{An object of class \code{GRanges} of length 47.}
\usage{
data(loci2path.demo)
}
\description{
Demo Data to show how to perform eQTL-geneset enrichment query.
}
\details{
\code{eqtl.set.list} A list of eQTLset objects; eQTL data are collected from 
GTeX.

\code{biocarta} A Geneset object; Geneset are from Broad Institute's MSigDB:
 curated gene set, category 'cp': BIOCARTA

\code{query.gr} A Genomic Region object; Query regions are from immunoBase,
 crohn's disease.
}
\examples{
data(loci2path.demo)
}
\keyword{datasets}
