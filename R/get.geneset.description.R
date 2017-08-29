#' Extract gene set description from geneSet object
#'
#' This function extracts the gene set description from geneSet object.
#' @param geneset A \code{geneSet} object
#' @param geneset.ids Character; the names of the gene sets.
#' @export
#' @examples
#' query.result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set.list=eset.list, gene.set=biocarta)
#' head(query.result$result.table)
#' query.result$cover.gene
get.geneset.description=function(geneset, geneset.ids){
  ix=match(geneset.ids, names(geneset@gene.set))
  res=geneset@description[ix]
  res
}
