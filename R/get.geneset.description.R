#' Extract gene set description from geneSet object
#'
#' This function extracts the gene set description from geneSet object.
#' @param geneset A \code{geneSet} object
#' @param geneset.ids Character; the names of the gene sets.
#' @export
#' @examples
#' get.geneset.description(biocarta, geneset.ids=result$name_pthw)
get.geneset.description=function(geneset, geneset.ids){
  ix=match(geneset.ids, names(geneset@gene.set))
  res=geneset@description[ix]
  res
}
