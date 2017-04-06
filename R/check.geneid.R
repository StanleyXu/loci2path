#' check compatibility of gene identifiers between eQTL set and gene set
#'
#' This function perform enrichment test between one eQTL set and one gene set
#' @param e.set an eqtlSet object; the eQTL set to be queried against
#' @param g.set an object of geneSet class; the gene set to be tested
#' @export
#' @examples
#' #to be added
check.geneid=function(e.set, g.set){
  eid=unique(e.set@gene)
  gid=unique(unlist(g.set@gene.set))
  res=c(length(eid), length(gid), length(intersect(eid, gid)))
  names(res)=c("eqtl.Set", "gene.Set", "common")
  res
}
