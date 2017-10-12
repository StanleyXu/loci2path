#' check compatibility of gene identifiers between eQTL set and gene set
#'
#' This function perform enrichment test between one eQTL set and one gene set
#' @param e.set an eqtlSet object; the eQTL set to be queried against
#' @param g.set an object of geneSet class; the gene set to be tested
#' @return a \code{data.frame} shows the number of genes from 
#'   (1) eqtl Set (2) gene Set (3) shared 
#' @export
#' @examples
#' check.geneid(eset.list$Skin, biocarta)
check.geneid <- function(e.set, g.set) {
    eid <- unique(eqtlGene(e.set))
    gid <- unique(unlist(geneSetList(g.set)))
    res <- c(length(eid), length(gid), length(intersect(eid, gid)))
    names(res) <- c("eqtl.Set", "gene.Set", "common")
    res
}
