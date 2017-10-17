#' @include allClasses.R allGenerics.R



## eqtlSet-class accessor

#' @param x An eqtlSet object
#' @return Object of class eqtlSet

#' @rdname eqtlSet-class
#' @aliases tissue,eqtlSet-method
setMethod("tissue", "eqtlSet", function(x) x@tissue)
#' @rdname eqtlSet-class
#' @aliases eqtlId,eqtlSet-method
setMethod("eqtlId", "eqtlSet", function(x) x@eqtlId)
#' @rdname eqtlSet-class
#' @aliases eqtlRange,eqtlSet-method
setMethod("eqtlRange", "eqtlSet", function(x) x@eqtlRange)
#' @rdname eqtlSet-class
#' @aliases eqtlGene,eqtlSet-method
setMethod("eqtlGene", "eqtlSet", function(x) x@gene)

## geneSet-class accessor

#' @param x An geneSet object
#' @return Object of class geneSet
#' @rdname geneSet-class
#' @aliases numGene,geneSet-method
setMethod("numGene", "geneSet", function(x) x@numGene)
#' @rdname geneSet-class
#' @aliases description,geneSet-method
setMethod("description", "geneSet", function(x) x@description)
#' @rdname geneSet-class
#' @aliases geneSetList,geneSet-method
setMethod("geneSetList", "geneSet", function(x) x@geneSetList)



## loci2pathResult-class accessor

#' @param x An geneSet object
#' @return Object of CLass loci2pathResult
#' @rdname loci2pathResult-class
#' @aliases resultTable,loci2pathResult-method
setMethod("resultTable", "loci2pathResult", function(x) x@resultTable)
#' @rdname loci2pathResult-class
#' @aliases coveredGene,loci2pathResult-method
setMethod("coveredGene", "loci2pathResult", function(x) x@coveredGene)
