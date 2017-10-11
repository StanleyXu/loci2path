## eqtlSet-class accessor

#' @include allClasses.R allGenerics.R

setMethod("tissue", "eqtlSet", function(x) x@tissue)
setMethod("eqtlId", "eqtlSet", function(x) x@snp.id)
setMethod("eqtlRange", "eqtlSet", function(x) x@snp.gr)
setMethod("eqtlGene", "eqtlSet", function(x) x@gene)

## geneSet-class accessor

setMethod("numGene", "geneSet", function(x) x@total.number.gene)
setMethod("description", "geneSet", function(x) x@description)
setMethod("geneSetList", "geneSet", function(x) x@gene.set)
