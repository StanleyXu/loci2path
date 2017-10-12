### all setGeneric functions

#' @include allClasses.R


#' @rdname eqtlSet-class
#' @export
setGeneric("tissue", function(x) standardGeneric("tissue"))


#' @rdname eqtlSet-class
#' @exportMethod eqtlId
setGeneric("eqtlId", function(x) standardGeneric("eqtlId"))

#' @rdname eqtlSet-class
#' @exportMethod eqtlRange
setGeneric("eqtlRange", function(x) standardGeneric("eqtlRange"))


#' @rdname eqtlSet-class
#' @exportMethod eqtlGene
setGeneric("eqtlGene", function(x) standardGeneric("eqtlGene"))


#' @rdname geneSet-class
#' @exportMethod numGene
setGeneric("numGene", function(x) standardGeneric("numGene"))


#' @rdname geneSet-class
#' @exportMethod description
setGeneric("description", function(x) standardGeneric("description"))


#' @rdname geneSet-class
#' @exportMethod geneSetList
setGeneric("geneSetList", function(x) standardGeneric("geneSetList"))

