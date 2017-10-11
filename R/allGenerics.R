#setGeneric("eqtlSet", function(object) standardGeneric("eqtlSet"))
#setGeneric("geneSet", function(object) standardGeneric("geneSet"))

#' Method tissue.
#' @name eqtlSet-class
#' @rdname eqtlSet-class
#' @exportMethod tissue
setGeneric("tissue", function(x) standardGeneric("tissue"))

#' Method eqtlId.
#' @name eqtlSet-class
#' @rdname eqtlSet-class
#' @exportMethod eqtlId
setGeneric("eqtlId", function(x) standardGeneric("eqtlId"))

#' Method eqtlRange.
#' @name eqtlSet-class
#' @rdname eqtlSet-class
#' @exportMethod eqtlRange
setGeneric("eqtlRange", function(x) standardGeneric("eqtlRange"))

#' Method eqtlGene.
#' @name eqtlSet-class
#' @rdname eqtlSet-class
#' @exportMethod eqtlGene
setGeneric("eqtlGene", function(x) standardGeneric("eqtlGene"))


#' Method numGene.
#' @name geneSet-class
#' @rdname geneSet-class
#' @exportMethod numGene
setGeneric("numGene", function(x) standardGeneric("numGene"))

#' Method description.
#' @name geneSet-class
#' @rdname geneSet-class
#' @exportMethod description
setGeneric("description", function(x) standardGeneric("description"))

#' Method geneSetList.
#' @name geneSet-class
#' @rdname geneSet-class
#' @exportMethod geneSetList
setGeneric("geneSetList", function(x) standardGeneric("geneSetList"))