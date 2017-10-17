#### all setGeneric functions
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


#' @rdname loci2pathResult-class
#' @exportMethod resultTable
setGeneric("resultTable", function(x) standardGeneric("resultTable"))

#' @rdname loci2pathResult-class
#' @exportMethod coveredGene
setGeneric("coveredGene", function(x) standardGeneric("coveredGene"))


#' @rdname getTissueDegree-methods
#' @exportMethod getTissueDegree
setGeneric("getTissueDegree", 
           function(res,...) standardGeneric("getTissueDegree"))

#' @rdname getMat-methods
#' @exportMethod getMat
setGeneric("getMat", function(res,...) standardGeneric("getMat"))

#' @rdname getPathDescription-methods
#' @exportMethod getPathDescription
setGeneric("getPathDescription", 
           function(res,...) standardGeneric("getPathDescription"))

#' @rdname getHeatmap-methods
#' @exportMethod getHeatmap
setGeneric("getHeatmap", function(res,...) standardGeneric("getHeatmap"))

#' @rdname getWordcloud-methods
#' @exportMethod getWordcloud
setGeneric("getWordcloud", function(res,...) standardGeneric("getWordcloud"))

#' @rdname getPval-methods
#' @exportMethod getPval
setGeneric("getPval", function(res,...) standardGeneric("getPval"))

#' @rdname query-methods
#' @exportMethod query
setGeneric("query", 
           function(query.gr, loci, path, ...) standardGeneric("query"))

