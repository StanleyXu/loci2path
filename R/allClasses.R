#' eqtlSet Class
#'
#' eqtlSet Class contains information for eqtl-gene association, gene identifier, position of SNPs, etc.
#'
#' \describe{
#'    \item{tissue}{character; name of the cell/tissue of the eQTL study.}
#'
#'    \item{snp.id}{character; name of the SNPs}
#'
#'    \item{snp.gr}{GenomicRanges; position of the SNPs}
#'
#'    \item{gene}{character; gene identifier}
#'  }
#' @name eqtlSet
#' @rdname eqtlSet
#' @exportClass eqtlSet

eqtlSet <- setClass("eqtlSet",
                    slot=c(
                      tissue="character",
                      snp.id="character",
                      snp.gr="GenomicRanges",
                      gene="character"))

#' geneSet Class
#'
#' geneSet Class contains information for names of gene sets and a list of gene sets
#'
#' \describe{
#'    \item{description}{character; additional information for gene sets, such as names, URLs, a short description, etc.}
#'
#'    \item{gene.set}{list; a list of gene sets; each member of the list is a vector containing a group of gene identifiers}
#'
#'  }
#' @name geneSet
#' @rdname geneSet
#' @exportClass geneSet

geneSet <- setClass("geneSet",
                    slot=c(
                      description="character",
                      gene.set="list" ))
