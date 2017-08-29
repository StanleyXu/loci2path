#' eqtlSet Class
#'
#' eqtlSet Class contains information for eqtl-gene association, gene identifier,
#' position of SNPs, etc.
#'
#' @slot tissue character; name of the cell/tissue of the eQTL study
#' @slot snp.id character; name of the SNPs
#' @slot snp.gr GenomicRanges; position of the SNPs
#' @slot gene character; gene identifier
#'
#' @name eqtlSet-class
#' @rdname eqtlSet-class
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
#' @slot total.number.gene numeric; the total number of all genes;
#'       This number is used in enrichment tests
#' @slot description vector of character; additional information for gene sets,
#'       such as names, URLs, a short description, etc.
#' @slot gene.set list;  a list of gene sets; each member of the list is a vector
#'       containing a group of gene identifiers
#' @name geneSet-class
#' @rdname geneSet-class
#' @exportClass geneSet
geneSet <- setClass("geneSet",
                    slot=c(
                      total.number.gene="numeric",
                      description="character",
                      gene.set="list" ))
