#' eqtlSet Class
#'
#' eqtlSet Class contains information for eqtl-gene association,
#'  gene identifier, position of SNPs, etc.
#'
#' @slot tissue character; name of the cell/tissue of the eQTL study
#' @slot snp.id character; name of the SNPs
#' @slot snp.gr GenomicRanges; position of the SNPs
#' @slot gene character; gene identifier
#'
#' @name eqtlSet-class
#' @rdname eqtlSet-class
#' @exportClass eqtlSet
#' @examples
#' brain.file = system.file("extdata", "eqtl/brain.gtex.txt", 
#'     package = "loci2path")
#' tab = read.table(brain.file, stringsAsFactors = FALSE, header = TRUE)
#' snp.gr = GRanges(seqnames = Rle(tab$snp.chr), 
#'     ranges=IRanges(start=tab$snp.pos, 
#'     width=1))
#' brain.eset = eqtlSet(tissue="brain",
#'     snp.id=tab$snp.id,
#'     snp.gr=snp.gr,
#'     gene=as.character(tab$gene.entrez.id))
eqtlSet <- setClass(
    "eqtlSet",
    slot = c(
        tissue = "character",
        snp.id = "character",
        snp.gr = "GenomicRanges",
        gene = "character"
    )
)


#' geneSet Class
#'
#' geneSet Class contains information for names of gene sets and a list of 
#' gene sets
#'
#' @slot total.number.gene numeric; the total number of all genes;
#'       This number is used in enrichment tests
#' @slot description vector of character; additional information for gene 
#'   sets, such as names, URLs, a short description, etc.
#' @slot gene.set list;  a list of gene sets; each member of the list is a 
#'   vector containing a group of gene identifiers
#' @name geneSet-class
#' @rdname geneSet-class
#' @exportClass geneSet
#' @examples
#' biocarta.link.file = system.file("extdata", 
#'     "geneSet/biocarta.txt", package = "loci2path")
#' biocarta.link = read.delim(biocarta.link.file, header = FALSE, 
#'     stringsAsFactors = FALSE)
#' biocarta.set.file = system.file("extdata", "geneSet/biocarta.set.txt", 
#'     package = "loci2path")
#' set.geneid = read.table(biocarta.set.file, stringsAsFactors = FALSE)
#' set.geneid = strsplit(set.geneid[,1], split=",")
#' names(set.geneid) = biocarta.link[,1]
#' biocarta = g eneSet(
#'     gene.set = set.geneid,
#'     description = biocarta.link[,2],
#'     total.number.gene = 31847)
geneSet <- setClass(
    "geneSet",
    slot = c(
        total.number.gene = "numeric",
        description = "character",
        gene.set = "list"
    )
)
