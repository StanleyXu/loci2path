#' eqtlSet Class
#'
#' eqtlSet Class contains information for eqtl-gene association,
#'  gene identifier, position of SNPs, etc.
#'
#' @slot tissue character; name of the cell/tissue of the eQTL study
#' @slot eqtlId character; name of the SNPs
#' @slot eqtlRange GenomicRanges; position of the SNPs
#' @slot gene character; gene identifier
#'
#' @name eqtlSet-class
#' @rdname eqtlSet-class
#' @import GenomicRanges
#' @importFrom methods new
#' @exportClass eqtlSet
#' @examples
#' require(GenomicRanges)
#' brain.file <- system.file("extdata", "eqtl/brain.gtex.txt", 
#'     package="loci2path")
#' tab <- read.table(brain.file, stringsAsFactors=FALSE, header=TRUE)
#' eqtlRange <- GRanges(seqnames=Rle(tab$snp.chr), 
#'     ranges=IRanges(start=tab$snp.pos, 
#'     width=1))
#' brain.eset <- eqtlSet(tissue="brain",
#'     eqtlId=tab$snp.id,
#'     eqtlRange=eqtlRange,
#'     gene=as.character(tab$gene.entrez.id))
#' tissue(brain.eset)
#' head(eqtlId(brain.eset))
#' eqtlRange(brain.eset)
#' head(eqtlGene(brain.eset))
#' @export eqtlSet
eqtlSet <- setClass(
    "eqtlSet",
    slot = c(
        tissue="character",
        eqtlId="character",
        eqtlRange="GenomicRanges",
        gene="character"
    )
)

#' geneSet Class
#'
#' geneSet Class contains information for names of gene sets and a list of 
#' gene sets
#'
#' @slot numGene numeric; the total number of all genes;
#'       This number is used in enrichment tests
#' @slot description vector of character; additional information for gene 
#'   sets, such as names, URLs, a short description, etc.
#' @slot geneSetList list;  a list of gene sets; each member is a 
#'   vector containing a group of gene identifiers
#' @name geneSet-class
#' @rdname geneSet-class
#' @exportClass geneSet
#' @examples
#' biocarta.link.file <- system.file("extdata", 
#'     "geneSet/biocarta.txt", package="loci2path")
#' biocarta.link <- read.delim(biocarta.link.file, header=FALSE, 
#'     stringsAsFactors=FALSE)
#' biocarta.set.file <- system.file("extdata", "geneSet/biocarta.set.txt", 
#'     package="loci2path")
#' set.geneid <- read.table(biocarta.set.file, stringsAsFactors=FALSE)
#' set.geneid <- strsplit(set.geneid[,1], split=",")
#' names(set.geneid) <- biocarta.link[,1]
#' biocarta <- geneSet(
#'     geneSetList=set.geneid,
#'     description=biocarta.link[,2],
#'     numGene=31847)
#' numGene(biocarta)
#' head(description(biocarta))
#' head(geneSetList(biocarta))
#' @export geneSet
geneSet <- setClass(
    "geneSet",
    slot=c(
        numGene="numeric",
        description="character",
        geneSetList="list"
    )
)




#' loci2pathResult Class
#'
#' loci2pathResult Class contains information for the query result from 
#' query function \code{query}. Result object contains a ranked
#' pathway table, and a vector of gene names that are associated with loci 
#' covered by query regions
#'
#' @slot resultTable data.frame; contains enrichment statistics, 
#' summary of eQTL and gene numbers, pathway names and gene names, etc.
#' @slot coveredGene list; each member is a vector of genes associated 
#' with one tissue, whose associating loci are covered by query regions
#' @name loci2pathResult-class
#' @rdname loci2pathResult-class
#' @exportClass loci2pathResult
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' result
#' resultTable(result) # a data.frame for enriched pathways
#' coveredGene(result)  
#' @export loci2pathResult
loci2pathResult <- setClass(
    "loci2pathResult",
    slot=c(
        resultTable="data.frame",
        coveredGene="list"
    )
)
