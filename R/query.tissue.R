#' Query enrichment in eQTL tissues through multiple eQTL sets.
#'
#' This function explore the tissue specificity for query regions. The user
#' need to specify
#' - (1) Query region
#' - (2) eQTL set list; this is usually more than one eQTL set.
#'       Only multiple eQTL set derived from different cells/tissues will
#'       show cell/tissue specificity
#' @param query.gr a GenomicRange object, representing query regions
#' @param eqtl.set.list a list of eqtlSet; each member should be an eqtlSet
#' object
#' @param N the total number of non-N nucleotides in the genome;
#' default N=2897310462 is for hg19
#' @return a data.frame showing the tissue enrichment of the query regions
#' by binomial test.
#' @importFrom stats p.adjust pbinom
#' @export
#' @examples
#' gr.tissue <- query.tissue(query.gr, eqtl.set.list=eset.list)
query.tissue <- function(query.gr, eqtl.set.list, N=2897310462) {
    #N=2897310462 is the non-N length of hg19;
    #consider change the background length if not hg19
    L <- sum(width(reduce(query.gr)))
    
    res.all <- NULL
    for (i in seq_len(length(eqtl.set.list))) {
        ee <- eqtl.set.list[[i]]
        t.gene <- length(unique(eqtlGene(ee)))
        p <- t.gene / N
        over <- findOverlaps(query.gr, eqtlRange(ee))
        hit <- unique(to(over))
        hit.gene <- length(unique(eqtlGene(ee)[hit]))
        pp <- pbinom(hit.gene, L, p, lower.tail=FALSE)
        res=c(t.gene, hit.gene, pp)
        res.all=rbind(res.all, res)
    }
    res.all <- as.data.frame(res.all)
    colnames(res.all) <- c("eQTL_gene_in_tissue", "eQTL_gene_in_query", "pval")
    rownames(res.all) <- names(eqtl.set.list)
    #adjust p-value
    res.all$padj <- p.adjust(res.all$pval)
    res.all <- res.all[order(res.all$padj), ]
    
    res.all
}
