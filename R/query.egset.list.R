#' Query enrichment in geneset through multiple eQTL sets.
#'
#' This is the main function for query. The user need to specify
#' - (1) Query region; if it is ordered, together with the Query score based
#' on which the regions are ranked
#' - (2) eQTL set list; this is usually more than one eQTL set.
#'       Only multiple eQTL set derived from different cells/tissues will
#'       show cell/tissue specificity
#' - (3) gene set; the gene sets that enrichment tests would be performed to.
#' @param query.gr a GenomicRange object, representing query regions
#' @param query.score optional, set to NULL if the regions are not ordered.
#' @param eqtl.set.list a list of eqtlSet; each member should be an eqtlSet
#' object
#' @param gene.set an object of geneSet class; the gene set to be tested
#' @param parallel bool; whether to enable parallel computing;
#' default is FALSE
#' @param verbose bool; whether to show more information during query;
#' default is FALSE
#' @return a list; \code{result.table} is the major result table showing
#' enrichment assessment;
#'  \code{cover.gene} is the list showing the genes from the eqtl Sets
#'  covered by the query region(s)
#' @export
#' @examples
#' result = uery.egset.list(query.gr=query.gr, query.score=NULL,
#'     eqtl.set.list=eset.list, gene.set=biocarta)
#' #enrichment result table
#' result$result.table
#' #all the genes associated with eQTLs covered by the query region
#' result$cover.gene
query.egset.list = function(query.gr,
                            query.score,
                            eqtl.set.list,
                            gene.set,
                            parallel = FALSE,
                            verbose = FALSE) {
    ts = names(eqtl.set.list)
    tl = length(ts)
    cat(paste0("Start query: ", tl, " eqtl Sets...\n"))
    
    if (parallel) {
        cat("Run in parallel mode...\n")
        res.list = bplapply(
            eqtl.set.list,
            FUN = function(x)
                query.egset(query.gr, query.score, x, gene.set, verbose),
            BPPARAM = MulticoreParam()
        )
        res = do.call(rbind, sapply(res.list, "[", 1))
        res = cbind(tissue = rep(ts, sapply(
            res.list,
            FUN = function(x)
                nrow(x$result.table)
        )),
        res)
        res.gene = sapply(res.list, "[", 2)
    } else{
        res = NULL
        res.gene = list()
        for (i in 1:tl) {
            cat(paste0(i, " of ", tl, ": ", ts[i], "...\n"))
            one.t = query.egset(
                query.gr = query.gr,
                query.score = query.score,
                eqtl.set = eqtl.set.list[[i]],
                gene.set = gene.set
            )
            if (nrow(one.t$result.table) > 0) {
                res = rbind(res,
                            data.frame(tissue = ts[i], one.t$result.table))
            }
            res.gene[[ts[i]]] = one.t$cover.gene
        }
    }
    
    #stop if nothing hit
    if (is.null(res)) {
        stop("There's no enrichment detected")
    }
    
    #reorder res
    if (sum(is.na(res$pval_lr)) == 0) {
        ## sort by lr p-value
        res = res[order(res$pval_lr),]
    } else{
        res = res[order(res$pval_fisher_gene),]
    }
    rownames(res) = NULL
    cat(paste0("\ndone!\n"))
    
    out = list(result.table = res, cover.gene = res.gene)
}
