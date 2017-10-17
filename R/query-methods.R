#' Query enrichment in geneset through multiple eQTL sets.
#'
#' This is the main function for loci2path query. Query can be made
#' on either pathway enrichment or tissue-specificity, depending on the input
#' Class. See \strong{Details} for more.
#' 
#' 
#' The user need to specify
#' \enumerate{
#'     \item Query region; 
#'     \item loci; one or more eQTL set; this is usually more than one eQTL 
#'     set. Only multiple eQTL set derived from different cells/tissues will
#'       show cell/tissue specificity. 
#'     \item path; pre-defined Pathways, or gene sets. the gene sets that 
#'  enrichment tests would be performed to.
#' }
#'
#' \code{loci} must be provided; \code{path} is optional. When \code{path} is
#' missing, the tissue-specificity query for the regions is performed.
#' 
#' The most common case for \code{loci} is an eQTL set list.
#' This function perform enrichment test between one eQTL set and a group of
#'  gene sets. Usually query are based on eQTL set list, rather than only one
#'  eQTL set. Several result exploring functions (\code{getMat}, 
#'  \code{getHeatmap}, \code{getPval}, etc...) are designed for query result 
#'  from eQTL set list and gene sets. The class \code{loci2pathResult} is also
#'  designed for eQTL set list query result only. The result returns a
#'  \code{loci2pathResult} only the class of \code{loci} is a list 
#'  of \code{eqtlSet}.
#'  
#'  If user input one eQTL set as argument \code{loci}, a simple
#'  list object will be returned for specific research purpose.
#'  
#' @param query.gr a GenomicRange object, representing query regions
#' @param loci a list of eqtlSet; each member should be an eqtlSet;
#' Or it can be a single eqtlSet.
#' @param path Pathways or geneSets to be tested for enrichment
#' @param query.score optional, set to NULL if the regions are not ordered.
#' If the query regions are ordered,  query.score is the quantity based
#' on which the regions are ordered
#' @param \dots additional params


#' @rdname query-methods
#' @aliases query,list-method
#' @param N the total number of non-N nucleotides in the genome;
#' default N=2897310462 is for hg19
#' @param verbose bool; whether to show eqtlSet/geneSet summary information;
#'  default is FALSE
#' @return a data.frame showing the tissue enrichment of the query regions
#' by binomial test.
#' @importFrom stats p.adjust pbinom
#' @importFrom S4Vectors to
#' @export
#' @examples
#' gr.tissue <- query(query.gr, eset.list)
setMethod("query", signature = c(query.gr="GenomicRanges", loci="list"), 
          function(query.gr, loci, N=2897310462) {
    #N=2897310462 is the non-N length of hg19;
    #consider change the background length if not hg19
    L <- sum(width(reduce(query.gr)))
    eqtl.set.list=loci
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
    colnames(res.all) <- c("eQTL_gene_in_tissue", 
                           "eQTL_gene_in_query", "pval")
    rownames(res.all) <- names(eqtl.set.list)
    #adjust p-value
    res.all$padj <- p.adjust(res.all$pval)
    res.all <- res.all[order(res.all$padj), ]
    
    res.all
})




















#' @rdname query-methods
#' @aliases query,eqtlSet,geneSet-method
#' @return a list; \code{result.table} is the major result table showing
#' enrichment assessment;
#'  \code{cover.gene} is the vector showing the genes from the eqtl Sets
#'  covered by the query region(s)
#' @importFrom data.table data.table .SD
#' @importFrom stats phyper
#' @export
#' @examples
#' #build one eqtlset
#' skin.eset <- eset.list$Skin
#' #query one egset
#' res.one <- query(query.gr, skin.eset, biocarta)
#' #enrichment result table
#' res.one$result.table
#' #all the genes associated with eQTLs covered by the query region
#' res.one$cover.gene
setMethod("query", signature = c(query.gr="GenomicRanges", 
                                 loci="eqtlSet",
                                 path="geneSet"), 
          function(query.gr,
                   loci,
                   path,
                   query.score=NULL,
                   verbose=FALSE){
    eqtl.set=loci
    gene.set=path
    ## check gene id compatibility
    comp <- check.geneid(eqtl.set, gene.set)
    if (comp[3] == 0) {
        warning("No genes in common between eQTL set and gene set!")
    }
    if (verbose) {
        ## report query region info
        cat(paste0(
            length(query.gr),
            " Query Region(s): ",
            sum(width(reduce(query.gr))),
            "bps in total\n"
        ))
        cat("--eQTL Set:\n")
        print(eqtl.set)
        cat("--gene Set:\n")
        print(gene.set)
    }
    
    
    ## get total snp number
    snp.gr <- eqtlRange(eqtl.set)
    snp.all <- length(snp.gr)
    ## get total gene number
    gene.all <- numGene(gene.set)
    ## overlapping between query regions and all snps
    over.all <- as.data.frame(findOverlaps(query.gr, snp.gr))
    snp.q <- length(unique(over.all$subjectHits))
    ## output overlapping gene list (new version)
    out.gene.list <- unique(eqtlGene(eqtl.set)[over.all$subjectHits])
    gene.q <- length(out.gene.list)
    
    gs <- geneSetList(gene.set)
    lgs <- length(gs)
    gene.j <- sapply(gs, length)
    #match hit genes with geneset
    ## find unique gene ids
    gg <- unique(eqtlGene(eqtl.set)[over.all$subjectHits])
    over.j.gene <- matrix(FALSE, nrow=length(gg), ncol=lgs)
    for (i in seq_len(ncol(over.j.gene))) {
        over.j.gene[, i] <- is.element(gg, gs[[i]])
    }
    ix <- match(eqtlGene(eqtl.set)[over.all$subjectHits], gg)
    over.j.snp <- over.j.gene[ix, , drop=FALSE]
    ## calculate snp.j: # eQTL snps associated with each pathway, match
    ## eqtl-gene-pathway,
    gg <- unique(eqtlGene(eqtl.set))  #find unique gene ids
    snp.j <- matrix(FALSE, nrow=length(gg), ncol=length(gs))
    for (i in seq_len(ncol(snp.j))) {
        snp.j[, i] <- is.element(gg, gs[[i]])
    }
    ix <- match(eqtlGene(eqtl.set), gg)
    tt <- data.table(snp.j)
    tt <- tt[ix, ]
    snp.j <- tt[, sapply(.SD, sum)]
    ## summary snp/gene hit
    snp.qj <- colSums(over.j.snp)
    gene.qj <- colSums(over.j.gene)
    ## get p-val
    pval.fisher.gene <- phyper(gene.qj - 1, gene.j, gene.all - gene.j, gene.q,
                               lower.tail=FALSE)
    pval.fisher.snp <- phyper(snp.qj - 1, snp.j, snp.all - snp.j, snp.q,
                              lower.tail=FALSE)
    ## add log ratio
    snp.log.ratio <- log((snp.qj / snp.q) / (snp.j / snp.all))
    gene.log.ratio <- log((gene.qj / gene.q) / (gene.j / gene.all))
    #gene hit
    #tt=apply(over.j.gene, 2, FUN=function(x) out.gene.list[x])
    tt <- apply(
        over.j.gene,
        2,
        FUN=function(x)
            if (length(out.gene.list[x])) {
                out.gene.list[x]
            } else{
                NA
            }
    )
    gene.hit <- sapply(
        tt,
        FUN=function(x)
            paste(x, collapse=";")
    )
    ##output
    res <- data.frame(
        names(gs),
        snp.j,
        # num of SNPs associated with gene
        rep(snp.all, lgs),
        # num of all GWAS snps in the tissue
        rep(snp.q, lgs),
        # num of eQTL overlapping with query region
        snp.qj,
        #num of eQTL associated with pathway j, overlap query
        snp.log.ratio,
        pval.lr=NA,#need update when glm implemented
        pval.fisher.snp,
        gene.j,
        # num of gene in current geneset
        rep(gene.q, lgs),
        # number of genes associated with SNPs overlap query
        gene.qj,
        # number of genes hit
        gene.hit,
        # display hit gene ids
        gene.log.ratio,
        # log ratio based on gene numbers
        pval.fisher.gene,
        
        stringsAsFactors=FALSE
    )
    colnames(res) <- c(
        "name_pthw",
        "eQTL_pthw",
        "eQTL_total_tissue",
        "eQTL_query",
        "eQTL_pthw_query",
        "log_ratio",
        "pval_lr",
        "pval_fisher",
        "num_gene_set",
        #gene.j,  num of gene in current geneset
        "num_gene_query",
        # gene.q, # number of genes- SNPs overlap query
        "num_gene_hit",
        #gene.jq, # number of genes hit
        "gene_hit",
        #gene.hit, # display hit gene ids
        "log_ratio_gene",
        #log.ratio.gene, # log ratio based on gene numbers
        "pval_fisher_gene" #pval.fisher.gene
    )
    ## remove NA
    num_gene_hit <- NULL #simply to surpress checking note; not used
    res <- subset(res, num_gene_hit > 0)
    
    out.res <- list(result.table=res, cover.gene=out.gene.list)
    out.res
})















#' @rdname query-methods
#' @aliases query,list,geneSet-method
#' @param parallel bool; whether to enable parallel computing;
#' default is FALSE

#' @return a \code{loci2pathResult} class object
#' @seealso loci2pathResult
#' @import BiocParallel
#' @export
#' @examples
#' result <- query(query.gr=query.gr,
#'     loci=eset.list, path=biocarta)
#' #enrichment result table
#' resultTable(result)
#' #all the genes associated with eQTLs covered by the query region
#' coveredGene(result)
setMethod("query", signature = c(query.gr="GenomicRanges", 
                                 loci="list",
                                 path="geneSet"), 
          function(query.gr,
                   loci,
                   path,
                   query.score=NULL,
                   parallel=FALSE,
                   verbose=FALSE) {
    eqtl.set.list=loci
    gene.set=path
    ts <- names(eqtl.set.list)
    tl <- length(ts)
    cat(paste0("Start query: ", tl, " eqtl Sets...\n"))
    
    if (parallel) {
        cat("Run in parallel mode...\n")
        res.list <- bplapply(
            eqtl.set.list,
            FUN=function(x)
                query(query.gr, loci=x, path=gene.set, 
                      query.score=query.score, verbose=verbose),
            BPPARAM=MulticoreParam()
        )
        res <- do.call(rbind, sapply(res.list, "[", 1))
        res <- cbind(tissue=rep(ts, 
                                sapply(res.list,
                                       FUN=function(x) nrow(x$result.table))),
                     res)
        res.gene <- sapply(res.list, "[", 2)
    } else{
        res <- NULL
        res.gene <- list()
        for (i in seq_len(tl)) {
            cat(paste0(i, " of ", tl, ": ", ts[i], "...\n"))
            one.t <- query(
                query.gr=query.gr,
                loci=eqtl.set.list[[i]],
                path=gene.set,
                query.score=query.score
            )
            if (nrow(one.t$result.table) > 0) {
                res <- rbind(res,
                             data.frame(tissue=ts[i], one.t$result.table))
            }
            res.gene[[ts[i]]] <- one.t$cover.gene
        }
    }
    
    #stop if nothing hit
    if (is.null(res)) {
        stop("There's no enrichment detected")
    }
    
    #reorder res
    if (sum(is.na(res$pval_lr)) == 0) {
        ## sort by lr p-value
        res <- res[order(res$pval_lr),]
    } else{
        res <- res[order(res$pval_fisher_gene),]
    }
    rownames(res) <- NULL
    cat(paste0("\ndone!\n"))
    
    out <- loci2pathResult(resultTable=res,
                           coveredGene=res.gene)
    out
})
