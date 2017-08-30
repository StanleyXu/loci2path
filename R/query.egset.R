#' Query enrichment in geneset through one eQTL set.
#'
#' This function perform enrichment test between one eQTL set and a group of
#'  gene sets
#' @param query.gr a GenomicRange object, representing query regions
#' @param query.score optional, set to NULL if the regions are not ordered.
#' @param eqtl.set an eqtlSet object; the eQTL set to be queried against
#' @param gene.set an object of geneSet class; the gene set to be tested
#' @param verbose bool; whether to show eqtlSet/geneSet summary information;
#'  default is FALSE
#' @return a list; \code{result.table} is the major result table showing
#' enrichment assessment;
#'  \code{cover.gene} is the vector showing the genes from the eqtl Sets
#'  covered by the query region(s)
#' @export
#' @examples
#' #build one eqtlset
#' skin.eset=eset.list$Skin
#' #query one egset
#' res.one=query.egset(
#'   query.gr=query.gr,
#'   query.score=NULL,
#'   eqtl.set=skin.eset,
#'   gene.set=biocarta
#' )
#' #enrichment result table
#' res.one$result.table
#' #all the genes associated with eQTLs covered by the query region
#' res.one$cover.gene
query.egset = function(query.gr,
                       query.score,
                       eqtl.set,
                       gene.set,
                       verbose = FALSE) {
    ## check gene id compatibility
    comp = check.geneid(eqtl.set, gene.set)
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
    snp.gr = eqtl.set@snp.gr
    snp.all = length(snp.gr)
    ## get total gene number
    gene.all = gene.set@total.number.gene
    ## overlapping between query regions and all snps
    over.all = as.data.frame(findOverlaps(query.gr, snp.gr))
    snp.q = length(unique(over.all$subjectHits))
    ## output overlapping gene list (new version)
    out.gene.list = unique(eqtl.set@gene[over.all$subjectHits])
    gene.q = length(out.gene.list)
    
    gs = gene.set@gene.set
    lgs = length(gs)
    gene.j = sapply(gs, length)
    #match hit genes with geneset
    ## find unique gene ids
    gg = unique(eqtl.set@gene[over.all$subjectHits])
    over.j.gene = matrix(FALSE, nrow = length(gg), ncol = lgs)
    for (i in 1:ncol(over.j.gene)) {
        over.j.gene[, i] = is.element(gg, gs[[i]])
    }
    ix = match(eqtl.set@gene[over.all$subjectHits], gg)
    over.j.snp = over.j.gene[ix, ]
    ## calculate snp.j: # eQTL snps associated with each pathway, match
    ## eqtl-gene-pathway,
    gg = unique(eqtl.set@gene)  #find unique gene ids
    snp.j = matrix(FALSE, nrow = length(gg), ncol = length(gs))
    for (i in 1:ncol(snp.j)) {
        snp.j[, i] = is.element(gg, gs[[i]])
    }
    ix = match(eqtl.set@gene, gg)
    tt = data.table(snp.j)
    tt = tt[ix, ]
    snp.j = tt[, sapply(.SD, sum)]
    ## summary snp/gene hit
    snp.qj = colSums(over.j.snp)
    gene.qj = colSums(over.j.gene)
    ## get p-val
    pval.fisher.gene = phyper(gene.qj - 1, gene.j, gene.all - gene.j, gene.q,
                              lower.tail = FALSE)
    pval.fisher.snp = phyper(snp.qj - 1, snp.j, snp.all - snp.j, snp.q,
                             lower.tail = FALSE)
    ## add log ratio
    snp.log.ratio = log((snp.qj / snp.q) / (snp.j / snp.all))
    gene.log.ratio = log((gene.qj / gene.q) / (gene.j / gene.all))
    #gene hit
    #tt=apply(over.j.gene, 2, FUN=function(x) out.gene.list[x])
    tt = apply(
        over.j.gene,
        2,
        FUN = function(x)
            if (length(out.gene.list[x])) {
                out.gene.list[x]
            } else{
                NA
            }
    )
    gene.hit = sapply(
        tt,
        FUN = function(x)
            paste(x, collapse = ";")
    )
    ##output
    res = data.frame(
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
        pval.lr = NA,
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
        
        stringsAsFactors = FALSE
    )
    colnames(res) = c(
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
    num_gene_hit = NULL #simply to surpress checking note; not used
    res = subset(res, num_gene_hit > 0)
    
    out.res = list(result.table = res, cover.gene = out.gene.list)
    out.res
}
