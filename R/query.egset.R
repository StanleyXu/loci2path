#' Query enrichment in geneset through one eQTL set.
#'
#' This function perform enrichment test between one eQTL set and a group of gene sets
#' @param query.gr a GenomicRange object, representing query regions
#' @param query.score optional, set to NULL if the regions are not ordered.
#' @param eqtl.set an eqtlSet object; the eQTL set to be queried against
#' @param gene.set an object of geneSet class; the gene set to be tested
#' @param parallel bool; whether to enable parallel computing;  default is F
#' @param verbose bool; whether to show eqtlSet/geneSet summary information; default is F
#' @export
#' @examples
#' #to be added
query.egset2=function(query.gr, query.score, eqtl.set, gene.set, parallel=F, verbose=F){
  ## check gene id compatibility
  comp=check.geneid(eqtl.set, gene.set)
  if(comp[3]==0){
    warning("No genes in common between eQTL set and gene set!")
  }
  if(verbose){
    ## report query region info
    cat(paste0(length(query.gr), " Query Region(s): ", sum(width(reduce(query.gr))),"bps in total\n"))
    cat("--eQTL Set:\n")
    print(eqtl.set)
    cat("--gene Set:\n")
    print(gene.set)
  }
  ## get total snp number
  snp.gr=eqtl.set@snp.gr
  n.snp.t=length(snp.gr)
  ## get total gene number
  gene.all=gene.set@total.number.gene
  ## overlapping between query regions and all snps
  over.all=findOverlaps(query.gr, snp.gr)
  q.all=length(unique(as.data.frame(over.all)$subjectHits))
  gene.q=length(unique(eqtl.set@gene[as.data.frame(over.all)$subjectHits]))
  #init
  res=NULL
  ## query one.set, loop across all sets
  if(parallel==T){
    res.list=bplapply(gene.set@gene.set, FUN=function(x){
      query.one.set(query.gr, query.score, eqtl.set, x, q.all, n.snp.t, gene.all, gene.q)})
  }else{
    res.list=lapply(gene.set@gene.set, FUN=function(x){
      query.one.set2(query.gr, query.score, eqtl.set, x, q.all, n.snp.t, gene.all, gene.q)})    
  }

  ## remove NA
  res.list=res.list[-which(is.na(res.list))]
  if(length(res.list)>0){
    tt=NULL
    for(i in 1:length(res.list)){
      tt=rbind(tt, res.list[[i]])
    }
    ## organize result
    res=data.frame(names(res.list), tt, stringsAsFactors = F)
    
    pval.fisher=phyper(res$q.j-1, res$j.all, n.snp.t-res$j.all, res$q.all, lower.tail=F)
    log.ratio=log((res$q.j/res$q.all)/(res$j.all/n.snp.t))
    pval.fisher.gene=phyper(res$gene.jq-1, res$gene.j, gene.all-res$gene.j, res$gene.q, lower.tail = F)
    log.ratio.gene=log((res$gene.jq/res$gene.q)/(res$gene.j/gene.all))
    
    res=cbind(res[,1:5], log.ratio, res[,6],pval.fisher, res[,7:10], log.ratio.gene, pval.fisher.gene)
    colnames(res)=c(
      "name_pthw",
      #"genes_pthw", #removed; duplicate with num_gene_set
      "eQTL_pthw",
      "eQTL_total_tissue",
      "eQTL_query",
      "eQTL_pthw_query",
      "log_ratio",
      "pval_lr",
      "pval_fisher",
      #"pval_hypergeom", #removed; redundant
      "num_gene_set", #gene.j,  num of gene in current geneset
      "num_gene_query", # gene.q, # number of genes associated with SNPs overlapped by query region
      "num_gene_hit", #gene.jq, # number of genes hit
      "gene_hit", #gene.hit, # display hit gene ids
      "log_ratio_gene", #log.ratio.gene, # log ratio based on gene numbers
      "pval_fisher_gene" #pval.fisher.gene
      #"pval_hypergeom_gene" #pval.hypergeom.gene #removed
    )
  }

  res
}
