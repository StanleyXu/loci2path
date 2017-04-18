#' Query enrichment in geneset through one eQTL set and one gene set
#'
#' This function perform enrichment test between one eQTL set and one gene set
#' @param query.gr a GenomicRange object, representing query regions
#' @param query.score optional, set to NULL if the regions are not ordered.
#' @param eqtl.set an eqtlSet object; the eQTL set to be queried against
#' @param set.j a character vector; a set of genes
#' @param q.all integer
#' @param n.snp.t integer
#' @export
#' @examples
#' #to be added
query.one.set=function(query.gr, query.score=NULL, eqtl.set, set.j, q.all, n.snp.t, gene.all, gene.q){
  snp.set.j.ix=which(eqtl.set@gene %in% set.j)
  snp.set.j.gr=eqtl.set@snp.gr[snp.set.j.ix] # gr of eQTL associated with genes within geneset j
  # intialize pval for all tests as 1
  res=rep(NA,16)

  # check overlapping relationship
  over.j=findOverlaps(query.gr, snp.set.j.gr)

  if(length(over.j) > 0){ # if there are overlapping between query region and geneset-j-associated eQTLs; otherwise return 1

    ## (1) Logistic Regression, use 0/1 as Y
    pval.lr=NA
    if(!(is.null(query.score))){

      o.q=as.data.frame(over.j)$queryHits;
      y=rep(0, length(query.gr)); y[o.q]=1
      x=query.score
      fit=gam(y~x, family="binomial")
      pval.lr=summary(fit)$parametric.anova$'Pr(>F)'[1]
    }

    ## (2) eQTL based enrichment tests

    ##prepare regions
    j.all=length(snp.set.j.gr) # all SNP in pathway j
    q.j=length(unique(as.data.frame(over.j)$subjectHits)) # query region overlapping with pathway-j-associated SNPs
    nq.j=j.all-q.j # all pathway-j-associated SNPs exclude SNPs overlapping with query
    q.nj=q.all-q.j #q.all
    nq.nj=n.snp.t-q.all-nq.j
    ##fisher exact test
    c.t=rbind(c(q.j, q.nj), c(nq.j, nq.nj))
    tt=fisher.test(c.t, alternative="greater")
    pval.fisher=tt$p.value
    log.ratio=log( tt$estimate )

    ##  hypergeom test
    pval.hypergeom=1-phyper(q.j, j.all, n.snp.t-j.all, q.all)

    ## (3) gene based enrichment tests
    #gene.all=31847## need allgene --**temp! remove this line
    gene.j=length(set.j)## total number of gene in geneset j
    #gene.q ## number of gene associated with eQTLs covered by query region
    gene.jq.id=unique(eqtl.set@gene[snp.set.j.ix[as.data.frame(over.j)$subjectHits]])## intersect
    gene.jq=length(gene.jq.id)
    gene.hit=paste(gene.jq.id, collapse = ";")
    gene.nq.j=gene.j-gene.jq
    gene.q.nj=gene.q-gene.jq
    gene.nq.nj=gene.all-gene.j-gene.q.nj
    ##fisher exact test
    c.t=rbind(c(gene.jq, gene.q.nj), c(gene.nq.j, gene.nq.nj))
    tt=fisher.test(c.t, alternative="greater")
    pval.fisher.gene=tt$p.value
    log.ratio.gene=log(tt$estimate)
    ##  hypergeom test
    pval.hypergeom.gene=phyper(gene.jq, gene.j, gene.all, gene.q, lower.tail = F)
    ## (4) combine results:
    res=data.frame(
      length(set.j), # num of genes
      j.all, # num of SNPs associated with gene
      n.snp.t, # num of all GWAS snps in the tissue
      q.all, # num of eQTL overlapping with query region
      q.j, #num of eQTL associated with pathway j, and overlapping with query region
      log.ratio,
      pval.lr,
      pval.fisher,
      pval.hypergeom,
      gene.j, # num of gene in current geneset
      gene.q, # number of genes associated with SNPs overlapped by query region
      gene.jq, # number of genes hit
      gene.hit, # display hit gene ids
      log.ratio.gene, # log ratio based on gene numbers
      pval.fisher.gene,
      pval.hypergeom.gene,

      stringsAsFactors = F
    )

  }#end of if

  res
}
