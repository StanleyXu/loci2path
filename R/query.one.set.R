
query.one.set=function(query.gr, query.score, eqtl.set, set.j, lr=T, q.all, n.snp.t){
  snp.set.j.ix=which(eqtl.set@gene %in% set.j)
  snp.set.j.gr=eqtl.set@snp.gr[snp.set.j.ix] # gr of eQTL associated with genes within geneset j
  # intialize pval for all tests as 1
  res=rep(NA,9)

  # check overlapping relationship
  over.j=findOverlaps(query.gr, snp.set.j.gr)

  if(length(over.j) > 0){ # if there are overlapping between query region and geneset-j-associated eQTLs; otherwise return 1
    j.all=length(snp.set.j.gr) # all SNP in pathway j
    # query region overlapping with all SNP in the tissue

    q.j=length(unique(as.data.frame(over.j)$subjectHits)) # query region overlapping with pathway-j-associated SNPs
    nq.j=j.all-q.j # all pathway-j-associated SNPs exclude SNPs overlapping with query
    q.nj=q.all-q.j #q.all
    nq.nj=n.snp.t-q.all-nq.j

    ## (1) Logistic Regression, use 0/1 as Y
    pval.lr=NA
    if(lr==T){

      o.q=as.data.frame(over.j)$queryHits;
      y=rep(0, length(query.gr)); y[o.q]=1
      x=query.score
      fit=gam(y~x, family="binomial")
      pval.lr=summary(fit)$parametric.anova$'Pr(>F)'[1]
    }

    ## (2) fisher exact test
    c.t=rbind(c(q.j, q.nj), c(nq.j, nq.nj))
    log.ratio=log( (q.j/q.nj) / (nq.j/nq.nj) )
    tt=fisher.test(c.t, alternative="greater")
    pval.fisher=tt$p.value

    ## (3) hypergeom test
    pval.hypergeom=1-phyper(q.j, j.all, n.snp.t-j.all, q.all)

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
      pval.hypergeom
    )

  }#end of if

  res
}
