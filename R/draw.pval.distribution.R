#' Extract tissue/geneset enrichment p-value distribution from query result
#'
#' This function extracts the enrichment p-value distribution from eQTL list query result.
#' P-values from different tissues/cell types are organized, and QQ-plot is generated against uniform distribution
#'
#' @param res query result from function query.egset.list()
#' @param test.method Choose which enrichment test should be used to retrive p-values from. Options include:"glm"(logistic regression),"fisher"(fisher exact test) and "hypergeom"(hypergeometric test)
#' @keywords result
#' @export
#' @examples
#' result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set.list=eset.list, gene.set=biocarta)$result.table
#' draw.pval.distribution(result, test.method="fisher")
draw.pval.distribution=function(res, test.method=c("glm","fisher")){
  # if(nrow(res)>1000){
  #   res=subset(res, genes_pthw>min.ptw.gene)
  # }
  test.method=match.arg(test.method)
  if(test.method=="glm"){
    pval=res$pval_lr
  }else if(test.method=="fisher"){
    pval=res$pval_fisher_gene
  }else{
    stop("please specify test method")
  }
  tissue=res$tissue
  qq=NULL
  for(ti in unique(tissue)){
    pp=pval[which(tissue==ti)]
    qq=rbind(qq,quantile(pp,seq(0,1, by=0.05), na.rm=T))
  }

  #qq=-log(qq)
  #mc=max(qq[-which(is.infinite(qq))]) + 100
  #qq[which(is.infinite(qq))]=mc

  colors=1:nrow(qq)
  plot(qq[1,], type='l', lty=2, xaxt='n', xlab="", ylab="p-value", col=colors[1])
  text(qq[1,], labels=1, cex = 0.7, col=colors[1])

  for(i in 2:nrow(qq)){
    points(qq[i,], type='l', lty=2, col=colors[i]); text(jitter(qq[i,],1), labels=i, cex = 0.7, col=colors[i])
  }
  legend("topleft", legend=unique(tissue), lty=2, col=colors, bty='n')
  # filter by quantile

}
