#' Extract tissue/geneset enrichment matrix from query result
#'
#' This function extracts the enrichment matrix from eQTL list query result. The rows of the matrixs are pathways; and the columns of the matrixs are tissues/cell lines of the eQTL sets. P-Values from enrichment tests are summarized in this matrix
#' @param res query result from function query.egset.list()
#' @param test.method Choose which enrichment test should be used to retrive p-values from. Options include:"glm"(logistic regression),"fisher"(fisher exact test) 
#' @param filter.quantile Filter option; choose the max quantile of all p-values being kept in the matrix; default is 0.5, which means p-values larger than median p-values are discarded
#' @param min.ptw.gene Filter option; minimum number of genes in a pathway; default is 30 (pathway with <30 genes are not included in the matrix)
#' @keywords result
#' @export
#' @examples
#' res.get.heat.mat(res, test.method="fisher")
res.get.heat.mat=function(res, test.method=c("glm","fisher"), filter.quantile=0.5, min.ptw.gene=30){
  # if(nrow(res)>1000){
  #   res=subset(res, genes_pthw>min.ptw.gene)
  # }
  test.method=match.arg(test.method)
  if(test.method=="glm"){
    pval=res$pval_lr
  }else if(test.method=="fisher"){
    pval=res$pval_fisher
  }
  tissue=res$tissue
  gene=res$name_pthw
  # filter by quantile
  ix=which(pval <= quantile(pval, filter.quantile))
  pval=pval[ix]
  tissue=tissue[ix]
  gene=gene[ix]
  uni.tissue=unique(tissue)
  uni.gene=unique(gene)
  mat=matrix(1,ncol=length(uni.tissue), nrow=length(uni.gene))
  rownames(mat)=uni.gene
  colnames(mat)=uni.tissue
  for(i in 1:length(pval)){
    mat[as.character(gene[i]), as.character(tissue[i])] = pval[i]
  }
  #reorder pthw by avg.pval, row
  if(length(uni.gene)>1){
    rs=rowSums(mat)
    mat=mat[order(rs),]
  }
  #reorder tissue by avg.pval, col
  cs=colSums(mat)
  mat=mat[,order(cs)]

  mat
}
