#' Extract tissue/geneset enrichment matrix from query result
#'
#' This function extracts the enrichment matrix from eQTL list query result. The rows of the matrixs are pathways; and the columns of the matrixs are tissues/cell lines of the eQTL sets. P-Values from enrichment tests are summarized in this matrix
#' @param res query result from function query.egset.list()
#' @param eqtl.set.list a list of eqtlSet; each member should be an eqtlSet object
#' @keywords result
#' @export
#' @examples
#' result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set.list=eqtl.set.list, gene.set=biocarta)
#' tissue.dgree=res.get.tissue.degree(result, eqtl.set.list)
#' head(tissue.degree$gene.tissue.map)
#' head(tissue.degree$gene.tissue.degree)
#' head(tissue.degree$mean.tissue.digree)
res.get.tissue.degree=function(res, eqtl.set.list){
  gt=list()
  tissues=names(eqtl.set.list)
  for(tt in tissues){
    eset=eqtl.set.list[[tt]]
    gene=sort(unique(eset@gene))
    gt[gene]=lapply(gt[gene], FUN=function(x) append(x, tt))
  }
  
  gt.length=sapply(gt, length)
  
  gg=res$gene_hit
  mean.tissue=rep(NA, nrow(res))
  for(i in 1:length(gg)){
    gl=strsplit(gg[i], split=";")[[1]]
    mean.tissue[i]=mean(gt.length[gl])
  }
  res=list(gene.tissue.map=gt, gene.tissue.degree=gt.length, mean.tissue.degree=mean.tissue)
  res
}
