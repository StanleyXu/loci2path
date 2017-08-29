#' Extract tissue degree from query result
#'
#' This function extracts the tissue degree from eQTL list query result for each pathway. 
#' 
#' @param res query result from function query.egset.list()
#' @param eqtl.set.list a list of eqtlSet; each member should be an eqtlSet object
#' @return a list; {gene.tissue.map} shows
#' @keywords result
#' @export
#' @examples
#' result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set.list=eset.list, gene.set=biocarta)$result.table
#' tissue.degree=res.get.tissue.degree(result, eset.list)
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
