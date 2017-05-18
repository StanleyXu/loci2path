#' Query enrichment in geneset through multiple eQTL sets.
#'
#' This is the main function for query. The user need to specify
#' - (1) Query region; if it is ordered, together with the Query score based on which the regions are ranked
#' - (2) eQTL set list; this is usually more than one eQTL set.
#'       Only multiple eQTL set derived from different cells/tissues will show cell/tissue specificity
#' - (3) gene set; the gene sets that enrichment tests would be performed to.
#' @param query.gr a GenomicRange object, representing query regions
#' @param query.score optional, set to NULL if the regions are not ordered.
#' @param eqtl.set.list a list of eqtlSet; each member should be an eqtlSet object
#' @param gene.set an object of geneSet class; the gene set to be tested
#' @param parallel bool; whether to enable parallel computing;  default is F
#' @export
#' @examples
#' #to be added
query.egset.list=function(query.gr, query.score, eqtl.set.list, gene.set, parallel=F){
  ts=names(eqtl.set.list)
  tl=length(ts)
  cat(paste0("Start query: ", tl, " eqtl Sets...\n"))
  res=NULL
  for(i in 1:tl){
    cat(paste0(i, " of ", tl, ": ",ts[i], "...\n"))
    one.t=query.egset(query.gr=query.gr, query.score=query.score,
                      eqtl.set=eqtl.set.list[[i]], gene.set=gene.set,
                      parallel=parallel)
    if(length(one.t)>0){
      res=rbind(res,data.frame(tissue=ts[i], one.t))
    }

  }
  #stop if nothing hit
  if(is.null(res)){
    message("There's no enrichment detected")
    return(res)
  }
  #reorder res
  if(sum(is.na(res$pval_lr))==0){## sort by lr p-value
    res=res[order(res$pval_lr),]
  }else{
    res=res[order(res$pval_fisher_gene),]
  }
  rownames(res)=NULL
  cat(paste0("\ndone!\n"))

  res
}
