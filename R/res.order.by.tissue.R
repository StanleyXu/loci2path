#' Order result by tissue-specificity
#'
#' This function re-rank the query result by the tissue specificity (from query.tissue())
#' @param result query result from function query.egset.list()
#' @param tissue.list a list of tissues to be re-ordered
#' @export
#' @examples
#' #to be added
res.order.by.tissue=function(result, tissue.list){
  res=NULL
  for(i in 1:length(tissue.list)){
    res=rbind(res, result[which(result$tissue==tissue.list[i]),])
  }
  res
}

