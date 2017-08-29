#' Generate heatmap of enrichment matrix
#'
#' This function generate the enrichment heatmap using pheatmap package.
#' @param mat query result matrix from function \code{get.heat.mat()}
#' @param main title of the heatmap, default is ""
#' @keywords heatmap
#' @export
#' @examples
#' result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set.list=eset.list, gene.set=biocarta)$result.table
#' mat=res.get.heat.mat(result, test.method = "fisher")
#' draw.heatmap(mat, main="enrichment between tissue/pathway")
draw.heatmap=function(mat, main=""){
  mat=-log(mat)
  if(sum(is.infinite(mat))>0){#when pval=0
    mc=max(mat[-which(is.infinite(mat))]) + 100
    mat[which(is.infinite(mat))]=mc
  }
  breaks=unique(quantile(mat,1:100/100))
  color = colorRampPalette(rev(brewer.pal(n=7, name ="RdYlBu")))(length(breaks)-1)
  pheatmap(mat,
              show_rownames=T, treeheight_row=0, treeheight_col=0,  cluster_rows=F, cluster_cols=F,
              #cellwidth = 10, cellheight = 4, fontsize_row = 4, fontsize_col=6,
              breaks=breaks,color=color,
              legend=F,
              main=main
 )
}
