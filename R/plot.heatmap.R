#' Generate heatmap of enrichment matrix
#'
#' This function generate the enrichment heatmap using pheatmap package.
#' @param mat query result matrix from function \code{get.heat.mat()}
#' @param main title of the heatmap, default is ""
#' @keywords heatmap
#' @export
#' @examples
#' result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set=eqtl.set.list, gene.set=biocarta)
#' mat=get.heat.mat(result, test.method = "fisher")
#' plot.heatmap(mat, main="enrichment between tissue/pathway")
plot.heatmap=function(mat, file, main="", silent=T){
  mat=-log(mat)
  mc=max(mat[-which(is.infinite(mat))]) + 100
  mat[which(is.infinite(mat))]=mc
  breaks=unique(quantile(mat,1:100/100))
  color = colorRampPalette(rev(brewer.pal(n=7, name ="RdYlBu")))(length(breaks)-1)
  #plot.height= floor((nrow(mat)+100)/20)
  #png(file, height=plot.height, width=6, res=200, unit="in")
  #par(mar=c(4,4,4,4))
  pheatmap(mat,
              show_rownames=T, treeheight_row=0, treeheight_col=0,  cluster_rows=F, cluster_cols=F,
              cellwidth = 10, cellheight = 4, fontsize_row = 4, fontsize_col=6,
              breaks=breaks,color=color,
              legend=F,
              main=main, filename = file
 )
  #grid.text("ylabel example", x=0.5, rot=90, gp=gpar(fontsize=10))
}
