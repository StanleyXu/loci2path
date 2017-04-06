#' Export heatmap to file
#'
#' This function exports the enrichment heatmap to a specific file using pheatmap package.
#' @param mat query result from function query.egset.list()
#' @param file path to the image file
#' @param main title of the heatmap, default is ""
#' @keywords heatmap
#' @export
#' @examples
#' export.heatmap(mat, file="result.heatmap.png", main="enrichment between tissue/pathway")
export.heatmap=function(mat, file, main=""){
  library(pheatmap)
  library(RColorBrewer)

  mat=-log(mat)
  mc=max(mat[-which(is.infinite(mat))]) + 100
  mat[which(is.infinite(mat))]=mc
  breaks=unique(quantile(mat,1:100/100))
  color = colorRampPalette(rev(brewer.pal(n=7, name ="RdYlBu")))(length(breaks)-1)
  #plot.height= floor((nrow(mat)+100)/20)
  #png(file, height=plot.height, width=6, res=200, unit="in")
  #par(mar=c(4,4,4,4))
  ph=pheatmap(mat,
              show_rownames=T, treeheight_row=0, treeheight_col=0,  cluster_rows=F, cluster_cols=F,
              cellwidth = 10, cellheight = 4, fontsize_row = 4, fontsize_col=6,
              breaks=breaks,color=color,
              legend=F,
              main=main, filename = file)
  #grid.text("ylabel example", x=0.5, rot=90, gp=gpar(fontsize=10))
}
