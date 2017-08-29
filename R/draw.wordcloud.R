#' Draw word cloud to file
#'
#' This function draw the enrichment heatmap to a specific file using pheatmap package.
#' @param result query result from function \code{query.egset.list()}
#' @param min.freq.tissue minimum frequency of tissue/cell to be plotted in the word cloud
#' @param min.freq.gset minimum frequency of geneset to be plotted in the word cloud
#' @param max.words maximum words to be generated
#' @keywords word cloud
#' @export
#' @examples
#' result=query.egset.list(query.gr=query.gr, query.score=NULL,
#'   eqtl.set.list=eset.list, gene.set=biocarta)$result.table
#' mat=res.get.heat.mat(result, test.method = "fisher")
#' draw.wordcloud(result, min.freq.tissue = 2, min.freq.gset = 1)
draw.wordcloud=function(result, min.freq.tissue=5, min.freq.gset=5, max.words=50){
  pthw=sort(table(as.character(result$name_pthw)), decreasing = T)
  tissue=sort(table(as.character(result$tissue)), decreasing = T)
  par(mfrow=c(1,2))
  wordcloud(words = names(pthw), freq = pthw, min.freq = min.freq.gset,
            max.words=max.words, random.order=FALSE, rot.per=0,
            colors=brewer.pal(8, "Dark2"), scale=c(0.8,0.5))
  wordcloud(words = names(tissue), freq = tissue, min.freq = min.freq.tissue,
            max.words=max.words, random.order=FALSE, rot.per=0,
            colors=brewer.pal(8, "Dark2"), scale=c(0.8,0.5))
  par(mfrow=c(1,1))
}

