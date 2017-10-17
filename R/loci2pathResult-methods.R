#' Extract tissue degree from query result
#'
#' This function extracts the tissue degree from eQTL list query result for
#' each pathway.
#'
#' @rdname getTissueDegree-methods
#' @aliases getTissueDegree,loci2pathResult-method
#' @param res query result from function query.egset.list()
#' @param loci a list of eqtlSet; each member should be an eqtlSet
#' object
#' @param \dots additional params
#' @return \item{gene.tissue.map}{shows mapping:gene<->tissue}
#' \item{gene.tissue.degree}{shows tissue degree for each gene}
#' \item{mean.tissue.degree}{shows the average tissue digree for each pathway
#'  in the result table}
#' @keywords result
#' @export
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' tissue.degree=getTissueDegree(result, eset.list)
#' head(tissue.degree$gene.tissue.map)
#' head(tissue.degree$gene.tissue.degree)
#' head(tissue.degree$mean.tissue.degree)
setMethod("getTissueDegree", "loci2pathResult", function(res, loci){
    gt <- list()
    eqtl.set.list <- loci
    tissues <- names(eqtl.set.list)
    res <- resultTable(res)
    for (tt in tissues) {
        eset <- eqtl.set.list[[tt]]
        gene <- sort(unique(eqtlGene(eset)))
        gt[gene] <- lapply(
            gt[gene],
            FUN=function(x)
                append(x, tt)
        )
    }
    
    gt.length <- sapply(gt, length)
    
    gg <- res$gene_hit
    mean.tissue <- rep(NA, nrow(res))
    for (i in 1:length(gg)) {
        gl <- strsplit(gg[i], split=";")[[1]]
        mean.tissue[i] <- mean(gt.length[gl])
    }
    res <- list(
        gene.tissue.map=gt,
        gene.tissue.degree=gt.length,
        mean.tissue.degree=mean.tissue
    )
    res 
})



#' Extract tissue/geneset enrichment matrix from query result
#'
#' This function extracts the enrichment matrix from eQTL list query result.
#' The rows of the matrixs are pathways; and the columns of the matrixs are
#' tissues/cell lines of the eQTL sets. P-Values from enrichment tests are
#' summarized in this matrix
#' @rdname getMat-methods
#' @aliases getMat,loci2pathResult-method
#' @param res query result from function query.egset.list()
#' @param test.method Choose which enrichment test should be used to retrive
#' p-values from. Options include:"gene"(default, gene-based fisher's exact 
#' test),"eqtl" (eqtl based fisher's exact test), "glm" (ordered query)
#' @param filter.quantile Filter option; choose the max quantile of all
#' p-values being kept in the matrix; default is 0.5, which means p-values
#' larger than median p-values are discarded
#' @param max.ptw.gene Filter option; minimum number of genes in a pathway;
#' default is 5000 (pathway with >5000 genes are not included in the matrix)
#' @param \dots additional params
#' @return p-value matrix collected from enrichment result table
#' @keywords result
#' @export
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' mat <- getMat(result, test.method="gene")
setMethod("getMat", "loci2pathResult", function(res,
                             test.method=c("gene", "eqtl", "glm"),
                             filter.quantile=0.5,
                             max.ptw.gene=5000) {
    res=resultTable(res)
    if(nrow(res)>max.ptw.gene){
      res <- res[which(res$num_gene_set <= max.ptw.gene), ]
      #res=subset(res, num_gene_set <= max.ptw.gene)
    }
    test.method <- match.arg(test.method)
    if (test.method == "gene") {
        pval <- res$pval_fisher_gene
    } else if (test.method == "eqtl") {
        pval <- res$pval_fisher
    } else if (test.method == "glm"){
        pval <- res$pval_lr
    }
    tissue <- res$tissue
    gene <- res$name_pthw
    # filter by quantile
    ix <- which(pval <= quantile(pval, filter.quantile))
    pval <- pval[ix]
    tissue <- tissue[ix]
    gene <- gene[ix]
    uni.tissue <- unique(tissue)
    uni.gene <- unique(gene)
    mat <- matrix(1,
                  ncol=length(uni.tissue),
                  nrow=length(uni.gene))
    rownames(mat) <- uni.gene
    colnames(mat) <- uni.tissue
    for (i in seq_len(length(pval))) {
        mat[as.character(gene[i]), as.character(tissue[i])] <- pval[i]
    }
    #reorder pthw by avg.pval, row
    if (length(uni.gene) > 1) {
        rs <- rowSums(mat)
        mat <- mat[order(rs), ]
    }
    #reorder tissue by avg.pval, col
    cs <- colSums(mat)
    mat <- mat[, order(cs)]
    
    mat
})



#' Extract description for enriched pathways from query result and 
#' geneSet object 
#'
#' This function extracts the pathway description from geneSet object.
#' @rdname getPathDescription-methods
#' @aliases getPathDescription,loci2pathResult-method
#' @param res query result from function query.egset.list()
#' @param geneset A \code{geneSet} object
#' @return a vector of gene set description from \code{geneSet} description 
#' slot
#' @param \dots additional params
#' @export
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' path.des <- getPathDescription(result, biocarta)
 
setMethod("getPathDescription", "loci2pathResult", function(res, geneset){
    res <- resultTable(res)
    geneset.ids <- res$name_pthw
    ix <- match(geneset.ids, names(geneSetList(geneset)))
    res <- geneset@description[ix]
  
    res  
})


#' Generate heatmap of enrichment matrix from query result
#'
#' This function generate the enrichment heatmap using pheatmap package.
#' @rdname getHeatmap-methods
#' @aliases getHeatmap,loci2pathResult-method
#' @param res query result from function query.egset.list()
#' @param main title of the heatmap, default is ""
#' @param test.method Choose which enrichment test should be used to retrive
#' p-values from. Options include:"gene"(default, gene-based fisher's exact 
#' test),"eqtl" (eqtl based fisher's exact test), "glm" (ordered query)
#' @param filter.quantile Filter option; choose the max quantile of all
#' p-values being kept in the matrix; default is 0.5, which means p-values
#' larger than median p-values are discarded
#' @param max.ptw.gene Filter option; minimum number of genes in a pathway;
#' default is 5000 (pathway with >5000 genes are not included in the matrix)
#' @param \dots additional params
#' @keywords heatmap
#' @return \item{pathways}{frequent pathways}
#' \item{tissues}{frequent tissues}
#' @import pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' getHeatmap(result)
setMethod("getHeatmap", "loci2pathResult", 
          function(res, main="",
                   test.method=c("gene", "eqtl", "glm"),
                   filter.quantile=0.5,
                   max.ptw.gene=5000){
    mat <- getMat(res, test.method, filter.quantile, max.ptw.gene)
    mat <- -log(mat)
    if (sum(is.infinite(mat)) > 0) {
        #when pval=0
        mc <- max(mat[-which(is.infinite(mat))]) + 100
        mat[which(is.infinite(mat))] <- mc
    }
    breaks <- unique(quantile(mat, 1:100 / 100))
    color <-
        colorRampPalette(
            rev(brewer.pal(n=7, name="RdYlBu")))(
                length(breaks) -1)
    res <- pheatmap(
        mat,
        show_rownames=TRUE,
        treeheight_row=0,
        treeheight_col=0,
        cluster_rows=FALSE,
        cluster_cols=FALSE,
        breaks=breaks,
        color=color,
        legend=FALSE,
        main=main
    )
    res 
})


#' Plot word cloud using frequent terms of pathways and genes
#'
#' This function draw the enrichment heatmap using 
#' \code{wordcloud} package.
#' @rdname getWordcloud-methods
#' @aliases getWordcloud,loci2pathResult-method
#' @param res query result from function query.egset.list()
#' @param min.freq.tissue minimum frequency of tissue/cell to be plotted in 
#' the word cloud
#' @param min.freq.gset minimum frequency of geneset to be plotted in the 
#' word cloud
#' @param max.words maximum words to be generated
#' @param \dots additional params
#' @return empty
#' @keywords word cloud
#' @import wordcloud
#' @export
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' getWordcloud(result, min.freq.tissue=2, min.freq.gset=1)
setMethod("getWordcloud", "loci2pathResult", 
          function(res, 
                   min.freq.tissue=5,
                   min.freq.gset=5,
                   max.words=50) {
    result <- resultTable(res)
    pthw <- sort(table(as.character(result$name_pthw)), decreasing=TRUE)
    tissue <- sort(table(as.character(result$tissue)), decreasing=TRUE)
    par(mfrow=c(1, 2))
    wordcloud(
        words=names(pthw),
        freq=pthw,
        min.freq=min.freq.gset,
        max.words=max.words,
        random.order=FALSE,
        rot.per=0,
        colors=brewer.pal(8, "Dark2"),
        scale=c(0.8, 0.5)
    )
    wordcloud(
        words=names(tissue),
        freq=tissue,
        min.freq=min.freq.tissue,
        max.words=max.words,
        random.order=FALSE,
        rot.per=0,
        colors=brewer.pal(8, "Dark2"),
        scale=c(0.8, 0.5)
    )
    par(mfrow=c(1, 1))
    res <- list(
        pathways=pthw,
        tissues=tissue
    )
    res
})


#' Extract tissue/geneset enrichment p-value distribution from query result
#'
#' This function extracts the enrichment p-value distribution from eQTL list 
#'  query result. P-values from different tissues/cell types are organized, 
#'  and QQ-plot is generated against uniform distribution
#' @rdname getPval-methods
#' @aliases getPval,loci2pathResult-method
#' @param res query result from function query.egset.list()
#' @param test.method Choose which enrichment test should be used to retrive
#' p-values from. Options include:"gene"(default, gene-based fisher's exact 
#' test),"eqtl" (eqtl based fisher's exact test), "glm" (ordered query)
#' @param \dots additional params
#' @keywords result
#' @return generate pval distribution plot
#' @importFrom stats quantile
#' @import graphics
#' @export
#' @examples
#' result <- query(query.gr=query.gr, 
#'   loci=eset.list, path=biocarta)
#' getPval(result, test.method="gene")
setMethod("getPval", "loci2pathResult", 
          function(res, 
                   test.method=c("gene", "eqtl", "glm")) {
    res <- resultTable(res)
    test.method <- match.arg(test.method)
    if (test.method == "gene") {
        pval <- res$pval_fisher_gene
    } else if (test.method == "eqtl") {
        pval <- res$pval_fisher
    } else if (test.method == "glm"){
        pval <- res$pval_lr
    }
    tissue <- res$tissue
    qq <- NULL
    qrowname <- NULL
    for (ti in unique(tissue)) {
        pp <- pval[which(tissue == ti)]
        qq <- rbind(qq, quantile(pp, seq(0, 1, by=0.05), na.rm=TRUE))
        qrowname <- c(qrowname, ti)
    }
    
    #qq=-log(qq)
    #mc=max(qq[-which(is.infinite(qq))]) + 100
    #qq[which(is.infinite(qq))]=mc
    
    rownames(qq) <- qrowname
    colors <- 1:nrow(qq)
    plot(
        qq[1, ],
        ylim=c(0, max(qq)*1.2),
        type='l',
        lty=2,
        xaxt='n',
        xlab="",
        ylab="p-value",
        col=colors[1]
    )
    text(qq[1, ],
         labels=1,
         cex=0.7,
         col=colors[1])
    
    for (i in 2:nrow(qq)) {
        points(
            qq[i, ],
            type='l',
            lty=2,
            col=colors[i])
        
        text(
            jitter(qq[i, ], 1),
            labels=i,
            cex=0.7,
            col=colors[i]
        )
    }
    legend(
        "topleft",
        legend=unique(tissue),
        lty=2,
        col=colors,
        bty='n'
    )
    # filter by quantile
   
    qq 
})

