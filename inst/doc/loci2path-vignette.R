## ----query_region--------------------------------------------------------
bed.file=system.file("extdata", "query/Psoriasis.BED", package = "loci2path")
query.bed=read.table(bed.file, header=F)
colnames(query.bed)=c("chr","start","end")
query.gr=makeGRangesFromDataFrame(query.bed)

## ----eset----------------------------------------------------------------
brain.file=system.file("extdata", "eqtl/brain.gtex.txt", package = "loci2path")
tab=read.table(brain.file, stringsAsFactors = F, header = T)
snp.gr=GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
brain.eset=eqtlSet(tissue="brain",
  snp.id=tab$snp.id,
  snp.gr=snp.gr,
  gene=as.character(tab$gene.entrez.id))
brain.eset

skin.file=system.file("extdata", "eqtl/skin.gtex.txt", package = "loci2path")
tab=read.table(skin.file, stringsAsFactors = F, header = T)
snp.gr=GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
skin.eset=eqtlSet(tissue="skin",
  snp.id=tab$snp.id,
  snp.gr=snp.gr,
  gene=as.character(tab$gene.entrez.id))
skin.eset

blood.file=system.file("extdata", "eqtl/blood.gtex.txt", package = "loci2path")
tab=read.table(blood.file, stringsAsFactors = F, header = T)
snp.gr=GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
blood.eset=eqtlSet(tissue="blood",
  snp.id=tab$snp.id,
  snp.gr=snp.gr,
  gene=as.character(tab$gene.entrez.id))
blood.eset

## ----esetlist------------------------------------------------------------
eset.list=list(Brain=brain.eset, Skin=skin.eset, Blood=blood.eset)
eset.list

## ------------------------------------------------------------------------
biocarta.link.file=system.file("extdata", "geneSet/biocarta.txt", package = "loci2path")
biocarta.set.file=system.file("extdata", "geneSet/biocarta.set.txt", package = "loci2path")

biocarta.link=read.delim(biocarta.link.file, header = F, stringsAsFactors = F)
set.geneid=read.table(biocarta.set.file, stringsAsFactors = F)
set.geneid=strsplit(set.geneid[,1], split=",")
names(set.geneid)=biocarta.link[,1]

head(biocarta.link)
head(set.geneid)

## ------------------------------------------------------------------------
#build geneSet
biocarta=geneSet(
  gene.set=set.geneid,
  description=biocarta.link[,2],
  total.number.gene=31847)
biocarta

## ----echo = FALSE--------------------------------------------------------

#show query region
query.gr
#show eqtl set
eqtl.set.list
#show gene set(biocarta)
biocarta

#query
result=query.egset.list(query.gr=query.gr, query.score=NULL, eqtl.set.list=eqtl.set.list, gene.set=biocarta)
head(result)
result.parallel=query.egset.list(query.gr=query.gr, query.score=NULL, eqtl.set=eqtl.set.list, gene.set=biocarta,
                                 parallel = T)
head(result.parallel)


#query tissue specificity
gr.tissue=query.tissue(query.gr, eqtl.set.list=eqtl.set.list)

#reorder result by tissue specificity
result2=res.order.by.tissue(result, tissue.list = rownames(gr.tissue))

#extract tissue/geneset matrix
mat=res.get.heat.mat(result, test.method = "fisher")

#plot heatmap
draw.heatmap(mat)

#plot word cloud
draw.wordcloud(result)

#plot p-value distribution of result
draw.pval.distribution(result, test.method="fisher")

#obtain geneset description from object
description=get.geneset.description(biocarta, geneset.ids=result$name_pthw)
head(description)

