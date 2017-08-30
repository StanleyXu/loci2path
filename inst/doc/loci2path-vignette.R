## ----query_region--------------------------------------------------------
require(GenomicRanges)
bed.file=system.file("extdata", "query/Psoriasis.BED", package = "loci2path")
query.bed=read.table(bed.file, header=FALSE)
colnames(query.bed)=c("chr","start","end")
query.gr=makeGRangesFromDataFrame(query.bed)

## ----eset----------------------------------------------------------------
library(loci2path)
brain.file=system.file("extdata", "eqtl/brain.gtex.txt", 
                       package = "loci2path")
tab=read.table(brain.file, stringsAsFactors = FALSE, header = TRUE)
snp.gr=GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
brain.eset=eqtlSet(tissue="brain",
  snp.id=tab$snp.id,
  snp.gr=snp.gr,
  gene=as.character(tab$gene.entrez.id))
brain.eset

skin.file=system.file("extdata", "eqtl/skin.gtex.txt", package = "loci2path")
tab=read.table(skin.file, stringsAsFactors = FALSE, header = TRUE)
snp.gr=GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
skin.eset=eqtlSet(tissue="skin",
  snp.id=tab$snp.id,
  snp.gr=snp.gr,
  gene=as.character(tab$gene.entrez.id))
skin.eset

blood.file=system.file("extdata", "eqtl/blood.gtex.txt", 
                       package = "loci2path")
tab=read.table(blood.file, stringsAsFactors = FALSE, header = TRUE)
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
biocarta.link.file=system.file("extdata", "geneSet/biocarta.txt", 
                               package = "loci2path")
biocarta.set.file=system.file("extdata", "geneSet/biocarta.set.txt", 
                              package = "loci2path")

biocarta.link=read.delim(biocarta.link.file, header = FALSE, 
                         stringsAsFactors = FALSE)
set.geneid=read.table(biocarta.set.file, stringsAsFactors = FALSE)
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

## ------------------------------------------------------------------------
#query from one eQTL set.
res.one=query.egset(
  query.gr=query.gr,
  query.score=NULL,
  eqtl.set=skin.eset, 
  gene.set=biocarta)

#enrichment result table
res.one$result.table

#all the genes associated with eQTLs covered by the query region
res.one$cover.gene

## ------------------------------------------------------------------------
#query from one eQTL set.
res.esetlist=query.egset.list(
  query.gr=query.gr, 
  query.score=NULL, 
  eqtl.set.list=eset.list, 
  gene.set=biocarta)  

#enrichment result table, tissue column added
res.esetlist$result.table

#all the genes associated with eQTLs covered by the query region; 
#names of the list are tissue names from eqtl set list
res.esetlist$cover.gene

## ------------------------------------------------------------------------
#query from one eQTL set.
res.paral = query.egset.list(
  query.gr = query.gr, 
  query.score = NULL, 
  eqtl.set.list = eset.list, 
  gene.set = biocarta, 
  parallel = TRUE)  
#should return the same result as res.esetlist

## ------------------------------------------------------------------------
result=res.esetlist$result.table

## ------------------------------------------------------------------------
#all the genes associated with eQTLs covered by the query region
res.one$cover.gene

#all the genes associated with eQTLs covered by the query region; 
#names of the list are tissue names from eqtl set list
res.esetlist$cover.gene

## ------------------------------------------------------------------------
tissue.degree=res.get.tissue.degree(
  result, 
  eset.list)

#check gene-tissue mapping for each gene
head(tissue.degree$gene.tissue.map)

#check degree for each gene
head(tissue.degree$gene.tissue.degree)

#average tissue degree for the input result table
tissue.degree$mean.tissue.degree

#add avg. tissue degree to existing result table
res.tissue=data.frame(res.esetlist$result.table,
                      t.degree=tissue.degree$mean.tissue.degree)

## ------------------------------------------------------------------------
#query tissue specificity
gr.tissue=query.tissue(query.gr, eqtl.set.list=eset.list)
gr.tissue

## ------------------------------------------------------------------------
#extract tissue-pathway matrix
mat=res.get.heat.mat(result, test.method = "fisher")

#plot heatmap
draw.heatmap(mat)

## ------------------------------------------------------------------------

#plot word cloud
draw.wordcloud(result)

## ------------------------------------------------------------------------
#plot p-value distribution of result
draw.pval.distribution(result, test.method="fisher")

## ------------------------------------------------------------------------
#obtain geneset description from object
description=get.geneset.description(biocarta, geneset.ids=result$name_pthw)
head(description)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

