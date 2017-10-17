---

# loci2path
  Annotating a given genomic locus or a set of genomic loci is an important
  yet challenging task. This is especially true for the non-coding part of the
  genome which is enormous yet poorly understood. Since gene set enrichment
  analyses have demonstrated to be effective approach to annotate a set of
  genes, this idea can be extended to explore the enrichment of functional
  elements or features in a set of genomic intervals to reveal potential
  functional connections. `loci2path` is a novel computational
  strategy that takes advantage of the newly emerged, genome-wide and
  tissue-specific expression quantitative trait loci (eQTL) information to
  help annotate a set of genomic intervals in terms of transcription
  regulation. By checking the presence or absence of millions of eQTLs in the
  set of genomic intervals of interest, loci2path build a bridge connecting
  genomic intervals to biological pathway or pre-defined biological-meaningful
  gene sets. Our method enjoys two key advantages over existing methods:
  first, we no longer rely on proximity to link a locus to a gene which has
  shown to be unreliable; second, eQTL allows us to provide the regulatory
  annotation under the context of specific tissue types which is important.
  
---


# Install

```
library(devtools)
install_github(stanleyxu/loci2path)
```
 * The Bioconductor package is under review process; Once released you can also install loci2path using `biocLite()` function from bioconductor.


# Query loci2path with genomic regions

## 1. Query regions

`loci2path` takes query regions in the format of `GenomicRanges`. Only the
Genomic Locations (chromosomes, start and end position) will be used. Strand
information and other metadata columns are ignored.
In the demo data, 47 regions associated with Psoriasis disease were downloaded
from **immunoBase.org** and used as demo query regions. 

```{r query_region}
require(GenomicRanges)
bed.file <- system.file("extdata", "query/Psoriasis.BED", package = "loci2path")
query.bed <- read.table(bed.file, header=FALSE)
colnames(query.bed) <- c("chr","start","end")
query.gr <- makeGRangesFromDataFrame(query.bed)
```

## 2. Prepare eQTL sets.

eQTL sets are entities recording 1-to-1 links between eQTL SNPs and genes.
eQTL set entity also contains the following information: tissue name for the
eQTL study, IDs and genomic ranges for the eQTL SNPs, IDs for the associated
genes.

eQTL set can be constructed manually by specifying the corresponding
information in each slot. 

eQTL set list is a list of multiple eQTL sets, usually collected from
different tissues.

Below is an example to construct customized eQTL set and eQTL set list using
demo data files. In the demo data folder, three eQTL sets downloaded from GTEx
project are included. Due to the large size, each eQTL dataset is down sampled
to 3000 records for demostration purpose. 

#### 2.1 construct eQTL set
```{r eset}
library(loci2path)
library(GenomicRanges)
brain.file <- system.file("extdata", "eqtl/brain.gtex.txt", 
                       package="loci2path")
tab <- read.table(brain.file, stringsAsFactors=FALSE, header=TRUE)
snp.gr <- GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
brain.eset <- eqtlSet(tissue="brain",
  eqtlId=tab$snp.id,
  eqtlRange=snp.gr,
  gene=as.character(tab$gene.entrez.id))
brain.eset

skin.file <- system.file("extdata", "eqtl/skin.gtex.txt", package="loci2path")
tab=read.table(skin.file, stringsAsFactors=FALSE, header=TRUE)
snp.gr <- GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
skin.eset <- eqtlSet(tissue="skin",
  eqtlId=tab$snp.id,
  eqtlRange=snp.gr,
  gene=as.character(tab$gene.entrez.id))
skin.eset

blood.file <- system.file("extdata", "eqtl/blood.gtex.txt", 
                       package="loci2path")
tab <- read.table(blood.file, stringsAsFactors=FALSE, header=TRUE)
snp.gr <- GRanges(seqnames=Rle(tab$snp.chr), 
  ranges=IRanges(start=tab$snp.pos, 
  width=1))
blood.eset <- eqtlSet(tissue="blood",
  eqtlId=tab$snp.id,
  eqtlRange=snp.gr,
  gene=as.character(tab$gene.entrez.id))
blood.eset
```

#### 2.2 construct eQTL set list

```{r esetlist}
eset.list <- list(Brain=brain.eset, Skin=skin.eset, Blood=blood.eset)
eset.list
```



## 3. Prepare gene set collection

A geneset collection contains a list of gene sets, with each gene set is
represented as a vector of member genes. A vector of description is also
provided as the metadata slot for each gene set. The total number of gene in
the geneset collection is also required to perform the enrichment test. In
this tutorial the BIOCARTA pathway collection was downloaded from MSigDB. 

```{r}
biocarta.link.file <- system.file("extdata", "geneSet/biocarta.txt", 
                               package="loci2path")
biocarta.set.file <- system.file("extdata", "geneSet/biocarta.set.txt", 
                              package="loci2path")

biocarta.link <- read.delim(biocarta.link.file, header=FALSE, 
                         stringsAsFactors=FALSE)
set.geneid <- read.table(biocarta.set.file, stringsAsFactors=FALSE)
set.geneid <- strsplit(set.geneid[,1], split=",")
names(set.geneid) <- biocarta.link[,1]

head(biocarta.link)
head(set.geneid)
```

In order to build gene set, we also need to know the total number of genes in
order to perform enrichment test. In this study, the total number of gene in
MSigDB pathway collection is 31,847(Liberzon et. al, 2015). 

```{r}
#build geneSet
biocarta <- geneSet(
    numGene=31847,
    description=biocarta.link[,2],
    geneSetList=set.geneid)
biocarta
```


## 4.Perform query

####  4.1 peroform query from one eQTL set
```{r}
res.one <- query(
  query.gr=query.gr,
  loci=skin.eset, 
  path=biocarta)
```
Enrichment result table

```
res.one$result.table
```

All the genes associated with eQTLs covered by the query region

```
res.one$cover.gene
```

####  4.2 peroform query from multiple eQTL sets 
```{r}
#query from one eQTL set.
res.esetlist <- query(
  query.gr=query.gr, 
  loci=eset.list, 
  path=biocarta)  
```

Enrichment result table

```
resultTable(res.esetlist)
```

All the genes associated with eQTLs covered by the query region; names of the list are tissue names from eqtl set list

```
coveredGene(res.esetlist)
```


#### 4.3 parallel query from multiple eQTL sets
```{r}
#query from one eQTL set.
res.paral <- query(
  query.gr=query.gr, 
  loci=eset.list, 
  path=biocarta, 
  parallel=TRUE) 
``` 
Should return the same result as res.esetlist



#### 4.4 obtain tissue enrichment for query regions
In this case, in the generic function `query`, only the argunment
`loci` need to be provided. This will trigger the query of tissue
specificity

```{r}
#query tissue specificity
gr.tissue <- query(query.gr, loci=eset.list)
gr.tissue
```




## 5. Explore query result




#### 5.1 obtain eQTL gene list

```{r}
#all the genes associated with eQTLs covered by the query region
res.one$cover.gene

#all the genes associated with eQTLs covered by the query region; 
#names of the list are tissue names from eqtl set list
coveredGene(res.esetlist)
```





#### 5.2 obtain average tissue degree for each pathway
```{r}
tissue.degree=getTissueDegree(res.esetlist, eset.list)

#check gene-tissue mapping for each gene
head(tissue.degree$gene.tissue.map)

#check degree for each gene
head(tissue.degree$gene.tissue.degree)

#average tissue degree for the input result table
tissue.degree$mean.tissue.degree

#add avg. tissue degree to existing result table
res.tissue=data.frame(resultTable(res.esetlist),
                      t.degree=tissue.degree$mean.tissue.degree)
```





#### 5.3 extract tissue-pathway heatmap
```{r}
#extract tissue-pathway matrix
mat <- getMat(res.esetlist, test.method = "gene")

#plot heatmap
heatmap.para <- getHeatmap(res.esetlist)
```


#### 5.4 extract word cloud from result
```{r}

#plot word cloud
wc <- getWordcloud(res.esetlist)
```

#### 5.5 plot p-value distribution of result
```{r}
#plot p-value distribution of result
pval <- getPval(res.esetlist, test.method="gene")
```


#### 5.6 obtain geneset description from object
```{r}
#obtain geneset description from object
description <- getPathDescription(res.esetlist, biocarta)
head(description)
```




