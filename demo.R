
library(devtools)
install_github("StanleyXu/egquery")
library(egquery)


#show query region
query.gr
#show eqtl set
eqtl.set.list
#show gene set(biocarta)
biocarta

#query
result=query.egset.list(query.gr=query.gr, query.score=NULL, eqtl.set=eqtl.set.list, gene.set=biocarta)
head(result)

#query tissue specificity
gr.tissue=query.tissue(query.gr, eqtl.set.list=eqtl.set.list)

#reorder result by tissue specificity
result2=order.by.tissue(result, tissue.list = rownames(gr.tissue))

#extract tissue/geneset matrix
mat=get.heat.mat(result, test.method = "fisher")

#plot heatmap
plot.heatmap(mat)

#plot word cloud
plot.wordcloud(result)

#plot p-value distribution of result
plot.pval.distribution(result, test.method="fisher")

#obtain geneset description from object
description=get.geneset.description(biocarta, geneset.ids=result$name_pthw)
head(description)

