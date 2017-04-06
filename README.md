# egquery
eQTL-geneset query with ordered regions

data(egquery.data)
#show query region
query.gr
#show eqtl set
eqtl.set.list
#show gene set(biocarta)
biocarta

#query
result=query.egset.list(query.gr=query.gr, query.score=NULL, eqtl.set=eqtl.set.list, gene.set=biocarta, lr=F)
head(result)
#extract tissue/geneset matrix
mat=get.heat.mat(result, test.method = "fisher")

#plot heatmap
export.heatmap(mat, file="test.png")
