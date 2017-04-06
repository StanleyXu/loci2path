eqtlSet <- setClass("eqtlSet",
                    slot=c(
                      tissue="character",
                      snp.id="character",
                      snp.gr="GenomicRanges",
                      gene="character"))

geneSet <- setClass("geneSet",
                    slot=c(
                      description="character",
                      gene.set="list" ))
