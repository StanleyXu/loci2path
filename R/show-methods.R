setMethod(
    "show",
    signature = "eqtlSet",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat(" eQTL collected from tissue:", tissue(object), "\n")
        cat(" number of eQTLs:", length(eqtlId(object)), "\n")
        cat(" number of associated genes:", length(unique(eqtlGene(object))), "\n")
    }
)
setMethod(
    "show",
    signature = "geneSet",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat(" Number of gene sets:", length(geneSetList(object)), "\n")
        if (length(geneSetList(object)) > 0) {
            rr = range(lapply(geneSetList(object), length))
            cat("   ", rr[1], "~", rr[2], " genes within sets\n")
        }
    }
)
